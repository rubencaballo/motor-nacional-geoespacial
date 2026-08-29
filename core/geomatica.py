#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pipeline geomático e hidrológico (D8) por zonas concéntricas (núcleo +
buffers). Calcula área 2D/3D corregida por pendiente, factor de relieve,
elevación, y red hidrográfica por acumulación de flujo D8 (pysheds).

REGLA DE ORO de este módulo (para no repetir los bugs de versiones previas
del script, ya documentados en el repo):
    Cada zona se recorta y se calcula de forma TOTALMENTE AISLADA -- su
    propio array de elevación, su propia máscara de válidos, su propio
    gradiente/pendiente. Nunca se combinan arrays de zonas distintas ni se
    usa la máscara de una zona para indexar el array de otra. Si algún día
    vuelve a salir un factor_relieve < 1.0, es una alarma real de bug --
    por eso NO se aplica ningún piso artificial (max(1.0, ...)) que pueda
    esconderlo.

Qué SÍ calcula (topografía/hidrología, solo con un DEM tipo SRTM):
    - Área 2D, área 3D (corregida por pendiente), factor de relieve.
    - Elevación mín/máx/promedio, rango altitudinal, pendiente promedio,
      % de área con pendiente crítica (>25°).
    - Red hidrográfica por acumulación de flujo D8, con relleno de
      depresiones -- sin esto, la acumulación da puntos dispersos en vez
      de cauces ramificados (ver commits previos donde se intentó
      aproximar con filtros locales y no funcionó).

Qué NO calcula (fuera de alcance de este módulo):
    - AGB (biomasa aérea), CO2, NDWI, GEDI -- eso sale de Earth Engine, no
      de un DEM. Si se necesitan esas métricas por zona, va en un módulo
      aparte que consulte GEE, no se mezcla aquí.

Limitación conocida, documentada -- no oculta:
    Este módulo carga el DEM completo del área solicitada en memoria de una
    sola vez (sin lectura por ventanas/teselas). Es apropiado para sitios
    de tamaño similar a un ANP compacta de hasta cientos de miles de
    hectáreas (probado hasta ~3,000 ha reales). Para ANPs muy grandes y/o
    muy dispersas (ej. con multipolígonos separados por decenas de km),
    el rectángulo envolvente usado para descargar el DEM puede ser mucho
    mayor que el área real y agotar memoria en equipos con poca RAM --
    en ese caso se requiere procesamiento por ventanas, no implementado
    todavía (ver plan de escalamiento del proyecto).
"""

import argparse
import os

import numpy as np

if not hasattr(np, "in1d"):
    np.in1d = np.isin

import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.features import rasterize
import pyproj
from shapely.geometry import Point
from shapely.ops import transform as shp_transform
from scipy.ndimage import gaussian_filter

from config import (
    ZONAS_ANALISIS_M, PERCENTIL_CAUCE_HIDROLOGIA, CARPETA_SRTM, log,
)


# ==============================================================================
# --- FUNCIONES DE BAJO NIVEL (puras: reciben datos, no leen archivos ni red) ---
# ==============================================================================
def recortar_y_calcular(dst_array, meta_utm, geom_poligono):
    """Recorta UN raster (ya en memoria, ya reproyectado) a UN polígono y
    calcula elevación/pendiente de forma completamente aislada. No debe
    recibir ni devolver nada que se mezcle con otro recorte -- ver REGLA DE
    ORO en el docstring del módulo."""
    with rasterio.io.MemoryFile() as memfile:
        with memfile.open(**meta_utm) as tmp:
            tmp.write(dst_array, 1)
            out_image, out_transform = mask(tmp, [geom_poligono], crop=True, nodata=np.nan, filled=True)

    z_raw = out_image[0].astype(np.float64)
    rasterized_mask = rasterize(
        [(geom_poligono, 1)], out_shape=z_raw.shape, transform=out_transform,
        fill=0, default_value=1, dtype=np.uint8,
    )
    valid_mask = (rasterized_mask == 1) & (~np.isnan(z_raw)) & (z_raw > 0)

    pw, ph = abs(out_transform[0]), abs(out_transform[4])

    # Limitación conocida: para polígonos angostos, una fracción alta de
    # píxeles válidos está cerca del borde, y el relleno con la media
    # introduce un sesgo leve ahí. No se corrige aquí -- documentado como
    # margen de incertidumbre (ver docstring del módulo).
    z_grad_base = np.where(valid_mask, z_raw, np.nanmean(z_raw[valid_mask]) if np.any(valid_mask) else 0.0)
    dy, gx = np.gradient(z_grad_base, ph, pw)
    slope_rad = np.arctan(np.sqrt(gx**2 + dy**2))
    slope_deg = np.degrees(slope_rad)

    return z_raw, valid_mask, slope_rad, slope_deg, out_transform, pw, ph


def metricas_de_zona(dst_array, meta_utm, geom_zona, buf_m):
    """Calcula el diccionario de métricas de UNA zona, de forma aislada.
    No tiene efectos secundarios (no imprime, no escribe archivos) para
    poder probarse en --demo sin red ni disco."""
    z, valid, slope_rad, slope_deg, _, pw, ph = recortar_y_calcular(dst_array, meta_utm, geom_zona)

    a2d_ha = geom_zona.area / 10000.0
    valid_pixels = z[valid]
    slope_valid_rad = slope_rad[valid]
    slope_valid_deg = slope_deg[valid]

    a3d_ha = np.nansum((pw * ph) / np.cos(slope_valid_rad)) / 10000.0
    factor_relieve = a3d_ha / a2d_ha if a2d_ha > 0 else float("nan")

    if factor_relieve < 1.0:
        log(f"ALERTA: factor_relieve={factor_relieve:.4f} < 1.0 es matemáticamente "
            f"imposible (1/cos(x) siempre es >=1). Hay un bug -- no usar este número.",
            nivel="WARN")

    return {
        "zona": "nucleo" if buf_m == 0 else f"buffer_{buf_m}m",
        "buffer_m": buf_m,
        "area_2d_ha": round(a2d_ha, 3),
        "area_3d_ha": round(a3d_ha, 3),
        "factor_relieve": round(factor_relieve, 4),
        "elev_min_m": round(float(np.nanmin(valid_pixels)), 1),
        "elev_max_m": round(float(np.nanmax(valid_pixels)), 1),
        "elev_promedio_m": round(float(np.nanmean(valid_pixels)), 1),
        "rango_altitudinal_dH": round(float(np.nanmax(valid_pixels) - np.nanmin(valid_pixels)), 1),
        "pendiente_promedio_deg": round(float(np.nanmean(slope_valid_deg)), 2),
        "pct_area_pendiente_critica": round(float(np.sum(slope_valid_deg > 25.0) / slope_valid_deg.size * 100), 2),
    }


def calcular_metricas_por_zona(geom_utm_nucleo, dst_array, meta_utm, zonas_m):
    """Calcula las métricas de todas las zonas (núcleo + buffers). Devuelve
    un DataFrame -- no escribe nada a disco (eso lo hace el llamador)."""
    filas = []
    for buf_m in zonas_m:
        geom_zona = geom_utm_nucleo.buffer(buf_m) if buf_m > 0 else geom_utm_nucleo
        log(f"Zona buffer={buf_m}m -- calculando de forma aislada...")
        fila = metricas_de_zona(dst_array, meta_utm, geom_zona, buf_m)
        filas.append(fila)
        log(f"  -> area_2d={fila['area_2d_ha']}ha  area_3d={fila['area_3d_ha']}ha  "
            f"factor_relieve={fila['factor_relieve']}  pendiente_prom={fila['pendiente_promedio_deg']}°")
    return pd.DataFrame(filas)


def calcular_hidrologia_d8(dst_array, meta_utm, geom_utm_nucleo, zonas_m, buf_visual_m,
                            utm_crs, percentil_cauce, carpeta_srtm, id_proyecto):
    """Corre D8 (pysheds) sobre la zona más grande y devuelve todo lo
    necesario para graficar: superficie suavizada, máscara de cauces, y a
    qué zona pertenece cada píxel de cauce. Requiere pysheds -- se importa
    aquí (no al inicio del módulo) para que el resto del módulo funcione
    aunque pysheds no esté instalado."""
    from pysheds.grid import Grid

    geom_visual = geom_utm_nucleo.buffer(buf_visual_m) if buf_visual_m > 0 else geom_utm_nucleo
    log(f"Hidrología D8 sobre la zona visual (buffer={buf_visual_m}m)...")

    z_viz, valid_viz, _, _, trans_viz, pw_v, ph_v = recortar_y_calcular(dst_array, meta_utm, geom_visual)
    z_filled_viz = np.where(valid_viz, z_viz, np.nanmean(z_viz[valid_viz]))

    meta_pysheds = {
        "driver": "GTiff", "height": z_filled_viz.shape[0], "width": z_filled_viz.shape[1],
        "count": 1, "dtype": "float32", "crs": utm_crs, "transform": trans_viz, "nodata": -9999,
    }
    temp_raster_path = os.path.join(carpeta_srtm, f"temp_pysheds_{id_proyecto}.tif")
    with rasterio.open(temp_raster_path, "w", **meta_pysheds) as dst:
        dst.write(np.where(valid_viz, z_filled_viz, -9999).astype(np.float32), 1)

    grid = Grid.from_raster(temp_raster_path)
    dem = grid.read_raster(temp_raster_path)
    pit_filled_dem = grid.fill_pits(dem)
    flooded_dem = grid.fill_depressions(pit_filled_dem)
    inflated_dem = grid.resolve_flats(flooded_dem)
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
    flow_dir = grid.flowdir(inflated_dem, dirmap=dirmap)
    acc = grid.accumulation(flow_dir, dirmap=dirmap)

    acc_arr = np.asarray(acc)
    threshold = np.nanpercentile(acc_arr[valid_viz & (acc_arr >= 0)], percentil_cauce) if np.any(valid_viz) else 10.0
    stream_mask = (acc_arr > threshold) & valid_viz
    log(f"  -> {stream_mask.sum()} píxeles clasificados como cauce (percentil {percentil_cauce})")

    zona_de_pixel = np.full(stream_mask.shape, "fuera", dtype=object)
    for buf_m in sorted(zonas_m):
        geom_zona = geom_utm_nucleo.buffer(buf_m) if buf_m > 0 else geom_utm_nucleo
        mascara_zona = rasterize(
            [(geom_zona, 1)], out_shape=stream_mask.shape, transform=trans_viz,
            fill=0, default_value=1, dtype=np.uint8,
        ).astype(bool)
        etiqueta = "nucleo" if buf_m == 0 else f"buffer_{buf_m}m"
        zona_de_pixel[mascara_zona & (zona_de_pixel == "fuera")] = etiqueta

    Z_smooth = gaussian_filter(z_filled_viz, sigma=1.0)
    Z_smooth[~valid_viz] = np.nan

    # Elevación SIN suavizar (para el hover -- el suavizado es solo para
    # que la superficie se vea bien visualmente, pero distorsiona el valor
    # real de elevación en un punto específico, sobre todo en picos/crestas
    # angostas. El hover debe mostrar el dato real, no el visual.
    Z_raw = np.where(valid_viz, z_viz, np.nan)

    return {
        "Z_smooth": Z_smooth, "Z_raw": Z_raw, "stream_mask": stream_mask, "zona_de_pixel": zona_de_pixel,
        "pw_v": pw_v, "ph_v": ph_v, "transform": trans_viz,
    }


def calcular_grid_latlon(transform, utm_crs, rows, cols):
    """Calcula lat/lon (WGS84) para cada píxel de la malla (centro de
    píxel), para poder mostrar coordenadas reales al pasar el mouse sobre
    el mapa 3D -- en vez de solo kilómetros locales sin referencia
    geográfica, que no sirven para ubicar el punto en un mapa real."""
    import pyproj
    proj_a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)
    cols_idx, rows_idx = np.meshgrid(np.arange(cols) + 0.5, np.arange(rows) + 0.5)
    xs, ys = transform * (cols_idx, rows_idx)
    lons, lats = proj_a_wgs84.transform(xs, ys)
    return lats, lons


def generar_mapa_3d(hidrologia, id_proyecto, html_path, subtitulo=None, utm_crs=None):
    """Genera el HTML interactivo con Plotly a partir del resultado de
    calcular_hidrologia_d8(). Import de plotly aquí mismo, mismo criterio
    que con pysheds. `subtitulo` es opcional y agnóstico de su origen (por
    ejemplo, core/carbono.py lo usa para añadir CO2e por zona) -- este
    módulo no sabe ni le importa de dónde viene el texto, solo lo agrega
    como segunda línea del título.

    `utm_crs`: si se da, se calcula lat/lon real por punto y se muestra en
    el hover (además de los km locales) -- así se puede ubicar cualquier
    punto del mapa de inmediato (copiar/pegar a Google Maps), sin tener
    que traducir manualmente coordenadas locales."""
    import plotly.graph_objects as go

    Z_smooth = hidrologia["Z_smooth"]
    Z_raw = hidrologia.get("Z_raw", Z_smooth)  # por compatibilidad si algún llamador viejo no lo trae
    stream_mask = hidrologia["stream_mask"]
    zona_de_pixel = hidrologia["zona_de_pixel"]
    pw_v, ph_v = hidrologia["pw_v"], hidrologia["ph_v"]
    transform = hidrologia.get("transform")

    rows, cols = Z_smooth.shape
    X, Y = np.meshgrid(np.arange(cols) * pw_v / 1000.0, np.flipud(np.arange(rows) * ph_v / 1000.0))

    lat_grid = lon_grid = None
    if utm_crs and transform is not None:
        lat_grid, lon_grid = calcular_grid_latlon(transform, utm_crs, rows, cols)

    fig = go.Figure()
    if lat_grid is not None:
        # customdata: [altitud REAL (sin suavizar), lat, lon] -- el %{z} del
        # trace sigue siendo la superficie suavizada (para que se vea bien),
        # pero el hover muestra la altitud real del pixel, no la suavizada.
        customdata_superficie = np.dstack([Z_raw, lat_grid, lon_grid])
        fig.add_trace(go.Surface(
            z=Z_smooth, x=X, y=Y, colorscale="Earth", connectgaps=False, name="Terreno 3D",
            customdata=customdata_superficie,
            hovertemplate="Altitud: %{customdata[0]:.0f} msnm<br>Lat: %{customdata[1]:.5f}<br>Lon: %{customdata[2]:.5f}<extra></extra>",
        ))
    else:
        fig.add_trace(go.Surface(z=Z_smooth, x=X, y=Y, colorscale="Earth", connectgaps=False, name="Terreno 3D"))

    colores_zona = {"nucleo": "red", "buffer_500m": "orange", "buffer_1000m": "gold"}
    river_y, river_x = np.where(stream_mask)
    if len(river_x) > 0:
        for etiqueta, color in colores_zona.items():
            sel = zona_de_pixel[river_y, river_x] == etiqueta
            if not np.any(sel):
                continue
            # BUG CORREGIDO: antes se indexaba X/Y/Z con "river_y_plotly"
            # (un segundo flip sobre river_y), cuando X y Y YA tienen el
            # flip codificado en su propia construcción (arriba). Usar
            # river_y_plotly aquí colocaba el punto en la fila espejo
            # equivocada -- posición Y y elevación de OTRO pixel del
            # raster, no del cauce real. river_y/river_x (los índices
            # reales del pixel) son los correctos para indexar X, Y, Z.
            altitudes_reales = Z_raw[river_y[sel], river_x[sel]]
            kwargs_extra = {}
            if lat_grid is not None:
                customdata_pts = np.column_stack([
                    altitudes_reales,
                    lat_grid[river_y[sel], river_x[sel]],
                    lon_grid[river_y[sel], river_x[sel]],
                ])
                kwargs_extra = dict(
                    customdata=customdata_pts,
                    hovertemplate=f"Cauce ({etiqueta})<br>Altitud: %{{customdata[0]:.0f}} msnm<br>"
                                  f"Lat: %{{customdata[1]:.5f}}<br>Lon: %{{customdata[2]:.5f}}<extra></extra>",
                )
            fig.add_trace(go.Scatter3d(
                x=X[river_y[sel], river_x[sel]],
                y=Y[river_y[sel], river_x[sel]],
                z=Z_smooth[river_y[sel], river_x[sel]] + 12,
                mode="markers", marker=dict(size=2.2, color=color, opacity=0.9),
                name=f"Cauce en {etiqueta}", **kwargs_extra,
            ))

    titulo = (f"Modelo Hidrológico D8 por zonas - {id_proyecto} "
              f"(núcleo=rojo, buffer 500m=naranja, buffer 1000m=dorado)")
    if subtitulo:
        titulo += f"<br><sub>{subtitulo}</sub>"

    fig.update_layout(
        title=titulo,
        scene=dict(xaxis_title="Este [km]", yaxis_title="Norte [km]", zaxis_title="Altitud [msnm]",
                   aspectmode="manual", aspectratio=dict(x=1, y=1, z=0.7)),
        autosize=True,
    )
    fig.write_html(html_path)
    return html_path


def cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm=None):
    """Descarga/lee el SRTM del sitio y lo reproyecta a UTM. Extraído de
    procesar_sitio_real() para poder reutilizarse desde otros módulos (ej.
    core/carbono.py) sin duplicar este bloque ni descargar el DEM dos
    veces -- si el .tif ya existe en carpeta_srtm, se reusa tal cual.
    Devuelve (geom_utm_nucleo, dst_array, meta_utm, utm_crs)."""
    import elevation

    carpeta_srtm = carpeta_srtm or CARPETA_SRTM
    os.makedirs(carpeta_srtm, exist_ok=True)

    if not os.path.exists(geojson_path):
        raise FileNotFoundError(f"No se encontró el archivo GeoJSON en: {geojson_path}")

    gdf = gpd.read_file(geojson_path)
    geom_wgs84 = gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all") else gdf.geometry.unary_union

    bounds = geom_wgs84.bounds
    margen_grados = 0.02 + (max(zonas_m) / 111000.0)
    west, south, east, north = (
        bounds[0] - margen_grados, bounds[1] - margen_grados,
        bounds[2] + margen_grados, bounds[3] + margen_grados,
    )
    id_proyecto_tif = os.path.splitext(os.path.basename(geojson_path))[0]
    tif_path = os.path.join(carpeta_srtm, f"srtm_{id_proyecto_tif}.tif")
    if not os.path.exists(tif_path):
        log("Descargando modelo de elevación SRTM...")
        elevation.clip(bounds=(west, south, east, north), output=os.path.abspath(tif_path))

    centroid = geom_wgs84.centroid
    utm_zone = int((centroid.x + 180) / 6) + 1
    utm_crs = f"EPSG:326{utm_zone}"
    log(f"Zona UTM: {utm_crs}")

    project = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True).transform
    geom_utm_nucleo = shp_transform(project, geom_wgs84)

    with rasterio.open(tif_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, utm_crs, src.width, src.height, *src.bounds
        )
        dst_array = np.empty((height, width), dtype=np.float32)
        reproject(
            source=rasterio.band(src, 1), destination=dst_array,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=transform, dst_crs=utm_crs, resampling=Resampling.bilinear,
        )

    meta_utm = {
        "driver": "GTiff", "height": height, "width": width,
        "count": 1, "dtype": "float32", "crs": utm_crs, "transform": transform,
    }
    return geom_utm_nucleo, dst_array, meta_utm, utm_crs


# ==============================================================================
# --- ORQUESTADOR CON DATOS REALES (SRTM real, requiere red la primera vez) ---
# ==============================================================================
def procesar_sitio_real(geojson_path, id_proyecto, zonas_m=None, percentil_cauce=None,
                         carpeta_salida=None, carpeta_srtm=None):
    """Pipeline completo con datos reales: descarga/lee SRTM, calcula
    métricas por zona, corre D8 y genera el mapa 3D. Devuelve
    (df_zonas, csv_path, html_path)."""
    zonas_m = zonas_m if zonas_m is not None else ZONAS_ANALISIS_M
    percentil_cauce = percentil_cauce if percentil_cauce is not None else PERCENTIL_CAUCE_HIDROLOGIA
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    carpeta_srtm = carpeta_srtm or CARPETA_SRTM
    os.makedirs(carpeta_salida, exist_ok=True)

    geom_utm_nucleo, dst_array, meta_utm, utm_crs = cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)

    df_zonas = calcular_metricas_por_zona(geom_utm_nucleo, dst_array, meta_utm, zonas_m)
    csv_path = os.path.join(carpeta_salida, f"features_por_zona_{id_proyecto.lower()}.csv")
    df_zonas.to_csv(csv_path, index=False)
    log(f"CSV por zona guardado en: {csv_path}")

    hidrologia = calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m),
        utm_crs, percentil_cauce, carpeta_srtm, id_proyecto,
    )
    html_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_3d_zonas.html")
    generar_mapa_3d(hidrologia, id_proyecto, html_path, utm_crs=utm_crs)
    log(f"Mapa 3D: {html_path}")

    return df_zonas, csv_path, html_path


# ==============================================================================
# --- MODO DEMO: DEM sintético en memoria, sin red, sin depender de archivos
#     del usuario -- mismo espíritu que core/motor.py --demo ---
# ==============================================================================
def _dem_sintetico(size=220, resolucion_m=30.0, semilla=42):
    """Genera un DEM sintético (una loma con una cañada), determinista, para
    poder probar la lógica de zonas + D8 sin descargar nada de internet."""
    rng = np.random.default_rng(semilla)
    y, x = np.mgrid[0:size, 0:size]
    cx, cy = size / 2, size / 2

    # Plano inclinado + un cono central (loma) + una cañada rasgada + ruido leve
    plano = 1000 + 0.4 * x
    loma = 300 * np.exp(-(((x - cx) ** 2 + (y - cy) ** 2) / (2 * (size / 4) ** 2)))
    canada = -150 * np.exp(-(((x - cx) ** 2) / (2 * (size / 12) ** 2)))
    ruido = rng.normal(0, 3, size=(size, size))
    Z = (plano + loma + canada + ruido).astype(np.float32)

    transform = rasterio.transform.from_origin(0, size * resolucion_m, resolucion_m, resolucion_m)
    meta_utm = {
        "driver": "GTiff", "height": size, "width": size,
        "count": 1, "dtype": "float32", "crs": "EPSG:32614", "transform": transform,
    }
    # Polígono núcleo: un círculo en el centro de la malla, radio moderado
    centro_x_m = (size / 2) * resolucion_m
    centro_y_m = (size / 2) * resolucion_m
    geom_utm_nucleo = Point(centro_x_m, centro_y_m).buffer(1200)  # ~452 ha aprox
    return Z, meta_utm, geom_utm_nucleo, "EPSG:32614"


def demo():
    """Corre el pipeline completo (zonas + D8 + mapa 3D) sobre un DEM
    sintético en memoria. No descarga nada, no necesita ningún geojson del
    usuario -- sirve para probar que la lógica sigue intacta después de
    cualquier cambio al módulo, igual que --demo en core/motor.py."""
    log("=== core.geomatica --demo (DEM sintético, sin red) ===")
    dst_array, meta_utm, geom_utm_nucleo, utm_crs = _dem_sintetico()
    zonas_m = [0, 300, 600]
    id_proyecto = "DEMO_SINTETICO"

    df_zonas = calcular_metricas_por_zona(geom_utm_nucleo, dst_array, meta_utm, zonas_m)
    print("\n--- Métricas por zona (demo) ---")
    print(df_zonas.to_string(index=False))

    try:
        carpeta_tmp = os.path.expanduser("~/resultados_demo_geomatica")
        os.makedirs(carpeta_tmp, exist_ok=True)
        hidrologia = calcular_hidrologia_d8(
            dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m),
            utm_crs, PERCENTIL_CAUCE_HIDROLOGIA, carpeta_tmp, id_proyecto,
        )
        html_path = os.path.join(carpeta_tmp, f"{id_proyecto.lower()}_3d_zonas.html")
        generar_mapa_3d(hidrologia, id_proyecto, html_path, utm_crs=utm_crs)
        log(f"Mapa 3D demo generado en: {html_path}")
    except ImportError as e:
        log(f"pysheds/plotly no instalado, se omite hidrología en el demo: {e}", nivel="WARN")

    return df_zonas


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(description="Módulo geomático e hidrológico (D8) por zonas -- Motor Nacional")
    ap.add_argument("--demo", action="store_true", help="Corre con un DEM sintético, sin red ni archivos del usuario")
    ap.add_argument("--geojson", type=str, help="Ruta al GeoJSON del polígono núcleo")
    ap.add_argument("--id-proyecto", type=str, help="Nombre identificador del sitio (para nombres de archivo)")
    ap.add_argument("--zonas", type=str, default=None,
                     help="Buffers en metros separados por coma, ej. '0,500,1000' (default: config.ZONAS_ANALISIS_M)")
    ap.add_argument("--percentil-cauce", type=float, default=None,
                     help="Percentil de acumulación de flujo para declarar cauce (default: config.PERCENTIL_CAUCE_HIDROLOGIA)")
    ap.add_argument("--carpeta-salida", type=str, default=None)
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not args.geojson or not args.id_proyecto:
        ap.error("--geojson y --id-proyecto son obligatorios fuera de --demo (o usa --demo para probar sin datos reales)")

    zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
    procesar_sitio_real(
        geojson_path=args.geojson,
        id_proyecto=args.id_proyecto,
        zonas_m=zonas_m,
        percentil_cauce=args.percentil_cauce,
        carpeta_salida=args.carpeta_salida,
    )


if __name__ == "__main__":
    main()

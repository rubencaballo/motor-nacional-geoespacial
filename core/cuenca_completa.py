#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Delimitación de cuenca hidrológica REAL entre un sitio de origen (ej.
Cofre de Perote) y un punto de salida río abajo (ej. la Cascada de Texolo),
para verificar con matemática -- no solo por el nombre de un río en algún
dataset -- si un área específica realmente alimenta hidrológicamente otro
punto específico, kilómetros más abajo.

POR QUÉ ESTE MÓDULO ES DISTINTO DE core/geomatica.py:
    geomatica.py corre D8 sobre un buffer CHICO alrededor de UN sitio (1-2km)
    y clasifica cauces por un umbral de percentil de acumulación de flujo --
    eso responde "¿por dónde escurre el agua DENTRO de este sitio?". Este
    módulo responde una pregunta distinta: "¿toda el agua que cae en esta
    zona A, termina saliendo por este punto B, kilómetros más abajo?" -- eso
    requiere un DEM que cubra TODO el corredor entre A y B (no solo un
    buffer chico de cada uno por separado), y una función distinta de
    pysheds: catchment() en vez de threshold de percentil.

CÓMO FUNCIONA LA DELIMITACIÓN DE CUENCA (catchment):
    1) Se calcula un bounding box que envuelve TANTO el polígono de origen
       como el punto de salida, con margen -- si el margen es insuficiente,
       la cuenca real puede "cortarse" en el borde del DEM y salir
       incompleta (ver MARGEN_CORREDOR_GRADOS en config.py).
    2) Se corre el mismo preprocesamiento D8 que geomatica.py
       (fill_pits -> fill_depressions -> resolve_flats -> flowdir), pero
       sobre esta malla grande del corredor, no sobre un buffer chico.
    3) El punto de salida dado (coordenadas aproximadas, ej. leídas de un
       mapa) casi nunca cae exacto sobre la celda de mayor acumulación de
       flujo -- por precisión de coordenadas o por la resolución del pixel
       (30m de SRTM). Si no se ajusta, catchment() puede delimitar una
       cuenca chueca, incompleta o vacía. Por eso se usa snap_to_mask()
       ANTES de catchment(): ajusta el punto dado a la celda de cauce real
       (alta acumulación de flujo) más cercana. La distancia de ese ajuste
       se reporta explícita -- si es grande (cientos de metros), puede
       significar que se ajustó al cauce equivocado (dos ríos cercanos) y
       hay que revisarlo a mano.
    4) catchment(x, y, fdir, ...) devuelve una máscara booleana: TODAS las
       celdas del DEM que drenan hacia ese punto de salida ajustado -- la
       cuenca completa, aguas arriba.
    5) Verificación: se compara esa máscara contra el polígono del sitio de
       origen (rasterizado sobre la misma malla). El % de traslape es la
       respuesta real a "¿esta zona alimenta ese punto?" -- no una
       suposición por el nombre del río en INEGI.

LIMITACIONES (documentadas, no escondidas):
    - Sigue dependiendo de SRTM de 30m -- la resolución fina de cauces
      angostos en cabeceras puede perderse, igual que en geomatica.py.
    - snap_to_mask() necesita un umbral de acumulación de flujo (percentil)
      para decidir qué cuenta como "cauce real" al buscar la celda más
      cercana -- un umbral mal elegido puede ajustar el punto a un cauce
      equivocado si hay dos ríos cercanos. Por eso SIEMPRE se reporta la
      distancia del ajuste, para que quede explícito qué tanto se movió el
      punto dado respecto al que se usó realmente.
    - Es hidrología de superficie pura, basada en pendiente -- no corrige
      por sumideros/karst ni por trasvases humanos (canales, acueductos,
      tomas de agua).
    - El corredor completo (decenas de km) es mucho más grande que el
      buffer de un solo sitio -- la descarga y el preprocesamiento tardan
      más y usan más memoria, aunque siguen siendo manejables en una
      computadora normal (del orden de 1 millón de celdas para ~30km a
      resolución SRTM).
"""

import argparse
import os

import numpy as np

if not hasattr(np, "in1d"):
    np.in1d = np.isin

import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.features import rasterize
import pyproj
from shapely.geometry import Point
from shapely.ops import transform as shp_transform

from config import (
    CARPETA_SRTM, MARGEN_CORREDOR_GRADOS, PERCENTIL_SNAP_CAUCE, TOLERANCIA_VALIDACION_HIDRO_M, log,
)


# ==============================================================================
# --- BBOX Y DEM DEL CORREDOR COMPLETO (no un buffer chico de un solo sitio) ---
# ==============================================================================
def calcular_bbox_corredor(geom_wgs84_origen, punto_salida_wgs84, margen_grados=None):
    """Calcula el bounding box (WGS84) que envuelve TANTO el polígono de
    origen como el punto de salida, con margen -- si el margen es chico, la
    cuenca real delimitada más adelante puede cortarse en el borde del DEM
    y salir incompleta. `punto_salida_wgs84` es (lon, lat), mismo criterio
    always_xy que el resto del proyecto."""
    margen_grados = margen_grados if margen_grados is not None else MARGEN_CORREDOR_GRADOS
    bounds_origen = geom_wgs84_origen.bounds  # (minx, miny, maxx, maxy) = (min_lon, min_lat, max_lon, max_lat)
    px, py = punto_salida_wgs84
    west = min(bounds_origen[0], px) - margen_grados
    south = min(bounds_origen[1], py) - margen_grados
    east = max(bounds_origen[2], px) + margen_grados
    north = max(bounds_origen[3], py) + margen_grados
    return west, south, east, north


def cargar_dem_corredor(bbox_wgs84, id_corredor, carpeta_srtm=None):
    """Descarga/lee el SRTM del CORREDOR completo (bbox grande, no un
    buffer chico de un solo sitio) y lo reproyecta a UTM. Mismo mecanismo
    que geomatica.cargar_dem_utm() (mismo elevation.clip(), mismo caché por
    archivo en disco) pero recibe el bbox directamente en vez de derivarlo
    de un solo geojson -- así puede cubrir dos sitios lejanos entre sí."""
    import elevation

    carpeta_srtm = carpeta_srtm or CARPETA_SRTM
    os.makedirs(carpeta_srtm, exist_ok=True)

    west, south, east, north = bbox_wgs84
    tif_path = os.path.join(carpeta_srtm, f"srtm_corredor_{id_corredor.lower()}.tif")
    if not os.path.exists(tif_path):
        ancho_km, alto_km = (east - west) * 111.0, (north - south) * 111.0
        log(f"Descargando SRTM del corredor completo (~{ancho_km:.0f}km x {alto_km:.0f}km)...")
        elevation.clip(bounds=(west, south, east, north), output=os.path.abspath(tif_path))
    else:
        log(f"Reusando SRTM del corredor ya descargado: {tif_path}")

    centroid_lon, centroid_lat = (west + east) / 2, (south + north) / 2
    utm_zone = int((centroid_lon + 180) / 6) + 1
    utm_crs = f"EPSG:326{utm_zone}"
    log(f"Zona UTM del corredor: {utm_crs}")

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
    log(f"DEM del corredor: {width}x{height} píxeles ({width * height:,} celdas)")
    return dst_array, meta_utm, utm_crs


# ==============================================================================
# --- DELIMITACIÓN DE CUENCA (pysheds catchment(), no threshold de percentil) ---
# ==============================================================================
def delimitar_cuenca(dst_array, meta_utm, punto_salida_utm, utm_crs, carpeta_srtm, id_corredor,
                      percentil_snap=None):
    """Corre D8 (pysheds) sobre el DEM COMPLETO del corredor y delimita la
    cuenca real que drena hacia `punto_salida_utm` (x, y en la proyección
    utm_crs) usando catchment(). A diferencia de
    geomatica.calcular_hidrologia_d8() (que clasifica cauces por un umbral
    de percentil dentro de un buffer chico), aquí el resultado es una
    cuenca COMPLETA delimitada matemáticamente desde un punto de salida
    específico -- responde "¿qué área entera drena hacia AQUÍ?", no "¿por
    dónde pasan los cauces en esta zona chica?".

    `percentil_snap`: ver snap_to_mask() en el docstring del módulo --
    ajusta el punto de salida dado a la celda de cauce real más cercana
    antes de delimitar. Devuelve un dict con la máscara de la cuenca, la
    elevación real, y el detalle del ajuste (snap) para que quede
    explícito qué tanto se movió el punto dado."""
    from pysheds.grid import Grid

    percentil_snap = percentil_snap if percentil_snap is not None else PERCENTIL_SNAP_CAUCE

    valid = ~np.isnan(dst_array) & (dst_array > -1000)  # descarta nodata / valores absurdos de borde
    z_filled = np.where(valid, dst_array, np.nanmean(dst_array[valid]) if np.any(valid) else 0.0)

    meta_pysheds = {
        "driver": "GTiff", "height": z_filled.shape[0], "width": z_filled.shape[1],
        "count": 1, "dtype": "float32", "crs": utm_crs, "transform": meta_utm["transform"], "nodata": -9999,
    }
    temp_raster_path = os.path.join(carpeta_srtm, f"temp_pysheds_corredor_{id_corredor}.tif")
    with rasterio.open(temp_raster_path, "w", **meta_pysheds) as dst:
        dst.write(np.where(valid, z_filled, -9999).astype(np.float32), 1)

    grid = Grid.from_raster(temp_raster_path)
    dem = grid.read_raster(temp_raster_path)
    pit_filled_dem = grid.fill_pits(dem)
    flooded_dem = grid.fill_depressions(pit_filled_dem)
    inflated_dem = grid.resolve_flats(flooded_dem)
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
    flow_dir = grid.flowdir(inflated_dem, dirmap=dirmap)
    acc = grid.accumulation(flow_dir, dirmap=dirmap)
    acc_arr = np.asarray(acc)

    umbral_snap = np.nanpercentile(acc_arr[valid & (acc_arr >= 0)], percentil_snap) if np.any(valid) else 10.0
    # OJO: snap_to_mask() exige que `mask` sea un pysheds.sview.Raster, no un
    # numpy.ndarray plano -- por eso la comparación se hace sobre `acc`
    # (todavía Raster) y no sobre `acc_arr` (ya convertido a numpy puro),
    # que es el que se usa para el resto de operaciones numéricas normales.
    mascara_cauces = (acc > umbral_snap) & valid

    x_dado, y_dado = punto_salida_utm
    x_snap, y_snap = grid.snap_to_mask(mascara_cauces, (x_dado, y_dado))
    distancia_snap_m = float(np.hypot(x_snap - x_dado, y_snap - y_dado))
    log(f"Punto de salida ajustado (snap) {distancia_snap_m:.0f}m al cauce real más cercano "
        f"(percentil {percentil_snap} de acumulación de flujo).")
    if distancia_snap_m > 500:
        log(f"OJO: el ajuste fue de {distancia_snap_m:.0f}m -- revisa que no haya caído en un cauce "
            f"equivocado (puede pasar si hay dos ríos cercanos, o si --percentil-snap está mal calibrado).",
            nivel="WARN")

    catchment = grid.catchment(x=x_snap, y=y_snap, fdir=flow_dir, dirmap=dirmap, xytype="coordinate")
    catchment_mask = np.asarray(catchment).astype(bool)
    log(f"  -> Cuenca delimitada: {catchment_mask.sum():,} celdas dentro (de {catchment_mask.size:,} totales del corredor).")

    # Cauces REALES dentro de la cuenca ya delimitada -- reusa la misma
    # `mascara_cauces` que ya se calculó para el snap (arriba), nada más
    # restringida al área de la cuenca. Esto es lo que se dibuja en
    # generar_mapa_3d_corredor() como "el río corriendo del polígono al
    # punto de salida" -- no se vuelve a correr pysheds para esto.
    stream_mask_corredor = np.asarray(mascara_cauces).astype(bool) & catchment_mask
    log(f"  -> {stream_mask_corredor.sum():,} píxeles de cauce dentro de la cuenca delimitada.")

    pw, ph = abs(meta_utm["transform"][0]), abs(meta_utm["transform"][4])

    return {
        "catchment_mask": catchment_mask,
        "stream_mask_corredor": stream_mask_corredor,
        "transform": meta_utm["transform"],
        "Z_raw": np.where(valid, dst_array, np.nan),
        "pw": pw, "ph": ph,
        "punto_dado_utm": (float(x_dado), float(y_dado)),
        "punto_snap_utm": (float(x_snap), float(y_snap)),
        "distancia_snap_m": distancia_snap_m,
        "umbral_snap_acumulacion": float(umbral_snap),
    }


def verificar_traslape(catchment_mask, transform, geom_utm_origen, pw, ph):
    """Compara el polígono de origen (UTM) contra la máscara de cuenca ya
    delimitada -- responde con números reales qué fracción del origen
    efectivamente escurre hacia el punto de salida, en vez de asumirlo por
    el nombre de un río en algún dataset."""
    mascara_origen = rasterize(
        [(geom_utm_origen, 1)], out_shape=catchment_mask.shape, transform=transform,
        fill=0, default_value=1, dtype=np.uint8,
    ).astype(bool)

    area_pixel_ha = (pw * ph) / 10000.0
    area_origen_total_ha = float(mascara_origen.sum() * area_pixel_ha)
    area_origen_dentro_ha = float((mascara_origen & catchment_mask).sum() * area_pixel_ha)
    area_cuenca_total_ha = float(catchment_mask.sum() * area_pixel_ha)
    pct_traslape = (area_origen_dentro_ha / area_origen_total_ha * 100.0) if area_origen_total_ha > 0 else 0.0

    return {
        "area_origen_total_ha": area_origen_total_ha,
        "area_origen_dentro_ha": area_origen_dentro_ha,
        "area_cuenca_total_ha": area_cuenca_total_ha,
        "pct_traslape": pct_traslape,
    }


def validar_cuenca_contra_inegi(cuenca, bbox_wgs84, utm_crs, shapefile_inegi_path, layer=None, tolerancia_m=None):
    """Compara los cauces D8 DENTRO de la cuenca ya delimitada
    (`stream_mask_corredor`) contra la Red Hidrográfica oficial de INEGI --
    mismo mecanismo que core/validacion_hidrologica.py (vectorizar a
    puntos + KD-tree + % dentro de tolerancia), REUTILIZADO aquí, no
    reescrito -- para no arriesgarse a que la lógica de comparación
    diverja entre los dos módulos.

    IMPORTANTE: valida solo los cauces DENTRO de la cuenca ya delimitada,
    no el corredor completo -- lo relevante aquí es si la porción de red
    que sí importa para esta conexión hidrológica específica (origen ->
    punto de salida) coincide con INEGI, no el corredor entero (que
    incluye zonas fuera de la cuenca, irrelevantes para esta pregunta)."""
    from core.validacion_hidrologica import _muestrear_lineas_a_puntos, validar_zona

    tolerancia_m = tolerancia_m if tolerancia_m is not None else TOLERANCIA_VALIDACION_HIDRO_M
    layer = layer or "corriente_ag_l"

    if not os.path.exists(shapefile_inegi_path):
        raise FileNotFoundError(f"No se encontró el shapefile de INEGI en: {shapefile_inegi_path}")

    log(f"Cargando red hidrográfica oficial de INEGI para el corredor (capa: {layer})...")
    gdf_inegi = gpd.read_file(shapefile_inegi_path, layer=layer, bbox=bbox_wgs84)
    log(f"  -> {len(gdf_inegi)} elementos cargados dentro del bbox del corredor (de todo el archivo, sin cargarlo completo)")
    gdf_inegi_utm = gdf_inegi.to_crs(utm_crs)
    puntos_ref = _muestrear_lineas_a_puntos(gdf_inegi_utm)
    log(f"  -> {len(puntos_ref)} puntos de referencia")

    transform = cuenca["transform"]
    river_y, river_x = np.where(cuenca["stream_mask_corredor"])
    if len(river_x) > 0:
        xs_utm, ys_utm = transform * (river_x + 0.5, river_y + 0.5)
        puntos_d8 = np.column_stack([xs_utm, ys_utm])
    else:
        puntos_d8 = np.empty((0, 2))

    resultado = validar_zona(puntos_d8, puntos_ref, tolerancia_m)
    log(f"  -> Cauces dentro de la cuenca vs INEGI: {resultado['pct_dentro_tolerancia']}% dentro de "
        f"{tolerancia_m}m (dist. promedio: {resultado['distancia_promedio_m']}m, "
        f"{resultado['n_puntos_d8']} píxeles de cauce comparados).")

    resultado["river_y"], resultado["river_x"] = river_y, river_x
    resultado["puntos_ref_utm"] = puntos_ref
    resultado["tolerancia_m"] = tolerancia_m
    return resultado


# ==============================================================================
# --- MAPA 2D EN PLANTA (vista cenital real, sin rotación -- ver conversación
#     sobre paralaje en geomatica.generar_mapa_3d: aquí no aplica porque no
#     hay cámara inclinada, es un mapa plano de verdad) ---
# ==============================================================================
def _utm_a_km_local(xs_utm, ys_utm, transform, rows, pw, ph):
    """Convierte coordenadas UTM (easting, northing) a las mismas 'km
    locales' que usa la malla del mapa (mismo criterio que
    geomatica.generar_mapa_3d: x_km = col*pw/1000, y_km = (rows-1-row)*ph/1000),
    para que un polígono o punto en UTM quede alineado exacto con el
    raster de fondo, sin tener que reconstruir la conversión a mano cada vez."""
    xs_utm = np.atleast_1d(np.asarray(xs_utm, dtype=np.float64))
    ys_utm = np.atleast_1d(np.asarray(ys_utm, dtype=np.float64))
    cols_f, rows_f = ~transform * (xs_utm, ys_utm)
    x_km = cols_f * pw / 1000.0
    y_km = (rows - 1 - rows_f) * ph / 1000.0
    return x_km, y_km


def generar_mapa_corredor(cuenca, geom_utm_origen, id_corredor, html_path, utm_crs=None, subtitulo=None):
    """Mapa 2D en planta de la cuenca completa delimitada: terreno de
    fondo, área de la cuenca resaltada, contorno del polígono de origen, y
    el punto de salida (dado y ya ajustado/snap). A propósito NO es un mapa
    3D como geomatica.generar_mapa_3d() -- a la escala de un corredor de
    decenas de km, una malla 3D sería pesada y, más importante, reintroduce
    el problema de paralaje de una vista rotada que ya discutimos. Un mapa
    2D real en planta no tiene ese problema."""
    import plotly.graph_objects as go
    from core.geomatica import calcular_grid_latlon

    Z_raw = cuenca["Z_raw"]
    catchment_mask = cuenca["catchment_mask"]
    transform = cuenca["transform"]
    pw, ph = cuenca["pw"], cuenca["ph"]
    rows, cols = Z_raw.shape

    x_km = np.arange(cols) * pw / 1000.0
    y_km = np.flipud(np.arange(rows) * ph / 1000.0)

    lat_grid = lon_grid = None
    if utm_crs:
        lat_grid, lon_grid = calcular_grid_latlon(transform, utm_crs, rows, cols)

    fig = go.Figure()

    customdata_terreno = np.dstack([Z_raw, lat_grid, lon_grid]) if lat_grid is not None else None
    fig.add_trace(go.Heatmap(
        z=Z_raw, x=x_km, y=y_km, colorscale="Earth", showscale=True,
        colorbar=dict(title="Altitud [msnm]", x=1.02), name="Terreno",
        customdata=customdata_terreno,
        hovertemplate=("Altitud: %{customdata[0]:.0f} msnm<br>Lat: %{customdata[1]:.5f}<br>"
                        "Lon: %{customdata[2]:.5f}<extra></extra>") if customdata_terreno is not None
                       else "Altitud: %{z:.0f} msnm<extra></extra>",
    ))

    fig.add_trace(go.Heatmap(
        z=np.where(catchment_mask, 1.0, np.nan), x=x_km, y=y_km,
        colorscale=[[0, "rgba(30,144,255,0.4)"], [1, "rgba(30,144,255,0.4)"]],
        showscale=False, hoverinfo="skip", name="Cuenca delimitada",
    ))

    xs_origen, ys_origen = geom_utm_origen.exterior.coords.xy
    x_km_origen, y_km_origen = _utm_a_km_local(xs_origen, ys_origen, transform, rows, pw, ph)
    fig.add_trace(go.Scatter(
        x=x_km_origen, y=y_km_origen, mode="lines", line=dict(color="red", width=3),
        name="Polígono de origen", hoverinfo="name",
    ))

    x_dado_km, y_dado_km = _utm_a_km_local(
        [cuenca["punto_dado_utm"][0]], [cuenca["punto_dado_utm"][1]], transform, rows, pw, ph)
    x_snap_km, y_snap_km = _utm_a_km_local(
        [cuenca["punto_snap_utm"][0]], [cuenca["punto_snap_utm"][1]], transform, rows, pw, ph)
    fig.add_trace(go.Scatter(
        x=x_dado_km, y=y_dado_km, mode="markers", marker=dict(size=10, color="black", symbol="x"),
        name="Punto de salida (dado)",
        hovertext=[f"Punto dado<br>Ajustado {cuenca['distancia_snap_m']:.0f}m al cauce real"], hoverinfo="text",
    ))
    fig.add_trace(go.Scatter(
        x=x_snap_km, y=y_snap_km, mode="markers",
        marker=dict(size=13, color="lime", symbol="star", line=dict(color="black", width=1)),
        name="Punto de salida (ajustado/snap)", hoverinfo="name",
    ))

    titulo = f"Cuenca hidrológica completa - {id_corredor}"
    if subtitulo:
        titulo += f"<br><sub>{subtitulo}</sub>"

    fig.update_layout(
        title=titulo, xaxis_title="Este [km]", yaxis_title="Norte [km]",
        yaxis=dict(scaleanchor="x", scaleratio=1),  # proporción real 1:1, sin distorsión
        autosize=True,
    )
    fig.write_html(html_path)
    return html_path


def generar_mapa_3d_corredor(cuenca, geom_utm_origen, id_corredor, html_path, utm_crs=None,
                              subtitulo=None, paso_downsample=None, validacion=None):
    """Mapa 3D del corredor completo -- misma idea que
    geomatica.generar_mapa_3d() (hover con altitud/lat/lon real vía
    customdata), pero adaptado a que aquí el área es 10-50x más grande que
    un buffer de un solo sitio:

    - El TERRENO de fondo se dibuja "downsampleado" (cada `paso_downsample`
      píxeles, no cada uno) -- una malla a resolución nativa de 30m sobre
      decenas de km tendría millones de vértices y sería pesada para que
      el navegador la renderice fluido. El downsample es solo visual/de
      geometría; el hover de CADA vértice sigue mostrando su altitud y
      lat/lon reales de ESE punto downsampleado (no un promedio ni un
      valor inventado).
    - Los CAUCES dentro de la cuenca (`stream_mask_corredor`, ya calculados
      en delimitar_cuenca()) se dibujan a resolución NATIVA completa (30m)
      -- son puntos dispersos, no una malla, así que renderizar miles de
      ellos es barato. Esto es lo que responde "¿por dónde corre el río
      del polígono al punto de salida, con qué coordenadas?".
    - Consecuencia cosmética conocida (documentada, no oculta): como el
      terreno de fondo es más grueso que los cauces, en zonas de relieve
      muy quebrado un punto de cauce puede verse ligeramente enterrado o
      flotando sobre la malla gruesa -- es una discrepancia visual entre
      dos resoluciones distintas del mismo DEM, no un dato incorrecto (el
      hover de cada cauce sigue mostrando SU altitud real, no la de la
      malla de fondo).

    `validacion`: resultado opcional de validar_cuenca_contra_inegi(). Si se
    pasa, los cauces YA NO se dibujan en un solo color (dodgerblue) -- se
    dividen en verde (dentro de tolerancia de la red oficial INEGI) / rojo
    (fuera), y se agrega la red oficial de INEGI en azul claro como
    referencia -- exactamente la misma gramática visual que
    validacion_hidrologica.generar_mapa_3d_validacion(), para que el usuario
    no tenga que aprender un lenguaje de colores distinto por módulo. Sin
    `validacion` (default), el comportamiento es idéntico al de antes: un
    solo trazo dodgerblue sin comparar contra INEGI (útil para --demo, donde
    no hay ningún shapefile de INEGI que cargar)."""
    import plotly.graph_objects as go
    from core.geomatica import calcular_grid_latlon

    Z_raw = cuenca["Z_raw"]
    catchment_mask = cuenca["catchment_mask"]
    stream_mask_corredor = cuenca["stream_mask_corredor"]
    transform = cuenca["transform"]
    pw, ph = cuenca["pw"], cuenca["ph"]
    rows, cols = Z_raw.shape

    paso = paso_downsample or max(1, round(max(rows, cols) / 350.0))
    log(f"Malla 3D del corredor: downsample cada {paso} píxeles para renderizado "
        f"({rows}x{cols} -> ~{rows // paso}x{cols // paso}) -- los cauces se dibujan a resolución nativa completa.")

    lat_full = lon_full = None
    if utm_crs:
        lat_full, lon_full = calcular_grid_latlon(transform, utm_crs, rows, cols)

    Z_mesh = Z_raw[::paso, ::paso]
    mrows, mcols = Z_mesh.shape
    x_km_mesh = np.arange(mcols) * (pw * paso) / 1000.0
    y_km_mesh = np.flipud(np.arange(mrows) * (ph * paso) / 1000.0)

    fig = go.Figure()

    customdata_terreno = None
    if lat_full is not None:
        customdata_terreno = np.dstack([Z_mesh, lat_full[::paso, ::paso], lon_full[::paso, ::paso]])
    fig.add_trace(go.Surface(
        z=Z_mesh, x=x_km_mesh, y=y_km_mesh, colorscale="Earth", opacity=0.93, name="Terreno",
        customdata=customdata_terreno, showscale=False,  # altitud ya se lee en Z y en el hover -- mismo
        # criterio que geomatica.generar_mapa_3d(): un colorbar de altitud aquí solo compite por
        # espacio con la leyenda (Polígono de origen, Cauces, Punto de salida, etc.)
        hovertemplate=("Altitud: %{customdata[0]:.0f} msnm<br>Lat: %{customdata[1]:.5f}<br>"
                        "Lon: %{customdata[2]:.5f}<extra></extra>") if customdata_terreno is not None
                       else "Altitud: %{z:.0f} msnm<extra></extra>",
    ))

    def _z_en_utm(x_utm, y_utm):
        col, row = ~transform * (x_utm, y_utm)
        col = int(np.clip(round(col), 0, cols - 1))
        row = int(np.clip(round(row), 0, rows - 1))
        return float(Z_raw[row, col]) if not np.isnan(Z_raw[row, col]) else 0.0

    river_y, river_x = np.where(stream_mask_corredor)
    if len(river_x) > 0:
        x_km_riv = river_x * pw / 1000.0
        y_km_riv = (rows - 1 - river_y) * ph / 1000.0
        offset_z = max(ph, pw) * 0.5  # chico, solo para que no se entierre visualmente en la malla gruesa
        z_riv = Z_raw[river_y, river_x] + offset_z

        # dentro_tolerancia (de validar_cuenca_contra_inegi) se calculó con
        # np.where(cuenca["stream_mask_corredor"]) -- el MISMO array, ya
        # congelado en el dict `cuenca` -- así que el orden de river_y/river_x
        # aquí es idéntico al que se usó para construir ese array booleano.
        dentro = validacion.get("dentro_tolerancia") if validacion is not None else None
        if dentro is not None and len(dentro) == len(river_x):
            for etiqueta, color, mascara in [
                ("Cauces dentro de tolerancia INEGI", "green", dentro),
                ("Cauces fuera de tolerancia INEGI", "red", ~dentro),
            ]:
                if not np.any(mascara):
                    continue
                kwargs_riv = {}
                if lat_full is not None:
                    customdata_riv = np.column_stack([
                        Z_raw[river_y[mascara], river_x[mascara]],
                        lat_full[river_y[mascara], river_x[mascara]],
                        lon_full[river_y[mascara], river_x[mascara]],
                    ])
                    kwargs_riv = dict(
                        customdata=customdata_riv,
                        hovertemplate=(f"{etiqueta}<br>Altitud: %{{customdata[0]:.0f}} msnm<br>"
                                        "Lat: %{customdata[1]:.5f}<br>Lon: %{customdata[2]:.5f}<extra></extra>"),
                    )
                fig.add_trace(go.Scatter3d(
                    x=x_km_riv[mascara], y=y_km_riv[mascara], z=z_riv[mascara], mode="markers",
                    marker=dict(size=2.3, color=color, opacity=0.9),
                    name=etiqueta, **kwargs_riv,
                ))
        else:
            kwargs_riv = {}
            if lat_full is not None:
                customdata_riv = np.column_stack([
                    Z_raw[river_y, river_x], lat_full[river_y, river_x], lon_full[river_y, river_x],
                ])
                kwargs_riv = dict(
                    customdata=customdata_riv,
                    hovertemplate=("Cauce<br>Altitud: %{customdata[0]:.0f} msnm<br>Lat: %{customdata[1]:.5f}<br>"
                                    "Lon: %{customdata[2]:.5f}<extra></extra>"),
                )
            fig.add_trace(go.Scatter3d(
                x=x_km_riv, y=y_km_riv, z=z_riv, mode="markers",
                marker=dict(size=2.3, color="dodgerblue", opacity=0.9),
                name="Cauces dentro de la cuenca", **kwargs_riv,
            ))

    # --- Red oficial de INEGI (referencia, azul claro) -- solo si se pasó `validacion` ---
    if validacion is not None and len(validacion.get("puntos_ref_utm", [])) > 0:
        puntos_ref = np.asarray(validacion["puntos_ref_utm"])
        xs_ref, ys_ref = puntos_ref[:, 0], puntos_ref[:, 1]
        x_km_ref, y_km_ref = _utm_a_km_local(xs_ref, ys_ref, transform, rows, pw, ph)
        # Filtra puntos de INEGI fuera del lienzo del corredor -- el bbox
        # usado para cargar INEGI puede extenderse un poco más allá del DEM
        # descargado (mismo margen, pero no exactamente el mismo rectángulo).
        dentro_lienzo = ((x_km_ref >= 0) & (x_km_ref <= cols * pw / 1000.0) &
                          (y_km_ref >= 0) & (y_km_ref <= rows * ph / 1000.0))
        if np.any(dentro_lienzo):
            xs_ref_v, ys_ref_v = xs_ref[dentro_lienzo], ys_ref[dentro_lienzo]
            x_km_ref_v, y_km_ref_v = x_km_ref[dentro_lienzo], y_km_ref[dentro_lienzo]
            z_ref_v = [_z_en_utm(x, y) + max(ph, pw) * 0.3 for x, y in zip(xs_ref_v, ys_ref_v)]
            kwargs_ref = {}
            if utm_crs:
                proj_a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)
                lon_ref, lat_ref = proj_a_wgs84.transform(xs_ref_v, ys_ref_v)
                kwargs_ref = dict(
                    customdata=np.column_stack([lat_ref, lon_ref]),
                    hovertemplate="Red oficial INEGI<br>Lat: %{customdata[0]:.5f}<br>Lon: %{customdata[1]:.5f}<extra></extra>",
                )
            fig.add_trace(go.Scatter3d(
                x=x_km_ref_v, y=y_km_ref_v, z=z_ref_v, mode="markers",
                marker=dict(size=1.6, color="deepskyblue", opacity=0.6),
                name="Red oficial INEGI (referencia)", **kwargs_ref,
            ))

    xs_o, ys_o = geom_utm_origen.exterior.coords.xy
    x_km_o, y_km_o = _utm_a_km_local(xs_o, ys_o, transform, rows, pw, ph)
    z_o = [_z_en_utm(x, y) + max(ph, pw) * 0.5 for x, y in zip(xs_o, ys_o)]
    fig.add_trace(go.Scatter3d(
        x=x_km_o, y=y_km_o, z=z_o, mode="lines", line=dict(color="red", width=5),
        name="Polígono de origen",
    ))

    x_snap_km, y_snap_km = _utm_a_km_local(
        [cuenca["punto_snap_utm"][0]], [cuenca["punto_snap_utm"][1]], transform, rows, pw, ph)
    z_snap = _z_en_utm(*cuenca["punto_snap_utm"]) + max(ph, pw) * 1.5
    fig.add_trace(go.Scatter3d(
        x=x_snap_km, y=y_snap_km, z=[z_snap], mode="markers+text",
        marker=dict(size=8, color="lime", symbol="diamond", line=dict(color="black", width=1)),
        text=["▼ Punto de salida"], textposition="top center",
        name="Punto de salida (ajustado/snap)",
    ))

    titulo = f"Cuenca hidrológica completa (3D) - {id_corredor}"
    if subtitulo:
        titulo += f"<br><sub>{subtitulo}</sub>"
    if validacion is not None and validacion.get("pct_dentro_tolerancia") is not None:
        titulo += (f"<br><sub>Cauces vs. Red oficial INEGI: {validacion['pct_dentro_tolerancia']}% dentro de "
                    f"{validacion['tolerancia_m']:.0f}m (dist. promedio: {validacion['distancia_promedio_m']}m) "
                    f"-- verde=dentro, rojo=fuera, azul claro=red oficial INEGI</sub>")

    fig.update_layout(
        title=titulo,
        scene=dict(xaxis_title="Este [km]", yaxis_title="Norte [km]", zaxis_title="Altitud [msnm]",
                   aspectmode="manual", aspectratio=dict(x=1, y=1, z=0.35),
                   camera=dict(eye=dict(x=1.35, y=-1.35, z=0.8), center=dict(x=0, y=0, z=-0.05))),
        # Mismo criterio que geomatica.generar_mapa_3d(): leyenda anclada arriba a la
        # izquierda (antes, sin posición fija, Plotly la mandaba a la esquina superior
        # derecha por default -- aquí no hay colorbar que le compita, pero sí queda
        # mejor separada del área del mapa en corredores muy anchos).
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.75)"),
        margin=dict(l=10, r=10, t=90, b=10),
        autosize=True,
    )
    fig.write_html(html_path)
    return html_path


# ==============================================================================
# --- ORQUESTADOR CON DATOS REALES ---
# ==============================================================================
def analizar_conexion_hidrologica(geojson_origen, punto_salida_latlon, id_corredor,
                                   carpeta_srtm=None, carpeta_salida=None,
                                   margen_grados=None, percentil_snap=None,
                                   shapefile_inegi=None, layer_inegi=None, tolerancia_m=None):
    """Pipeline completo: bbox del corredor -> DEM grande -> D8 sobre el
    corredor -> snap + catchment() desde el punto de salida -> traslape
    contra el polígono de origen -> (opcional) validación contra la Red
    Hidrográfica oficial de INEGI -> CSV + mapa 2D + mapa 3D.
    `punto_salida_latlon` es (lat, lon) -- el orden natural en que se copian
    coordenadas de Google Earth/Maps, se convierte internamente a (lon, lat)
    para pyproj.

    `shapefile_inegi`: ruta opcional al cnit50k.gpkg (u otro shapefile) de
    INEGI. Si se pasa, los cauces D8 DENTRO de la cuenca delimitada se
    comparan contra la red oficial (validar_cuenca_contra_inegi()) y el
    resultado se agrega al CSV y se dibuja en el mapa 3D (verde/rojo +
    referencia azul) -- exactamente igual que ya se hace por sitio en
    core/validacion_hidrologica.py. Sin este parámetro (default None), el
    comportamiento es idéntico al de antes: no se toca INEGI para nada.
    Devuelve (df_resumen, csv_path, html_path, html_path_3d)."""
    carpeta_srtm = carpeta_srtm or CARPETA_SRTM
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_corredor.lower()}")
    os.makedirs(carpeta_salida, exist_ok=True)
    percentil_snap = percentil_snap if percentil_snap is not None else PERCENTIL_SNAP_CAUCE

    if not os.path.exists(geojson_origen):
        raise FileNotFoundError(f"No se encontró el GeoJSON de origen en: {geojson_origen}")

    gdf_origen = gpd.read_file(geojson_origen)
    geom_wgs84_origen = (gdf_origen.geometry.union_all() if hasattr(gdf_origen.geometry, "union_all")
                          else gdf_origen.geometry.unary_union)

    lat_salida, lon_salida = punto_salida_latlon
    punto_salida_wgs84 = (lon_salida, lat_salida)

    bbox = calcular_bbox_corredor(geom_wgs84_origen, punto_salida_wgs84, margen_grados)
    log(f"=== Analizando conexión hidrológica: {id_corredor} ===")
    log(f"Corredor bbox: {tuple(round(b, 4) for b in bbox)} "
        f"(~{(bbox[2] - bbox[0]) * 111:.0f}km x {(bbox[3] - bbox[1]) * 111:.0f}km aprox)")

    dst_array, meta_utm, utm_crs = cargar_dem_corredor(bbox, id_corredor, carpeta_srtm)

    proyector = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True).transform
    geom_utm_origen = shp_transform(proyector, geom_wgs84_origen)
    x_utm_salida, y_utm_salida = proyector(lon_salida, lat_salida)

    cuenca = delimitar_cuenca(
        dst_array, meta_utm, (x_utm_salida, y_utm_salida), utm_crs, carpeta_srtm, id_corredor, percentil_snap,
    )
    traslape = verificar_traslape(cuenca["catchment_mask"], cuenca["transform"], geom_utm_origen,
                                   cuenca["pw"], cuenca["ph"])

    log(f"Traslape: {traslape['area_origen_dentro_ha']:.1f} ha de {traslape['area_origen_total_ha']:.1f} ha "
        f"del origen caen DENTRO de la cuenca delimitada ({traslape['pct_traslape']:.1f}%).")

    # --- Validación opcional contra la Red Hidrográfica oficial de INEGI ---
    # NO se usa ningún dato de INEGI en ningún paso anterior de este pipeline
    # (el DEM es SRTM, los cauces son del D8 propio) -- esto SOLO compara,
    # después de ya tener la cuenca delimitada, qué tan cerca cae de la red
    # oficial. Si no se pasa --shapefile-inegi, este bloque no se ejecuta y
    # el comportamiento es idéntico al de antes.
    validacion = None
    if shapefile_inegi:
        log("Validando cauces D8 de la cuenca contra la Red Hidrográfica oficial de INEGI...")
        validacion = validar_cuenca_contra_inegi(
            cuenca, bbox, utm_crs, shapefile_inegi, layer=layer_inegi, tolerancia_m=tolerancia_m,
        )

    fila = {
        "id_corredor": id_corredor,
        "geojson_origen": geojson_origen,
        "lat_salida_dado": lat_salida, "lon_salida_dado": lon_salida,
        "distancia_snap_m": round(cuenca["distancia_snap_m"], 1),
        "area_origen_total_ha": round(traslape["area_origen_total_ha"], 2),
        "area_origen_dentro_cuenca_ha": round(traslape["area_origen_dentro_ha"], 2),
        "pct_traslape": round(traslape["pct_traslape"], 2),
        "area_cuenca_total_ha": round(traslape["area_cuenca_total_ha"], 2),
    }
    if validacion is not None:
        fila.update({
            "n_puntos_cauce_validados_inegi": validacion["n_puntos_d8"],
            "pct_cauces_dentro_tolerancia_inegi": validacion["pct_dentro_tolerancia"],
            "distancia_promedio_inegi_m": validacion["distancia_promedio_m"],
            "distancia_mediana_inegi_m": validacion["distancia_mediana_m"],
            "tolerancia_inegi_m": validacion["tolerancia_m"],
        })
    df = pd.DataFrame([fila])
    csv_path = os.path.join(carpeta_salida, f"conexion_hidrologica_{id_corredor.lower()}.csv")
    df.to_csv(csv_path, index=False)
    log(f"CSV de conexión hidrológica guardado en: {csv_path}")

    html_path = html_path_3d = None
    try:
        subtitulo = (f"{traslape['pct_traslape']:.1f}% del área de origen escurre hacia el punto de salida dado "
                     f"({traslape['area_origen_dentro_ha']:.0f} de {traslape['area_origen_total_ha']:.0f} ha) "
                     f"-- ajuste del punto de salida: {cuenca['distancia_snap_m']:.0f}m")
        if validacion is not None and validacion.get("pct_dentro_tolerancia") is not None:
            subtitulo += (f" -- cauces vs INEGI: {validacion['pct_dentro_tolerancia']}% dentro de "
                          f"{validacion['tolerancia_m']:.0f}m")
        html_path = os.path.join(carpeta_salida, f"{id_corredor.lower()}_cuenca_completa.html")
        generar_mapa_corredor(cuenca, geom_utm_origen, id_corredor, html_path, utm_crs=utm_crs, subtitulo=subtitulo)
        log(f"Mapa 2D de la cuenca completa: {html_path}")

        html_path_3d = os.path.join(carpeta_salida, f"{id_corredor.lower()}_cuenca_completa_3d.html")
        generar_mapa_3d_corredor(cuenca, geom_utm_origen, id_corredor, html_path_3d, utm_crs=utm_crs,
                                  subtitulo=subtitulo, validacion=validacion)
        log(f"Mapa 3D de la cuenca completa (con cauces y hover lat/lon/altitud): {html_path_3d}")
    except ImportError as e:
        log(f"plotly no instalado, se omiten los mapas: {e}", nivel="WARN")

    return df, csv_path, html_path, html_path_3d


# ==============================================================================
# --- MODO DEMO: DEM sintético en memoria, sin red -- mismo espíritu que
#     --demo en geomatica.py y carbono.py ---
# ==============================================================================
def _dem_sintetico_corredor(size=300, resolucion_m=100.0, semilla=7):
    """DEM sintético de un 'corredor' completo (mucho más grande que el de
    geomatica._dem_sintetico, que es de un solo sitio): una ladera larga
    que baja de una esquina alta (análoga a Cofre de Perote, el origen)
    hacia la esquina opuesta más baja (análoga a la Cascada de Texolo, el
    punto de salida), con una cañada diagonal que concentra el flujo --
    determinista, sin red."""
    rng = np.random.default_rng(semilla)
    y, x = np.mgrid[0:size, 0:size]

    plano_inclinado = 4200 - 3.2 * (x + y)  # alto en (0,0), bajo en (size,size)
    canada = -200 * np.exp(-(((x - y) ** 2) / (2 * (size / 6) ** 2)))  # concentra flujo en la diagonal
    ruido = rng.normal(0, 4, size=(size, size))
    Z = np.clip(plano_inclinado + canada + ruido, 800, None).astype(np.float32)

    transform = rasterio.transform.from_origin(0, size * resolucion_m, resolucion_m, resolucion_m)
    meta_utm = {
        "driver": "GTiff", "height": size, "width": size,
        "count": 1, "dtype": "float32", "crs": "EPSG:32614", "transform": transform,
    }
    # Origen: círculo cerca de la esquina alta (donde "nace" el agua)
    origen_x_m, origen_y_m = size * 0.12 * resolucion_m, (size - size * 0.12) * resolucion_m
    geom_utm_origen = Point(origen_x_m, origen_y_m).buffer(900)
    # Punto de salida: cerca de la esquina baja, sobre la cañada modelada
    salida_utm = (size * 0.88 * resolucion_m, (size - size * 0.88) * resolucion_m)
    return Z, meta_utm, geom_utm_origen, salida_utm, "EPSG:32614"


def demo():
    """Corre el pipeline completo (delimitación de cuenca + traslape + mapa)
    sobre un DEM sintético en memoria. No descarga nada, no necesita ningún
    geojson del usuario -- sirve para confirmar que la lógica sigue intacta
    después de cualquier cambio a este módulo, igual que --demo en
    geomatica.py y carbono.py."""
    log("=== core.cuenca_completa --demo (DEM sintético, sin red) ===")
    dst_array, meta_utm, geom_utm_origen, punto_salida_utm, utm_crs = _dem_sintetico_corredor()
    id_corredor = "DEMO_CORREDOR"

    try:
        carpeta_tmp = os.path.expanduser("~/resultados_demo_cuenca_completa")
        os.makedirs(carpeta_tmp, exist_ok=True)
        cuenca = delimitar_cuenca(dst_array, meta_utm, punto_salida_utm, utm_crs, carpeta_tmp, id_corredor)
        traslape = verificar_traslape(cuenca["catchment_mask"], cuenca["transform"], geom_utm_origen,
                                       cuenca["pw"], cuenca["ph"])

        print("\n--- Traslape origen vs cuenca delimitada (demo sintético) ---")
        for k, v in traslape.items():
            print(f"  {k}: {v:.2f}")
        print(f"  distancia_snap_m: {cuenca['distancia_snap_m']:.1f}")

        subtitulo = f"{traslape['pct_traslape']:.1f}% del origen escurre hacia el punto de salida (demo sintético)"
        html_path = os.path.join(carpeta_tmp, f"{id_corredor.lower()}_cuenca_completa.html")
        generar_mapa_corredor(cuenca, geom_utm_origen, id_corredor, html_path, utm_crs=utm_crs, subtitulo=subtitulo)
        log(f"Mapa 2D demo generado en: {html_path}")

        html_path_3d = os.path.join(carpeta_tmp, f"{id_corredor.lower()}_cuenca_completa_3d.html")
        generar_mapa_3d_corredor(cuenca, geom_utm_origen, id_corredor, html_path_3d, utm_crs=utm_crs, subtitulo=subtitulo)
        log(f"Mapa 3D demo generado en: {html_path_3d}")
        return traslape
    except ImportError as e:
        log(f"pysheds/plotly no instalado, se omite la demo completa: {e}", nivel="WARN")
        return None


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Delimita la cuenca hidrológica real entre un sitio de origen y un punto de salida río "
                    "abajo (ej. ¿el Cofre de Perote realmente alimenta la Cascada de Texolo?) -- Motor Nacional"
    )
    ap.add_argument("--demo", action="store_true", help="Corre con un DEM sintético, sin red ni archivos del usuario")
    ap.add_argument("--geojson-origen", type=str, help="GeoJSON del polígono de origen (ej. Cofre de Perote)")
    ap.add_argument("--punto-salida", type=str,
                     help="Coordenadas 'lat,lon' del punto de salida río abajo (ej. la Cascada de Texolo)")
    ap.add_argument("--id-corredor", type=str, default=None, help="Nombre identificador (para nombres de archivo)")
    ap.add_argument("--margen-grados", type=float, default=None,
                     help=f"Margen extra alrededor del bbox origen+salida (default: config.MARGEN_CORREDOR_GRADOS={MARGEN_CORREDOR_GRADOS})")
    ap.add_argument("--percentil-snap", type=float, default=None,
                     help=f"Percentil de acumulación de flujo para el snap del punto de salida (default: config.PERCENTIL_SNAP_CAUCE={PERCENTIL_SNAP_CAUCE})")
    ap.add_argument("--carpeta-srtm", type=str, default=None)
    ap.add_argument("--carpeta-salida", type=str, default=None)
    ap.add_argument("--shapefile-inegi", type=str, default=None,
                     help="Ruta al shapefile/gpkg de la Red Hidrográfica de INEGI (ej. cnit50k.gpkg) -- si se "
                          "pasa, valida los cauces D8 DENTRO de la cuenca delimitada contra la red oficial "
                          "(igual que core/validacion_hidrologica.py, pero solo para la porción de cauces que "
                          "importa para esta conexión origen->salida). Sin esto, no se toca ningún dato de INEGI.")
    ap.add_argument("--layer", type=str, default="corriente_ag_l",
                     help="Capa dentro del GeoPackage de INEGI a usar (default: 'corriente_ag_l', la red LINEAL "
                          "de cauces del CNIT50k -- no 'corriente_ag_a', que es polígono para ríos anchos)")
    ap.add_argument("--tolerancia-m", type=float, default=None,
                     help=f"Tolerancia de distancia en metros para la validación contra INEGI "
                          f"(default: config.TOLERANCIA_VALIDACION_HIDRO_M={TOLERANCIA_VALIDACION_HIDRO_M})")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not args.geojson_origen or not args.punto_salida:
        ap.error("--geojson-origen y --punto-salida son obligatorios fuera de --demo "
                 "(ej. --punto-salida '19.4213,-97.0295')")

    try:
        lat_str, lon_str = args.punto_salida.split(",")
        punto_salida_latlon = (float(lat_str.strip()), float(lon_str.strip()))
    except ValueError:
        ap.error("--punto-salida debe ser 'lat,lon', ej. '19.4213,-97.0295'")
        return

    id_corredor = args.id_corredor or os.path.splitext(os.path.basename(args.geojson_origen))[0]

    analizar_conexion_hidrologica(
        geojson_origen=args.geojson_origen, punto_salida_latlon=punto_salida_latlon, id_corredor=id_corredor,
        carpeta_srtm=args.carpeta_srtm, carpeta_salida=args.carpeta_salida,
        margen_grados=args.margen_grados, percentil_snap=args.percentil_snap,
        shapefile_inegi=args.shapefile_inegi, layer_inegi=args.layer, tolerancia_m=args.tolerancia_m,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Valida la red de cauces D8 (generada por core/geomatica.py) contra la Red
Hidrográfica Digital oficial de INEGI (escala 1:50,000), calculando qué
porcentaje de los cauces modelados cae dentro de una tolerancia de distancia
respecto a la red oficial.

POR QUÉ ESTE MÓDULO ESTÁ SEPARADO DE core/geomatica.py Y core/carbono.py:
    Mismo principio que ya aplicamos con carbono.py: cada módulo habla con
    UNA fuente de datos externa (geomatica.py con SRTM local, carbono.py con
    Earth Engine, este módulo con un shapefile de INEGI). Si algo falla, se
    sabe de inmediato cuál pieza fue -- no se mezclan fuentes en un mismo
    archivo.

FUENTE DE DATOS Y SU LIMITACIÓN -- documentada, no escondida:
    Red Hidrográfica Digital de México, escala 1:50,000, edición 2.0
    (INEGI). Descarga manual desde el portal de INEGI (Conjunto Nacional
    de Información Topográfica a Escala 1:50,000):
        https://gaia.inegi.org.mx/app/geo2/cnit50k/
    No hay una API pública de descarga directa -- hay que bajar el
    shapefile del estado/área de interés desde el portal a mano.

    - Escala 1:50,000 significa que la red oficial NO captura arroyos
      efímeros/intermitentes muy pequeños que el modelo D8 sí puede marcar
      (el D8 corre sobre SRTM de 30m, más fino que la digitalización
      1:50,000). Una discrepancia baja NO necesariamente significa que el
      modelo esté mal -- puede significar que detectó un cauce real que
      INEGI no digitalizó a esa escala.
    - Este módulo solo usa la GEOMETRÍA del shapefile de INEGI (las
      líneas), no depende de nombres de columnas de atributos -- para no
      repetir el riesgo de asumir un nombre de campo que resulte distinto
      en la práctica.

Qué SÍ calcula:
    - % de píxeles de cauce D8 que caen dentro de una tolerancia de
      distancia (default 30m, igual a la resolución del SRTM) respecto a
      la línea más cercana de la red oficial, por zona.
    - Distancia promedio y mediana de los cauces D8 a la red oficial.

Qué NO calcula:
    - Ninguna corrección automática del modelo D8 -- esto es solo una
      métrica de validación, no ajusta el algoritmo de acumulación de flujo.
"""

import argparse
import os

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString
from scipy.spatial import cKDTree

from config import TOLERANCIA_VALIDACION_HIDRO_M, ZONAS_ANALISIS_M, PERCENTIL_CAUCE_HIDROLOGIA, log


# ==============================================================================
# --- FUNCIONES PURAS ---
# ==============================================================================
def vectorizar_cauces_d8(stream_mask, zona_de_pixel, transform):
    """Convierte la máscara raster de cauces D8 (booleana) a puntos
    vectoriales en coordenadas del CRS del transform (normalmente UTM,
    metros). Devuelve un GeoDataFrame con columnas x, y, zona, geometry."""
    filas_y, filas_x = np.where(stream_mask)
    if len(filas_x) == 0:
        return gpd.GeoDataFrame(columns=["zona", "geometry"])

    xs, ys = [], []
    for row, col in zip(filas_y, filas_x):
        x, y = transform * (col + 0.5, row + 0.5)  # centro del píxel
        xs.append(x)
        ys.append(y)

    zonas = [zona_de_pixel[row, col] for row, col in zip(filas_y, filas_x)]
    gdf = gpd.GeoDataFrame(
        {"zona": zonas, "x": xs, "y": ys},
        geometry=[Point(x, y) for x, y in zip(xs, ys)],
    )
    return gdf


def _muestrear_lineas_a_puntos(gdf_lineas_utm, paso_m=15.0):
    """Convierte las líneas de la red oficial a una nube densa de puntos
    (muestreados cada `paso_m` metros a lo largo de cada línea), para poder
    usar un árbol de vecinos más cercanos (KD-tree) en vez de calcular
    distancia línea-a-línea, que es mucho más lento sobre redes grandes."""
    puntos = []
    for geom in gdf_lineas_utm.geometry:
        if geom is None or geom.is_empty:
            continue
        lineas = [geom] if isinstance(geom, LineString) else list(geom.geoms)
        for linea in lineas:
            largo = linea.length
            if largo == 0:
                puntos.append(linea.coords[0])
                continue
            n_muestras = max(2, int(largo / paso_m) + 1)
            for i in range(n_muestras):
                d = (i / (n_muestras - 1)) * largo
                p = linea.interpolate(d)
                puntos.append((p.x, p.y))
    return np.array(puntos)


def validar_zona(puntos_d8_xy, puntos_referencia_xy, tolerancia_m):
    """Calcula, para un conjunto de puntos D8, la distancia al punto de
    referencia (INEGI) más cercano, y el % que cae dentro de la
    tolerancia. Usa un KD-tree -- rápido incluso con miles de puntos.
    También devuelve 'distancias_m' (array por punto) y
    'dentro_tolerancia' (booleano por punto), para poder colorear un mapa
    3D punto a punto -- no solo el agregado."""
    if len(puntos_d8_xy) == 0:
        return {"n_puntos_d8": 0, "pct_dentro_tolerancia": None,
                "distancia_promedio_m": None, "distancia_mediana_m": None,
                "distancias_m": np.array([]), "dentro_tolerancia": np.array([], dtype=bool)}
    if len(puntos_referencia_xy) == 0:
        log("Sin puntos de referencia INEGI en esta zona -- no se puede validar.", nivel="WARN")
        return {"n_puntos_d8": len(puntos_d8_xy), "pct_dentro_tolerancia": None,
                "distancia_promedio_m": None, "distancia_mediana_m": None,
                "distancias_m": np.full(len(puntos_d8_xy), np.nan), "dentro_tolerancia": np.zeros(len(puntos_d8_xy), dtype=bool)}

    arbol = cKDTree(puntos_referencia_xy)
    distancias, _ = arbol.query(puntos_d8_xy, k=1)
    dentro = distancias <= tolerancia_m

    return {
        "n_puntos_d8": len(puntos_d8_xy),
        "pct_dentro_tolerancia": round(float(np.mean(dentro) * 100), 1),
        "distancia_promedio_m": round(float(np.mean(distancias)), 1),
        "distancia_mediana_m": round(float(np.median(distancias)), 1),
        "distancias_m": distancias,
        "dentro_tolerancia": dentro,
    }


def _validar_por_zona_con_total(zona_por_punto, xy_puntos_d8, puntos_ref, zonas_m, tolerancia_m):
    """Corre validar_zona() por cada zona (nucleo + buffers) y agrega una
    fila TOTAL al final. Las zonas de zona_por_punto ya son EXCLUSIVAS por
    diseño -- cada píxel de cauce D8 pertenece a una sola zona (ver
    core/geomatica.py, calcular_hidrologia_d8: itera los buffers de menor a
    mayor y solo asigna los píxeles que siguen marcados "fuera"), así que
    sumar/juntar zonas aquí no duplica ningún punto.

    El TOTAL se calcula corriendo validar_zona() sobre la UNIÓN real de
    todos los puntos D8 -- nunca como promedio de los % por zona. Promediar
    los porcentajes sería el mismo error de fondo ya corregido en
    core/carbono_perdida.py (generar_balance_stock_vs_perdida): si las
    zonas tienen tamaños de muestra distintos, un promedio simple no
    representa el total real. Ej. con los datos reales de Cofre de Perote
    (nucleo 1806 pts al 41.7%, buffer_500m 1154 pts al 56.2%, buffer_1000m
    821 pts al 54.7%), el promedio simple de los tres % da 50.9%, pero el
    TOTAL correcto (ponderado por cuántos puntos tiene cada zona,
    equivalente a recalcularlo sobre todos los puntos juntos) da 48.9% --
    una diferencia real, no redondeo.

    Devuelve una lista de dicts (uno por zona + uno TOTAL al final), lista
    para pd.DataFrame(...). Usado por validar_sitio_real() y
    generar_mapa_3d_validacion() -- antes cada una tenía su propio bucle
    casi idéntico para armar el CSV por zona, con el riesgo de que un
    cambio futuro se aplicara a una copia y no a la otra."""
    zona_por_punto = np.asarray(zona_por_punto)
    filas = []
    for etiqueta in ["nucleo"] + [f"buffer_{b}m" for b in sorted(zonas_m) if b > 0]:
        sel = zona_por_punto == etiqueta
        sub_puntos = xy_puntos_d8[sel] if len(xy_puntos_d8) else np.empty((0, 2))
        r = validar_zona(sub_puntos, puntos_ref, tolerancia_m)
        r["zona"] = etiqueta
        r["tolerancia_m"] = tolerancia_m
        filas.append(r)
        log(f"  Zona {etiqueta}: {r['n_puntos_d8']} píxeles de cauce D8, "
            f"{r['pct_dentro_tolerancia']}% dentro de {tolerancia_m}m de la red oficial "
            f"(dist. promedio: {r['distancia_promedio_m']}m)")

    r_total = validar_zona(xy_puntos_d8, puntos_ref, tolerancia_m)
    r_total["zona"] = "TOTAL (todas las zonas, sin traslape)"
    r_total["tolerancia_m"] = tolerancia_m
    filas.append(r_total)
    log(f"  TOTAL: {r_total['n_puntos_d8']} píxeles de cauce D8, "
        f"{r_total['pct_dentro_tolerancia']}% dentro de {tolerancia_m}m de la red oficial "
        f"(dist. promedio: {r_total['distancia_promedio_m']}m)")
    return filas


# ==============================================================================
# --- ORQUESTADOR CON DATOS REALES ---
# ==============================================================================
def validar_sitio_real(geojson_path, shapefile_inegi_path, id_proyecto, zonas_m=None,
                        tolerancia_m=None, percentil_cauce=None, carpeta_salida=None, carpeta_srtm=None,
                        layer=None):
    """Corre geomatica.py para obtener los cauces D8 (reusa el SRTM ya
    descargado si existe), carga la red oficial de INEGI (filtrada por
    bbox del sitio -- no carga el archivo completo si es un GeoPackage
    nacional grande), y calcula el índice de validación por zona (+ una
    fila TOTAL final, calculada sobre la unión real de los puntos D8 de
    todas las zonas -- ver _validar_por_zona_con_total). Devuelve
    (df_validacion, csv_path).

    `layer`: nombre de la capa dentro del GeoPackage/shapefile, si aplica
    (ej. 'corriente_ag_l' del CNIT50k de INEGI, que es la red LINEAL de
    cauces -- no usar 'corriente_ag_a', que es la versión en polígono para
    ríos anchos, no sirve para esta comparación punto-a-línea)."""
    from core import geomatica
    import config as cfg

    zonas_m = zonas_m if zonas_m is not None else ZONAS_ANALISIS_M
    tolerancia_m = tolerancia_m if tolerancia_m is not None else TOLERANCIA_VALIDACION_HIDRO_M
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    os.makedirs(carpeta_salida, exist_ok=True)

    if not os.path.exists(shapefile_inegi_path):
        raise FileNotFoundError(f"No se encontró el shapefile de INEGI en: {shapefile_inegi_path}")

    log("Cargando terreno y recalculando cauces D8 (reusa SRTM si ya está en caché)...")
    geom_utm_nucleo, dst_array, meta_utm, utm_crs = geomatica.cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        percentil_cauce or cfg.PERCENTIL_CAUCE_HIDROLOGIA,
        carpeta_srtm or cfg.CARPETA_SRTM, id_proyecto,
    )

    log(f"Cargando red hidrográfica oficial de INEGI: {shapefile_inegi_path}"
        + (f" (capa: {layer})" if layer else ""))

    # Filtrar por bbox del sitio en vez de cargar el archivo completo --
    # importante si es un GeoPackage nacional grande (ej. cnit50k.gpkg,
    # ~800MB con todo México). El bbox se calcula del propio geojson del
    # sitio, con margen generoso para cubrir el buffer más grande.
    gdf_sitio = gpd.read_file(geojson_path)
    geom_sitio_wgs84 = gdf_sitio.geometry.union_all() if hasattr(gdf_sitio.geometry, "union_all") else gdf_sitio.geometry.unary_union
    b = geom_sitio_wgs84.bounds
    margen_grados = 0.05 + (max(zonas_m) / 111000.0)
    bbox = (b[0] - margen_grados, b[1] - margen_grados, b[2] + margen_grados, b[3] + margen_grados)

    gdf_inegi = gpd.read_file(shapefile_inegi_path, layer=layer, bbox=bbox)
    log(f"  -> {len(gdf_inegi)} elementos cargados dentro del bbox del sitio (de todo el archivo, sin cargarlo completo)")
    gdf_inegi_utm = gdf_inegi.to_crs(utm_crs)
    log(f"  -> reproyectados a {utm_crs}")

    log("Muestreando la red oficial a puntos (para comparación por vecino más cercano)...")
    puntos_ref = _muestrear_lineas_a_puntos(gdf_inegi_utm)
    log(f"  -> {len(puntos_ref)} puntos de referencia")

    gdf_d8 = vectorizar_cauces_d8(hidrologia["stream_mask"], hidrologia["zona_de_pixel"],
                                   hidrologia["transform"])

    xy_todos = np.column_stack([gdf_d8["x"].values, gdf_d8["y"].values]) if len(gdf_d8) > 0 else np.empty((0, 2))
    filas = _validar_por_zona_con_total(gdf_d8["zona"].values, xy_todos, puntos_ref, zonas_m, tolerancia_m)

    df = pd.DataFrame(filas)[["zona", "n_puntos_d8", "tolerancia_m", "pct_dentro_tolerancia",
                               "distancia_promedio_m", "distancia_mediana_m"]]
    csv_path = os.path.join(carpeta_salida, f"validacion_hidrologica_{id_proyecto.lower()}_p{int(percentil_cauce or cfg.PERCENTIL_CAUCE_HIDROLOGIA)}.csv")
    df.to_csv(csv_path, index=False)
    log(f"CSV de validación guardado en: {csv_path}")
    return df, csv_path


# ==============================================================================
# --- MAPA 3D: D8 (verde=dentro tolerancia, rojo=fuera) + red oficial INEGI ---
# ==============================================================================
def generar_mapa_3d_validacion(geojson_path, shapefile_inegi_path, id_proyecto, zonas_m=None,
                                tolerancia_m=None, percentil_cauce=None, carpeta_salida=None,
                                carpeta_srtm=None, layer="corriente_ag_l"):
    """Genera un HTML interactivo con el terreno 3D, los cauces D8
    coloreados por si validaron o no (verde/rojo), y la red oficial de
    INEGI superpuesta en azul -- para poder ver a ojo, no solo con el
    número, si el D8 se equivoca o si la red oficial simplemente no
    digitalizó ese cauce a esa escala. Reusa geomatica.py para el
    terreno/D8 (mismo SRTM en caché, no descarga nada de nuevo)."""
    from core import geomatica
    import plotly.graph_objects as go
    import config as cfg

    zonas_m = zonas_m if zonas_m is not None else ZONAS_ANALISIS_M
    tolerancia_m = tolerancia_m if tolerancia_m is not None else TOLERANCIA_VALIDACION_HIDRO_M
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    carpeta_srtm = carpeta_srtm or cfg.CARPETA_SRTM
    os.makedirs(carpeta_salida, exist_ok=True)

    log("Cargando terreno y cauces D8 (reusa SRTM en caché)...")
    geom_utm_nucleo, dst_array, meta_utm, utm_crs = geomatica.cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        percentil_cauce or cfg.PERCENTIL_CAUCE_HIDROLOGIA, carpeta_srtm, id_proyecto,
    )
    Z_smooth = hidrologia["Z_smooth"]
    Z_raw = hidrologia.get("Z_raw", Z_smooth)
    stream_mask = hidrologia["stream_mask"]
    transform = hidrologia["transform"]
    pw_v, ph_v = hidrologia["pw_v"], hidrologia["ph_v"]
    rows, cols = Z_smooth.shape

    log(f"Cargando red hidrográfica oficial de INEGI (capa: {layer})...")
    gdf_sitio = gpd.read_file(geojson_path)
    geom_sitio_wgs84 = gdf_sitio.geometry.union_all() if hasattr(gdf_sitio.geometry, "union_all") else gdf_sitio.geometry.unary_union
    b = geom_sitio_wgs84.bounds
    margen_grados = 0.05 + (max(zonas_m) / 111000.0)
    bbox = (b[0] - margen_grados, b[1] - margen_grados, b[2] + margen_grados, b[3] + margen_grados)
    gdf_inegi_utm = gpd.read_file(shapefile_inegi_path, layer=layer, bbox=bbox).to_crs(utm_crs)
    puntos_ref = _muestrear_lineas_a_puntos(gdf_inegi_utm)
    log(f"  -> {len(puntos_ref)} puntos de referencia")

    # --- Puntos D8: coordenadas UTM (para medir distancia) y fila/col (para graficar) ---
    river_y, river_x = np.where(stream_mask)
    xs_utm, ys_utm = [], []
    for row, col in zip(river_y, river_x):
        x, y = transform * (col + 0.5, row + 0.5)
        xs_utm.append(x)
        ys_utm.append(y)
    puntos_d8_utm = np.column_stack([xs_utm, ys_utm]) if len(xs_utm) > 0 else np.empty((0, 2))
    resultado = validar_zona(puntos_d8_utm, puntos_ref, tolerancia_m)
    dentro = resultado["dentro_tolerancia"]

    # --- CSV por zona (mismo desglose que validar_sitio_real, para que
    #     --mapa-3d también deje el CSV -- antes solo generaba el HTML y
    #     el CSV se perdía si no se corría el otro comando por separado. ---
    zona_de_pixel = hidrologia["zona_de_pixel"]
    zona_por_punto_d8 = zona_de_pixel[river_y, river_x] if len(river_y) else np.array([], dtype=object)
    filas_csv = _validar_por_zona_con_total(zona_por_punto_d8, puntos_d8_utm, puntos_ref, zonas_m, tolerancia_m)
    df_zonas = pd.DataFrame(filas_csv)[["zona", "n_puntos_d8", "tolerancia_m", "pct_dentro_tolerancia",
                                         "distancia_promedio_m", "distancia_mediana_m"]]
    pct_str = int(percentil_cauce or cfg.PERCENTIL_CAUCE_HIDROLOGIA)
    csv_path = os.path.join(carpeta_salida, f"validacion_hidrologica_{id_proyecto.lower()}_p{pct_str}.csv")
    df_zonas.to_csv(csv_path, index=False)
    log(f"CSV de validación guardado en: {csv_path}")

    X, Y = np.meshgrid(np.arange(cols) * pw_v / 1000.0, np.flipud(np.arange(rows) * ph_v / 1000.0))

    # --- lat/lon real para el hover, en vez de solo km locales sin referencia ---
    import pyproj
    proj_a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)
    lat_grid, lon_grid = geomatica.calcular_grid_latlon(transform, utm_crs, rows, cols)

    fig = go.Figure()
    customdata_superficie = np.dstack([Z_raw, lat_grid, lon_grid])
    fig.add_trace(go.Surface(
        z=Z_smooth, x=X, y=Y, colorscale="Earth", connectgaps=False, opacity=0.85, name="Terreno 3D",
        customdata=customdata_superficie, showscale=False,  # altitud ya se lee en Z y en el hover --
        # mismo criterio que geomatica.generar_mapa_3d(): sin esto, el colorbar de altitud
        # choca con la leyenda (dentro/fuera de tolerancia, red oficial INEGI).
        hovertemplate="Altitud: %{customdata[0]:.0f} msnm<br>Lat: %{customdata[1]:.5f}<br>Lon: %{customdata[2]:.5f}<extra></extra>",
    ))

    for etiqueta, color, mascara in [
        ("D8 dentro de tolerancia", "green", dentro if len(dentro) else np.array([], dtype=bool)),
        ("D8 fuera de tolerancia", "red", ~dentro if len(dentro) else np.array([], dtype=bool)),
    ]:
        if not np.any(mascara):
            continue
        lons_pts, lats_pts = proj_a_wgs84.transform(puntos_d8_utm[mascara, 0], puntos_d8_utm[mascara, 1])
        # BUG CORREGIDO: se indexaba con "river_y_plotly" (doble flip), que
        # colocaba el punto en la fila espejo equivocada -- posición Y y
        # elevación de OTRO pixel del raster. river_y/river_x (índices
        # reales) son los correctos, porque X,Y ya tienen el flip
        # codificado en su propia construcción arriba.
        altitudes_reales = Z_raw[river_y[mascara], river_x[mascara]]
        fig.add_trace(go.Scatter3d(
            x=X[river_y[mascara], river_x[mascara]],
            y=Y[river_y[mascara], river_x[mascara]],
            z=Z_smooth[river_y[mascara], river_x[mascara]] + 12,
            mode="markers", marker=dict(size=2.5, color=color, opacity=0.9), name=etiqueta,
            customdata=np.column_stack([altitudes_reales, lats_pts, lons_pts]),
            hovertemplate=f"{etiqueta}<br>Altitud: %{{customdata[0]:.0f}} msnm<br>"
                          f"Lat: %{{customdata[1]:.5f}}<br>Lon: %{{customdata[2]:.5f}}<extra></extra>",
        ))

    # --- Red oficial INEGI: convertir de UTM a fila/col del raster (~transform inversa) ---
    inv = ~transform
    ref_x_km, ref_y_km, ref_z, ref_lat, ref_lon = [], [], [], [], []
    for x, y in puntos_ref:
        col_f, row_f = inv * (x, y)
        row_i, col_i = int(round(row_f)), int(round(col_f))
        if 0 <= row_i < rows and 0 <= col_i < cols and not np.isnan(Z_smooth[row_i, col_i]):
            ref_x_km.append(col_f * pw_v / 1000.0)
            ref_y_km.append((rows - 1 - row_f) * ph_v / 1000.0)
            ref_z.append(Z_smooth[row_i, col_i] + 8)
            lon_p, lat_p = proj_a_wgs84.transform(x, y)
            ref_lat.append(lat_p)
            ref_lon.append(lon_p)
    if ref_x_km:
        fig.add_trace(go.Scatter3d(
            x=ref_x_km, y=ref_y_km, z=ref_z,
            mode="markers", marker=dict(size=1.6, color="deepskyblue", opacity=0.6),
            name="Red oficial INEGI (referencia)",
            customdata=np.column_stack([ref_lat, ref_lon]),
            hovertemplate="Red oficial INEGI<br>Lat: %{customdata[0]:.5f}<br>Lon: %{customdata[1]:.5f}<extra></extra>",
        ))

    # Aviso honesto en el propio mapa, no solo en el docstring del módulo (mismo criterio ya
    # aplicado en deforestacion.py y validacion_incendios.py para el aviso de 723 vs 771 ha):
    # un % bajo aquí no necesariamente significa que el modelo D8 esté mal -- INEGI 1:50,000 no
    # digitaliza cauces efímeros/pequeños que el D8 (SRTM 30m, más fino) sí puede detectar. Antes
    # esta explicación solo vivía en el docstring del módulo -- nadie que abriera el mapa la veía.
    # `resultado` ya es el TOTAL sobre TODOS los puntos D8 juntos (calculado arriba para colorear
    # verde/rojo) -- se reusa aquí en vez de recalcularlo, para no correr KD-tree dos veces.
    def _fmt(v, decimales=0):
        return f"{v:.{decimales}f}" if v is not None else "N/D"

    pct_total = resultado["pct_dentro_tolerancia"]
    n_total = resultado["n_puntos_d8"]
    promedio_total = resultado["distancia_promedio_m"]
    mediana_total = resultado["distancia_mediana_m"]

    # Líneas CORTAS a propósito (mismo problema ya resuelto en deforestacion.py: un título de
    # Plotly no hace word-wrap solo, una línea larga se corta en el borde del gráfico) Y dentro
    # del MISMO bloque <sub>...</sub> que el resto del subtítulo (fuente más chica -- si se deja
    # una línea fuera de <sub>, Plotly la dibuja en el tamaño grande del título principal, que
    # ocupa más ancho Y más alto por línea de lo que asume calcular_margen_top_titulo()).
    subtitulo_extra = (
        f"<br>TOTAL: {_fmt(pct_total, 1)}% de {n_total} píxeles D8 dentro de {tolerancia_m:.0f}m "
        f"(mediana {_fmt(mediana_total)}m, promedio {_fmt(promedio_total)}m)."
        f"<br>% bajo no implica error del D8: INEGI (1:50,000) no digitaliza cauces pequeños "
        f"que el D8 (SRTM 30m) sí detecta."
    )
    if promedio_total is not None and mediana_total is not None and mediana_total > 0 \
            and promedio_total > 2 * mediana_total:
        # Mediana << promedio: la mayoría de los puntos SÍ está cerca (la mediana lo refleja),
        # pero un grupo minoritario está muy lejos de cualquier cauce oficial y jala el promedio
        # hacia arriba -- vale la pena decirlo porque de otro modo el promedio solo se lee como
        # "el modelo está uniformemente mal", cuando el patrón real es distinto (posible zona sin
        # digitalizar en INEGI, no un error parejo del D8 en todo el sitio).
        subtitulo_extra += (
            f"<br>Mediana ({_fmt(mediana_total)}m) << promedio ({_fmt(promedio_total)}m): la "
            "mayoría SÍ está cerca, una minoría está muy lejos (posible zona sin digitalizar)."
        )

    titulo = (f"Validación hidrológica D8 vs INEGI - {id_proyecto}<br>"
              f"<sub>Verde=dentro de {tolerancia_m:.0f}m de la red oficial | Rojo=fuera | "
              f"Azul=red oficial INEGI (capa {layer})" + subtitulo_extra + "</sub>")
    margin_t = geomatica.calcular_margen_top_titulo(titulo)

    fig.update_layout(
        title=titulo,
        scene=dict(xaxis_title="Este [km]", yaxis_title="Norte [km]", zaxis_title="Altitud [msnm]",
                   aspectmode="manual", aspectratio=dict(x=1, y=1, z=0.7),
                   camera=dict(eye=dict(x=1.35, y=-1.35, z=0.8), center=dict(x=0, y=0, z=-0.05))),
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.75)"),
        margin=dict(l=10, r=10, t=margin_t, b=10),
        autosize=True,
    )
    html_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_3d_validacion_p{int(percentil_cauce or cfg.PERCENTIL_CAUCE_HIDROLOGIA)}.html")
    fig.write_html(html_path, div_id="mapa3d", post_script=geomatica.CLICK_ABRE_GOOGLE_EARTH_JS)
    log(f"Mapa 3D de validación: {html_path}")
    return html_path


# ==============================================================================
# --- MODO DEMO: sin archivos externos, geometrías sintéticas deterministas ---
# ==============================================================================
def demo():
    """Prueba la lógica de comparación (vectorización + KD-tree + %
    dentro de tolerancia) con puntos D8 y una línea de referencia
    sintéticos, sin depender de ningún shapefile ni de geomatica.py."""
    log("=== core.validacion_hidrologica --demo (geometrías sintéticas) ===")

    # Línea de referencia sintética: un "río" recto de 2km
    linea_sintetica = LineString([(0, 0), (2000, 500)])
    gdf_ref = gpd.GeoDataFrame(geometry=[linea_sintetica])
    puntos_ref = _muestrear_lineas_a_puntos(gdf_ref)

    # Puntos D8 sintéticos: unos cerca de la línea (deberían validar bien),
    # otros lejos (a propósito, para probar que el % SÍ baja cuando corresponde)
    rng = np.random.default_rng(42)
    cerca = np.column_stack([
        np.linspace(0, 2000, 40),
        np.linspace(0, 500, 40) + rng.normal(0, 5, 40),  # ruido de ±5m alrededor de la línea
    ])
    lejos = np.column_stack([
        rng.uniform(0, 2000, 10),
        rng.uniform(800, 1200, 10),  # deliberadamente separados de la línea
    ])
    puntos_d8 = np.vstack([cerca, lejos])

    tolerancia_m = TOLERANCIA_VALIDACION_HIDRO_M
    resultado = validar_zona(puntos_d8, puntos_ref, tolerancia_m)
    print(f"\n--- Validación D8 vs referencia (demo, {len(cerca)} puntos cercanos + "
          f"{len(lejos)} deliberadamente lejanos) ---")
    for k, v in resultado.items():
        if k in ("distancias_m", "dentro_tolerancia"):
            continue  # arrays por punto -- no imprimir completos, solo para uso interno (ej. mapa 3D)
        print(f"  {k}: {v}")
    print(f"\n  (Se espera ~{len(cerca)}/{len(puntos_d8)} = {len(cerca)/len(puntos_d8)*100:.0f}% "
          f"dentro de tolerancia, ya que 'lejos' se generó a propósito fuera de rango)")
    return resultado


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(description="Validación de cauces D8 contra la Red Hidrográfica INEGI -- Motor Nacional")
    ap.add_argument("--demo", action="store_true", help="Corre con geometrías sintéticas, sin archivos externos")
    ap.add_argument("--geojson", type=str, help="Ruta al GeoJSON del polígono núcleo")
    ap.add_argument("--shapefile-inegi", type=str, help="Ruta al shapefile/gpkg de la Red Hidrográfica de INEGI")
    ap.add_argument("--layer", type=str, default="corriente_ag_l",
                     help="Capa dentro del GeoPackage a usar (default: 'corriente_ag_l', la red LINEAL de cauces del CNIT50k)")
    ap.add_argument("--id-proyecto", type=str, help="Nombre identificador del sitio")
    ap.add_argument("--zonas", type=str, default=None, help="Buffers en metros separados por coma, ej. '0,500,1000'")
    ap.add_argument("--tolerancia-m", type=float, default=None,
                     help=f"Tolerancia de distancia en metros (default: config.TOLERANCIA_VALIDACION_HIDRO_M={TOLERANCIA_VALIDACION_HIDRO_M})")
    ap.add_argument("--percentil-cauce", type=float, default=None,
                     help=f"Percentil de acumulación de flujo D8 para declarar cauce (default: "
                          f"config.PERCENTIL_CAUCE_HIDROLOGIA={PERCENTIL_CAUCE_HIDROLOGIA}). "
                          "Subirlo exige más acumulación de flujo para marcar un cauce -- útil para probar si el ruido "
                          "en zonas planas se debe a un umbral demasiado permisivo.")
    ap.add_argument("--carpeta-salida", type=str, default=None)
    ap.add_argument("--mapa-3d", action="store_true",
                     help="Genera el mapa 3D (D8 verde/rojo + red INEGI en azul) en vez de solo el CSV")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not args.geojson or not args.shapefile_inegi or not args.id_proyecto:
        ap.error("--geojson, --shapefile-inegi y --id-proyecto son obligatorios fuera de --demo")

    zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None

    if args.mapa_3d:
        generar_mapa_3d_validacion(
            geojson_path=args.geojson, shapefile_inegi_path=args.shapefile_inegi, id_proyecto=args.id_proyecto,
            zonas_m=zonas_m, tolerancia_m=args.tolerancia_m, percentil_cauce=args.percentil_cauce,
            carpeta_salida=args.carpeta_salida, layer=args.layer,
        )
        return

    validar_sitio_real(
        geojson_path=args.geojson, shapefile_inegi_path=args.shapefile_inegi, id_proyecto=args.id_proyecto,
        zonas_m=zonas_m, tolerancia_m=args.tolerancia_m, percentil_cauce=args.percentil_cauce,
        carpeta_salida=args.carpeta_salida, layer=args.layer,
    )


if __name__ == "__main__":
    main()

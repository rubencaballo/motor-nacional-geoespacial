#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Deforestación por año (Hansen Global Forest Change) + tendencia de
índices espectrales (NDVI/NDMI/NDWI) por año, para CUALQUIER polígono --
una ANP completa (por zona, núcleo+buffers, igual que geomatica.py/
carbono.py) o la parcela chica de un productor (una sola geometría, sin
zonas) -- INDEPENDIENTE del análisis de riesgo EUDR (core/motor.py usa un
corte binario 2020 para EUDR; este módulo da la serie completa por año,
para cualquier rango, para cualquier polígono).

POR QUÉ ESTE MÓDULO ESTÁ SEPARADO DE core/carbono.py:
    Mismo principio que ya aplicamos en todo el proyecto: un módulo, una
    fuente de datos por responsabilidad. carbono.py habla con Earth Engine
    para BIOMASA (ESA CCI + GEDI); este módulo habla con Earth Engine para
    PÉRDIDA DE COBERTURA (Hansen) y para índices espectrales derivados de
    Sentinel-2 (NDVI/NDMI/NDWI). Si algo falla, se sabe de inmediato cuál
    pieza fue.

TRES ÍNDICES ESPECTRALES POR AÑO -- documentados, no confundidos entre sí:
    - NDVI (B8, B4): salud/vigor de la vegetación.
    - NDMI (B11, B8): estrés hídrico / contenido de humedad de la
      vegetación -- en un script previo de este proyecto (irdcloudramsar.py)
      esto se llamaba por error "NBDI"; es la MISMA fórmula, se corrige el
      nombre aquí porque no es un índice de agua superficial.
    - NDWI (B3, B8, fórmula de McFeeters): agua superficial REAL -- cuerpos
      de agua visibles (humedales, ríos, charcas). Nuevo en este módulo, no
      existía en el script anterior.
    Los tres salen de LA MISMA imagen Sentinel-2 mediana por año, así que
    calcular los tres no cuesta descargas extra.

DATASET HANSEN -- documentado, no escondido:
    UMD/hansen/global_forest_change_2025_v1_13 (config.HANSEN_DATASET,
    mismo dataset que usa el resto del proyecto para EUDR -- consistencia
    interna). Banda 'lossyear': 0 = sin pérdida, 1-N = año de pérdida
    (2000+valor). Cubre desde 2001 -- MUCHO antes que Sentinel-2
    (2016+ para NDVI/NDMI/NDWI confiables, ver DEFORESTACION_ANIO_MIN_SENTINEL2
    en config.py). Por eso la pérdida Hansen y los índices espectrales
    pueden tener rangos de años distintos dentro del mismo reporte -- se
    reporta explícito, nunca se rellena un año sin Sentinel-2 con un valor
    inventado.

CAPA 3D DE "MANCHAS" DE DEFORESTACIÓN (lo que hace este módulo distinto de
solo un CSV):
    Hansen (30m) y el SRTM que ya usa geomatica.py (30m) son resoluciones
    compatibles -- eso permite descargar la banda 'lossyear' recortada al
    sitio y REALINEARLA pixel a pixel a la MISMA malla de terreno que ya
    construye geomatica.py (misma transform, mismo shape), en vez de
    inventar una malla nueva. El resultado: las manchas de deforestación se
    pintan literalmente sobre el terreno 3D ya existente, coloreadas por
    año (amarillo=más antiguo, rojo oscuro=más reciente, misma paleta que
    ya se usaba en irdcloudramsar.py). Ver construir_capa_deforestacion_3d()
    y geomatica.generar_mapa_3d(capas_extra=...).

    OJO: el remuestreo de 'lossyear' a la malla del terreno usa SIEMPRE
    Resampling.nearest, nunca bilinear -- es un código categórico (año),
    no una magnitud continua como la elevación. Ver el docstring de
    _descargar_lossyear_alineado() para el detalle.

Qué NO calcula:
    - Ninguna atribución de causa (agrícola, incendio, tala legal/ilegal) --
      Hansen detecta pérdida de cobertura, no la causa.
    - Ninguna corrección por nubes persistentes en Sentinel-2 más allá del
      filtro de % de nubosidad por escena -- si un año no tiene ninguna
      imagen limpia, se reporta None, nunca se inventa un promedio.
"""

import argparse
import os
from datetime import datetime

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.warp import reproject, Resampling
import requests
from shapely.geometry import mapping

from config import (
    HANSEN_DATASET, DEFORESTACION_ANIO_INICIO_DEFAULT, DEFORESTACION_ANIO_MIN_SENTINEL2,
    DEFORESTACION_NUBOSIDAD_MAX_PCT, ZONAS_ANALISIS_M, CARPETA_SRTM, PERCENTIL_CAUCE_HIDROLOGIA, log,
)


# ==============================================================================
# --- GEOMETRÍA (UTM) -- copia local a propósito, mismo criterio que
#     core/carbono.py: cada módulo aislado, no se importa de geomatica.py/
#     carbono.py (si algo truena, se sabe de inmediato qué pieza fue) ---
# ==============================================================================
def _reproyectores_utm(geom_wgs84_nucleo):
    import pyproj
    centroid = geom_wgs84_nucleo.centroid
    utm_zone = int((centroid.x + 180) / 6) + 1
    utm_crs = f"EPSG:326{utm_zone}"
    a_utm = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True).transform
    a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True).transform
    return a_utm, a_wgs84, utm_crs


def geom_zona_precisa(geom_wgs84_nucleo, buf_m, a_utm, a_wgs84):
    """Buffer EXACTO en metros (vía UTM), igual que geomatica.py/carbono.py
    -- no una aproximación en grados. Devuelve (geometria_wgs84_zona, area_ha)."""
    from shapely.ops import transform as shp_transform

    geom_utm_nucleo = shp_transform(a_utm, geom_wgs84_nucleo)
    geom_utm_zona = geom_utm_nucleo.buffer(buf_m) if buf_m > 0 else geom_utm_nucleo
    area_ha = geom_utm_zona.area / 10000.0
    geom_wgs84_zona = shp_transform(a_wgs84, geom_utm_zona)
    return geom_wgs84_zona, area_ha


# ==============================================================================
# --- HANSEN: PÉRDIDA POR AÑO PARA UN POLÍGONO (agregada -- hectáreas/año) ---
# ==============================================================================
def calcular_deforestacion_anual_gee(ee, geom_wgs84_geojson, anio_inicio=None, anio_fin=None, dataset=None):
    """Pérdida de cobertura Hansen por año, para UN polígono (dict geojson,
    WGS84). Devuelve área total, pérdida total y el historial año->hectáreas
    (incluyendo años sin pérdida como 0.0 -- nunca se omiten, para que una
    serie temporal graficada no tenga huecos silenciosos)."""
    dataset = dataset or HANSEN_DATASET
    anio_inicio = anio_inicio or DEFORESTACION_ANIO_INICIO_DEFAULT
    anio_fin = anio_fin or (datetime.now().year - 1)
    if anio_fin < anio_inicio:
        raise ValueError(f"anio_fin ({anio_fin}) no puede ser menor que anio_inicio ({anio_inicio})")

    aoi = ee.Geometry(geom_wgs84_geojson)
    area_ha = ee.Number(aoi.area(1)).divide(10000).getInfo()

    hansen = ee.Image(dataset).select("lossyear")
    codigo_inicio, codigo_fin = anio_inicio - 2000, anio_fin - 2000
    en_rango = hansen.gte(codigo_inicio).And(hansen.lte(codigo_fin))
    hansen_en_rango = hansen.updateMask(en_rango)

    stats = hansen_en_rango.reduceRegion(
        reducer=ee.Reducer.frequencyHistogram(), geometry=aoi, scale=30, maxPixels=1e9, bestEffort=True,
    ).getInfo()
    hist_raw = stats.get("lossyear") or {}

    # 900 m2/pixel (30x30m de Hansen) -> hectareas
    historial_ha = {2000 + int(k): (v * 900 / 10000) for k, v in hist_raw.items()}
    for anio in range(anio_inicio, anio_fin + 1):
        historial_ha.setdefault(anio, 0.0)

    perdida_total_ha = sum(historial_ha.values())
    perdida_total_pct = (perdida_total_ha / area_ha * 100) if area_ha > 0 else 0.0

    return {
        "area_ha": area_ha, "perdida_total_ha": perdida_total_ha, "perdida_total_pct": perdida_total_pct,
        "historial_ha": historial_ha, "dataset": dataset, "anio_inicio": anio_inicio, "anio_fin": anio_fin,
    }


# ==============================================================================
# --- SENTINEL-2: NDVI/NDMI/NDWI POR AÑO PARA UN POLÍGONO ---
# ==============================================================================
def calcular_indices_anuales_gee(ee, geom_wgs84_geojson, anio_inicio=None, anio_fin=None,
                                  nubosidad_max_pct=None, incluir_ndvi=True, incluir_ndmi=True, incluir_ndwi=True):
    """NDVI/NDMI/NDWI promedio por año (mediana Sentinel-2 SR anual, filtro
    de nubes por escena). Años anteriores a DEFORESTACION_ANIO_MIN_SENTINEL2
    se reportan como None (Sentinel-2 no tiene cobertura confiable ahí) --
    NUNCA se inventa un valor. Si un año SÍ está en rango pero no hay
    ninguna imagen con nubosidad aceptable, también se reporta None,
    explícito en 'n_imagenes_s2'=0."""
    anio_inicio = anio_inicio or DEFORESTACION_ANIO_INICIO_DEFAULT
    anio_fin = anio_fin or (datetime.now().year - 1)
    nubosidad_max_pct = nubosidad_max_pct if nubosidad_max_pct is not None else DEFORESTACION_NUBOSIDAD_MAX_PCT
    aoi = ee.Geometry(geom_wgs84_geojson)

    anio_inicio_s2 = max(anio_inicio, DEFORESTACION_ANIO_MIN_SENTINEL2)
    if anio_inicio_s2 > anio_inicio:
        log(f"Sentinel-2 no tiene cobertura confiable antes de {DEFORESTACION_ANIO_MIN_SENTINEL2} -- "
            f"los años {anio_inicio}-{anio_inicio_s2 - 1} se reportan sin NDVI/NDMI/NDWI (None); la pérdida "
            f"Hansen de esos años SÍ se calcula normal (Hansen cubre desde 2001).", nivel="WARN")

    resultados = {}
    for anio in range(anio_inicio_s2, anio_fin + 1):
        s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
              .filterBounds(aoi).filterDate(f"{anio}-01-01", f"{anio}-12-31")
              .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", nubosidad_max_pct))
              .select(["B3", "B4", "B8", "B11"]))
        n_imagenes = s2.size().getInfo()
        if n_imagenes == 0:
            log(f"  Año {anio}: sin imágenes Sentinel-2 con <{nubosidad_max_pct}% de nubes -- "
                f"se reporta None, no se inventa.", nivel="WARN")
            resultados[anio] = {"ndvi": None, "ndmi": None, "ndwi": None, "n_imagenes_s2": 0}
            continue

        compuesto = s2.median()
        bandas = []
        if incluir_ndvi:
            bandas.append(compuesto.normalizedDifference(["B8", "B4"]).rename("ndvi"))
        if incluir_ndmi:
            bandas.append(compuesto.normalizedDifference(["B11", "B8"]).rename("ndmi"))
        if incluir_ndwi:
            bandas.append(compuesto.normalizedDifference(["B3", "B8"]).rename("ndwi"))

        stats = ee.Image.cat(bandas).reduceRegion(
            reducer=ee.Reducer.mean(), geometry=aoi, scale=10, maxPixels=1e9, bestEffort=True,
        ).getInfo()

        resultados[anio] = {
            "ndvi": stats.get("ndvi") if incluir_ndvi else None,
            "ndmi": stats.get("ndmi") if incluir_ndmi else None,
            "ndwi": stats.get("ndwi") if incluir_ndwi else None,
            "n_imagenes_s2": n_imagenes,
        }
        log(f"  Año {anio}: NDVI={resultados[anio]['ndvi']} NDMI={resultados[anio]['ndmi']} "
            f"NDWI={resultados[anio]['ndwi']} ({n_imagenes} imágenes S2)")

    for anio in range(anio_inicio, anio_inicio_s2):
        resultados[anio] = {"ndvi": None, "ndmi": None, "ndwi": None, "n_imagenes_s2": 0}

    return resultados


# ==============================================================================
# --- ORQUESTADOR TABULAR (por zona -- núcleo + buffers, o una sola geometría) ---
# ==============================================================================
def calcular_deforestacion_por_zona_real(ee, geom_wgs84_nucleo, zonas_m, anio_inicio=None, anio_fin=None,
                                          incluir_indices=True, nubosidad_max_pct=None):
    """Corre Hansen + (opcional) NDVI/NDMI/NDWI para cada zona (núcleo +
    buffers, mismo patrón que carbono.calcular_carbono_por_zona_real). Para
    una parcela chica de un productor, se llama con zonas_m=[0] -- una sola
    zona, sin buffers. Devuelve (df_resumen, df_historial):
      - df_resumen: una fila por zona (área, pérdida total, % del área).
      - df_historial: FORMATO LARGO, una fila por (zona, año) -- se eligió
        largo en vez de ancho porque el rango de años puede variar y un CSV
        ancho con una columna por año sería irregular/frágil de leer."""
    a_utm, a_wgs84, _ = _reproyectores_utm(geom_wgs84_nucleo)
    filas_resumen, filas_historial = [], []

    for buf_m in zonas_m:
        etiqueta = "nucleo" if buf_m == 0 else f"buffer_{buf_m}m"
        geom_zona_wgs84, _ = geom_zona_precisa(geom_wgs84_nucleo, buf_m, a_utm, a_wgs84)

        log(f"Zona {etiqueta} -- calculando pérdida Hansen por año...")
        defo = calcular_deforestacion_anual_gee(ee, mapping(geom_zona_wgs84), anio_inicio, anio_fin)
        log(f"  -> {defo['perdida_total_ha']:.2f} ha perdidas de {defo['area_ha']:.2f} ha totales "
            f"({defo['perdida_total_pct']:.2f}%) entre {defo['anio_inicio']}-{defo['anio_fin']}")

        indices_por_anio = {}
        if incluir_indices:
            log(f"Zona {etiqueta} -- calculando NDVI/NDMI/NDWI por año...")
            indices_por_anio = calcular_indices_anuales_gee(
                ee, mapping(geom_zona_wgs84), defo["anio_inicio"], defo["anio_fin"], nubosidad_max_pct,
            )

        filas_resumen.append({
            "zona": etiqueta, "buffer_m": buf_m, "area_ha": round(defo["area_ha"], 3),
            "anio_inicio": defo["anio_inicio"], "anio_fin": defo["anio_fin"],
            "perdida_total_ha": round(defo["perdida_total_ha"], 3),
            "perdida_total_pct": round(defo["perdida_total_pct"], 3),
            "dataset_hansen": defo["dataset"],
        })

        for anio in range(defo["anio_inicio"], defo["anio_fin"] + 1):
            idx = indices_por_anio.get(anio, {})
            filas_historial.append({
                "zona": etiqueta, "buffer_m": buf_m, "anio": anio,
                "perdida_ha": round(defo["historial_ha"].get(anio, 0.0), 4),
                "ndvi": round(idx["ndvi"], 4) if idx.get("ndvi") is not None else None,
                "ndmi": round(idx["ndmi"], 4) if idx.get("ndmi") is not None else None,
                "ndwi": round(idx["ndwi"], 4) if idx.get("ndwi") is not None else None,
                "n_imagenes_s2": idx.get("n_imagenes_s2", 0),
            })

    return pd.DataFrame(filas_resumen), pd.DataFrame(filas_historial)


def procesar_sitio_real(geojson_path, id_proyecto, zonas_m=None, anio_inicio=None, anio_fin=None,
                         carpeta_salida=None, proyecto_gee=None, incluir_indices=True, nubosidad_max_pct=None):
    """Pipeline tabular completo: inicializa Earth Engine, calcula
    Hansen+índices por zona, guarda los DOS CSV (resumen y historial largo).
    Para una parcela sin zonas, pasa zonas_m=[0]. Devuelve
    (df_resumen, df_historial, csv_resumen_path, csv_historial_path)."""
    import ee
    try:
        if proyecto_gee:
            ee.Initialize(project=proyecto_gee)
        else:
            ee.Initialize()
        log("Earth Engine inicializado.", nivel="OK")
    except Exception as e:
        raise RuntimeError(f"No se pudo inicializar Earth Engine: {e}")

    zonas_m = zonas_m if zonas_m is not None else ZONAS_ANALISIS_M
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    os.makedirs(carpeta_salida, exist_ok=True)

    if not os.path.exists(geojson_path):
        raise FileNotFoundError(f"No se encontró el archivo GeoJSON en: {geojson_path}")

    gdf = gpd.read_file(geojson_path)
    geom_wgs84 = gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all") else gdf.geometry.unary_union

    df_resumen, df_historial = calcular_deforestacion_por_zona_real(
        ee, geom_wgs84, zonas_m, anio_inicio=anio_inicio, anio_fin=anio_fin,
        incluir_indices=incluir_indices, nubosidad_max_pct=nubosidad_max_pct,
    )

    csv_resumen = os.path.join(carpeta_salida, f"deforestacion_resumen_{id_proyecto.lower()}.csv")
    csv_historial = os.path.join(carpeta_salida, f"deforestacion_historial_anual_{id_proyecto.lower()}.csv")
    df_resumen.to_csv(csv_resumen, index=False)
    df_historial.to_csv(csv_historial, index=False)
    log(f"CSV resumen por zona: {csv_resumen}")
    log(f"CSV historial anual (largo, zona+año): {csv_historial}")
    generar_resumen_no_traslapado(df_historial, id_proyecto, carpeta_salida)

    return df_resumen, df_historial, csv_resumen, csv_historial


def generar_resumen_no_traslapado(df_historial, id_proyecto, carpeta_salida):
    """Anexo de lectura rápida a deforestacion_historial_anual_*.csv -- NO lo
    reemplaza. `perdida_ha` ahí es ACUMULATIVA por diseño (buffer_1000m
    incluye completo a buffer_500m, que incluye completo al núcleo -- misma
    convención en todo el proyecto, ver core/carbono.py y
    core/validacion_incendios.py) -- así que buffer_500m y buffer_1000m
    siempre se van a parecer entre sí (cada uno ya trae adentro lo del
    anterior), lo cual a simple vista parece traslape/muñeca rusa aunque sea
    matemáticamente correcto.

    Esta función arma un CSV aparte con `perdida_ha` por ANILLO EXCLUSIVO
    (restando cada buffer del siguiente más grande -- válido porque cada
    buffer_Xm = nucleo.buffer(X), así que geométricamente uno SIEMPRE
    contiene completo al anterior, sin excepción) -- estas franjas nunca se
    traslapan, así que SÍ se pueden sumar directo, con una fila TOTAL por
    año que ya es esa suma real. NDVI/NDMI/NDWI se quedan fuera de este
    resumen a propósito: son promedios de zona, no una cantidad que se
    pueda restar entre dos buffers para obtener el promedio real del
    anillo -- incluirlos aquí sería inventar un número, no reportarlo."""
    faltantes = [c for c in ["zona", "buffer_m", "anio", "perdida_ha"] if c not in df_historial.columns]
    if faltantes:
        log(f"generar_resumen_no_traslapado: faltan columnas {faltantes} -- no se genera el resumen.", nivel="WARN")
        return None

    buffers_ordenados = sorted(df_historial["buffer_m"].unique())
    nombre_anillo = {}
    for i, buf_m in enumerate(buffers_ordenados):
        if i == 0:
            nombre_anillo[buf_m] = "nucleo (0m)" if buf_m == 0 else f"anillo_0-{buf_m}m"
        else:
            nombre_anillo[buf_m] = f"anillo_{buffers_ordenados[i - 1]}-{buf_m}m"

    pivot = df_historial.pivot_table(index="anio", columns="buffer_m", values="perdida_ha", aggfunc="first")
    filas = []
    for anio in sorted(pivot.index):
        acumulado_prev, total_anio = 0.0, 0.0
        for buf_m in buffers_ordenados:
            valor = pivot.loc[anio, buf_m]
            valor = 0.0 if pd.isna(valor) else float(valor)
            valor_anillo = round(valor - acumulado_prev, 4)
            acumulado_prev = valor
            total_anio += valor_anillo
            filas.append({"anillo": nombre_anillo[buf_m], "anio": anio, "perdida_ha": valor_anillo})
        filas.append({"anillo": "TOTAL (suma sin traslape)", "anio": anio, "perdida_ha": round(total_anio, 4)})

    df_out = pd.DataFrame(filas).sort_values(["anio", "anillo"]).reset_index(drop=True)

    # --- GRAN TOTAL: una fila por anillo (+ el total general) sumando TODOS los años del historial ---
    # "anio" queda como texto en estas filas (ej. "TOTAL 2010-2025") -- no es un año real, es el
    # periodo completo. El resto de la columna sigue siendo años enteros, es intencional: en un CSV
    # no hay problema en mezclar ambos, y así queda claro a simple vista cuál fila es cuál.
    anio_min, anio_max = df_historial["anio"].min(), df_historial["anio"].max()
    etiqueta_periodo = f"TOTAL {anio_min}-{anio_max}"
    filas_gran_total = [
        {"anillo": anillo, "anio": etiqueta_periodo,
         "perdida_ha": round(df_out.loc[df_out["anillo"] == anillo, "perdida_ha"].sum(), 4)}
        for anillo in list(nombre_anillo.values()) + ["TOTAL (suma sin traslape)"]
    ]
    df_out = pd.concat([df_out, pd.DataFrame(filas_gran_total)], ignore_index=True)

    nombre_out = f"deforestacion_resumen_sin_traslape_{id_proyecto.lower()}.csv"
    csv_path = os.path.join(carpeta_salida, nombre_out)
    df_out.to_csv(csv_path, index=False)
    log(f"Resumen SIN traslape (perdida_ha por anillo exclusivo, sumable) guardado en: {csv_path}", nivel="OK")
    return csv_path


# ==============================================================================
# --- MAPA 3D: MANCHAS DE DEFORESTACIÓN POR AÑO ENCIMA DEL TERRENO ---
# ==============================================================================
def _descargar_lossyear_alineado(ee, geom_wgs84_visual, transform_ref, shape_ref, utm_crs, dataset=None,
                                  carpeta_tmp=None):
    """Descarga la banda 'lossyear' de Hansen recortada al polígono visual,
    y la REALINEA exacto a la misma malla (transform_ref/shape_ref) que ya
    usa el terreno 3D de geomatica.py -- así cada mancha de deforestación
    cae físicamente sobre el pixel de terreno correcto, sin desplazamiento
    ni tener que reconstruir una malla nueva.

    OJO CRÍTICO: el remuestreo usa Resampling.nearest, NUNCA bilinear.
    'lossyear' es un CÓDIGO categórico (0=sin pérdida, 1=2001 ... 24=2024),
    no una magnitud continua como la elevación -- interpolar entre el
    código 5 y el código 20 daría un 'año 12.5' sin significado, y peor:
    mezclaría 'sin pérdida' (0) con años reales justo en el borde de cada
    mancha. geomatica.py SÍ usa bilinear para el DEM porque ahí la
    elevación sí es continua -- ese criterio NO aplica aquí."""
    import tempfile

    dataset = dataset or HANSEN_DATASET
    carpeta_tmp = carpeta_tmp or tempfile.gettempdir()
    os.makedirs(carpeta_tmp, exist_ok=True)

    aoi = ee.Geometry(geom_wgs84_visual)
    img = ee.Image(dataset).select("lossyear").clip(aoi)
    url = img.getDownloadURL({"region": aoi, "scale": 30, "crs": "EPSG:4326", "format": "GEO_TIFF"})
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    tif_crudo = os.path.join(carpeta_tmp, "temp_hansen_lossyear_crudo.tif")
    with open(tif_crudo, "wb") as f:
        f.write(r.content)

    rows, cols = shape_ref
    lossyear_alineado = np.zeros((rows, cols), dtype=np.uint8)
    with rasterio.open(tif_crudo) as src:
        reproject(
            source=rasterio.band(src, 1), destination=lossyear_alineado,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=transform_ref, dst_crs=utm_crs,
            resampling=Resampling.nearest,
        )
    os.remove(tif_crudo)
    return lossyear_alineado


def construir_capa_deforestacion_3d(lossyear_alineado, hidrologia, anio_inicio=None, anio_fin=None, utm_crs=None):
    """Construye el trazo go.Scatter3d de las manchas de deforestación,
    coloreado por año (amarillo=más antiguo, rojo oscuro=más reciente --
    misma paleta que ya se usaba en irdcloudramsar.py), listo para pasarse
    como `capas_extra` a geomatica.generar_mapa_3d(). `lossyear_alineado`
    debe venir YA alineado a la malla de `hidrologia` (ver
    _descargar_lossyear_alineado()). Devuelve None (sin trazo) si no hay
    pérdida en el rango de años dentro de la zona visual -- eso es una
    buena noticia, no un error, y no se agrega una capa vacía al mapa."""
    import plotly.graph_objects as go
    from core.geomatica import calcular_grid_latlon

    anio_inicio = anio_inicio or DEFORESTACION_ANIO_INICIO_DEFAULT
    anio_fin = anio_fin or (datetime.now().year - 1)
    codigo_inicio, codigo_fin = anio_inicio - 2000, anio_fin - 2000

    Z_raw = hidrologia["Z_raw"]
    Z_smooth = hidrologia.get("Z_smooth", Z_raw)
    pw_v, ph_v = hidrologia["pw_v"], hidrologia["ph_v"]
    transform = hidrologia["transform"]
    rows, cols = Z_raw.shape

    mascara_perdida = (lossyear_alineado >= codigo_inicio) & (lossyear_alineado <= codigo_fin)
    # OJO: cargar_dem_utm() recorta el SRTM al BOUNDING BOX (rectángulo) que envuelve la zona
    # visual, no al círculo/buffer real -- así que lossyear_alineado trae también las esquinas
    # del rectángulo, FUERA del buffer circular. Sin este filtro, esas esquinas (con pérdida
    # Hansen real pero fuera del área que en realidad se está analizando) se cuelan en el mapa
    # y en el total de hectáreas -- así se detectó un caso real: el subtítulo decía 1316 ha
    # cuando el buffer_1000m real (el mismo universo que usa el historial anual por zona) eran
    # 771 ha. zona_de_pixel ya distingue esto ("fuera" = fuera de cualquier zona real).
    if "zona_de_pixel" in hidrologia:
        mascara_perdida &= (hidrologia["zona_de_pixel"] != "fuera")
    fy, fx = np.where(mascara_perdida)
    if len(fx) > 0:
        validos = ~np.isnan(Z_raw[fy, fx]) & ~np.isnan(Z_smooth[fy, fx])
        fy, fx = fy[validos], fx[validos]

    if len(fx) == 0:
        log("Sin pérdida Hansen detectada en el rango de años pedido dentro de la zona visual -- "
            "no se agrega capa de deforestación al mapa 3D (buena noticia, no un error).")
        return None

    anios_reales = 2000 + lossyear_alineado[fy, fx].astype(int)
    x_km = fx * pw_v / 1000.0
    y_km = (rows - 1 - fy) * ph_v / 1000.0
    # offset chico hacia arriba, mismo criterio que las demás capas de puntos
    # de este proyecto (cauces, puntos más altos): para que no se entierren
    # visualmente en la malla del terreno.
    z_km = Z_smooth[fy, fx] + max(pw_v, ph_v) * 0.6

    customdata_cols = [Z_raw[fy, fx], anios_reales]
    hovertemplate = "Deforestación<br>Año detectado: %{customdata[1]}<br>Altitud: %{customdata[0]:.0f} msnm"
    if utm_crs:
        lat_grid, lon_grid = calcular_grid_latlon(transform, utm_crs, rows, cols)
        customdata_cols += [lat_grid[fy, fx], lon_grid[fy, fx]]
        hovertemplate += "<br>Lat: %{customdata[2]:.5f}<br>Lon: %{customdata[3]:.5f}"
    hovertemplate += "<extra></extra>"
    customdata = np.column_stack(customdata_cols)

    # Ticks explícitos del colorbar: el tickmode automático de Plotly elige un
    # paso "bonito" (ej. cada 2 años) desde anio_inicio -- si anio_fin cae en
    # un año que ese paso no toca (ej. rango par-a-impar 2010-2025), el último
    # año NUNCA aparece impreso en la barra, aunque el punto sí esté pintado
    # con el color extremo correcto (cmax). Visualmente parece que el dato
    # "no cuadra" con la escala cuando en realidad el hover es el valor real
    # y la escala solo se saltó ese número. Se fuerza aquí que anio_inicio y
    # anio_fin siempre queden etiquetados, para no dejarle esa duda a quien
    # vea el mapa sin poder preguntar.
    tickvals = list(range(anio_inicio, anio_fin, 2)) + [anio_fin]

    capa = go.Scatter3d(
        x=x_km, y=y_km, z=z_km, mode="markers",
        marker=dict(
            size=2.6, color=anios_reales,
            colorscale=[[0.0, "#FFFF00"], [0.33, "#FFA500"], [0.66, "#FF4500"], [1.0, "#8B0000"]],
            cmin=anio_inicio, cmax=anio_fin, opacity=0.85,
            colorbar=dict(title="Año de<br>pérdida", x=1.15, tickvals=tickvals),
        ),
        customdata=customdata, hovertemplate=hovertemplate,
        name=f"Deforestación Hansen {anio_inicio}-{anio_fin}",
    )
    log(f"Capa de deforestación 3D: {len(fx)} píxeles de pérdida entre {anio_inicio}-{anio_fin} "
        f"(~{len(fx) * pw_v * ph_v / 10000:.2f} ha en la zona visual).")
    return capa


def generar_mapa_3d_deforestacion(geojson_path, id_proyecto, zonas_m=None, anio_inicio=None, anio_fin=None,
                                   percentil_cauce=None, carpeta_salida=None, carpeta_srtm=None, proyecto_gee=None):
    """Pipeline del mapa 3D con las manchas de deforestación por año encima
    del MISMO terreno que ya usa core/geomatica.py -- reusa su SRTM (si ya
    está en caché, no vuelve a descargar) y su hidrología D8 (cauces +
    puntos más altos), y solo agrega la capa nueva de pérdida Hansen
    alineada pixel a pixel a esa misma malla. Devuelve la ruta del HTML."""
    from core import geomatica
    import pyproj
    from shapely.ops import transform as shp_transform
    import ee

    try:
        if proyecto_gee:
            ee.Initialize(project=proyecto_gee)
        else:
            ee.Initialize()
        log("Earth Engine inicializado.", nivel="OK")
    except Exception as e:
        raise RuntimeError(f"No se pudo inicializar Earth Engine: {e}")

    zonas_m = zonas_m if zonas_m is not None else ZONAS_ANALISIS_M
    percentil_cauce = percentil_cauce if percentil_cauce is not None else PERCENTIL_CAUCE_HIDROLOGIA
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    carpeta_srtm = carpeta_srtm or CARPETA_SRTM
    os.makedirs(carpeta_salida, exist_ok=True)
    anio_inicio = anio_inicio or DEFORESTACION_ANIO_INICIO_DEFAULT
    anio_fin = anio_fin or (datetime.now().year - 1)

    geom_utm_nucleo, dst_array, meta_utm, utm_crs = geomatica.cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        percentil_cauce, carpeta_srtm, id_proyecto,
    )

    a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True).transform
    geom_visual_utm = geom_utm_nucleo.buffer(max(zonas_m))
    geom_visual_wgs84 = shp_transform(a_wgs84, geom_visual_utm)

    log(f"Descargando Hansen lossyear alineado a la malla del terreno ({anio_inicio}-{anio_fin})...")
    lossyear_alineado = _descargar_lossyear_alineado(
        ee, mapping(geom_visual_wgs84), hidrologia["transform"], hidrologia["Z_raw"].shape, utm_crs,
        carpeta_tmp=carpeta_srtm,
    )

    capa = construir_capa_deforestacion_3d(lossyear_alineado, hidrologia, anio_inicio, anio_fin, utm_crs=utm_crs)
    capas_extra = [capa] if capa is not None else []

    codigo_inicio, codigo_fin = anio_inicio - 2000, anio_fin - 2000
    mascara_total = (lossyear_alineado >= codigo_inicio) & (lossyear_alineado <= codigo_fin)
    # Mismo filtro que construir_capa_deforestacion_3d: sin esto, el total cuenta también las
    # esquinas del bounding box rectangular (fuera del buffer circular real) -- ver ese docstring.
    mascara_total &= (hidrologia["zona_de_pixel"] != "fuera")
    ha_perdidas_visual = mascara_total.sum() * hidrologia["pw_v"] * hidrologia["ph_v"] / 10000.0
    subtitulo = (f"Deforestación Hansen {anio_inicio}-{anio_fin}: {ha_perdidas_visual:.2f} ha perdidas dentro "
                 f"del área visual (buffer {max(zonas_m)}m) -- amarillo=más antiguo, rojo oscuro=más reciente")

    # Aviso honesto de que esta cifra es específica de ESTE mapa (conteo de píxeles ya remuestreados a la
    # malla del terreno, ver docstring de generar_desglose_anual_visual) y casi seguro no va a coincidir con
    # el total "oficial" de la plataforma (deforestacion_resumen_sin_traslape_*.csv, consulta vectorial directa
    # en Earth Engine, la que alimenta core/carbono_perdida.py) -- típicamente unos pocos puntos porcentuales
    # de diferencia por el remuestreo, no un error. Si ese CSV ya existe (se corrió procesar_sitio_real() antes
    # o después), se cita su cifra real en vez de dejar la duda -- si no existe todavía, el aviso queda
    # genérico, sin inventar un número.
    ruta_resumen_oficial = os.path.join(carpeta_salida, f"deforestacion_resumen_sin_traslape_{id_proyecto.lower()}.csv")
    if os.path.exists(ruta_resumen_oficial):
        try:
            df_oficial = pd.read_csv(ruta_resumen_oficial)
            fila_oficial = df_oficial[(df_oficial["anillo"] == "TOTAL (suma sin traslape)")
                                       & df_oficial["anio"].astype(str).str.startswith("TOTAL ")]
            if not fila_oficial.empty:
                ha_oficial = fila_oficial["perdida_ha"].iloc[0]
                periodo_oficial = fila_oficial["anio"].iloc[0]
                # Título de Plotly: NO hace word-wrap solo, una línea larga se corta en el borde del gráfico
                # ("aquí se pierden los números", mismo problema ya resuelto antes para otros mapas) -- por eso
                # aquí se parte en líneas cortas con <br> en vez de una sola frase larga.
                subtitulo += (f"<br>Esta cifra es del conteo de píxeles de este mapa (no el total oficial)."
                              f"<br>Total oficial de la plataforma ({periodo_oficial}): {ha_oficial:,.1f} ha "
                              "-- ver mapa de CO2e liberado.")
        except Exception as e:
            log(f"No se pudo leer {ruta_resumen_oficial} para citar el total oficial en el subtítulo: {e}",
                nivel="WARN")
    else:
        subtitulo += ("<br>Esta cifra es del conteo de píxeles de este mapa, no la medición oficial de la "
                       "plataforma."
                       "<br>La oficial sale de deforestacion_resumen_sin_traslape_*.csv (consulta directa en "
                       "Earth Engine).")

    generar_desglose_anual_visual(lossyear_alineado, hidrologia, anio_inicio, anio_fin, id_proyecto, carpeta_salida)

    html_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_3d_deforestacion.html")
    titulo_base = f"{id_proyecto} -- Modelo de terreno 3D + deforestación Hansen {anio_inicio}-{anio_fin}"
    geomatica.generar_mapa_3d(
        hidrologia, id_proyecto, html_path, subtitulo=subtitulo, utm_crs=utm_crs, capas_extra=capas_extra,
        titulo_base=titulo_base,
    )
    log(f"Mapa 3D con deforestación por año: {html_path}")
    return html_path


def generar_desglose_anual_visual(lossyear_alineado, hidrologia, anio_inicio, anio_fin, id_proyecto, carpeta_salida):
    """Desglose EXACTO, año por año y por anillo de píxel EXCLUSIVO
    (hidrologia['zona_de_pixel'] -- nucleo/buffer_500m/buffer_1000m, SIN
    traslape, ver core/geomatica.py: cada píxel se asigna a un solo anillo,
    el más chico que lo contiene -- misma convención que anillo_exclusivo
    en generar_resumen_no_traslapado, aunque aquí reutiliza los nombres
    'buffer_Xm' porque así los etiquetó geomatica.py, no 'anillo_...').

    Es el desglose de LO QUE ESTE MISMO MAPA 3D ESTÁ PINTANDO, píxel a
    píxel, sobre la malla ya descargada y alineada -- por diseño, casi
    seguro NO va a sumar exacto igual al total de deforestacion_resumen_*.csv
    (ese usa una consulta de polígono precisa a Earth Engine, no esta malla
    rasterizada localmente) -- la diferencia típica observada en este
    proyecto es de unos pocos puntos porcentuales (ver docstring de
    generar_mapa_3d_deforestacion). No es un error de conteo doble: son dos
    pipelines independientes sobre el mismo dataset Hansen, y este es el
    único de los dos que corresponde EXACTO a los colores que se ven en el
    mapa 3D."""
    pw_ha = hidrologia["pw_v"] * hidrologia["ph_v"] / 10000.0
    zona_px = hidrologia["zona_de_pixel"]
    anillos_presentes = sorted(a for a in np.unique(zona_px) if a != "fuera")

    filas = []
    for anio in range(anio_inicio, anio_fin + 1):
        codigo = anio - 2000
        mascara_anio = (lossyear_alineado == codigo)
        total_anio = 0.0
        for anillo in anillos_presentes:
            ha = round(float((mascara_anio & (zona_px == anillo)).sum() * pw_ha), 4)
            filas.append({"anio": anio, "anillo": anillo, "perdida_ha_visual": ha})
            total_anio += ha
        filas.append({"anio": anio, "anillo": "TOTAL (area visual, sin traslape)", "perdida_ha_visual": round(total_anio, 4)})

    df_out = pd.DataFrame(filas)
    etiqueta_periodo = f"TOTAL {anio_inicio}-{anio_fin}"
    filas_gran_total = [
        {"anio": etiqueta_periodo, "anillo": a,
         "perdida_ha_visual": round(df_out.loc[df_out["anillo"] == a, "perdida_ha_visual"].sum(), 4)}
        for a in list(anillos_presentes) + ["TOTAL (area visual, sin traslape)"]
    ]
    df_out = pd.concat([df_out, pd.DataFrame(filas_gran_total)], ignore_index=True)

    csv_path = os.path.join(carpeta_salida, f"deforestacion_desglose_anual_visual_{id_proyecto.lower()}.csv")
    df_out.to_csv(csv_path, index=False)
    log(f"Desglose anual EXACTO del mapa 3D (píxel a píxel, por anillo sin traslape) guardado en: {csv_path}",
        nivel="OK")
    return csv_path


# ==============================================================================
# --- MODO DEMO: sin Earth Engine, sin red -- mismo espíritu que --demo en
#     el resto del proyecto ---
# ==============================================================================
def _lossyear_sintetico(shape, anio_inicio, anio_fin, semilla=11):
    """Array 'lossyear' sintético (determinista, sin red): algunas manchas
    circulares de 'pérdida' repartidas en distintos años, para poder probar
    la alineación + el mapa 3D sin depender de Earth Engine."""
    rng = np.random.default_rng(semilla)
    rows, cols = shape
    lossyear = np.zeros(shape, dtype=np.uint8)
    anios_disponibles = list(range(anio_inicio, anio_fin + 1))
    n_manchas = min(4, len(anios_disponibles))
    anios_elegidos = rng.choice(anios_disponibles, size=n_manchas, replace=False)
    yy, xx = np.mgrid[0:rows, 0:cols]
    for anio in anios_elegidos:
        cy, cx = rng.integers(0, rows), rng.integers(0, cols)
        radio = rng.integers(4, 10)
        mancha = ((yy - cy) ** 2 + (xx - cx) ** 2) <= radio ** 2
        lossyear[mancha] = int(anio) - 2000
    return lossyear


def demo():
    """Prueba TODA la lógica que no depende de red/Earth Engine: la parte
    tabular usa valores sintéticos directos (igual que core.carbono.demo(),
    que tampoco llama a GEE en su demo), y la capa 3D reusa el DEM
    sintético de geomatica._dem_sintetico() + un array 'lossyear' sintético
    propio -- así se prueba la alineación pixel a pixel y el mapa 3D
    completo sin descargar nada ni tocar Earth Engine."""
    log("=== core.deforestacion --demo (sin Earth Engine, valores sintéticos) ===")
    anio_inicio, anio_fin = 2018, 2023

    rng = np.random.default_rng(5)
    perdida_sintetica_ha = {anio: round(float(rng.uniform(0, 3.0)), 3) for anio in range(anio_inicio, anio_fin + 1)}
    ndvi_sintetico = {anio: round(float(0.75 - 0.02 * (anio - anio_inicio) + rng.normal(0, 0.02)), 3)
                       for anio in range(anio_inicio, anio_fin + 1)}
    print("\n--- Pérdida Hansen por año (demo, sintético) ---")
    for anio, ha in perdida_sintetica_ha.items():
        print(f"  {anio}: {ha} ha  |  NDVI promedio sintético: {ndvi_sintetico[anio]}")
    print(f"  TOTAL sintético: {sum(perdida_sintetica_ha.values()):.3f} ha")

    try:
        from core import geomatica
        dst_array, meta_utm, geom_utm_nucleo, utm_crs = geomatica._dem_sintetico()
        zonas_m = [0, 300, 600]
        id_proyecto = "DEMO_DEFORESTACION"
        carpeta_tmp = os.path.expanduser("~/resultados_demo_deforestacion")
        os.makedirs(carpeta_tmp, exist_ok=True)

        hidrologia = geomatica.calcular_hidrologia_d8(
            dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
            22, carpeta_tmp, id_proyecto,
        )
        lossyear_sint = _lossyear_sintetico(hidrologia["Z_raw"].shape, anio_inicio, anio_fin)
        capa = construir_capa_deforestacion_3d(lossyear_sint, hidrologia, anio_inicio, anio_fin, utm_crs=utm_crs)
        capas_extra = [capa] if capa is not None else []

        html_path = os.path.join(carpeta_tmp, f"{id_proyecto.lower()}_3d_deforestacion.html")
        subtitulo = f"Deforestación sintética {anio_inicio}-{anio_fin} (demo -- NO son datos reales de Hansen)"
        titulo_base = f"{id_proyecto} -- Modelo de terreno 3D + deforestación (DEMO sintético)"
        geomatica.generar_mapa_3d(
            hidrologia, id_proyecto, html_path, subtitulo=subtitulo, utm_crs=utm_crs, capas_extra=capas_extra,
            titulo_base=titulo_base,
        )
        log(f"Mapa 3D demo (con manchas de deforestación sintéticas) generado en: {html_path}")
    except ImportError as e:
        log(f"pysheds/plotly no instalado, se omite la parte 3D del demo: {e}", nivel="WARN")

    return perdida_sintetica_ha


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Deforestación Hansen por año + tendencia NDVI/NDMI/NDWI, para cualquier polígono "
                    "(ANP o parcela) -- Motor Nacional"
    )
    ap.add_argument("--demo", action="store_true", help="Corre con valores sintéticos, sin Earth Engine ni red")
    ap.add_argument("--geojson", type=str, help="Ruta al GeoJSON del polígono (ANP completa o parcela de un productor)")
    ap.add_argument("--id-proyecto", type=str, help="Nombre identificador (para nombres de archivo)")
    ap.add_argument("--zonas", type=str, default=None,
                     help="Buffers en metros separados por coma, ej. '0,500,1000'. Usa '0' solo (una parcela sin buffers)")
    ap.add_argument("--anio-inicio", type=int, default=None,
                     help=f"default: config.DEFORESTACION_ANIO_INICIO_DEFAULT={DEFORESTACION_ANIO_INICIO_DEFAULT}")
    ap.add_argument("--anio-fin", type=int, default=None, help="default: año actual - 1")
    ap.add_argument("--proyecto-gee", type=str, default=None, help="ID de proyecto de Google Cloud para ee.Initialize(project=...)")
    ap.add_argument("--sin-indices", action="store_true", help="No calcular NDVI/NDMI/NDWI, solo Hansen (más rápido)")
    ap.add_argument("--nubosidad-max", type=float, default=None,
                     help=f"default: config.DEFORESTACION_NUBOSIDAD_MAX_PCT={DEFORESTACION_NUBOSIDAD_MAX_PCT}")
    ap.add_argument("--percentil-cauce", type=float, default=None)
    ap.add_argument("--carpeta-srtm", type=str, default=None)
    ap.add_argument("--carpeta-salida", type=str, default=None)
    ap.add_argument("--mapa-3d", action="store_true",
                     help="Genera el mapa 3D con las manchas de deforestación por año encima del terreno "
                          "(requiere pysheds/plotly; descarga el SRTM del sitio si no está en caché)")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not args.geojson or not args.id_proyecto:
        ap.error("--geojson y --id-proyecto son obligatorios fuera de --demo")

    zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
    carpeta_salida = args.carpeta_salida or os.path.expanduser(f"~/resultados_{args.id_proyecto.lower()}")

    if args.mapa_3d:
        generar_mapa_3d_deforestacion(
            geojson_path=args.geojson, id_proyecto=args.id_proyecto, zonas_m=zonas_m,
            anio_inicio=args.anio_inicio, anio_fin=args.anio_fin, percentil_cauce=args.percentil_cauce,
            carpeta_salida=carpeta_salida, carpeta_srtm=args.carpeta_srtm, proyecto_gee=args.proyecto_gee,
        )
        return

    procesar_sitio_real(
        geojson_path=args.geojson, id_proyecto=args.id_proyecto, zonas_m=zonas_m,
        anio_inicio=args.anio_inicio, anio_fin=args.anio_fin, carpeta_salida=carpeta_salida,
        proyecto_gee=args.proyecto_gee, incluir_indices=not args.sin_indices, nubosidad_max_pct=args.nubosidad_max,
    )


if __name__ == "__main__":
    main()

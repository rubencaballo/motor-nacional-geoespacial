#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Estimación de carbono/CO2e por zona a partir de DOS fuentes satelitales de
biomasa (Earth Engine), por polígono: ESA CCI Above-Ground Biomass v6.0
(satelital, continuo) y GEDI L4A (mediciones LiDAR directas por huella,
dispersas) -- reportadas lado a lado en la misma fila del CSV para poder
comparar "lo local" (GEDI) contra "lo satelital" (ESA CCI).

POR QUÉ ESTE MÓDULO ESTÁ SEPARADO DE core/geomatica.py:
    geomatica.py trabaja con SRTM (elevación/pendiente) y no habla con Earth
    Engine. Este módulo sí habla con Earth Engine, igual que
    core/gee_backend.py -- mezclar ambas fuentes en un solo archivo dificulta
    saber qué falló si algo truena (¿fue el DEM? ¿fue GEE?), y ya nos costó
    varios bugs reales mezclar cosas que debían estar aisladas (ver REGLA DE
    ORO en geomatica.py). Este módulo se mantiene independiente por la misma
    razón.

DATASET 1 -- ESA CCI (satelital, cobertura completa) -- documentadas, no escondidas:
    ESA/CCI/Above_Ground_Biomass/V6_0 (Santoro & Cartus, 2025), catálogo
    OFICIAL de Earth Engine (developers.google.com/earth-engine/datasets),
    100m de resolución. Años disponibles: 2007, 2010, 2015-2022. Se accede
    como imagen individual por año: 'ESA/CCI/Above_Ground_Biomass/V6_0/AAAA'.

    - OJO CON LAS UNIDADES -- esto cambió respecto a la versión anterior de
      este módulo (que usaba NASA/ORNL, ya en Mg de carbono/ha): este
      dataset da la banda 'agb' en Mg de BIOMASA seca/ha (peso seco de
      tallo, corteza, ramas -- no en carbono). Por eso aquí la conversión a
      CO2e es en DOS pasos, no uno:
        1) carbono (Mg C/ha) = biomasa (Mg/ha) x FRACCION_CARBONO_BIOMASA (~0.47)
        2) CO2e (t) = area_ha x carbono (Mg C/ha) x FACTOR_C_A_CO2 (44/12)
      Aplicar solo el paso 2 sobre la biomasa directamente (sin el 0.47)
      DUPLICARIA el CO2e reportado -- ese es exactamente el tipo de error
      de unidades que ya nos costó tiempo detectar antes en este proyecto.
    - La banda de incertidumbre se llama 'agb_sd' en este dataset (no
      'agb_uncertainty' como en NASA/ORNL) -- nombre distinto, hay que
      usarlo tal cual o la consulta a Earth Engine falla.
    - Resolución 100m (más fina que los 300m del dataset anterior).
    - El año se fija explícitamente (config.CARBONO_ANIO, default 2022 --
      el más reciente disponible en agosto 2026). Sigue siendo una foto de
      un año concreto, no una serie temporal ni el "año en curso".

DATASET 2 -- GEDI L4A (local, huellas LiDAR dispersas) -- documentadas, no escondidas:
    LARSE/GEDI/GEDI04_A_002_MONTHLY, catálogo OFICIAL de Earth Engine.
    Bandas 'agbd' (Mg/ha, biomasa -- MISMA unidad que ESA CCI, se convierte
    con la misma fórmula de dos pasos) y 'agbd_se' (error estándar).
    ~25m de resolución nominal por huella.

    - A DIFERENCIA de ESA CCI, esto NO es un raster continuo -- son
      disparos LiDAR reales a lo largo de la órbita del ISS (mission
      2019-2023). Un predio chico puede tener varias huellas, pocas, o
      NINGUNA. Cuando no hay huellas dentro del polígono, se reporta
      'gedi_n_muestras=0' y los valores como None -- nunca se inventa un
      número ni se hace pasar por un cero real.
    - 'gedi_n_muestras' indica cuántos pixeles de 25m con huella real
      cayeron en el polígono -- un promedio de 1 huella no tiene la misma
      confiabilidad que uno de 40; siempre revisar este número antes de
      interpretar el AGBD de GEDI como representativo de toda la zona.
    - Cobertura: 51.6°N a 51.6°S (México completo está dentro de rango).
    - No es más "correcto" que ESA CCI por ser LiDAR directo -- cada huella
      SÍ mide estructura vertical real, pero son puntos dispersos que
      pueden no representar la variabilidad de todo el polígono. La
      comparación entre ambos es justo el punto: si concuerdan, da más
      confianza; si difieren mucho, hay que investigar por qué (¿el predio
      tiene alta variabilidad interna? ¿las pocas huellas GEDI cayeron en
      un parche atípico?).

Qué SÍ calcula:
    - Biomasa aérea promedio (Mg/ha) y su incertidumbre, por zona, de
      AMBAS fuentes en la misma fila.
    - Carbono equivalente (Mg C/ha) y CO2 equivalente total (t) por zona,
      con la incertidumbre propagada, para cada fuente por separado.

Qué NO calcula:
    - Biomasa subterránea, necromasa, carbono del suelo -- el dataset solo
      cubre biomasa aérea viva.
    - Ninguna proyección a futuro ni tendencia -- es una foto de un año fijo.
"""

import argparse
import os

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, mapping

from config import (
    CARBONO_DATASET_GEE, CARBONO_ANIO, FRACCION_CARBONO_BIOMASA, FACTOR_C_A_CO2,
    ZONAS_ANALISIS_M, PERCENTIL_CAUCE_HIDROLOGIA, CARPETA_SRTM,
    GEDI_DATASET_L4A, GEDI_ESCALA_M, log,
)


# ==============================================================================
# --- FUNCIÓN PURA: recibe `ee` ya inicializado, no lo importa ni lo inicializa
#     ella misma -- mismo patrón que core/gee_backend.py ---
# ==============================================================================
def medir_biomasa_zona_gee(ee, geom_wgs84_geojson, dataset=None, anio=None):
    """Consulta el dataset real de biomasa aérea en GEE para UN polígono (en
    coordenadas WGS84, como dict geojson) y devuelve biomasa promedio +
    incertidumbre, en Mg/ha (biomasa seca, NO carbono todavía). La
    conversión a carbono/CO2e se hace aparte, en biomasa_a_co2e(), para
    que quede claro en qué paso se aplica cada factor."""
    dataset = dataset or CARBONO_DATASET_GEE
    anio = anio or CARBONO_ANIO
    aoi = ee.Geometry(geom_wgs84_geojson)
    img = ee.Image(f"{dataset}/{anio}")

    stats = img.select(["agb", "agb_sd"]).reduceRegion(
        reducer=ee.Reducer.mean(), geometry=aoi, scale=100, maxPixels=1e9, bestEffort=True,
    ).getInfo()

    agb_mgha = stats.get("agb")
    agb_sd_mgha = stats.get("agb_sd")
    if agb_mgha is None:
        log(f"Sin datos de biomasa para esta zona en {dataset}/{anio} (fuera de cobertura, agua, o sin bosque).", nivel="WARN")
    return {"agb_mgha": agb_mgha, "agb_sd_mgha": agb_sd_mgha}


def medir_biomasa_gedi_zona(ee, geom_wgs84_geojson, dataset=None, escala_m=None):
    """Consulta biomasa GEDI L4A para UN polígono -- a diferencia de ESA CCI
    (un raster continuo, con valor en cada pixel siempre), GEDI son
    disparos LiDAR reales a lo largo de la órbita: en un predio chico
    puede haber varias huellas, pocas, o NINGUNA. Si no hay huellas dentro
    del polígono, se reporta explícitamente como None -- NUNCA se rellena
    con un valor inventado ni se hace pasar por un cero real.

    Devuelve también 'gedi_n_muestras': cuántos pixeles de 25m con datos
    reales de huella cayeron dentro del polígono -- para que quede claro
    qué tan confiable es el promedio (1 huella no es lo mismo que 40)."""
    dataset = dataset or GEDI_DATASET_L4A
    escala_m = escala_m or GEDI_ESCALA_M
    aoi = ee.Geometry(geom_wgs84_geojson)

    coleccion = ee.ImageCollection(dataset).filterBounds(aoi).select(["agbd", "agbd_se"])
    img_media = coleccion.mean()

    reducer = ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True)
    stats = img_media.reduceRegion(
        reducer=reducer, geometry=aoi, scale=escala_m, maxPixels=1e9, bestEffort=True,
    ).getInfo()

    agbd_mgha = stats.get("agbd_mean")
    agbd_se_mgha = stats.get("agbd_se_mean")
    n_muestras = stats.get("agbd_count") or 0

    if agbd_mgha is None or n_muestras == 0:
        log(f"Sin huellas GEDI dentro de esta zona ({dataset}) -- el predio es chico o la "
            f"órbita nunca pasó exactamente por aquí. No se inventa un valor.", nivel="WARN")
        return {"agbd_mgha": None, "agbd_se_mgha": None, "gedi_n_muestras": 0}

    return {"agbd_mgha": agbd_mgha, "agbd_se_mgha": agbd_se_mgha, "gedi_n_muestras": int(n_muestras)}


def biomasa_a_co2e(area_ha, agb_mgha, agb_sd_mgha=None):
    """Convierte biomasa aérea (Mg/ha) a carbono (Mg C/ha) y a CO2
    equivalente total (t), con incertidumbre propagada. Único lugar del
    módulo donde se aplican FRACCION_CARBONO_BIOMASA y FACTOR_C_A_CO2 --
    para no repetir (y arriesgarse a aplicar mal) la conversión en más de
    un sitio."""
    if agb_mgha is None:
        return {"carbono_mgc_ha": None, "co2e_t": None, "co2e_incertidumbre_t": None}

    carbono_mgc_ha = agb_mgha * FRACCION_CARBONO_BIOMASA
    co2e_t = area_ha * carbono_mgc_ha * FACTOR_C_A_CO2

    carbono_sd_mgc_ha = (agb_sd_mgha * FRACCION_CARBONO_BIOMASA) if agb_sd_mgha is not None else None
    co2e_unc_t = (area_ha * carbono_sd_mgc_ha * FACTOR_C_A_CO2) if carbono_sd_mgc_ha is not None else None

    return {"carbono_mgc_ha": carbono_mgc_ha, "co2e_t": co2e_t, "co2e_incertidumbre_t": co2e_unc_t}


# ==============================================================================
# --- GEOMETRÍA EXACTA POR ZONA (UTM, igual que geomatica.py) ---
# ==============================================================================
def _reproyectores_utm(geom_wgs84_nucleo):
    """Arma los transformadores WGS84<->UTM para el sitio, usando la misma
    zona UTM que geomatica.py calcularía (por el centroide del núcleo)."""
    import pyproj
    centroid = geom_wgs84_nucleo.centroid
    utm_zone = int((centroid.x + 180) / 6) + 1
    utm_crs = f"EPSG:326{utm_zone}"
    a_utm = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True).transform
    a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True).transform
    return a_utm, a_wgs84


def geom_zona_precisa(geom_wgs84_nucleo, buf_m, a_utm, a_wgs84):
    """Buffer EXACTO en metros (vía UTM), igual que geomatica.py -- no una
    aproximación en grados. Devuelve (geometria_wgs84_para_GEE, area_ha_exacta).
    Con esto, el área de cada zona coincide con la que calcula geomatica.py."""
    from shapely.ops import transform as shp_transform

    geom_utm_nucleo = shp_transform(a_utm, geom_wgs84_nucleo)
    geom_utm_zona = geom_utm_nucleo.buffer(buf_m) if buf_m > 0 else geom_utm_nucleo
    area_ha_exacta = geom_utm_zona.area / 10000.0
    geom_wgs84_zona = shp_transform(a_wgs84, geom_utm_zona)
    return geom_wgs84_zona, area_ha_exacta


def _area_ha_wgs84_aprox(geom_wgs84):
    """Área aproximada en hectáreas, reproyectando a un UTM estimado por el
    centroide. Se usa solo en --demo (donde no hay geomatica.py de por medio)."""
    import pyproj
    from shapely.ops import transform as shp_transform

    centroid = geom_wgs84.centroid
    utm_zone = int((centroid.x + 180) / 6) + 1
    utm_crs = f"EPSG:326{utm_zone}"
    project = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True).transform
    geom_utm = shp_transform(project, geom_wgs84)
    return geom_utm.area / 10000.0


# ==============================================================================
# --- ORQUESTADOR CON DATOS REALES (requiere ee.Initialize() ya corrido) ---
# ==============================================================================
def calcular_carbono_por_zona_real(ee, geom_wgs84_nucleo, zonas_m, anio=None, incluir_gedi=True):
    """Calcula biomasa/carbono/CO2e por zona (núcleo + buffers), cada zona
    aislada -- mismo principio que geomatica.calcular_metricas_por_zona().
    El buffer se hace en UTM (metros exactos), igual que geomatica.py.

    Si incluir_gedi=True (default), además consulta GEDI L4A (mediciones
    LiDAR directas por huella) y agrega sus columnas al lado de las de
    ESA CCI en la MISMA fila -- para comparar "lo local" (GEDI, real pero
    disperso) contra "lo satelital" (ESA CCI, modelado pero continuo) sin
    tener que cruzar dos archivos distintos."""
    anio = anio or CARBONO_ANIO
    a_utm, a_wgs84 = _reproyectores_utm(geom_wgs84_nucleo)
    filas = []
    for buf_m in zonas_m:
        geom_zona, area_ha = geom_zona_precisa(geom_wgs84_nucleo, buf_m, a_utm, a_wgs84)

        log(f"Zona buffer={buf_m}m -- consultando biomasa en Earth Engine ({CARBONO_DATASET_GEE}/{anio})...")
        medida = medir_biomasa_zona_gee(ee, mapping(geom_zona), anio=anio)
        conv = biomasa_a_co2e(area_ha, medida["agb_mgha"], medida["agb_sd_mgha"])

        fila = {
            "zona": "nucleo" if buf_m == 0 else f"buffer_{buf_m}m",
            "buffer_m": buf_m,
            "area_ha": round(area_ha, 3),
            "agb_mgha": round(medida["agb_mgha"], 2) if medida["agb_mgha"] is not None else None,
            "agb_incertidumbre_mgha": round(medida["agb_sd_mgha"], 2) if medida["agb_sd_mgha"] is not None else None,
            "carbono_mgc_ha": round(conv["carbono_mgc_ha"], 2) if conv["carbono_mgc_ha"] is not None else None,
            "co2e_t": round(conv["co2e_t"], 1) if conv["co2e_t"] is not None else None,
            "co2e_incertidumbre_t": round(conv["co2e_incertidumbre_t"], 1) if conv["co2e_incertidumbre_t"] is not None else None,
            "dataset": CARBONO_DATASET_GEE,
            "anio_dataset": anio,
        }

        if incluir_gedi:
            log(f"Zona buffer={buf_m}m -- consultando biomasa GEDI L4A (huellas locales)...")
            medida_gedi = medir_biomasa_gedi_zona(ee, mapping(geom_zona))
            conv_gedi = biomasa_a_co2e(area_ha, medida_gedi["agbd_mgha"], medida_gedi["agbd_se_mgha"])
            fila.update({
                "gedi_agbd_mgha": round(medida_gedi["agbd_mgha"], 2) if medida_gedi["agbd_mgha"] is not None else None,
                "gedi_agbd_incertidumbre_mgha": round(medida_gedi["agbd_se_mgha"], 2) if medida_gedi["agbd_se_mgha"] is not None else None,
                "gedi_n_muestras": medida_gedi["gedi_n_muestras"],
                "gedi_carbono_mgc_ha": round(conv_gedi["carbono_mgc_ha"], 2) if conv_gedi["carbono_mgc_ha"] is not None else None,
                "gedi_co2e_t": round(conv_gedi["co2e_t"], 1) if conv_gedi["co2e_t"] is not None else None,
                "gedi_co2e_incertidumbre_t": round(conv_gedi["co2e_incertidumbre_t"], 1) if conv_gedi["co2e_incertidumbre_t"] is not None else None,
                "gedi_dataset": GEDI_DATASET_L4A,
            })
            log(f"  -> GEDI: {medida_gedi['gedi_n_muestras']} huellas, "
                f"AGBD={fila.get('gedi_agbd_mgha')}±{fila.get('gedi_agbd_incertidumbre_mgha')} Mg/ha")

        filas.append(fila)
        log(f"  -> ESA CCI: AGB={fila['agb_mgha']}±{fila['agb_incertidumbre_mgha']} Mg/ha  "
            f"carbono={fila['carbono_mgc_ha']} MgC/ha  CO2e={fila['co2e_t']}±{fila['co2e_incertidumbre_t']} t")

    return pd.DataFrame(filas)


def procesar_sitio_real(geojson_path, id_proyecto, zonas_m=None, carpeta_salida=None, proyecto_gee=None,
                         anio=None, incluir_gedi=True):
    """Pipeline completo con datos reales: inicializa Earth Engine, calcula
    biomasa/carbono/CO2e por zona (ESA CCI + GEDI si incluir_gedi=True),
    guarda CSV. Devuelve (df, csv_path)."""
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

    df = calcular_carbono_por_zona_real(ee, geom_wgs84, zonas_m, anio=anio, incluir_gedi=incluir_gedi)
    csv_path = os.path.join(carpeta_salida, f"carbono_por_zona_{id_proyecto.lower()}.csv")
    df.to_csv(csv_path, index=False)
    log(f"CSV de carbono por zona guardado en: {csv_path}")
    return df, csv_path


# ==============================================================================
# --- COMBINAR CON geomatica.py (solo lectura de CSVs, sin tocar su código) ---
# ==============================================================================
def _agregar_co2e_incremental(df, col_co2e, col_incertidumbre, col_prefijo):
    """Calcula el CO2e por anillo (sin traslape, sí sumable) para UNA fuente
    de datos (columna col_co2e), agregando col_prefijo+'_incremental_t' y
    col_prefijo+'_incremental_incertidumbre_t' al DataFrame. df debe venir
    ordenado por buffer_m ascendente (núcleo primero). Si col_co2e no existe
    en el DataFrame (ej. se corrió con --sin-gedi), no hace nada -- no
    truena, simplemente no agrega esas columnas."""
    if col_co2e not in df.columns:
        return df

    prev, unc_prev = 0.0, 0.0
    incrementales, incertidumbres_inc = [], []
    for _, fila in df.iterrows():
        valor = fila.get(col_co2e)
        unc = fila.get(col_incertidumbre)
        if pd.isna(valor):
            incrementales.append(None)
            incertidumbres_inc.append(None)
            continue
        inc = valor - prev
        unc_inc = (unc**2 + unc_prev**2) ** 0.5 if pd.notna(unc) else None
        incrementales.append(round(inc, 1))
        incertidumbres_inc.append(round(unc_inc, 1) if unc_inc is not None else None)
        prev, unc_prev = valor, (unc if pd.notna(unc) else 0.0)

    df[f"{col_prefijo}_incremental_t"] = incrementales
    df[f"{col_prefijo}_incremental_incertidumbre_t"] = incertidumbres_inc
    return df


def combinar_con_geomatica(carpeta_resultados, id_proyecto):
    """Lee el CSV de geomatica.py (features_por_zona_*) y el de este módulo
    (carbono_por_zona_*), y arma un CSV combinado por zona -- terreno +
    carbono juntos. No modifica ni importa código de geomatica.py, solo
    junta sus salidas por la columna 'zona'.

    IMPORTANTE -- ZONAS ACUMULATIVAS, NO SUMABLES ENTRE SÍ:
    Cada zona (nucleo, buffer_500m, buffer_1000m) es ACUMULADA: el
    buffer_500m ya incluye el área del núcleo, y el buffer_1000m ya
    incluye el área del buffer_500m. Sumar las columnas 'co2e_t' o
    'gedi_co2e_t' directamente triplicaría el carbono del núcleo. Por eso
    esta función agrega 'co2e_incremental_t' Y 'gedi_co2e_incremental_t'
    -- el carbono SOLO del anillo que agrega cada zona sobre la anterior,
    para AMBAS fuentes (ESA CCI y GEDI) por separado, que sí es correcto
    sumar dentro de cada fuente (no mezclar incrementales de una fuente
    con totales de la otra).

    El cálculo del incremental es una resta directa de masas totales, lo
    cual es válido incluso si la densidad de carbono no es uniforme dentro
    de cada zona -- porque co2e_t ya es masa total, no densidad. La
    incertidumbre del incremental es una aproximación conservadora (raíz
    de la suma de cuadrados, asumiendo independencia) -- no se corrigió
    por la correlación espacial real entre zonas superpuestas."""
    path_geo = os.path.join(carpeta_resultados, f"features_por_zona_{id_proyecto.lower()}.csv")
    path_carbono = os.path.join(carpeta_resultados, f"carbono_por_zona_{id_proyecto.lower()}.csv")

    if not os.path.exists(path_geo) or not os.path.exists(path_carbono):
        log(f"Falta uno de los CSVs (geomatica: {os.path.exists(path_geo)}, "
            f"carbono: {os.path.exists(path_carbono)}) -- corre ambos módulos antes de combinar.", nivel="WARN")
        return None

    df_geo = pd.read_csv(path_geo)
    df_carbono = pd.read_csv(path_carbono)
    df_combinado = df_geo.merge(df_carbono, on=["zona", "buffer_m"], suffixes=("", "_carbono"))
    df_combinado = df_combinado.sort_values("buffer_m").reset_index(drop=True)

    df_combinado = _agregar_co2e_incremental(df_combinado, "co2e_t", "co2e_incertidumbre_t", "co2e")
    df_combinado = _agregar_co2e_incremental(df_combinado, "gedi_co2e_t", "gedi_co2e_incertidumbre_t", "gedi_co2e")

    zonas_m_usadas = sorted(df_combinado["buffer_m"].dropna().unique().tolist())
    df_combinado = _agregar_fila_total_carbono(df_combinado, zonas_m_usadas)

    out_path = os.path.join(carpeta_resultados, f"resumen_terreno_y_carbono_{id_proyecto.lower()}.csv")
    df_combinado.to_csv(out_path, index=False)
    log(f"Resumen combinado (terreno + carbono, con CO2e incremental por anillo, ESA CCI + GEDI, "
        f"+ fila TOTAL) guardado en: {out_path}")
    return df_combinado


def _agregar_fila_total_carbono(df_combinado, zonas_m):
    """Agrega UNA fila 'TOTAL' al final de resumen_terreno_y_carbono con el
    gran total de CO2e (anillo exclusivo, sí sumable) para ESA CCI y GEDI
    -- para no tener que sumar a mano las filas por zona ni ir a buscarlo
    al subtítulo del mapa 3D (que ya lo calcula, pero solo como texto).

    El resto de columnas (terreno, AGB/ha por zona, etc.) se dejan vacías
    a propósito: no existe un 'total' que tenga sentido para ellas (sumar
    elev_promedio_m entre zonas, por ejemplo, no significa nada). Solo se
    llenan las columnas donde sumar SÍ es correcto -- los *_incremental_t
    (ya en anillo exclusivo, por diseño de _agregar_co2e_incremental) y
    gedi_n_muestras (conteo, sí sumable)."""
    fila_total = {col: None for col in df_combinado.columns}
    fila_total["zona"] = f"TOTAL (anillo exclusivo, 0-{max(zonas_m)}m, sí sumable)"
    fila_total["buffer_m"] = max(zonas_m)

    if "co2e_incremental_t" in df_combinado.columns:
        vals = df_combinado["co2e_incremental_t"].dropna()
        fila_total["co2e_incremental_t"] = round(vals.sum(), 1) if len(vals) else None
    if "co2e_incremental_incertidumbre_t" in df_combinado.columns:
        unc = df_combinado["co2e_incremental_incertidumbre_t"].dropna()
        fila_total["co2e_incremental_incertidumbre_t"] = round(float(np.sqrt((unc ** 2).sum())), 1) if len(unc) else None
    if "gedi_co2e_incremental_t" in df_combinado.columns:
        gedi_vals = df_combinado["gedi_co2e_incremental_t"].dropna()
        fila_total["gedi_co2e_incremental_t"] = round(gedi_vals.sum(), 1) if len(gedi_vals) else None
    if "gedi_co2e_incremental_incertidumbre_t" in df_combinado.columns:
        unc_g = df_combinado["gedi_co2e_incremental_incertidumbre_t"].dropna()
        fila_total["gedi_co2e_incremental_incertidumbre_t"] = round(float(np.sqrt((unc_g ** 2).sum())), 1) if len(unc_g) else None
    if "gedi_n_muestras" in df_combinado.columns:
        n = df_combinado["gedi_n_muestras"].dropna()
        fila_total["gedi_n_muestras"] = int(n.sum()) if len(n) else None
    for col_pasada in ("dataset", "anio_dataset", "gedi_dataset"):
        if col_pasada in df_combinado.columns and not df_combinado[col_pasada].dropna().empty:
            fila_total[col_pasada] = df_combinado[col_pasada].dropna().iloc[0]

    return pd.concat([df_combinado, pd.DataFrame([fila_total])], ignore_index=True)


def construir_capas_carbono_3d(hidrologia, df_zonas_reales, incluir_gedi=True):
    """Arma capas_extra para geomatica.generar_mapa_3d(): el ANILLO visual
    de cada zona (núcleo/buffer_500m/buffer_1000m, dibujado como borde de
    puntos sobre el propio terreno -- mismo estilo que 'Cauce en <zona>',
    ya usado en este mapa). El CO2e de cada zona se puede leer al pasar el
    mouse sobre su anillo (hover), pero YA NO se escribe como texto
    flotando en la escena 3D -- se probó así (ver versión anterior de esta
    función) y en la práctica el texto queda cortado por los ejes, tapado
    por el tooltip nativo de Plotly, o ilegible según el ángulo de cámara
    ("aquí se pierden los números", feedback real del usuario). Los
    números viven ahora en las tarjetas HTML alrededor del mapa (ver
    generar_mapa_3d_con_carbono()) -- ahí SIEMPRE son legibles sin
    importar cómo esté rotado el modelo. El texto en 3D world-space no es
    un lugar confiable para números que alguien necesita leer con
    certeza; el hover y las tarjetas sí lo son.

    Reto pedido: "que cada anillo... me dijera cuánto carbono almacena" --
    el anillo visual queda en la escena 3D (posición espacial, que sí
    importa aquí); el número en sí vive en la tarjeta HTML de al lado
    (legibilidad, que también importa, y en 3D no se puede garantizar).

    NO recalcula nada: el número del hover es el mismo
    co2e_incremental_t / gedi_co2e_incremental_t que ya trae
    resumen_terreno_y_carbono_*.csv (anillo exclusivo, el mismo que ya se
    verificó en balance_stock_vs_perdida y el mismo que usan las tarjetas
    HTML) -- una sola fuente de verdad, nunca dos cálculos por separado.

    La detección de bordes en sí (dónde dibujar el anillo) vive en
    geomatica.construir_anillos_visuales_3d() -- esta función solo arma
    el texto de hover con el CO2e de cada zona y se lo pasa; así
    core/carbono_perdida.py puede dibujar SUS anillos (con hover de CO2e
    LIBERADO en vez de almacenado) sin duplicar la detección de bordes.

    hidrologia: salida de geomatica.calcular_hidrologia_d8() -- usa
    'zona_de_pixel', que YA viene en anillo exclusivo por construcción
    (ver ese docstring).
    df_zonas_reales: resumen_terreno_y_carbono_*.csv ya cargado, SIN la
    fila TOTAL (pásale df_combinado filtrado, o el DataFrame que arma
    combinar_con_geomatica() antes de agregar esa fila)."""
    from core import geomatica

    hover_por_zona = {}
    for _, fila in df_zonas_reales.iterrows():
        zona = fila.get("zona")
        if not isinstance(zona, str):
            continue
        co2e = fila.get("co2e_incremental_t")
        if pd.isna(co2e):
            continue
        co2e_unc = fila.get("co2e_incremental_incertidumbre_t")
        gedi = fila.get("gedi_co2e_incremental_t") if incluir_gedi else None
        gedi_unc = fila.get("gedi_co2e_incremental_incertidumbre_t") if incluir_gedi else None

        hover_txt = f"{zona} -- CO2e almacenado (ESA CCI, anillo exclusivo): {co2e:,.0f}"
        hover_txt += f" ± {co2e_unc:,.0f} t" if pd.notna(co2e_unc) else " t"
        if gedi is not None and pd.notna(gedi):
            hover_txt += f"<br>{zona} -- CO2e almacenado (GEDI L4A, anillo exclusivo): {gedi:,.0f}"
            hover_txt += f" ± {gedi_unc:,.0f} t" if pd.notna(gedi_unc) else " t"
        hover_por_zona[zona] = hover_txt

    return geomatica.construir_anillos_visuales_3d(hidrologia, hover_por_zona=hover_por_zona)


def generar_mapa_3d_con_carbono(geojson_path, id_proyecto, zonas_m=None, percentil_cauce=None,
                                 carpeta_salida=None, carpeta_srtm=None):
    """Regenera el mapa 3D de geomatica.py, pero con el CO2e por zona
    agregado al título. Reutiliza geomatica.cargar_dem_utm() (que reusa el
    .tif de SRTM ya descargado -- no vuelve a bajar nada) y
    geomatica.calcular_hidrologia_d8()/generar_mapa_3d(), pasando el texto
    de carbono como 'subtitulo'. Requiere que ya exista el CSV combinado
    (corre combinar_con_geomatica() antes, o este mismo lo intenta)."""
    from core import geomatica

    zonas_m = zonas_m if zonas_m is not None else ZONAS_ANALISIS_M
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    carpeta_srtm = carpeta_srtm or CARPETA_SRTM

    df_combinado = combinar_con_geomatica(carpeta_salida, id_proyecto)
    if df_combinado is None:
        log("No se pudo generar el mapa con carbono: falta correr geomatica.py y/o carbono.py primero.", nivel="WARN")
        return None

    # df_combinado (desde esta sesión) ya trae, al final, una fila TOTAL
    # (ver _agregar_fila_total_carbono) -- para armar las tarjetas por
    # zona y el total aquí abajo hay que trabajar SOLO sobre las zonas
    # reales (nucleo/buffer_...), o esa fila TOTAL se sumaría una segunda
    # vez sobre sí misma (bug real: dejaba el total al doble del valor
    # correcto -- detectado y corregido en esta misma sesión).
    df_zonas_reales = df_combinado[~df_combinado["zona"].astype(str).str.startswith("TOTAL")]
    incluir_gedi = "gedi_co2e_incremental_t" in df_zonas_reales.columns
    anio_dataset = df_zonas_reales["anio_dataset"].iloc[0] if "anio_dataset" in df_zonas_reales.columns else "?"

    geom_utm_nucleo, dst_array, meta_utm, utm_crs = geomatica.cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        percentil_cauce or PERCENTIL_CAUCE_HIDROLOGIA,
        carpeta_srtm, id_proyecto,
    )
    capas_extra = construir_capas_carbono_3d(hidrologia, df_zonas_reales, incluir_gedi=incluir_gedi)

    # Título del propio Plotly: corto y sin repetir el nombre del
    # proyecto (eso ya lo dice el <h1> de la página) -- los números por
    # zona ya NO viven aquí (ver docstring de construir_capas_carbono_3d:
    # el texto en 3D world-space se pierde según el ángulo de cámara,
    # feedback real del usuario -- "aquí se pierden los números"). Las
    # cifras completas viven en las tarjetas HTML de arriba, siempre
    # legibles sin importar cómo esté rotado el modelo -- aquí solo queda
    # una pista de interacción.
    titulo_base = f"{id_proyecto} -- Modelo de terreno 3D + carbono (CO2e) por zona"
    titulo_interno = "Vista 3D interactiva"
    subtitulo_interno = "rota, acerca, y pasa el mouse sobre los anillos de color para ver su CO2e"
    fig = geomatica.generar_mapa_3d(hidrologia, id_proyecto, html_path=None, subtitulo=subtitulo_interno,
                                     utm_crs=utm_crs, titulo_base=titulo_interno, capas_extra=capas_extra,
                                     devolver_fig=True)

    html_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_3d_zonas_con_carbono.html")
    html_completo = _construir_html_carbono_con_tarjetas(fig, df_zonas_reales, id_proyecto, titulo_base,
                                                          anio_dataset, max(zonas_m), incluir_gedi)
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_completo)
    log(f"Mapa 3D con carbono (con tarjetas): {html_path}")
    return html_path


# ==============================================================================
# --- TARJETAS HTML: los números de CO2e, fuera de la escena 3D -----------------
# ==============================================================================
# El diseño de la tarjeta y de la página (CSS, layout) vive en
# core/reportes_html.py -- compartido con core/carbono_perdida.py (CO2e
# LIBERADO), para que un mismo ajuste de diseño no haya que repetirlo en
# dos archivos. Aquí solo queda la lógica de QUÉ dato va en cada tarjeta,
# específica de "CO2e almacenado".
def _construir_html_carbono_con_tarjetas(fig, df_zonas_reales, id_proyecto, titulo_base, anio_dataset,
                                          buffer_max_m, incluir_gedi):
    """Envuelve el go.Figure (ya armado por geomatica.generar_mapa_3d con
    devolver_fig=True) en una página propia: encabezado + una tarjeta por
    zona con su CO2e (grande, legible, sin importar cómo esté rotado el
    mapa) + el total + el mapa 3D debajo. Reemplaza el título/subtítulo
    largo de Plotly (donde los números "se perdían", feedback real del
    usuario) sin tocar geomatica.py ni ningún otro módulo que lo use --
    esto es exclusivo de este mapa."""
    from core import carbono_perdida, reportes_html

    # Hectáreas por anillo EXCLUSIVO (no acumuladas) -- mismo criterio que
    # el CO2e de cada tarjeta, para no mezclar una cifra acumulada con una
    # exclusiva en la misma tarjeta (esa mezcla sí sería muñeca rusa). Se
    # reusa el cálculo ya verificado de carbono_perdida.py en vez de
    # restar aquí por segunda vez -- una sola fuente de verdad para
    # "cuánta área tiene cada anillo".
    anillos_area, zona_de_buffer, _ = carbono_perdida._anillo_exclusivo_de_carbono(df_zonas_reales)
    areas_por_zona = {zona_de_buffer[buf_m]: datos["area_ha"] for buf_m, datos in anillos_area.items()}

    tarjetas = []
    total_co2e = total_unc_sq = total_gedi = total_gedi_unc_sq = total_huellas = total_area = 0.0
    hay_gedi_total = incluir_gedi

    for _, fila in df_zonas_reales.iterrows():
        zona = fila.get("zona")
        co2e = fila.get("co2e_incremental_t")
        if pd.isna(co2e):
            continue
        co2e_unc = fila.get("co2e_incremental_incertidumbre_t")
        gedi = fila.get("gedi_co2e_incremental_t") if incluir_gedi else None
        gedi_unc = fila.get("gedi_co2e_incremental_incertidumbre_t") if incluir_gedi else None
        n_huellas = fila.get("gedi_n_muestras") if incluir_gedi else None
        area_ha = areas_por_zona.get(zona)

        lineas_secundarias = []
        if gedi is not None and pd.notna(gedi):
            unc_txt = f" ± {gedi_unc:,.0f} t" if pd.notna(gedi_unc) else ""
            huellas_txt = f" · {int(n_huellas):,} huellas" if pd.notna(n_huellas) else ""
            lineas_secundarias.append(
                reportes_html.linea_secundaria_html(f"{gedi:,.0f} t CO2e{unc_txt}", f"GEDI L4A{huellas_txt}"))
        nota = f"± {co2e_unc:,.0f} t · ESA CCI Biomass" if pd.notna(co2e_unc) else "ESA CCI Biomass"

        tarjetas.append(reportes_html.tarjeta_html(
            zona, reportes_html.COLORES_ZONA_HEX.get(zona, "#666"), f"{co2e:,.0f}", "t CO2e",
            nota_principal=nota, lineas_secundarias=lineas_secundarias, area_ha=area_ha,
        ))
        total_co2e += co2e
        total_unc_sq += (co2e_unc ** 2) if pd.notna(co2e_unc) else 0.0
        if area_ha is not None:
            total_area += area_ha
        if gedi is not None and pd.notna(gedi):
            total_gedi += gedi
            total_gedi_unc_sq += (gedi_unc ** 2) if pd.notna(gedi_unc) else 0.0
            total_huellas += n_huellas if pd.notna(n_huellas) else 0.0
        else:
            hay_gedi_total = False

    lineas_total = []
    if hay_gedi_total:
        lineas_total.append(reportes_html.linea_secundaria_html(
            f"{total_gedi:,.0f} t CO2e ± {total_gedi_unc_sq ** 0.5:,.0f} t", f"GEDI L4A · {int(total_huellas):,} huellas"))
    tarjeta_total = reportes_html.tarjeta_html(
        "TOTAL", reportes_html.COLOR_TOTAL_HEX, f"{total_co2e:,.0f}", "t CO2e",
        nota_principal=f"± {total_unc_sq ** 0.5:,.0f} t · ESA CCI Biomass", lineas_secundarias=lineas_total,
        nombre_mostrado=f"Total (0-{buffer_max_m} m, anillos sumados)", area_ha=total_area, es_total=True,
    )

    div_mapa = fig.to_html(full_html=False, include_plotlyjs=True, config={"displaylogo": False})
    subtitulo = (f"CO2e almacenado por zona -- anillo exclusivo, sí sumable entre tarjetas (dataset {anio_dataset}, "
                 f"ESA CCI Above-Ground Biomass{' + GEDI L4A' if incluir_gedi else ''})")
    nota_pie = ("Los anillos de colores sobre el modelo marcan el límite de cada zona (rojo=núcleo, "
                "naranja=buffer 500m, dorado=buffer 1000m) -- pasa el mouse sobre el anillo para ver su CO2e "
                "también ahí. Incertidumbre de cada tarjeta: desviación estándar reportada por el dataset "
                "satelital, no un intervalo de confianza. GEDI y ESA CCI son dos estimaciones independientes -- "
                "no se promedian entre sí.")
    return reportes_html.pagina_html_con_tarjetas(
        f"{id_proyecto} -- Carbono por zona", titulo_base, subtitulo,
        "".join(tarjetas) + tarjeta_total, div_mapa, nota_pie,
    )


# ==============================================================================
# --- MODO DEMO: sin Earth Engine, valores sintéticos deterministas ---
# ==============================================================================
def demo():
    """Corre la lógica de cálculo (área + conversión biomasa->carbono->CO2e)
    sobre un círculo sintético y valores de biomasa inventados pero
    deterministas, sin tocar Earth Engine ni red. Incluye columnas GEDI
    sintéticas también, para probar que la comparación local/satelital
    no rompe nada aguas abajo (combinar_con_geomatica, mapas 3D)."""
    log("=== core.carbono --demo (sin Earth Engine, valores sintéticos) ===")
    centro = Point(-96.96, 19.40)  # coordenadas de referencia, no reales de ningún sitio
    zonas_m = [0, 300, 600]

    # Valores de BIOMASA (Mg/ha) sintéticos pero deterministas -- no vienen de GEE.
    # Rango orientativo para bosque mesófilo/secundario tropical (no es una medición).
    biomasa_sintetica = {0: (180.0, 25.0), 300: (155.0, 22.0), 600: (135.0, 19.0)}
    # GEDI sintético: deliberadamente distinto a ESA CCI (para simular la
    # discrepancia real que se espera entre "local" y "satelital"), y con
    # un caso de "sin huellas" en la zona más grande, a propósito -- para
    # probar que el código maneja bien la ausencia de datos, no solo el
    # caso feliz.
    gedi_sintetico = {0: (165.0, 30.0, 12), 300: (140.0, 35.0, 4), 600: (None, None, 0)}

    filas = []
    for buf_m in zonas_m:
        buf_grados = buf_m / 111000.0
        geom_zona = centro.buffer(0.02 + buf_grados)
        area_ha = _area_ha_wgs84_aprox(geom_zona)
        agb_mgha, agb_sd = biomasa_sintetica[buf_m]
        conv = biomasa_a_co2e(area_ha, agb_mgha, agb_sd)

        gedi_agbd, gedi_se, gedi_n = gedi_sintetico[buf_m]
        conv_gedi = biomasa_a_co2e(area_ha, gedi_agbd, gedi_se)

        fila = {
            "zona": "nucleo" if buf_m == 0 else f"buffer_{buf_m}m",
            "buffer_m": buf_m,
            "area_ha": round(area_ha, 3),
            "agb_mgha": agb_mgha,
            "agb_incertidumbre_mgha": agb_sd,
            "carbono_mgc_ha": round(conv["carbono_mgc_ha"], 2),
            "co2e_t": round(conv["co2e_t"], 1),
            "co2e_incertidumbre_t": round(conv["co2e_incertidumbre_t"], 1),
            "gedi_agbd_mgha": gedi_agbd,
            "gedi_agbd_incertidumbre_mgha": gedi_se,
            "gedi_n_muestras": gedi_n,
            "gedi_carbono_mgc_ha": round(conv_gedi["carbono_mgc_ha"], 2) if conv_gedi["carbono_mgc_ha"] is not None else None,
            "gedi_co2e_t": round(conv_gedi["co2e_t"], 1) if conv_gedi["co2e_t"] is not None else None,
            "gedi_co2e_incertidumbre_t": round(conv_gedi["co2e_incertidumbre_t"], 1) if conv_gedi["co2e_incertidumbre_t"] is not None else None,
        }
        filas.append(fila)

    df = pd.DataFrame(filas)
    print("\n--- Biomasa/Carbono/CO2e por zona (demo, valores sintéticos) ---")
    print(df.to_string(index=False))
    return df


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(description="Estimación de biomasa/carbono/CO2e por zona vía Earth Engine -- Motor Nacional")
    ap.add_argument("--demo", action="store_true", help="Corre con valores sintéticos, sin Earth Engine ni red")
    ap.add_argument("--geojson", type=str, help="Ruta al GeoJSON del polígono núcleo")
    ap.add_argument("--id-proyecto", type=str, help="Nombre identificador del sitio (para nombres de archivo)")
    ap.add_argument("--zonas", type=str, default=None, help="Buffers en metros separados por coma, ej. '0,500,1000'")
    ap.add_argument("--anio", type=int, default=None, help=f"Año de ESA CCI Biomass a usar (default: config.CARBONO_ANIO={CARBONO_ANIO})")
    ap.add_argument("--proyecto-gee", type=str, default=None, help="ID de proyecto de Google Cloud para ee.Initialize(project=...)")
    ap.add_argument("--sin-gedi", action="store_true", help="No consultar GEDI L4A (solo ESA CCI) -- por default SÍ se consultan ambos")
    ap.add_argument("--carpeta-salida", type=str, default=None)
    ap.add_argument("--combinar", action="store_true",
                     help="Solo combina CSVs ya existentes de geomatica.py + carbono.py en un resumen, sin consultar GEE")
    ap.add_argument("--mapa-3d", action="store_true",
                     help="Regenera el mapa 3D de geomatica.py con el CO2e por zona en el título (requiere haber corrido ambos módulos antes)")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not args.id_proyecto:
        ap.error("--id-proyecto es obligatorio")

    carpeta_salida = args.carpeta_salida or os.path.expanduser(f"~/resultados_{args.id_proyecto.lower()}")

    if args.mapa_3d:
        if not args.geojson:
            ap.error("--mapa-3d requiere --geojson (para poder regenerar la malla)")
        zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
        generar_mapa_3d_con_carbono(args.geojson, args.id_proyecto, zonas_m=zonas_m, carpeta_salida=carpeta_salida)
        return

    if args.combinar:
        combinar_con_geomatica(carpeta_salida, args.id_proyecto)
        return

    if not args.geojson:
        ap.error("--geojson es obligatorio fuera de --demo/--combinar/--mapa-3d")

    zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
    procesar_sitio_real(
        geojson_path=args.geojson, id_proyecto=args.id_proyecto,
        zonas_m=zonas_m, carpeta_salida=carpeta_salida, proyecto_gee=args.proyecto_gee, anio=args.anio,
        incluir_gedi=not args.sin_gedi,
    )


if __name__ == "__main__":
    main()

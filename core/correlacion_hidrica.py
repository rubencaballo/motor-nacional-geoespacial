#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prueba si el estrés hídrico (NDMI ya calculado por core/deforestacion.py,
o precipitación real vía CHIRPS) de un año dado PRECEDE la pérdida Hansen de
años posteriores -- en vez de solo mirar si coinciden en el mismo año.

POR QUÉ SEPARADO DE core/deforestacion.py (aunque reusa su NDMI):
    Nació de una pregunta real sobre Texolo (Ramsar 1601): "¿hay correlación
    entre la pérdida de 2025 y el estrés hídrico?" -- al probarla con el NDMI
    del MISMO año (lag=0), el resultado fue una correlación casi nula
    (r² entre 0.002 y 0.07 según la zona) Y el año con peor NDMI (2017) tuvo
    casi cero pérdida ese mismo año -- el dato contradice la hipótesis tal
    como se planteó inicialmente.

    Además hay un problema metodológico de fondo con lag=0: si una parte
    del bosque YA se taló ese año, esos píxeles quedan como suelo
    desnudo/vegetación rala, lo cual TAMBIÉN baja el NDMI promedio de la
    zona -- así que una correlación en lag=0 no distingue "había estrés
    hídrico antes y por eso lo cortaron/murió" de "lo cortaron y por eso
    bajó el NDMI ese mismo año". Por diseño, este módulo prueba lag>=1 (el
    estrés de un año, contra la pérdida de años SIGUIENTES) para evitar ese
    artefacto -- lag=0 se deja disponible solo como referencia/comparación,
    nunca como la lectura principal.

DOS FUENTES DE ESTRÉS HÍDRICO, no intercambiables:
    - NDMI (de core/deforestacion.py, YA calculado, sin costo GEE extra
      aquí): estrés hídrico de la VEGETACIÓN -- mezcla clima real con
      estructura/especie del dosel, y con el artefacto de arriba.
    - Precipitación CHIRPS (config.CHIRPS_DATASET): milímetros de lluvia
      reales, independientes de si hay o no bosque en el pixel -- más
      lento/caro de descargar (una reducción por año/zona), pero es la
      señal climática más directa, sin el sesgo de "el pixel ya no tiene
      vegetación que perder".

MUESTRA CHICA, SIEMPRE ADVERTIDO: con 8-10 años de historial (lo típico en
    este proyecto, limitado por Sentinel-2 desde 2016), cualquier
    coeficiente de correlación es estadísticamente débil. Este módulo NUNCA
    presenta un r como "comprobado" o "concluyente" sin reportar también n
    (pares de datos) al lado, y marca explícito 'muestra_chica'=True si
    n < config.CORRELACION_HIDRICA_MIN_ANIOS_PARA_R.

Qué NO hace:
    - No calcula un índice de sequía estandarizado (SPI/SPEI) -- solo
      milímetros totales anuales y su desviación simple contra el promedio
      del periodo. Para un análisis climático riguroso de verdad, un SPI/
      SPEI real sería más correcto; esto es un primer proxy, documentado
      como tal, no una verdad definitiva.
    - No prueba causalidad -- una correlación con rezago, aunque evita el
      artefacto de lag=0, sigue siendo correlación: no prueba que el
      estrés hídrico CAUSÓ la tala (podría haber otro factor de fondo,
      como accesibilidad, precio de la madera, o presión inmobiliaria, que
      también varíe año con año junto con el clima)."""

import argparse
import os
from datetime import datetime

import numpy as np
import pandas as pd

from config import (
    CHIRPS_DATASET, CORRELACION_HIDRICA_LAGS_DEFAULT, CORRELACION_HIDRICA_MIN_ANIOS_PARA_R,
    ZONAS_ANALISIS_M, DEFORESTACION_ANIO_INICIO_DEFAULT, log,
)


# ==============================================================================
# --- PARTE PURA: correlación con rezago, sin Earth Engine --------------------
# ==============================================================================
def analizar_correlacion_lag(anios, valores_clima, valores_perdida, lags=None):
    """Para cada lag en `lags` (años de rezago), correlaciona
    valores_clima[año] contra valores_perdida[año + lag] -- SOLO usando los
    pares de años donde ambos existen (ni None ni NaN). Pura, sin red.

    anios: lista de años (int), misma longitud/orden que valores_clima Y
    valores_perdida (un valor de clima y uno de pérdida por año -- el
    desfase lo aplica esta función, no se lo des ya corrido).

    Devuelve una lista de dicts, uno por lag: {lag, n_pares, r, r2,
    muestra_chica, anios_usados}. r/r2 son None si n_pares < 2 (no se
    puede calcular una correlación con menos de 2 puntos)."""
    lags = lags if lags is not None else CORRELACION_HIDRICA_LAGS_DEFAULT
    anios = list(anios)
    idx_por_anio = {a: i for i, a in enumerate(anios)}

    resultados = []
    for lag in lags:
        pares_clima, pares_perdida, anios_usados = [], [], []
        for anio in anios:
            anio_perdida = anio + lag
            if anio_perdida not in idx_por_anio:
                continue
            c = valores_clima[idx_por_anio[anio]]
            p = valores_perdida[idx_por_anio[anio_perdida]]
            if c is None or p is None:
                continue
            if isinstance(c, float) and np.isnan(c):
                continue
            if isinstance(p, float) and np.isnan(p):
                continue
            pares_clima.append(c)
            pares_perdida.append(p)
            anios_usados.append(anio)

        n = len(pares_clima)
        if n < 2:
            r = r2 = None
        else:
            r = float(np.corrcoef(pares_clima, pares_perdida)[0, 1])
            r2 = r ** 2
        resultados.append({
            "lag": lag, "n_pares": n, "r": round(r, 4) if r is not None else None,
            "r2": round(r2, 4) if r2 is not None else None,
            "muestra_chica": n < CORRELACION_HIDRICA_MIN_ANIOS_PARA_R,
            "anios_usados": anios_usados,
        })
    return resultados


def _perdida_exclusiva_por_zona(df_historial):
    """MUÑECA RUSA, arreglo compartido: `perdida_ha` en df_historial es
    ACUMULATIVA por zona (buffer_1000m incluye completo a buffer_500m, que
    incluye completo al núcleo -- misma convención que core/carbono.py y
    core/validacion_incendios.py). Correlacionar clima contra esa columna
    directa hace que buffer_500m y buffer_1000m salgan artificialmente
    parecidos entre sí (comparten los mismos píxeles de pérdida adentro) --
    dos "coincidencias" que en realidad son la misma señal contada dos
    veces, no dos confirmaciones independientes.

    Esta función resta cada buffer del siguiente más grande (idéntica
    lógica a core.deforestacion.generar_resumen_no_traslapado, válida
    porque cada buffer_Xm = nucleo.buffer(X), así que uno SIEMPRE contiene
    completo al anterior) y devuelve:
    - exclusiva: {buffer_m: {anio: perdida_ha del anillo exclusivo}}
    - zona_de_buffer: {buffer_m: etiqueta cruda usada en df_historial/CHIRPS, ej. 'buffer_500m'}
    - nombre_anillo: {buffer_m: etiqueta legible sin traslape, ej. 'anillo_0-500m'}"""
    buffers_ordenados = sorted(df_historial["buffer_m"].unique())
    zona_de_buffer = {buf_m: ("nucleo" if buf_m == 0 else f"buffer_{buf_m}m") for buf_m in buffers_ordenados}
    nombre_anillo = {}
    for i, buf_m in enumerate(buffers_ordenados):
        if i == 0:
            nombre_anillo[buf_m] = "nucleo (0m)" if buf_m == 0 else f"anillo_0-{buf_m}m"
        else:
            nombre_anillo[buf_m] = f"anillo_{buffers_ordenados[i - 1]}-{buf_m}m"

    pivot = df_historial.pivot_table(index="anio", columns="buffer_m", values="perdida_ha", aggfunc="first")
    exclusiva = {buf_m: {} for buf_m in buffers_ordenados}
    for anio in sorted(pivot.index):
        acumulado_prev = 0.0
        for buf_m in buffers_ordenados:
            valor = pivot.loc[anio, buf_m]
            valor = 0.0 if pd.isna(valor) else float(valor)
            exclusiva[buf_m][anio] = round(valor - acumulado_prev, 4)
            acumulado_prev = valor
    return exclusiva, zona_de_buffer, nombre_anillo


def analizar_correlacion_historial(df_historial, campo_clima="ndmi", lags=None):
    """Igual que analizar_correlacion_lag, pero recibe directo el
    df_historial de core/deforestacion.py (columnas zona/buffer_m/anio/
    perdida_ha/ndvi/ndmi/ndwi) y corre el análisis POR ZONA. No gasta
    cupo GEE -- reusa el NDVI/NDMI/NDWI que ya se calculó ahí.

    OJO -- MUÑECA RUSA: el lado de PÉRDIDA se convierte SIEMPRE a anillo
    exclusivo (sin traslape) vía _perdida_exclusiva_por_zona antes de
    correlacionar -- nunca se usa perdida_ha cruda/acumulativa aquí. El
    lado de CLIMA (NDMI/NDWI) se queda como promedio del disco acumulativo
    (0 a X m) porque un promedio no se puede "restar" entre zonas anidadas
    sin inventar un número (ver docstring del módulo). Esa asimetría
    (clima=disco, pérdida=anillo) queda documentada en las columnas
    geometria_clima/geometria_perdida del CSV de salida, para que nunca se
    lea como si fueran la misma geometría."""
    faltantes = [c for c in ["zona", "buffer_m", "anio", "perdida_ha"] if c not in df_historial.columns]
    if faltantes:
        log(f"analizar_correlacion_historial: faltan columnas {faltantes} en df_historial -- no se puede "
            "quitar la muñeca rusa (anillo exclusivo). Abortando.", nivel="ERROR")
        return None

    exclusiva, zona_de_buffer, nombre_anillo = _perdida_exclusiva_por_zona(df_historial)
    buffer_de_zona = {etiqueta: buf_m for buf_m, etiqueta in zona_de_buffer.items()}

    filas = []
    for zona, grupo in df_historial.groupby("zona"):
        if zona not in buffer_de_zona:
            log(f"  zona '{zona}' no coincide con ningún buffer_m conocido -- se omite.", nivel="WARN")
            continue
        buf_m = buffer_de_zona[zona]
        grupo = grupo.sort_values("anio")
        anios = grupo["anio"].tolist()
        clima = grupo[campo_clima].tolist()
        perdida = [exclusiva[buf_m].get(a, 0.0) for a in anios]
        for r in analizar_correlacion_lag(anios, clima, perdida, lags=lags):
            r["zona"] = nombre_anillo[buf_m]
            r["campo_clima"] = campo_clima
            r["geometria_clima"] = "nucleo" if buf_m == 0 else f"disco acumulativo 0-{buf_m}m"
            r["geometria_perdida"] = "anillo exclusivo (sin traslape)"
            filas.append(r)

    df_out = pd.DataFrame(filas)[
        ["zona", "campo_clima", "lag", "n_pares", "r", "r2", "muestra_chica",
         "geometria_clima", "geometria_perdida", "anios_usados"]
    ]
    for _, row in df_out.iterrows():
        etiqueta_lag = "mismo año (ojo: puede ser artefacto, ver docstring)" if row["lag"] == 0 else f"{row['lag']} año(s) antes"
        aviso = " -- MUESTRA CHICA, no concluyente" if row["muestra_chica"] else ""
        r_txt = f"r={row['r']}, r²={row['r2']}" if row["r"] is not None else "sin suficientes pares"
        log(f"  {row['zona']} | {row['campo_clima']} {etiqueta_lag} vs pérdida (anillo exclusivo) | n={row['n_pares']} | {r_txt}{aviso}")
    return df_out


# ==============================================================================
# --- CHIRPS: PRECIPITACIÓN ANUAL REAL PARA UN POLÍGONO -----------------------
# ==============================================================================
def calcular_precipitacion_anual_gee(ee, geom_wgs84_geojson, anio_inicio, anio_fin, dataset=None):
    """Precipitación TOTAL anual (mm, suma de CHIRPS diario) promediada
    sobre el polígono, año por año. Devuelve mm por año, el promedio del
    periodo, y la anomalía (mm y %) de cada año contra ese promedio -- NO
    es un SPI/SPEI real, es una desviación simple contra el promedio del
    propio periodo analizado (ver docstring del módulo). Si CHIRPS no
    tiene dato para un año (raro, pero posible en el borde de cobertura),
    se reporta None explícito, nunca se inventa."""
    dataset = dataset or CHIRPS_DATASET
    aoi = ee.Geometry(geom_wgs84_geojson)
    coleccion = ee.ImageCollection(dataset).select("precipitation")

    mm_por_anio = {}
    for anio in range(anio_inicio, anio_fin + 1):
        anual = coleccion.filterDate(f"{anio}-01-01", f"{anio}-12-31").sum()
        stats = anual.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=aoi, scale=5000, maxPixels=1e9, bestEffort=True,
        ).getInfo()
        mm = stats.get("precipitation")
        mm_por_anio[anio] = mm
        if mm is not None:
            log(f"  Año {anio}: {mm:.1f} mm totales (CHIRPS)")
        else:
            log(f"  Año {anio}: sin dato CHIRPS -- se reporta None, no se inventa.", nivel="WARN")

    valores_validos = [v for v in mm_por_anio.values() if v is not None]
    promedio_periodo = float(np.mean(valores_validos)) if valores_validos else None

    anomalia_mm, anomalia_pct = {}, {}
    for anio, mm in mm_por_anio.items():
        if mm is None or not promedio_periodo:
            anomalia_mm[anio] = anomalia_pct[anio] = None
        else:
            anomalia_mm[anio] = round(mm - promedio_periodo, 1)
            anomalia_pct[anio] = round((mm - promedio_periodo) / promedio_periodo * 100, 1)

    return {
        "mm_por_anio": mm_por_anio, "promedio_periodo_mm": promedio_periodo,
        "anomalia_mm": anomalia_mm, "anomalia_pct": anomalia_pct, "dataset": dataset,
    }


# ==============================================================================
# --- ORQUESTADOR CON DATOS REALES ---------------------------------------------
# ==============================================================================
def procesar_correlacion_hidrica_real(geojson_path, id_proyecto, zonas_m=None, anio_inicio=None, anio_fin=None,
                                       historial_csv_existente=None, lags=None, carpeta_salida=None,
                                       proyecto_gee=None, incluir_chirps=True, precipitacion_csv_existente=None):
    """Pipeline completo, en dos partes independientes:
    (1) SIEMPRE corre GRATIS la correlación contra el NDMI ya calculado, si
        le pasas `historial_csv_existente` (el CSV de core/deforestacion.py)
        -- cero cupo GEE, es solo aritmética sobre datos que ya tienes.
    (2) La precipitación real, de dos formas posibles:
        (a) `precipitacion_csv_existente`: reusa un precipitacion_anual_*.csv
            que ya descargaste en una corrida anterior -- CERO cupo GEE
            nuevo, solo recalcula la correlación (útil para volver a correr
            con el fix de anillo exclusivo sin pagar CHIRPS otra vez).
        (b) si no se da eso y `incluir_chirps=True`, descarga CHIRPS real
            por zona/año -- esto SÍ gasta cupo GEE (una reducción por
            año/zona, barato comparado con Sentinel-2).
    En ambos casos (1) y (2), el lado de PÉRDIDA se convierte a anillo
    exclusivo (sin traslape) antes de correlacionar -- ver docstring de
    _perdida_exclusiva_por_zona. Sin esto, buffer_500m y buffer_1000m salen
    artificialmente parecidos entre sí porque comparten los mismos píxeles
    de pérdida adentro (muñeca rusa)."""
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    os.makedirs(carpeta_salida, exist_ok=True)
    zonas_m = zonas_m if zonas_m is not None else ZONAS_ANALISIS_M

    df_historial = None
    resultados_ndmi = None
    exclusiva = zona_de_buffer = nombre_anillo = None
    if historial_csv_existente and os.path.exists(historial_csv_existente):
        log(f"Analizando correlación con NDMI ya calculado (sin costo GEE): {historial_csv_existente}")
        df_historial = pd.read_csv(historial_csv_existente)
        resultados_ndmi = analizar_correlacion_historial(df_historial, campo_clima="ndmi", lags=lags)
        if resultados_ndmi is not None:
            csv_ndmi = os.path.join(carpeta_salida, f"correlacion_ndmi_perdida_{id_proyecto.lower()}.csv")
            resultados_ndmi.to_csv(csv_ndmi, index=False)
            log(f"CSV correlación NDMI (con rezago, anillo exclusivo) guardado en: {csv_ndmi}", nivel="OK")
            exclusiva, zona_de_buffer, nombre_anillo = _perdida_exclusiva_por_zona(df_historial)
    else:
        log("No se dio --historial-csv-existente (o no se encontró) -- se omite la correlación NDMI gratis.",
            nivel="WARN")

    def _correlacionar_precip_por_zona(filas_precip_por_zona):
        """filas_precip_por_zona: {etiqueta ('nucleo'/'buffer_Xm'): [(anio, anomalia_pct), ...]}
        -- ya usa exclusiva/zona_de_buffer/nombre_anillo del scope exterior."""
        buffer_de_zona = {etiqueta: buf_m for buf_m, etiqueta in zona_de_buffer.items()}
        filas_corr = []
        for etiqueta, pares in filas_precip_por_zona.items():
            buf_m = buffer_de_zona.get(etiqueta)
            if buf_m is None:
                log(f"  zona '{etiqueta}' en el CSV de precipitación no coincide con ningún buffer_m del "
                    "historial -- se omite esa correlación (no se puede quitar la muñeca rusa sin saber "
                    "a qué anillo corresponde).", nivel="WARN")
                continue
            pares = sorted(pares, key=lambda t: t[0])
            anios_p = [a for a, _ in pares]
            anomalia_vals = [v for _, v in pares]
            perdida_vals = [exclusiva[buf_m].get(a, 0.0) for a in anios_p]
            for r in analizar_correlacion_lag(anios_p, anomalia_vals, perdida_vals, lags=lags):
                r["zona"] = nombre_anillo[buf_m]
                r["campo_clima"] = "precip_anomalia_pct"
                r["geometria_clima"] = "nucleo" if buf_m == 0 else f"disco acumulativo 0-{buf_m}m"
                r["geometria_perdida"] = "anillo exclusivo (sin traslape)"
                filas_corr.append(r)
        return filas_corr

    resultados_chirps = None

    if precipitacion_csv_existente and os.path.exists(precipitacion_csv_existente):
        if exclusiva is None:
            log("--precipitacion-csv-existente se dio pero no hay --historial-csv-existente válido -- "
                "no se puede calcular anillo exclusivo, se omite la correlación de precipitación.", nivel="ERROR")
        else:
            log(f"Recalculando correlación de precipitación con CHIRPS ya descargado (CERO cupo GEE nuevo, "
                f"con el fix de anillo exclusivo): {precipitacion_csv_existente}")
            df_precip_prev = pd.read_csv(precipitacion_csv_existente)
            filas_precip_por_zona = {}
            for etiqueta, grupo in df_precip_prev.groupby("zona"):
                filas_precip_por_zona[etiqueta] = list(zip(grupo["anio"], grupo["anomalia_pct"]))
            filas_corr = _correlacionar_precip_por_zona(filas_precip_por_zona)
            if filas_corr:
                resultados_chirps = pd.DataFrame(filas_corr)[
                    ["zona", "campo_clima", "lag", "n_pares", "r", "r2", "muestra_chica",
                     "geometria_clima", "geometria_perdida", "anios_usados"]
                ]
                csv_corr_chirps = os.path.join(
                    carpeta_salida, f"correlacion_precipitacion_perdida_{id_proyecto.lower()}.csv"
                )
                resultados_chirps.to_csv(csv_corr_chirps, index=False)
                log(f"CSV correlación precipitación-pérdida (recalculado, anillo exclusivo) guardado en: "
                    f"{csv_corr_chirps}", nivel="OK")

    elif incluir_chirps:
        import ee
        import geopandas as gpd
        from shapely.geometry import mapping as shp_mapping
        from core.deforestacion import _reproyectores_utm, geom_zona_precisa

        try:
            if proyecto_gee:
                ee.Initialize(project=proyecto_gee)
            else:
                ee.Initialize()
            log("Earth Engine inicializado.", nivel="OK")
        except Exception as e:
            raise RuntimeError(f"No se pudo inicializar Earth Engine: {e}")

        if not geojson_path or not os.path.exists(geojson_path):
            raise FileNotFoundError(
                f"--sin-chirps no está activo pero no se encontró el geojson ({geojson_path}) -- "
                "hace falta para descargar precipitación por zona."
            )

        gdf = gpd.read_file(geojson_path)
        geom_wgs84 = gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all") else gdf.geometry.unary_union
        a_utm, a_wgs84, _ = _reproyectores_utm(geom_wgs84)

        anio_inicio_f = anio_inicio or DEFORESTACION_ANIO_INICIO_DEFAULT
        anio_fin_f = anio_fin or (datetime.now().year - 1)

        filas_precip, filas_precip_por_zona = [], {}
        for buf_m in zonas_m:
            etiqueta = "nucleo" if buf_m == 0 else f"buffer_{buf_m}m"
            geom_zona_wgs84, _ = geom_zona_precisa(geom_wgs84, buf_m, a_utm, a_wgs84)
            log(f"Zona {etiqueta} -- descargando precipitación CHIRPS anual ({anio_inicio_f}-{anio_fin_f})...")
            precip = calcular_precipitacion_anual_gee(ee, shp_mapping(geom_zona_wgs84), anio_inicio_f, anio_fin_f)
            filas_precip_por_zona[etiqueta] = []
            for anio, mm in precip["mm_por_anio"].items():
                filas_precip.append({
                    "zona": etiqueta, "anio": anio, "precip_mm": mm,
                    "anomalia_mm": precip["anomalia_mm"][anio], "anomalia_pct": precip["anomalia_pct"][anio],
                })
                filas_precip_por_zona[etiqueta].append((anio, precip["anomalia_pct"][anio]))

        df_precip = pd.DataFrame(filas_precip)
        csv_precip = os.path.join(carpeta_salida, f"precipitacion_anual_{id_proyecto.lower()}.csv")
        df_precip.to_csv(csv_precip, index=False)
        log(f"CSV precipitación CHIRPS por zona/año guardado en: {csv_precip}", nivel="OK")

        if exclusiva is not None:
            filas_corr = _correlacionar_precip_por_zona(filas_precip_por_zona)
            if filas_corr:
                resultados_chirps = pd.DataFrame(filas_corr)[
                    ["zona", "campo_clima", "lag", "n_pares", "r", "r2", "muestra_chica",
                     "geometria_clima", "geometria_perdida", "anios_usados"]
                ]
                csv_corr_chirps = os.path.join(
                    carpeta_salida, f"correlacion_precipitacion_perdida_{id_proyecto.lower()}.csv"
                )
                resultados_chirps.to_csv(csv_corr_chirps, index=False)
                log(f"CSV correlación precipitación-pérdida guardado en: {csv_corr_chirps}", nivel="OK")
        else:
            log("No hay --historial-csv-existente válido -- se descargó precipitación pero se omite la "
                "correlación (no se puede quitar la muñeca rusa sin el historial de pérdida).", nivel="WARN")

    return resultados_ndmi, resultados_chirps


# ==============================================================================
# --- MODO DEMO: sin Earth Engine, sin red -- mismo espíritu que --demo en
#     el resto del proyecto ---
# ==============================================================================
def demo():
    """Prueba la parte pura (analizar_correlacion_lag) con una serie
    sintética construida A PROPÓSITO para que lag=1 tenga señal real
    (r alto) y lag=0 no -- para verificar que la función distingue
    correctamente entre rezagos, sin red ni Earth Engine."""
    log("=== core.correlacion_hidrica --demo (sin Earth Engine, valores sintéticos) ===")
    anios = list(range(2016, 2026))
    rng = np.random.default_rng(3)
    estres = rng.uniform(0, 1, size=len(anios))  # "estrés" sintético, 0=sin estrés, 1=máximo
    perdida = np.zeros(len(anios))
    for i in range(len(anios) - 1):
        perdida[i + 1] = max(0.0, estres[i] * 5 + rng.normal(0, 0.3))  # perdida del año SIGUIENTE depende del estres

    resultados = analizar_correlacion_lag(anios, list(estres), list(perdida), lags=(0, 1, 2))
    print("\n--- Correlación estrés-clima vs pérdida, por rezago (demo, sintético) ---")
    for r in resultados:
        print(f"  lag={r['lag']}: n={r['n_pares']}, r={r['r']}, r2={r['r2']}, muestra_chica={r['muestra_chica']}")
    print("\n  (Se espera lag=1 con r alto -- así se construyó la serie sintética a propósito -- "
          "y lag=0 mucho más débil, confirmando que la función distingue el rezago real)")
    return resultados


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Correlación (con rezago) entre estrés hídrico (NDMI o precipitación CHIRPS) y "
                    "pérdida Hansen -- Motor Nacional"
    )
    ap.add_argument("--demo", action="store_true", help="Corre con valores sintéticos, sin Earth Engine ni red")
    ap.add_argument("--geojson", type=str, help="Ruta al GeoJSON (obligatorio salvo con --sin-chirps)")
    ap.add_argument("--id-proyecto", type=str)
    ap.add_argument("--zonas", type=str, default=None, help="Buffers en metros separados por coma, ej. '0,500,1000'")
    ap.add_argument("--anio-inicio", type=int, default=None)
    ap.add_argument("--anio-fin", type=int, default=None)
    ap.add_argument("--historial-csv-existente", type=str, default=None,
                     help="CSV de core/deforestacion.py (deforestacion_historial_anual_*.csv) -- si se da, "
                          "corre GRATIS la correlación contra el NDMI ya calculado, sin tocar GEE")
    ap.add_argument("--lags", type=str, default=None, help="Años de rezago separados por coma, ej. '0,1,2'")
    ap.add_argument("--sin-chirps", action="store_true", help="No descargar precipitación CHIRPS (solo NDMI)")
    ap.add_argument("--precipitacion-csv-existente", type=str, default=None,
                     help="CSV precipitacion_anual_*.csv de una corrida previa -- si se da, recalcula la "
                          "correlación de precipitación CERO cupo GEE nuevo (ej. para aplicar el fix de "
                          "anillo exclusivo sin volver a pagar CHIRPS). Tiene prioridad sobre --sin-chirps.")
    ap.add_argument("--proyecto-gee", type=str, default=None)
    ap.add_argument("--carpeta-salida", type=str, default=None)
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not args.id_proyecto:
        ap.error("--id-proyecto es obligatorio fuera de --demo")

    zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
    lags = tuple(int(x) for x in args.lags.split(",")) if args.lags else None

    procesar_correlacion_hidrica_real(
        geojson_path=args.geojson, id_proyecto=args.id_proyecto, zonas_m=zonas_m,
        anio_inicio=args.anio_inicio, anio_fin=args.anio_fin,
        historial_csv_existente=args.historial_csv_existente, lags=lags,
        carpeta_salida=args.carpeta_salida, proyecto_gee=args.proyecto_gee,
        incluir_chirps=not args.sin_chirps,
        precipitacion_csv_existente=args.precipitacion_csv_existente,
    )


if __name__ == "__main__":
    main()

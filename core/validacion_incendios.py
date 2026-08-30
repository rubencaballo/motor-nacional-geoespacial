#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Valida los píxeles de pérdida Hansen de UN año específico contra la señal
real de quemado (dNBR -- Normalized Burn Ratio diferencial, Sentinel-2
pre/post evento), para distinguir un evento real (incendio) de un posible
artefacto del dataset Hansen cerca del treeline/zona alpina.

DE DÓNDE SALIÓ ESTE MÓDULO:
    Corriendo core/deforestacion.py sobre Cofre de Perote apareció una mancha
    grande y concentrada de pérdida Hansen 2025 a ~3899msnm, cerca de la
    cumbre. Antes de asumir nada, se buscó evidencia real: hubo un incendio
    forestal documentado en el Cofre de Perote en abril de 2025 (~200-400 ha,
    Plan DN-III-E activado). La correlación espacial/temporal es fuerte, pero
    Hansen NO distingue causa -- solo dice "aquí hubo pérdida de cobertura".
    Este módulo cierra esa brecha: compara, píxel a píxel, la pérdida Hansen
    de un año contra la firma espectral real de quemado (dNBR), en vez de
    quedarse en la inferencia por fecha/ubicación. Se construyó GENÉRICO --
    sirve para cualquier ANP y cualquier año, no solo para este caso.

POR QUÉ ESTE MÓDULO ESTÁ SEPARADO DE core/deforestacion.py:
    Mismo principio que en todo el proyecto: una fuente de datos por
    responsabilidad. deforestacion.py habla con Earth Engine para PÉRDIDA
    (Hansen) y tendencia de índices anuales (NDVI/NDMI/NDWI); este módulo
    habla con Earth Engine para un evento puntual (par de compuestos
    Sentinel-2 pre/post una fecha) y para la comparación estadística
    resultante. Si algo falla, se sabe de inmediato cuál pieza fue.

QUÉ ES dNBR -- documentado, no asumido:
    NBR = (B8 - B12) / (B8 + B12) -- banda NIR (B8) y SWIR2 (B12) de
    Sentinel-2. Vegetación sana tiene NBR alto; suelo/vegetación quemada
    tiene NBR bajo (la quema destruye la estructura celular que refleja NIR
    y aumenta la reflectancia SWIR). dNBR = NBR_pre - NBR_post: un valor
    POSITIVO y alto significa que el NBR bajó después del evento -- señal de
    quemado. Un valor cercano a 0 o negativo significa sin cambio o incluso
    reverdecimiento.

    Umbrales de severidad (UMBRALES_SEVERIDAD_DNBR en config.py): referencia
    de Key & Benson (USGS FIREMON), calibrados originalmente para Landsat.
    Se aplican aquí a Sentinel-2 como GUÍA APROXIMADA de severidad, no como
    umbral exacto validado para este sensor -- documentado, no una verdad
    absoluta.

    NOTA DE RESOLUCIÓN: B8 es nativo de 10m, B12 es nativo de 20m -- Earth
    Engine remuestrea internamente al combinarlas en normalizedDifference(),
    igual que ya hace este proyecto con B11 (20m) en el NDMI de
    deforestacion.py. Es una simplificación aceptada, no oculta.

VENTANA PRE/POST -- por qué días, no años:
    A diferencia de deforestacion.py (que compone Sentinel-2 por AÑO
    completo para una tendencia), aquí el compuesto "pre" y "post" son
    ventanas CORTAS (default 45 días, config.VALIDACION_INCENDIO_VENTANA_
    DIAS_PRE/POST_DEFAULT) alrededor de una fecha de evento puntual -- un
    compuesto anual completo diluiría/promediaría la señal de quemado con
    meses de vegetación sana antes y después.

QUÉ SÍ CALCULA:
    - Para los píxeles con lossyear==año_evento dentro de la zona: qué % cae
      por encima del umbral dNBR de "quemado" (confirmado), separado de los
      que no alcanzan ese umbral o no tienen dato dNBR limpio (sin confirmar
      -- nunca se asume "no fue incendio", se reporta explícito como
      "sin confirmar").
    - dNBR promedio en los píxeles Hansen del año vs. dNBR promedio en el
      resto de la zona visual (control) -- para contexto.
    - Desglose por clase de severidad (ver arriba).
    - Capa 3D opcional: verde=confirmado por dNBR, gris=sin confirmar,
      superpuesta al MISMO terreno/malla que ya usa geomatica.py.

QUÉ NO CALCULA:
    - Ninguna atribución de causa MÁS ALLÁ de "hubo/no hubo señal espectral
      de quemado" -- dNBR alto es consistente con incendio, pero también con
      otros disturbios severos (tala rasa reciente, por ejemplo). Para
      atribución de causa definitiva se requiere evidencia adicional (ej.
      reportes oficiales, como se hizo para el caso de Cofre de Perote).
    - Ninguna corrección automática de Hansen -- esto es solo una métrica de
      validación cruzada, igual que core/validacion_hidrologica.py con la
      Red Hidrográfica INEGI.
"""

import argparse
import os
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling
import requests
from shapely.geometry import mapping
from shapely.ops import transform as shp_transform

from config import (
    ZONAS_ANALISIS_M, PERCENTIL_CAUCE_HIDROLOGIA, CARPETA_SRTM, DEFORESTACION_ANIO_MIN_SENTINEL2,
    VALIDACION_INCENDIO_VENTANA_DIAS_PRE_DEFAULT, VALIDACION_INCENDIO_VENTANA_DIAS_POST_DEFAULT,
    VALIDACION_INCENDIO_NUBOSIDAD_MAX_PCT, VALIDACION_INCENDIO_UMBRAL_DNBR_QUEMADO,
    VALIDACION_INCENDIO_HA_MIN_PARA_EVALUAR, VALIDACION_INCENDIO_SCREENING_MESES_PRE,
    VALIDACION_INCENDIO_SCREENING_MESES_POST, UMBRALES_SEVERIDAD_DNBR, log,
)


# ==============================================================================
# --- FUNCIONES PURAS (sin red, sin Earth Engine -- testeables en --demo) ---
# ==============================================================================
def clasificar_severidad_dnbr(dnbr_array, umbrales=None):
    """Clasifica cada píxel de un array dNBR en una categoría de severidad
    (ver UMBRALES_SEVERIDAD_DNBR en config.py). Píxeles NaN se etiquetan
    'sin_dato', nunca se les asigna una categoría inventada."""
    umbrales = umbrales or UMBRALES_SEVERIDAD_DNBR
    etiquetas = np.full(dnbr_array.shape, "sin_dato", dtype=object)
    for nombre, (lo, hi) in umbrales.items():
        m = (dnbr_array >= lo) & (dnbr_array < hi) & ~np.isnan(dnbr_array)
        etiquetas[m] = nombre
    return etiquetas


def validar_perdida_contra_incendio(lossyear_alineado, dnbr_alineado, anio_evento, umbral_dnbr_quemado=None):
    """Compara, píxel a píxel, la pérdida Hansen de UN año (lossyear_alineado
    == anio_evento-2000) contra el dNBR real en esos mismos píxeles. Ambos
    arrays deben venir YA alineados a la misma malla (mismo shape/transform
    -- ver _descargar_lossyear_alineado() de core/deforestacion.py y
    _descargar_dnbr_alineado() de este módulo). Devuelve un dict con las
    estadísticas y una máscara booleana ('mascara_confirmados_dnbr') lista
    para colorear un mapa 3D punto a punto -- no solo el agregado."""
    umbral_dnbr_quemado = (umbral_dnbr_quemado if umbral_dnbr_quemado is not None
                            else VALIDACION_INCENDIO_UMBRAL_DNBR_QUEMADO)
    codigo_evento = anio_evento - 2000
    mascara_hansen = (lossyear_alineado == codigo_evento)
    n_pixeles_hansen = int(mascara_hansen.sum())

    if n_pixeles_hansen == 0:
        log(f"Sin píxeles Hansen con pérdida en {anio_evento} dentro de la zona visual -- nada que validar.",
            nivel="WARN")
        return {
            "anio_evento": anio_evento, "n_pixeles_hansen_ese_anio": 0, "n_pixeles_sin_dato_dnbr": 0,
            "n_pixeles_confirmados_dnbr": 0, "pct_confirmado_por_dnbr": None,
            "dnbr_promedio_pixeles_hansen": None, "dnbr_promedio_zona_control": None,
            "umbral_dnbr_usado": umbral_dnbr_quemado, "desglose_severidad_pixeles": {},
            "mascara_confirmados_dnbr": np.zeros_like(mascara_hansen, dtype=bool),
        }

    mascara_sin_dato = mascara_hansen & np.isnan(dnbr_alineado)
    mascara_hansen_con_dato = mascara_hansen & ~np.isnan(dnbr_alineado)
    n_sin_dato = int(mascara_sin_dato.sum())
    if n_sin_dato > 0:
        log(f"  {n_sin_dato} de {n_pixeles_hansen} píxeles Hansen de {anio_evento} SIN dato dNBR "
            f"(fuera del área cubierta por Sentinel-2 pre/post, o sin imagen limpia en la ventana) -- "
            f"se excluyen del % de confirmación, nunca se inventan.", nivel="WARN")

    dnbr_en_hansen = dnbr_alineado[mascara_hansen_con_dato]
    mascara_confirmados = mascara_hansen_con_dato & (dnbr_alineado >= umbral_dnbr_quemado)
    n_confirmados = int(mascara_confirmados.sum())
    n_con_dato = int(mascara_hansen_con_dato.sum())
    pct_confirmado = round(float(n_confirmados / n_con_dato * 100), 1) if n_con_dato > 0 else None

    mascara_control = ~mascara_hansen & ~np.isnan(dnbr_alineado)
    dnbr_control = dnbr_alineado[mascara_control]

    etiquetas = clasificar_severidad_dnbr(dnbr_alineado)
    desglose = {nombre: int(np.sum((etiquetas == nombre) & mascara_hansen_con_dato))
                for nombre in UMBRALES_SEVERIDAD_DNBR.keys()}

    return {
        "anio_evento": anio_evento,
        "n_pixeles_hansen_ese_anio": n_pixeles_hansen,
        "n_pixeles_sin_dato_dnbr": n_sin_dato,
        "n_pixeles_confirmados_dnbr": n_confirmados,
        "pct_confirmado_por_dnbr": pct_confirmado,
        "dnbr_promedio_pixeles_hansen": round(float(np.mean(dnbr_en_hansen)), 4) if len(dnbr_en_hansen) else None,
        "dnbr_promedio_zona_control": round(float(np.mean(dnbr_control)), 4) if len(dnbr_control) else None,
        "umbral_dnbr_usado": umbral_dnbr_quemado,
        "desglose_severidad_pixeles": desglose,
        "mascara_confirmados_dnbr": mascara_confirmados,
    }


# ==============================================================================
# --- EARTH ENGINE: dNBR PRE/POST ALINEADO A LA MALLA DEL TERRENO ---
# ==============================================================================
def _descargar_dnbr_alineado(ee, geom_wgs84_visual, fecha_evento, transform_ref, shape_ref, utm_crs,
                              ventana_dias_pre=None, ventana_dias_post=None, nubosidad_max_pct=None,
                              carpeta_tmp=None):
    """Construye compuestos Sentinel-2 (mediana) ANTES y DESPUÉS de
    `fecha_evento` (str 'YYYY-MM-DD'), calcula dNBR = NBR_pre - NBR_post, lo
    descarga vía getDownloadURL y lo REALINEA a la misma malla (transform_ref/
    shape_ref) que ya usa el terreno 3D de geomatica.py -- mismo patrón que
    _descargar_lossyear_alineado() de core/deforestacion.py, pero con
    Resampling.bilinear en vez de nearest: dNBR es una magnitud CONTINUA
    (como la elevación), no un código categórico como 'lossyear' -- ahí SÍ
    aplica bilinear.

    Devuelve un dict, nunca solo el array -- para reportar explícito cuántas
    imágenes Sentinel-2 entraron a cada compuesto (si alguno da 0, no se
    inventa un dNBR, se devuelve dnbr_alineado=None y el llamador decide)."""
    import tempfile

    ventana_dias_pre = ventana_dias_pre if ventana_dias_pre is not None else VALIDACION_INCENDIO_VENTANA_DIAS_PRE_DEFAULT
    ventana_dias_post = ventana_dias_post if ventana_dias_post is not None else VALIDACION_INCENDIO_VENTANA_DIAS_POST_DEFAULT
    nubosidad_max_pct = nubosidad_max_pct if nubosidad_max_pct is not None else VALIDACION_INCENDIO_NUBOSIDAD_MAX_PCT
    carpeta_tmp = carpeta_tmp or tempfile.gettempdir()
    os.makedirs(carpeta_tmp, exist_ok=True)

    fecha_dt = datetime.strptime(fecha_evento, "%Y-%m-%d")
    pre_ini = (fecha_dt - timedelta(days=ventana_dias_pre)).strftime("%Y-%m-%d")
    pre_fin = fecha_evento
    post_ini = fecha_evento
    post_fin = (fecha_dt + timedelta(days=ventana_dias_post)).strftime("%Y-%m-%d")

    aoi = ee.Geometry(geom_wgs84_visual)

    def _compuesto_nbr(ini, fin, etiqueta):
        col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
               .filterBounds(aoi).filterDate(ini, fin)
               .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", nubosidad_max_pct))
               .select(["B8", "B12"]))
        n = col.size().getInfo()
        log(f"  Compuesto {etiqueta} ({ini} a {fin}): {n} imágenes Sentinel-2 con <{nubosidad_max_pct}% nubes")
        if n == 0:
            return None, 0
        nbr = col.median().normalizedDifference(["B8", "B12"]).rename(f"nbr_{etiqueta}")
        return nbr, n

    nbr_pre, n_pre = _compuesto_nbr(pre_ini, pre_fin, "pre")
    nbr_post, n_post = _compuesto_nbr(post_ini, post_fin, "post")
    rangos = {"pre_rango": (pre_ini, pre_fin), "post_rango": (post_ini, post_fin)}

    if nbr_pre is None or nbr_post is None:
        log("Sin suficientes imágenes Sentinel-2 limpias en la ventana pre y/o post -- no se puede calcular "
            "dNBR. No se inventa un resultado; prueba con --ventana-dias-pre/--ventana-dias-post más grandes "
            "o --nubosidad-max más alto.", nivel="WARN")
        return {"dnbr_alineado": None, "n_imagenes_pre": n_pre, "n_imagenes_post": n_post, **rangos}

    dnbr = nbr_pre.subtract(nbr_post).rename("dnbr").clip(aoi)
    url = dnbr.getDownloadURL({"region": aoi, "scale": 30, "crs": "EPSG:4326", "format": "GEO_TIFF"})
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    tif_crudo = os.path.join(carpeta_tmp, "temp_dnbr_crudo.tif")
    with open(tif_crudo, "wb") as f:
        f.write(r.content)

    rows, cols = shape_ref
    # NaN, NO cero -- a diferencia de lossyear_alineado (donde 0 = "sin
    # pérdida" es un valor válido de relleno), aquí 0.0 SÍ es un valor real
    # de dNBR (sin cambio de NBR). Rellenar con 0 confundiría "sin dato" con
    # "confirmado que no se quemó" -- se excluye explícitamente con NaN.
    dnbr_alineado = np.full((rows, cols), np.nan, dtype=np.float32)
    with rasterio.open(tif_crudo) as src:
        reproject(
            source=rasterio.band(src, 1), destination=dnbr_alineado,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=transform_ref, dst_crs=utm_crs,
            resampling=Resampling.bilinear,  # dNBR es continuo, igual que el DEM -- aquí SÍ aplica bilinear
            src_nodata=src.nodata,
        )
    os.remove(tif_crudo)

    return {"dnbr_alineado": dnbr_alineado, "n_imagenes_pre": n_pre, "n_imagenes_post": n_post, **rangos}


def _descargar_dnbr_screening_alineado(ee, geom_wgs84_visual, anio, transform_ref, shape_ref, utm_crs,
                                        meses_pre=None, meses_post=None, nubosidad_max_pct=None,
                                        carpeta_tmp=None):
    """Igual que _descargar_dnbr_alineado(), pero SIN necesitar una fecha de
    evento externa (útil para años donde no tenemos una noticia/fecha
    confirmada de incendio -- la mayoría). En vez de una ventana pre/post
    alrededor de un evento puntual, compara dos ventanas FIJAS dentro del
    MISMO año (default: ene-abr vs sep-dic, ver config.
    VALIDACION_INCENDIO_SCREENING_MESES_PRE/POST) -- temporada seca/alta de
    incendios en el centro de México vs. fuera de temporada.

    OJO -- esto es un SCREENING, no una validación de evento: si el año SÍ
    tuvo un incendio real pero en, por ejemplo, noviembre, esta comparación
    (ene-abr vs sep-dic) puede no capturarlo bien, y al revés, un año con
    tala/plaga fuerte en esos mismos meses puede dar un dNBR alto sin haber
    sido fuego. Por eso el resultado de este método se etiqueta distinto
    ('screening_anual') del resultado de _descargar_dnbr_alineado()
    ('evento_confirmado') en procesar_historial_incendios_real() -- nunca se
    mezclan como si tuvieran la misma certeza."""
    meses_pre = meses_pre or VALIDACION_INCENDIO_SCREENING_MESES_PRE
    meses_post = meses_post or VALIDACION_INCENDIO_SCREENING_MESES_POST
    pre_ini = f"{anio}-{meses_pre[0]:02d}-01"
    pre_fin = f"{anio}-{meses_pre[1]:02d}-28"
    post_ini = f"{anio}-{meses_post[0]:02d}-01"
    post_fin = f"{anio}-{meses_post[1]:02d}-31"

    nubosidad_max_pct = nubosidad_max_pct if nubosidad_max_pct is not None else VALIDACION_INCENDIO_NUBOSIDAD_MAX_PCT
    import tempfile
    carpeta_tmp = carpeta_tmp or tempfile.gettempdir()
    os.makedirs(carpeta_tmp, exist_ok=True)

    aoi = ee.Geometry(geom_wgs84_visual)

    def _compuesto_nbr(ini, fin, etiqueta):
        col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
               .filterBounds(aoi).filterDate(ini, fin)
               .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", nubosidad_max_pct))
               .select(["B8", "B12"]))
        n = col.size().getInfo()
        log(f"    Compuesto {etiqueta} {anio} ({ini} a {fin}): {n} imágenes S2 con <{nubosidad_max_pct}% nubes")
        if n == 0:
            return None, 0
        return col.median().normalizedDifference(["B8", "B12"]).rename(f"nbr_{etiqueta}"), n

    nbr_pre, n_pre = _compuesto_nbr(pre_ini, pre_fin, "pre")
    nbr_post, n_post = _compuesto_nbr(post_ini, post_fin, "post")
    rangos = {"pre_rango": (pre_ini, pre_fin), "post_rango": (post_ini, post_fin)}
    if nbr_pre is None or nbr_post is None:
        log(f"  Año {anio}: sin suficientes imágenes S2 limpias para el screening -- se omite este año, "
            f"no se inventa.", nivel="WARN")
        return {"dnbr_alineado": None, "n_imagenes_pre": n_pre, "n_imagenes_post": n_post, **rangos}

    dnbr = nbr_pre.subtract(nbr_post).rename("dnbr").clip(aoi)
    url = dnbr.getDownloadURL({"region": aoi, "scale": 30, "crs": "EPSG:4326", "format": "GEO_TIFF"})
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    tif_crudo = os.path.join(carpeta_tmp, f"temp_dnbr_screening_{anio}.tif")
    with open(tif_crudo, "wb") as f:
        f.write(r.content)

    rows, cols = shape_ref
    dnbr_alineado = np.full((rows, cols), np.nan, dtype=np.float32)  # NaN, no cero -- mismo criterio
    # que _descargar_dnbr_alineado(): 0.0 dNBR es un valor real, no "sin dato".
    with rasterio.open(tif_crudo) as src:
        reproject(
            source=rasterio.band(src, 1), destination=dnbr_alineado,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=transform_ref, dst_crs=utm_crs,
            resampling=Resampling.bilinear, src_nodata=src.nodata,
        )
    os.remove(tif_crudo)

    return {"dnbr_alineado": dnbr_alineado, "n_imagenes_pre": n_pre, "n_imagenes_post": n_post, **rangos}


# ==============================================================================
# --- MAPA 3D: PÉRDIDA CONFIRMADA (verde) vs SIN CONFIRMAR (gris) POR dNBR ---
# ==============================================================================
def construir_capa_validacion_incendio_3d(lossyear_alineado, dnbr_alineado, hidrologia, resultado_validacion,
                                           utm_crs=None):
    """Construye los trazos go.Scatter3d (uno para confirmados, otro para
    sin confirmar) listos para pasarse como `capas_extra` a
    geomatica.generar_mapa_3d() -- mismo patrón que
    core/deforestacion.construir_capa_deforestacion_3d() y el split
    verde/rojo de core/validacion_hidrologica.py. Los píxeles 'sin dato
    dNBR' se agrupan dentro de 'sin confirmar' (honesto: no se pudieron
    confirmar, no se afirma que no se quemaron)."""
    import plotly.graph_objects as go
    from core.geomatica import calcular_grid_latlon

    anio_evento = resultado_validacion["anio_evento"]
    mascara_confirmados = resultado_validacion["mascara_confirmados_dnbr"]
    codigo_evento = anio_evento - 2000
    mascara_hansen = (lossyear_alineado == codigo_evento)
    mascara_no_confirmados = mascara_hansen & ~mascara_confirmados

    Z_raw = hidrologia["Z_raw"]
    Z_smooth = hidrologia.get("Z_smooth", Z_raw)
    pw_v, ph_v = hidrologia["pw_v"], hidrologia["ph_v"]
    transform = hidrologia["transform"]
    rows, cols = Z_raw.shape

    lat_grid = lon_grid = None
    if utm_crs:
        lat_grid, lon_grid = calcular_grid_latlon(transform, utm_crs, rows, cols)

    capas = []
    for etiqueta, color, mascara in [
        (f"Pérdida {anio_evento} confirmada por dNBR (incendio probable)", "lime", mascara_confirmados),
        (f"Pérdida {anio_evento} SIN confirmar por dNBR", "gray", mascara_no_confirmados),
    ]:
        fy, fx = np.where(mascara)
        if len(fx) > 0:
            validos = ~np.isnan(Z_raw[fy, fx]) & ~np.isnan(Z_smooth[fy, fx])
            fy, fx = fy[validos], fx[validos]
        if len(fx) == 0:
            continue

        x_km = fx * pw_v / 1000.0
        y_km = (rows - 1 - fy) * ph_v / 1000.0
        z_km = Z_smooth[fy, fx] + max(pw_v, ph_v) * 0.6

        dnbr_pts = dnbr_alineado[fy, fx]
        customdata_cols = [Z_raw[fy, fx], dnbr_pts]
        hovertemplate = f"{etiqueta}<br>Altitud: %{{customdata[0]:.0f}} msnm<br>dNBR: %{{customdata[1]:.3f}}"
        if lat_grid is not None:
            customdata_cols += [lat_grid[fy, fx], lon_grid[fy, fx]]
            hovertemplate += "<br>Lat: %{customdata[2]:.5f}<br>Lon: %{customdata[3]:.5f}"
        hovertemplate += "<extra></extra>"
        customdata = np.column_stack(customdata_cols)

        capas.append(go.Scatter3d(
            x=x_km, y=y_km, z=z_km, mode="markers",
            marker=dict(size=2.8, color=color, opacity=0.85),
            customdata=customdata, hovertemplate=hovertemplate, name=etiqueta,
        ))
        log(f"Capa '{etiqueta}': {len(fx)} píxeles.")

    return capas


# ==============================================================================
# --- ORQUESTADOR CON DATOS REALES ---
# ==============================================================================
def procesar_validacion_incendio_real(geojson_path, id_proyecto, fecha_evento, anio_hansen, zonas_m=None,
                                       ventana_dias_pre=None, ventana_dias_post=None, nubosidad_max_pct=None,
                                       umbral_dnbr_quemado=None, percentil_cauce=None, carpeta_salida=None,
                                       carpeta_srtm=None, proyecto_gee=None, mapa_3d=False):
    """Pipeline completo: reusa el terreno/D8 y el lossyear de Hansen (SRTM
    en caché, sin volver a descargar si ya existe -- misma lógica que
    core/deforestacion.py), descarga dNBR alineado a la misma malla, valida
    píxel a píxel, guarda el CSV y (opcional) genera el mapa 3D. Devuelve
    (resultado_dict, csv_path, html_path_o_None)."""
    from core import geomatica
    from core.deforestacion import _descargar_lossyear_alineado
    import pyproj
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

    log("Cargando terreno y cauces D8 (reusa SRTM en caché)...")
    geom_utm_nucleo, dst_array, meta_utm, utm_crs = geomatica.cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        percentil_cauce, carpeta_srtm, id_proyecto,
    )

    a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True).transform
    geom_visual_utm = geom_utm_nucleo.buffer(max(zonas_m))
    geom_visual_geojson = mapping(shp_transform(a_wgs84, geom_visual_utm))

    log("Descargando Hansen lossyear alineado a la malla del terreno...")
    lossyear_alineado = _descargar_lossyear_alineado(
        ee, geom_visual_geojson, hidrologia["transform"], hidrologia["Z_raw"].shape, utm_crs,
        carpeta_tmp=carpeta_srtm,
    )
    # OJO: cargar_dem_utm() recorta el SRTM al BOUNDING BOX (rectángulo) que envuelve la zona
    # visual, no al círculo/buffer real de zonas_m -- mismo caso que se detectó y corrigió en
    # core/deforestacion.generar_mapa_3d_deforestacion() (ahí el subtítulo llegó a decir 1316 ha
    # cuando el buffer real eran 771 ha). Sin este filtro, las esquinas del rectángulo (fuera del
    # buffer circular) se cuelan en "n_pixeles_hansen_ese_anio" -- inflando el total del subtítulo
    # y de paso metiendo puntos grises/verdes fuera del volcán en el mapa 3D. 255 = código inválido,
    # nunca coincide con ningún año real (mismo truco que usa procesar_historial_incendios_real).
    lossyear_alineado = np.where(hidrologia["zona_de_pixel"] != "fuera", lossyear_alineado, 255)

    log(f"Descargando dNBR (Sentinel-2 pre/post {fecha_evento}) alineado a la malla del terreno...")
    dnbr_resultado = _descargar_dnbr_alineado(
        ee, geom_visual_geojson, fecha_evento, hidrologia["transform"], hidrologia["Z_raw"].shape, utm_crs,
        ventana_dias_pre=ventana_dias_pre, ventana_dias_post=ventana_dias_post,
        nubosidad_max_pct=nubosidad_max_pct, carpeta_tmp=carpeta_srtm,
    )
    if dnbr_resultado["dnbr_alineado"] is None:
        raise RuntimeError(
            "No se pudo calcular dNBR -- sin suficientes imágenes Sentinel-2 limpias en la ventana pre/post "
            f"({dnbr_resultado['pre_rango']} / {dnbr_resultado['post_rango']}). Prueba con "
            "--ventana-dias-pre/--ventana-dias-post más grandes o --nubosidad-max más alto."
        )
    dnbr_alineado = dnbr_resultado["dnbr_alineado"]

    resultado = validar_perdida_contra_incendio(lossyear_alineado, dnbr_alineado, anio_hansen, umbral_dnbr_quemado)
    pixel_area_ha = hidrologia["pw_v"] * hidrologia["ph_v"] / 10000.0
    ha_hansen_total_log = resultado["n_pixeles_hansen_ese_anio"] * pixel_area_ha
    ha_confirmadas_log = resultado["n_pixeles_confirmados_dnbr"] * pixel_area_ha
    log(f"Validación {anio_hansen}: {resultado['n_pixeles_hansen_ese_anio']} píxeles Hansen "
        f"(~{ha_hansen_total_log:.2f} ha) -- {resultado['pct_confirmado_por_dnbr']}% confirmados por "
        f"dNBR>={resultado['umbral_dnbr_usado']} = ~{ha_confirmadas_log:.2f} ha CONFIRMADAS por incendio/"
        f"quema real (vs. ~{ha_hansen_total_log - ha_confirmadas_log:.2f} ha sin confirmar) "
        f"(dNBR prom. zona Hansen: {resultado['dnbr_promedio_pixeles_hansen']}, "
        f"dNBR prom. control/resto: {resultado['dnbr_promedio_zona_control']})")

    ha_hansen_total = round(resultado["n_pixeles_hansen_ese_anio"] * pixel_area_ha, 3)
    ha_confirmadas_dnbr = round(resultado["n_pixeles_confirmados_dnbr"] * pixel_area_ha, 3)
    ha_sin_confirmar_dnbr = round(ha_hansen_total - ha_confirmadas_dnbr, 3)

    fila_csv = {
        "id_proyecto": id_proyecto, "anio_evento": anio_hansen, "fecha_evento": fecha_evento,
        "n_pixeles_hansen_ese_anio": resultado["n_pixeles_hansen_ese_anio"],
        "ha_hansen_ese_anio": ha_hansen_total,
        "n_pixeles_sin_dato_dnbr": resultado["n_pixeles_sin_dato_dnbr"],
        "n_pixeles_confirmados_dnbr": resultado["n_pixeles_confirmados_dnbr"],
        "ha_confirmadas_dnbr": ha_confirmadas_dnbr,
        "ha_sin_confirmar_dnbr": ha_sin_confirmar_dnbr,
        "pct_confirmado_por_dnbr": resultado["pct_confirmado_por_dnbr"],
        "dnbr_promedio_pixeles_hansen": resultado["dnbr_promedio_pixeles_hansen"],
        "dnbr_promedio_zona_control": resultado["dnbr_promedio_zona_control"],
        "umbral_dnbr_usado": resultado["umbral_dnbr_usado"],
        "n_imagenes_s2_pre": dnbr_resultado["n_imagenes_pre"], "n_imagenes_s2_post": dnbr_resultado["n_imagenes_post"],
        "rango_pre": f"{dnbr_resultado['pre_rango'][0]} a {dnbr_resultado['pre_rango'][1]}",
        "rango_post": f"{dnbr_resultado['post_rango'][0]} a {dnbr_resultado['post_rango'][1]}",
    }
    for nombre, n in resultado["desglose_severidad_pixeles"].items():
        fila_csv[f"ha_severidad_{nombre}"] = round(n * pixel_area_ha, 3)

    df = pd.DataFrame([fila_csv])
    csv_path = os.path.join(carpeta_salida, f"validacion_incendio_{id_proyecto.lower()}_{anio_hansen}.csv")
    df.to_csv(csv_path, index=False)
    log(f"CSV de validación de incendio guardado en: {csv_path}")

    html_path = None
    if mapa_3d:
        capas_extra = construir_capa_validacion_incendio_3d(
            lossyear_alineado, dnbr_alineado, hidrologia, resultado, utm_crs=utm_crs,
        )
        subtitulo = (f"Validación incendio {fecha_evento} vs pérdida Hansen {anio_hansen}: "
                     f"{ha_confirmadas_dnbr:.2f} ha confirmadas por dNBR ({resultado['pct_confirmado_por_dnbr']}% "
                     f"de {ha_hansen_total:.2f} ha totales de pérdida Hansen {anio_hansen}) -- "
                     f"verde=confirmado, gris=sin confirmar ({ha_sin_confirmar_dnbr:.2f} ha)")
        html_path = os.path.join(
            carpeta_salida, f"{id_proyecto.lower()}_3d_validacion_incendio_{anio_hansen}.html",
        )
        geomatica.generar_mapa_3d(
            hidrologia, id_proyecto, html_path, subtitulo=subtitulo, utm_crs=utm_crs, capas_extra=capas_extra,
        )
        log(f"Mapa 3D de validación de incendio: {html_path}")

    return resultado, csv_path, html_path


def _causa_y_ha(r, metodo, pixel_area_ha):
    """Traduce el resultado de validar_perdida_contra_incendio() a una
    etiqueta de causa_probable + hectáreas confirmadas -- lógica compartida
    entre el cálculo ACUMULATIVO y el de ANILLO EXCLUSIVO en
    procesar_historial_incendios_real(), para no repetir el mismo if/elif
    dos veces con riesgo de que se desincronicen."""
    ha_confirmadas = (r["n_pixeles_confirmados_dnbr"] * pixel_area_ha) if r["n_pixeles_hansen_ese_anio"] else 0.0
    if r["n_pixeles_hansen_ese_anio"] == 0:
        causa = "perdida_dispersa_sin_evaluar"  # esta zona/anillo en particular no tenía pérdida ese año
    elif r["pct_confirmado_por_dnbr"] is None:
        causa = "sin_evaluar_sin_dato_dnbr"
    elif r["pct_confirmado_por_dnbr"] >= 50.0:
        causa = "confirmado_incendio" if metodo == "evento_confirmado" else "posible_incendio_screening"
    else:
        causa = "no_confirmado_por_dnbr"  # pérdida Hansen real, pero dNBR no da señal de quema
    return causa, ha_confirmadas


# ==============================================================================
# --- ORQUESTADOR MULTI-AÑO: causa probable por (zona, año), consolidado ---
# ==============================================================================
def procesar_historial_incendios_real(geojson_path, id_proyecto, zonas_m=None, anio_inicio=None,
                                       anio_fin=None, eventos_confirmados=None, ha_min_para_evaluar=None,
                                       nubosidad_max_pct=None, umbral_dnbr_quemado=None, percentil_cauce=None,
                                       historial_csv_existente=None, carpeta_salida=None, carpeta_srtm=None,
                                       proyecto_gee=None):
    """Corre la validación de causa (incendio confirmado / screening / sin
    evaluar) para VARIOS años a la vez, por zona, y arma UN solo CSV
    consolidado -- opcionalmente fusionado con el historial de
    core/deforestacion.py (perdida_ha, ndvi, ndmi, ndwi) si le pasas
    `historial_csv_existente`.

    DISEÑO PENSADO PARA NO GASTAR CUPO GEE DE MÁS (ver docstring del módulo):
    - El lossyear de Hansen se descarga UNA sola vez para todo el rango de
      años (igual que ya hace core/deforestacion.py) -- es gratis calcular
      la pérdida por (zona, año) a partir de ahí, sin ninguna llamada nueva.
    - El dNBR (lo caro: 2 compuestos Sentinel-2 + una descarga) se calcula
      UNA vez POR AÑO -- no una vez por (zona, año) -- y se reusa para las 3
      zonas enmascarando el mismo array. 3 zonas × 15 años sería 45
      descargas; así son máximo `anio_fin - anio_inicio + 1`.
    - Un año se SALTA por completo (cero llamadas a GEE) si la pérdida Hansen
      de ese año, en la zona más grande, es menor a `ha_min_para_evaluar`
      (default: config.VALIDACION_INCENDIO_HA_MIN_PARA_EVALUAR) -- unas
      pocas hectáreas dispersas en miles de hectáreas de zona no van a dar
      una señal de dNBR limpia, así que evaluarlas es gastar cupo sin
      ganar certeza real. Esos años quedan como 'perdida_dispersa_sin_evaluar',
      NUNCA se les inventa una causa.
    - Años sin cobertura Sentinel-2 (antes de config.DEFORESTACION_ANIO_MIN_
      SENTINEL2 + 1) se saltan automáticamente, mismo criterio.
    - ANTES de gastar ningún cupo, se imprime el plan completo (qué años se
      van a evaluar de verdad y cuáles se omiten y por qué) -- así puedes
      leerlo y cancelar (Ctrl+C) si el número de años a evaluar te parece
      alto, antes de que se haga ninguna descarga.

    `eventos_confirmados`: dict opcional {año: 'YYYY-MM-DD'} para años donde
    SÍ tienes una fecha real (ej. {2025: '2025-04-17'}) -- esos años usan el
    método de evento (ventana ajustada a la fecha real, más preciso) en vez
    del screening genérico ene-abr/sep-dic."""
    from core import geomatica
    from core.deforestacion import _descargar_lossyear_alineado
    import pyproj
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
    ha_min_para_evaluar = ha_min_para_evaluar if ha_min_para_evaluar is not None else VALIDACION_INCENDIO_HA_MIN_PARA_EVALUAR
    eventos_confirmados = eventos_confirmados or {}
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    carpeta_srtm = carpeta_srtm or CARPETA_SRTM
    os.makedirs(carpeta_salida, exist_ok=True)
    anio_fin = anio_fin or (datetime.now().year - 1)
    anio_inicio = anio_inicio or DEFORESTACION_ANIO_MIN_SENTINEL2

    log("Cargando terreno y cauces D8 (reusa SRTM en caché)...")
    geom_utm_nucleo, dst_array, meta_utm, utm_crs = geomatica.cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        percentil_cauce, carpeta_srtm, id_proyecto,
    )
    zona_de_pixel = hidrologia["zona_de_pixel"]  # OJO: son ANILLOS EXCLUSIVOS (nucleo / 0-500m / 500-1000m)
    # -- sirven para colorear cauces por anillo en geomatica.py, pero NO para esto: core/deforestacion.py
    # define cada "zona" como el buffer ACUMULATIVO completo (buffer_500m = TODO dentro de 500m, incluyendo
    # el núcleo; buffer_1000m = TODO dentro de 1000m, incluyendo núcleo + el anillo de 500m). Para que
    # perdida_ha (ya en el historial) y ha_confirmadas_dnbr (lo que calculamos aquí) hablen del MISMO
    # universo de píxeles al fusionarse, hay que rasterizar el buffer acumulativo de cada zona por separado
    # -- NO reusar zona_de_pixel para esto (ver bug real: buffer_1000m salía con ~4ha confirmadas en vez de
    # ~380ha porque solo contaba el anillo 500-1000m, no el disco completo de 0-1000m).
    from rasterio.features import rasterize
    mascaras_cumulativas = {}
    for buf_m in sorted(zonas_m):
        geom_zona = geom_utm_nucleo.buffer(buf_m) if buf_m > 0 else geom_utm_nucleo
        etiqueta = "nucleo" if buf_m == 0 else f"buffer_{buf_m}m"
        mascaras_cumulativas[etiqueta] = rasterize(
            [(geom_zona, 1)], out_shape=hidrologia["Z_raw"].shape, transform=hidrologia["transform"],
            fill=0, default_value=1, dtype=np.uint8,
        ).astype(bool)

    a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True).transform
    geom_visual_geojson = mapping(shp_transform(a_wgs84, geom_utm_nucleo.buffer(max(zonas_m))))

    log("Descargando Hansen lossyear alineado (una sola vez para todo el rango de años)...")
    lossyear_alineado = _descargar_lossyear_alineado(
        ee, geom_visual_geojson, hidrologia["transform"], hidrologia["Z_raw"].shape, utm_crs,
        carpeta_tmp=carpeta_srtm,
    )
    pixel_area_ha = hidrologia["pw_v"] * hidrologia["ph_v"] / 10000.0
    etiquetas_zona = ["nucleo"] + [f"buffer_{b}m" for b in sorted(zonas_m) if b > 0]

    # --- Plan (gratis, sin GEE): decidir qué años se evalúan de verdad ---
    anio_min_sentinel = DEFORESTACION_ANIO_MIN_SENTINEL2 + 1  # +1: el primer año con Sentinel-2 ya
    # completo suele venir con muy pocas imágenes (ver caso real: 2016 tuvo 0-1 imagen) -- no confiable
    plan_evaluar, plan_omitir = [], []
    for anio in range(anio_inicio, anio_fin + 1):
        if anio < anio_min_sentinel:
            plan_omitir.append((anio, f"sin cobertura Sentinel-2 confiable (< {anio_min_sentinel})"))
            continue
        codigo = anio - 2000
        ha_zona_mayor = float((lossyear_alineado[zona_de_pixel != "fuera"] == codigo).sum()) * pixel_area_ha
        if ha_zona_mayor < ha_min_para_evaluar:
            plan_omitir.append((anio, f"pérdida dispersa ({ha_zona_mayor:.2f} ha < {ha_min_para_evaluar} ha)"))
            continue
        metodo = "evento_confirmado" if anio in eventos_confirmados else "screening_anual"
        plan_evaluar.append((anio, metodo, ha_zona_mayor))

    log(f"=== PLAN (sin gastar cupo GEE todavía) ===")
    log(f"  Años a EVALUAR con dNBR ({len(plan_evaluar)}): " +
        ", ".join(f"{a}[{m}, ~{h:.1f}ha]" for a, m, h in plan_evaluar))
    log(f"  Años OMITIDOS ({len(plan_omitir)}): " + ", ".join(f"{a}({r})" for a, r in plan_omitir))
    log(f"  Esto son {len(plan_evaluar)} descargas de dNBR (2 compuestos S2 + 1 descarga cada una) -- "
        f"si te parece demasiado, cancela ahora (Ctrl+C) y ajusta --ha-min-evaluar o el rango de años.")

    # --- Ejecutar el plan (aquí sí se gasta cupo GEE, un año a la vez) ---
    filas = []
    for anio, metodo, _ in plan_evaluar:
        log(f"Año {anio} ({metodo})...")
        if metodo == "evento_confirmado":
            dnbr_resultado = _descargar_dnbr_alineado(
                ee, geom_visual_geojson, eventos_confirmados[anio], hidrologia["transform"],
                hidrologia["Z_raw"].shape, utm_crs, nubosidad_max_pct=nubosidad_max_pct, carpeta_tmp=carpeta_srtm,
            )
        else:
            dnbr_resultado = _descargar_dnbr_screening_alineado(
                ee, geom_visual_geojson, anio, hidrologia["transform"], hidrologia["Z_raw"].shape, utm_crs,
                nubosidad_max_pct=nubosidad_max_pct, carpeta_tmp=carpeta_srtm,
            )

        # perdida_ha_anillo_exclusivo: gratis (ya tenemos lossyear_alineado), y NO depende de si hubo
        # o no imágenes S2 limpias -- se calcula siempre, para las dos ramas de abajo.
        codigo_anio = anio - 2000
        perdida_anillo_ha = {
            etiqueta: float(((zona_de_pixel == etiqueta) & (lossyear_alineado == codigo_anio)).sum()) * pixel_area_ha
            for etiqueta in etiquetas_zona
        }

        if dnbr_resultado["dnbr_alineado"] is None:
            for etiqueta in etiquetas_zona:
                filas.append({
                    "zona": etiqueta, "anio": anio, "metodo_validacion": metodo,
                    "causa_probable": "sin_evaluar_sin_imagenes_s2_limpias",
                    "pct_confirmado_por_dnbr": None, "ha_confirmadas_dnbr": None,
                    "perdida_ha_anillo_exclusivo": round(perdida_anillo_ha[etiqueta], 3),
                    "causa_probable_anillo_exclusivo": "sin_evaluar_sin_imagenes_s2_limpias",
                    "pct_confirmado_por_dnbr_anillo_exclusivo": None, "ha_confirmadas_dnbr_anillo_exclusivo": None,
                    "n_imagenes_s2_pre": dnbr_resultado["n_imagenes_pre"],
                    "n_imagenes_s2_post": dnbr_resultado["n_imagenes_post"],
                })
            continue

        dnbr_alineado = dnbr_resultado["dnbr_alineado"]
        for etiqueta in etiquetas_zona:
            # (A) ACUMULATIVO -- "todo lo que está a menos de X metros del núcleo", el mismo universo
            #     que perdida_ha en el historial de deforestacion.py. Buffer_1000m SIEMPRE incluye lo
            #     del núcleo -- por diseño, no se puede sumar ha_confirmadas_dnbr de las 3 filas, cada
            #     una ya es un total independiente para su propio radio.
            lossyear_cumulativo = np.where(mascaras_cumulativas[etiqueta], lossyear_alineado, 255)
            r_cum = validar_perdida_contra_incendio(lossyear_cumulativo, dnbr_alineado, anio, umbral_dnbr_quemado)
            causa_cum, ha_cum = _causa_y_ha(r_cum, metodo, pixel_area_ha)

            # (B) ANILLO EXCLUSIVO -- "solo lo que está en ESTA franja, sin repetir lo de las franjas
            #     más chicas". Estas SÍ se pueden sumar: nucleo + buffer_500m_exclusivo +
            #     buffer_1000m_exclusivo = el total real quemado en toda la zona visual, una vez cada ha.
            lossyear_anillo = np.where(zona_de_pixel == etiqueta, lossyear_alineado, 255)
            r_ani = validar_perdida_contra_incendio(lossyear_anillo, dnbr_alineado, anio, umbral_dnbr_quemado)
            causa_ani, ha_ani = _causa_y_ha(r_ani, metodo, pixel_area_ha)

            filas.append({
                "zona": etiqueta, "anio": anio, "metodo_validacion": metodo,
                "causa_probable": causa_cum, "pct_confirmado_por_dnbr": r_cum["pct_confirmado_por_dnbr"],
                "ha_confirmadas_dnbr": round(ha_cum, 3),
                "perdida_ha_anillo_exclusivo": round(perdida_anillo_ha[etiqueta], 3),
                "causa_probable_anillo_exclusivo": causa_ani,
                "pct_confirmado_por_dnbr_anillo_exclusivo": r_ani["pct_confirmado_por_dnbr"],
                "ha_confirmadas_dnbr_anillo_exclusivo": round(ha_ani, 3),
                "n_imagenes_s2_pre": dnbr_resultado["n_imagenes_pre"],
                "n_imagenes_s2_post": dnbr_resultado["n_imagenes_post"],
            })

        total_confirmado_exclusivo = sum(
            f["ha_confirmadas_dnbr_anillo_exclusivo"] for f in filas
            if f["anio"] == anio and f["ha_confirmadas_dnbr_anillo_exclusivo"] is not None
        )
        log(f"  Año {anio}: total REAL confirmado por incendio en toda la zona visual (suma de los 3 "
            f"anillos exclusivos, sin repetir ninguna ha) = {total_confirmado_exclusivo:.2f} ha")

    for anio, _razon in plan_omitir:
        codigo_anio = anio - 2000
        for etiqueta in etiquetas_zona:
            perdida_anillo = float(((zona_de_pixel == etiqueta) & (lossyear_alineado == codigo_anio)).sum()) * pixel_area_ha
            filas.append({
                "zona": etiqueta, "anio": anio, "metodo_validacion": "no_evaluado", "causa_probable": _razon,
                "pct_confirmado_por_dnbr": None, "ha_confirmadas_dnbr": None,
                "perdida_ha_anillo_exclusivo": round(perdida_anillo, 3),
                "causa_probable_anillo_exclusivo": _razon,
                "pct_confirmado_por_dnbr_anillo_exclusivo": None, "ha_confirmadas_dnbr_anillo_exclusivo": None,
                "n_imagenes_s2_pre": None, "n_imagenes_s2_post": None,
            })

    df_causa = pd.DataFrame(filas).sort_values(["zona", "anio"]).reset_index(drop=True)

    if historial_csv_existente and os.path.exists(historial_csv_existente):
        log(f"Fusionando con el historial existente: {historial_csv_existente}")
        df_historial = pd.read_csv(historial_csv_existente)
        df_final = df_historial.merge(df_causa, on=["zona", "anio"], how="left")
        nombre_out = f"historial_completo_{id_proyecto.lower()}.csv"
    else:
        log("No se dio --historial-csv-existente (o no se encontró) -- se guarda solo la tabla de causa.")
        df_final = df_causa
        nombre_out = f"causa_incendios_{id_proyecto.lower()}.csv"

    csv_path = os.path.join(carpeta_salida, nombre_out)
    df_final.to_csv(csv_path, index=False)
    log(f"CSV consolidado guardado en: {csv_path}")
    generar_resumen_no_traslapado(df_final, id_proyecto, carpeta_salida)
    return df_final, csv_path


# Nombres de anillo NO acumulativos, para el resumen de lectura rápida (ver
# generar_resumen_no_traslapado abajo). Si algún día se agregan más buffers
# en ZONAS_ANALISIS_M, agrega aquí su traducción o el código cae al nombre
# original (zona_XXXm) en vez de reventar.
_NOMBRE_ANILLO_LEGIBLE = {
    "nucleo": "nucleo (0m)",
    "buffer_500m": "anillo_0-500m",
    "buffer_1000m": "anillo_500-1000m",
}


def generar_resumen_no_traslapado(df_final, id_proyecto, carpeta_salida):
    """Anexo de lectura rápida al CSV completo -- NO lo reemplaza.

    El CSV completo (historial_completo_*.csv / causa_incendios_*.csv) trae,
    a propósito, DOS familias de columnas para la misma fila de 'zona': las
    acumulativas (perdida_ha, ha_confirmadas_dnbr, causa_probable -- el
    universo completo desde el centro hasta ese radio, la convención de
    core/deforestacion.py usada en TODA la plataforma) y las de anillo
    exclusivo (*_anillo_exclusivo -- solo la franja, sin repetir lo de las
    zonas más chicas). Como ambas comparten el mismo texto de zona
    ('nucleo'/'buffer_500m'/'buffer_1000m'), a simple vista las acumulativas
    parecen 'muñeca rusa' (números parecidos entre sí porque cada zona
    grande contiene completa a la chica) -- son correctas, pero confunden
    si solo se lee esa columna sin el contexto.

    Esta función arma un CSV aparte con SOLO las columnas de anillo
    exclusivo, renombradas para que ninguna se pueda confundir con un
    buffer acumulativo (nucleo / anillo_0-500m / anillo_500-1000m), más una
    fila TOTAL por año que es la suma real de las tres franjas -- pensado
    para que alguien sin el contexto técnico de arriba lo pueda sumar
    directo y le dé el número correcto."""
    cols_necesarias = ["zona", "anio", "metodo_validacion", "perdida_ha_anillo_exclusivo",
                        "causa_probable_anillo_exclusivo", "pct_confirmado_por_dnbr_anillo_exclusivo",
                        "ha_confirmadas_dnbr_anillo_exclusivo"]
    faltantes = [c for c in cols_necesarias if c not in df_final.columns]
    if faltantes:
        log(f"generar_resumen_no_traslapado: faltan columnas {faltantes} -- no se genera el resumen.", nivel="WARN")
        return None

    df = df_final[cols_necesarias].copy().rename(columns={
        "zona": "anillo",
        "perdida_ha_anillo_exclusivo": "perdida_ha",
        "causa_probable_anillo_exclusivo": "causa_probable",
        "pct_confirmado_por_dnbr_anillo_exclusivo": "pct_confirmado_por_dnbr",
        "ha_confirmadas_dnbr_anillo_exclusivo": "ha_confirmadas_dnbr",
    })
    df["anillo"] = df["anillo"].map(lambda z: _NOMBRE_ANILLO_LEGIBLE.get(z, z))

    filas_total = []
    for anio, grupo in df.groupby("anio"):
        metodos = grupo["metodo_validacion"].dropna().unique()
        perdida_vals = grupo["perdida_ha"].dropna()
        if len(metodos) == 0 and len(perdida_vals) == 0:
            # Año que ni siquiera se tocó en esta corrida (ej. fuera del rango
            # --anio-inicio/--anio-fin usado) -- NO es lo mismo que "no_evaluado"
            # (que significa: sí se revisó y se decidió omitir por umbral/cobertura).
            # Aquí no sabemos nada del año, así que no se inventa un total en 0.
            metodo_anio, perdida_total = "sin_datos_en_esta_corrida", None
        else:
            metodo_anio = metodos[0] if len(metodos) == 1 else ("mixto" if len(metodos) > 1 else "no_evaluado")
            perdida_total = round(perdida_vals.sum(), 3) if len(perdida_vals) else None
        ha_conf = grupo["ha_confirmadas_dnbr"].dropna()
        filas_total.append({
            "anillo": "TOTAL 0-1000m (suma sin traslape)", "anio": anio, "metodo_validacion": metodo_anio,
            "perdida_ha": perdida_total, "causa_probable": "",
            "pct_confirmado_por_dnbr": None,
            "ha_confirmadas_dnbr": round(ha_conf.sum(), 3) if len(ha_conf) else None,
        })
    df = pd.concat([df, pd.DataFrame(filas_total)], ignore_index=True)
    df = df.sort_values(["anio", "anillo"]).reset_index(drop=True)

    nombre_out = f"resumen_incendios_sin_traslape_{id_proyecto.lower()}.csv"
    csv_path = os.path.join(carpeta_salida, nombre_out)
    df.to_csv(csv_path, index=False)
    log(f"Resumen SIN traslape (listo para sumar directo, sin muñeca rusa) guardado en: {csv_path}", nivel="OK")
    return csv_path


# ==============================================================================
# --- MODO DEMO: sin Earth Engine, sin red -- mismo espíritu que --demo en
#     el resto del proyecto ---
# ==============================================================================
def demo():
    """Prueba TODA la lógica que no depende de red/Earth Engine: la parte
    pura (clasificar_severidad_dnbr + validar_perdida_contra_incendio) con
    dos manchas sintéticas -- una con dNBR alto a propósito ('incendio
    real'), otra que se queda en el ruido de fondo ('posible artefacto') --
    y la capa 3D reusa geomatica._dem_sintetico(), igual que
    core/deforestacion.demo()."""
    log("=== core.validacion_incendios --demo (sin Earth Engine, valores sintéticos) ===")
    anio_evento = 2025
    codigo_evento = anio_evento - 2000
    rng = np.random.default_rng(7)

    # --- Parte pura, malla chica 50x50, sin geomatica ---
    dnbr_sint = rng.normal(0.05, 0.05, size=(50, 50)).astype(np.float32)  # ruido de fondo ~sin quemar
    lossyear_sint = np.zeros((50, 50), dtype=np.uint8)
    yy, xx = np.mgrid[0:50, 0:50]

    mancha_real = ((yy - 15) ** 2 + (xx - 15) ** 2) <= 8 ** 2
    lossyear_sint[mancha_real] = codigo_evento
    dnbr_sint[mancha_real] = rng.uniform(0.35, 0.7, size=int(mancha_real.sum()))  # severidad moderada-alta

    mancha_artefacto = ((yy - 35) ** 2 + (xx - 35) ** 2) <= 6 ** 2
    lossyear_sint[mancha_artefacto] = codigo_evento
    # dNBR se deja en el ruido de fondo -- simula "sin señal real de quemado"

    resultado = validar_perdida_contra_incendio(lossyear_sint, dnbr_sint, anio_evento)
    print("\n--- Validación incendio (demo, sintético) ---")
    for k, v in resultado.items():
        if k == "mascara_confirmados_dnbr":
            continue
        print(f"  {k}: {v}")
    print(f"\n  (Se espera ~{int(mancha_real.sum())} confirmados de "
          f"{int(mancha_real.sum() + mancha_artefacto.sum())} totales -- la mancha 1 tiene dNBR alto a "
          f"propósito, la mancha 2 se dejó en el ruido de fondo)")

    # --- Parte 3D: reusa geomatica._dem_sintetico(), igual que deforestacion.demo() ---
    try:
        from core import geomatica
        dst_array, meta_utm, geom_utm_nucleo, utm_crs = geomatica._dem_sintetico()
        zonas_m = [0, 300, 600]
        id_proyecto = "DEMO_VALIDACION_INCENDIO"
        carpeta_tmp = os.path.expanduser("~/resultados_demo_validacion_incendio")
        os.makedirs(carpeta_tmp, exist_ok=True)

        hidrologia = geomatica.calcular_hidrologia_d8(
            dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
            22, carpeta_tmp, id_proyecto,
        )
        rows, cols = hidrologia["Z_raw"].shape
        lossyear_grid = np.zeros((rows, cols), dtype=np.uint8)
        dnbr_grid = rng.normal(0.05, 0.05, size=(rows, cols)).astype(np.float32)
        cy1, cx1 = int(rows * 0.3), int(cols * 0.3)
        cy2, cx2 = int(rows * 0.7), int(cols * 0.7)
        yy2, xx2 = np.mgrid[0:rows, 0:cols]
        m1 = ((yy2 - cy1) ** 2 + (xx2 - cx1) ** 2) <= (min(rows, cols) * 0.1) ** 2
        m2 = ((yy2 - cy2) ** 2 + (xx2 - cx2) ** 2) <= (min(rows, cols) * 0.08) ** 2
        lossyear_grid[m1] = codigo_evento
        lossyear_grid[m2] = codigo_evento
        dnbr_grid[m1] = rng.uniform(0.35, 0.7, size=int(m1.sum()))

        resultado_grid = validar_perdida_contra_incendio(lossyear_grid, dnbr_grid, anio_evento)
        capas_extra = construir_capa_validacion_incendio_3d(
            lossyear_grid, dnbr_grid, hidrologia, resultado_grid, utm_crs=utm_crs,
        )

        html_path = os.path.join(carpeta_tmp, f"{id_proyecto.lower()}_3d_validacion_incendio.html")
        subtitulo = (f"Validación incendio (demo, sintético) {anio_evento}: "
                     f"{resultado_grid['pct_confirmado_por_dnbr']}% confirmados por dNBR "
                     f"(verde=confirmado, gris=sin confirmar -- NO son datos reales)")
        geomatica.generar_mapa_3d(
            hidrologia, id_proyecto, html_path, subtitulo=subtitulo, utm_crs=utm_crs, capas_extra=capas_extra,
        )
        log(f"Mapa 3D demo (validación de incendio) generado en: {html_path}")
    except ImportError as e:
        log(f"pysheds/plotly no instalado, se omite la parte 3D del demo: {e}", nivel="WARN")

    return resultado


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Valida píxeles de pérdida Hansen de UN año contra la señal real de quemado (dNBR, "
                    "Sentinel-2 pre/post evento) -- Motor Nacional"
    )
    ap.add_argument("--demo", action="store_true", help="Corre con valores sintéticos, sin Earth Engine ni red")
    ap.add_argument("--geojson", type=str, help="Ruta al GeoJSON del polígono núcleo")
    ap.add_argument("--id-proyecto", type=str, help="Nombre identificador (para nombres de archivo)")
    ap.add_argument("--fecha-evento", type=str, help="Fecha aproximada del evento (incendio), formato YYYY-MM-DD")
    ap.add_argument("--anio-hansen", type=int, help="Año de pérdida Hansen a validar (ej. 2025)")
    ap.add_argument("--zonas", type=str, default=None, help="Buffers en metros separados por coma, ej. '0,500,1000'")
    ap.add_argument("--ventana-dias-pre", type=int, default=None,
                     help=f"default: config.VALIDACION_INCENDIO_VENTANA_DIAS_PRE_DEFAULT="
                          f"{VALIDACION_INCENDIO_VENTANA_DIAS_PRE_DEFAULT}")
    ap.add_argument("--ventana-dias-post", type=int, default=None,
                     help=f"default: config.VALIDACION_INCENDIO_VENTANA_DIAS_POST_DEFAULT="
                          f"{VALIDACION_INCENDIO_VENTANA_DIAS_POST_DEFAULT}")
    ap.add_argument("--nubosidad-max", type=float, default=None,
                     help=f"default: config.VALIDACION_INCENDIO_NUBOSIDAD_MAX_PCT="
                          f"{VALIDACION_INCENDIO_NUBOSIDAD_MAX_PCT}")
    ap.add_argument("--umbral-dnbr", type=float, default=None,
                     help=f"default: config.VALIDACION_INCENDIO_UMBRAL_DNBR_QUEMADO="
                          f"{VALIDACION_INCENDIO_UMBRAL_DNBR_QUEMADO}")
    ap.add_argument("--percentil-cauce", type=float, default=None)
    ap.add_argument("--carpeta-srtm", type=str, default=None)
    ap.add_argument("--carpeta-salida", type=str, default=None)
    ap.add_argument("--proyecto-gee", type=str, default=None, help="ID de proyecto de Google Cloud para ee.Initialize(project=...)")
    ap.add_argument("--mapa-3d", action="store_true",
                     help="Genera también el mapa 3D (verde=confirmado por dNBR, gris=sin confirmar)")
    ap.add_argument("--historial", action="store_true",
                     help="Modo multi-año: evalúa causa probable (incendio/screening/sin evaluar) para un "
                          "rango de años a la vez, por zona, y guarda UN CSV consolidado (ver "
                          "procesar_historial_incendios_real). Usa --anio-inicio/--anio-fin en vez de "
                          "--fecha-evento/--anio-hansen.")
    ap.add_argument("--anio-inicio", type=int, default=None,
                     help="Solo con --historial: primer año a evaluar (default: config.DEFORESTACION_ANIO_MIN_SENTINEL2)")
    ap.add_argument("--anio-fin", type=int, default=None, help="Solo con --historial: default: año actual - 1")
    ap.add_argument("--eventos-confirmados", type=str, default=None,
                     help="Solo con --historial: años con fecha real conocida, formato 'anio:YYYY-MM-DD' "
                          "separados por coma, ej. '2025:2025-04-17'. Esos años usan el método de evento "
                          "(más preciso) en vez del screening genérico.")
    ap.add_argument("--ha-min-evaluar", type=float, default=None,
                     help=f"Solo con --historial: default: config.VALIDACION_INCENDIO_HA_MIN_PARA_EVALUAR="
                          f"{VALIDACION_INCENDIO_HA_MIN_PARA_EVALUAR}")
    ap.add_argument("--historial-csv-existente", type=str, default=None,
                     help="Solo con --historial: ruta al CSV de core.deforestacion (deforestacion_historial_anual_*.csv) "
                          "para fusionar perdida_ha/ndvi/ndmi/ndwi + causa en UN solo archivo")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if args.historial:
        if not args.geojson or not args.id_proyecto:
            ap.error("--geojson y --id-proyecto son obligatorios con --historial")
        zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
        eventos_confirmados = {}
        if args.eventos_confirmados:
            for par in args.eventos_confirmados.split(","):
                anio_str, fecha_str = par.split(":", 1)
                eventos_confirmados[int(anio_str.strip())] = fecha_str.strip()
        procesar_historial_incendios_real(
            geojson_path=args.geojson, id_proyecto=args.id_proyecto, zonas_m=zonas_m,
            anio_inicio=args.anio_inicio, anio_fin=args.anio_fin, eventos_confirmados=eventos_confirmados,
            ha_min_para_evaluar=args.ha_min_evaluar, nubosidad_max_pct=args.nubosidad_max,
            umbral_dnbr_quemado=args.umbral_dnbr, percentil_cauce=args.percentil_cauce,
            historial_csv_existente=args.historial_csv_existente, carpeta_salida=args.carpeta_salida,
            carpeta_srtm=args.carpeta_srtm, proyecto_gee=args.proyecto_gee,
        )
        return

    if not args.geojson or not args.id_proyecto or not args.fecha_evento or not args.anio_hansen:
        ap.error("--geojson, --id-proyecto, --fecha-evento y --anio-hansen son obligatorios fuera de --demo/--historial")

    zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
    procesar_validacion_incendio_real(
        geojson_path=args.geojson, id_proyecto=args.id_proyecto, fecha_evento=args.fecha_evento,
        anio_hansen=args.anio_hansen, zonas_m=zonas_m, ventana_dias_pre=args.ventana_dias_pre,
        ventana_dias_post=args.ventana_dias_post, nubosidad_max_pct=args.nubosidad_max,
        umbral_dnbr_quemado=args.umbral_dnbr, percentil_cauce=args.percentil_cauce,
        carpeta_salida=args.carpeta_salida, carpeta_srtm=args.carpeta_srtm, proyecto_gee=args.proyecto_gee,
        mapa_3d=args.mapa_3d,
    )


if __name__ == "__main__":
    main()

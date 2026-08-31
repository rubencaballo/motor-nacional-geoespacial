#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prueba si el incendio/pérdida de un año dado tiene un efecto hidrológico
DETECTABLE laderas abajo -- en celdas que NUNCA se quemaron ni se talaron,
pero que según el modelo D8 (core/geomatica.py) reciben el agua que escurre
desde la cicatriz.

DE DÓNDE SALIÓ ESTE MÓDULO:
    core/correlacion_hidrica.py (pérdida Hansen -> agua superficial JRC, por
    anillo concéntrico) no encontró ninguna correlación robusta en Cofre de
    Perote, ni siquiera para el incendio de 2025 (ver RUNBOOK.md, Paso 7c).
    Una explicación honesta -- no la única posible -- es que un anillo
    concéntrico promedia agua sobre TODA un área alrededor del núcleo, mezclando
    laderas que reciben el escurrimiento de una zona quemada con laderas que
    no tienen nada que ver: si el efecto real es direccional (agua abajo de
    la cicatriz, no en todas direcciones), un anillo lo diluye hasta volverlo
    invisible. Este módulo usa la dirección de flujo D8 (que
    core/geomatica.py ya calculaba internamente para llegar a la red de
    cauces, pero no exponía -- ver el campo 'flow_dir' que se le agregó) para
    construir el corredor EXACTO de celdas laderas abajo de la cicatriz, y
    compara SOLO esas celdas (nunca quemadas ellas mismas) contra un grupo
    control de terreno similar (mismo rango de elevación/pendiente) que no
    recibe ese flujo.

QUÉ SÍ CALCULA:
    - El "corredor descendente": celdas alcanzadas siguiendo la dirección de
      flujo D8 desde cada píxel de la cicatriz confirmada por dNBR (ver
      core/validacion_incendios.py), EXCLUYENDO la cicatriz misma -- el
      corredor es lo que está aguas abajo, no la quema en sí.
    - Un grupo control pareado por elevación y pendiente (mismos bins que el
      corredor), tomado de celdas que NO reciben ese flujo y que tampoco
      tienen pérdida Hansen propia en los últimos ~10 años (para no comparar
      contra terreno que en realidad también está disturbado por otra razón).
    - NDMI (Sentinel-2, humedad de la vegetación) ANTES del evento y DESPUÉS
      -- la ventana "después" es la MISMA ventana calendario, un año más
      tarde (no inmediatamente después del incendio), para no confundir el
      ciclo estacional normal (seco/lluvias) con un efecto real. Ver el
      porqué en `procesar_corredor_descendente`.
    - ΔNDMI = post - pre, por píxel, comparado corredor vs. control: un
      resumen global (media/mediana/n) y un PERFIL POR DISTANCIA (pasos D8
      desde la cicatriz) -- si el efecto es real, se esperaría más fuerte
      cerca de la cicatriz y desvanecerse con la distancia; si el corredor
      se comporta igual que el control a toda distancia, esa es la evidencia
      más fuerte posible de que este método no detecta un efecto.

QUÉ NO CALCULA / LIMITACIONES -- dichas de frente, no escondidas:
    - n=1: esto es UN caso de estudio (un incendio, un sitio), no una
      correlación repetida entre varios eventos como sí permite
      correlacion_hidrica.py (11-21 años de historial). Un resultado
      positivo O negativo aquí no se puede generalizar a "el fuego SIEMPRE/
      NUNCA afecta el agua laderas abajo" -- es evidencia de un solo evento.
    - PSEUDO-REPLICACIÓN: los píxeles vecinos NO son observaciones
      independientes (la humedad de un píxel se parece mucho a la de su
      vecino inmediato) -- cualquier prueba estadística píxel a píxel
      (Mann-Whitney, t de Welch) sobre miles de píxeles autocorrelacionados
      da un p-valor OPTIMISTA (parece más significativo de lo que en
      realidad es). Se reporta igual, con esta advertencia al lado siempre
      -- nunca como prueba concluyente por sí sola. El perfil por distancia
      es más informativo que el p-valor global: un efecto real debería
      desvanecerse lejos de la cicatriz; ruido no tiene por qué hacerlo.
    - Solo cubre el historial post-evento disponible al momento de correr
      esto -- típicamente un puñado de meses/un año, no varios ciclos
      completos.
    - Los píxeles DIRECTAMENTE quemados NO entran ni al corredor ni al
      control -- su NDMI cae por pérdida de vegetación misma (la quema), eso
      NO es evidencia de un efecto hidrológico, es la quema en sí. Se
      reportan aparte, solo de referencia.
    - El D8 usado viene de un DEM ESTÁTICO (SRTM, ~2000) -- no cambia con el
      incendio. Es una aproximación razonable de hacia dónde escurre el
      agua, no una medición real de escurrimiento post-incendio (que
      típicamente aumenta por hidrofobicidad del suelo quemado -- efecto
      conocido en la literatura, pero no medido directamente aquí).
    - Requiere que core/geomatica.calcular_hidrologia_d8() devuelva la clave
      'flow_dir' -- si se usa una copia vieja de ese módulo que no la trae,
      este módulo falla explícito con RuntimeError, nunca intenta adivinar.
"""

import argparse
import os
from collections import deque
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling
import requests
from scipy import stats as sstats

from config import (
    ZONAS_ANALISIS_M, PERCENTIL_CAUCE_HIDROLOGIA, CARPETA_SRTM, DEFORESTACION_NUBOSIDAD_MAX_PCT,
    VALIDACION_INCENDIO_UMBRAL_DNBR_QUEMADO, log,
)

# Códigos ESRI D8 -> (delta_fila, delta_columna), fila creciendo hacia abajo (convención numpy,
# fila 0 = arriba de la malla) -- EXACTO el mismo dirmap = (64, 128, 1, 2, 4, 8, 16, 32) que
# core/geomatica.calcular_hidrologia_d8() le pasa a pysheds internamente:
#     32  64 128
#     16   X   1
#      8   4   2
# Si algún día ese dirmap cambia en geomatica.py, este diccionario debe cambiar junto -- por eso
# queda documentado aquí en vez de "obvio"/implícito.
_CODIGO_A_DELTA_D8 = {
    1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1),
    16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1),
}


# ==============================================================================
# --- FUNCIONES PURAS (sin red, sin Earth Engine -- testeables en --demo) ------
# ==============================================================================
def trazar_corredor_descendente(flow_dir, mascara_semilla, mascara_valida, max_pasos=300):
    """BFS laderas-abajo: desde cada píxel True en `mascara_semilla` (ej. la
    cicatriz quemada), sigue `flow_dir` (códigos ESRI, ver _CODIGO_A_DELTA_D8)
    paso a paso hasta salir de la malla, salir de `mascara_valida`, llegar a
    un código sin dirección válida (sumidero/nodata), o alcanzar una celda ya
    visitada (por esta u otra semilla -- evita reprocesar y protege contra
    cualquier ciclo, aunque un D8 bien resuelto -- pit-filled + resolve_flats,
    como hace geomatica.py -- no debería tener ciclos). `max_pasos` es un
    tope de seguridad adicional, no la forma normal de terminar un trazo.

    Devuelve (mascara_corredor, distancia): `mascara_corredor` es
    `mascara_semilla`-EXCLUSIVA (el corredor es lo que está aguas abajo, no
    la semilla misma) y `distancia` es un array int32 con el número de pasos
    D8 desde la semilla más cercana a cada celda alcanzada (-1 donde no se
    alcanzó nada)."""
    rows, cols = flow_dir.shape
    visitado = np.zeros(flow_dir.shape, dtype=bool)
    distancia = np.full(flow_dir.shape, -1, dtype=np.int32)

    cola = deque()
    ys, xs = np.where(mascara_semilla)
    for y, x in zip(ys, xs):
        visitado[y, x] = True
        distancia[y, x] = 0
        cola.append((int(y), int(x), 0))

    pasos_truncados = 0
    while cola:
        y, x, d = cola.popleft()
        if d >= max_pasos:
            pasos_truncados += 1
            continue
        delta = _CODIGO_A_DELTA_D8.get(int(flow_dir[y, x]))
        if delta is None:
            continue  # sin dirección válida (sumidero/borde/nodata) -- el flujo termina aquí
        ny, nx = y + delta[0], x + delta[1]
        if not (0 <= ny < rows and 0 <= nx < cols):
            continue  # sale de la malla
        if not mascara_valida[ny, nx]:
            continue  # sale del área con dato de terreno válido
        if visitado[ny, nx]:
            continue  # ya alcanzada -- no se reprocesa
        visitado[ny, nx] = True
        distancia[ny, nx] = d + 1
        cola.append((ny, nx, d + 1))

    if pasos_truncados > 0:
        log(f"  Trazo de corredor: {pasos_truncados} rutas truncadas a los {max_pasos} pasos máximos "
            f"(protección contra rutas absurdamente largas, no debería activarse en terreno normal).",
            nivel="WARN")

    mascara_corredor = visitado & ~mascara_semilla
    return mascara_corredor, distancia


def seleccionar_grupo_control(mascara_corredor, mascara_excluir, elevacion, pendiente,
                               n_bins_elev=5, n_bins_pendiente=3, factor_control=1, semilla=42):
    """Grupo control pareado: para cada combinación (bin de elevación, bin de
    pendiente) presente en el corredor, toma al azar (sin reemplazo) hasta
    `factor_control` veces esa cantidad de píxeles de terreno con la MISMA
    combinación de bins, EXCLUYENDO `mascara_excluir` (típicamente cicatriz +
    corredor + cualquier pérdida Hansen reciente propia -- terreno que ya
    está disturbado por otra razón no es un buen "control sin disturbio").

    Los bins se definen por los CUANTILES de elevación/pendiente DEL
    CORREDOR (no de toda la zona) -- así el control queda pareado al rango
    real que ocupa el corredor, no a la distribución de todo el terreno
    (que puede incluir elevaciones/pendientes que el corredor ni toca).

    Devuelve (mascara_control, resumen_bins) -- `resumen_bins` es una lista
    de dicts con cuántos píxeles se pidieron/encontraron/tomaron por bin,
    para poder reportar explícito cuando un bin se queda con control
    incompleto (nunca se inventan píxeles para completar la cuota)."""
    rng = np.random.default_rng(semilla)
    valido = ~np.isnan(elevacion) & ~np.isnan(pendiente)
    candidatos = valido & ~mascara_corredor & ~mascara_excluir

    if not np.any(mascara_corredor) or not np.any(candidatos):
        return np.zeros_like(mascara_corredor, dtype=bool), []

    elev_corredor = elevacion[mascara_corredor]
    pend_corredor = pendiente[mascara_corredor]
    bordes_elev = np.unique(np.quantile(elev_corredor, np.linspace(0, 1, n_bins_elev + 1)))
    bordes_pend = np.unique(np.quantile(pend_corredor, np.linspace(0, 1, n_bins_pendiente + 1)))

    bin_elev_corredor = np.digitize(elev_corredor, bordes_elev[1:-1])
    bin_pend_corredor = np.digitize(pend_corredor, bordes_pend[1:-1])
    bin_elev_todo = np.digitize(elevacion, bordes_elev[1:-1])
    bin_pend_todo = np.digitize(pendiente, bordes_pend[1:-1])

    conteo_necesario = {}
    for be, bp in zip(bin_elev_corredor, bin_pend_corredor):
        conteo_necesario[(be, bp)] = conteo_necesario.get((be, bp), 0) + 1

    mascara_control = np.zeros_like(mascara_corredor, dtype=bool)
    resumen_bins = []
    for (be, bp), n_necesario in sorted(conteo_necesario.items()):
        n_pedido = n_necesario * factor_control
        candidatos_bin = candidatos & (bin_elev_todo == be) & (bin_pend_todo == bp)
        ys, xs = np.where(candidatos_bin)
        n_disponible = len(ys)
        n_tomar = min(n_pedido, n_disponible)
        resumen_bins.append({
            "bin_elevacion": int(be), "bin_pendiente": int(bp),
            "n_corredor": int(n_necesario), "n_pedido_control": int(n_pedido),
            "n_disponible_control": int(n_disponible), "n_tomado_control": int(n_tomar),
        })
        if n_tomar == 0:
            continue
        idx = rng.choice(n_disponible, size=n_tomar, replace=False)
        mascara_control[ys[idx], xs[idx]] = True

    return mascara_control, resumen_bins


def comparar_ndmi_corredor_control(ndmi_pre, ndmi_post, mascara_corredor, mascara_control, mascara_quemado,
                                    distancia_corredor=None, n_bins_distancia=4):
    """ΔNDMI = post - pre, por píxel, resumido por grupo (corredor/control/
    cicatriz -- la cicatriz solo de referencia, ver docstring del módulo) y,
    si se da `distancia_corredor` (de trazar_corredor_descendente), un
    perfil de ΔNDMI del corredor por bandas de distancia D8 desde la
    cicatriz. La prueba corredor-vs-control (Mann-Whitney + Welch t) SIEMPRE
    se reporta junto con el aviso de pseudo-replicación (ver docstring del
    módulo) -- nunca sola."""
    delta = ndmi_post - ndmi_pre

    def _stats_grupo(mascara, etiqueta):
        vals = delta[mascara & ~np.isnan(delta)]
        if len(vals) == 0:
            return {"grupo": etiqueta, "n_pixeles": 0, "delta_ndmi_media": None,
                     "delta_ndmi_mediana": None, "delta_ndmi_std": None}
        return {
            "grupo": etiqueta, "n_pixeles": int(len(vals)),
            "delta_ndmi_media": round(float(np.mean(vals)), 5),
            "delta_ndmi_mediana": round(float(np.median(vals)), 5),
            "delta_ndmi_std": round(float(np.std(vals)), 5),
        }

    resumen_global = [
        _stats_grupo(mascara_corredor, "corredor_descendente"),
        _stats_grupo(mascara_control, "control_terreno_similar"),
        _stats_grupo(mascara_quemado, "cicatriz_quemada_referencia_no_es_parte_de_la_prueba"),
    ]

    vals_corredor = delta[mascara_corredor & ~np.isnan(delta)]
    vals_control = delta[mascara_control & ~np.isnan(delta)]
    prueba = {"n_corredor": int(len(vals_corredor)), "n_control": int(len(vals_control))}
    if len(vals_corredor) >= 2 and len(vals_control) >= 2:
        u_stat, p_mw = sstats.mannwhitneyu(vals_corredor, vals_control, alternative="two-sided")
        t_stat, p_t = sstats.ttest_ind(vals_corredor, vals_control, equal_var=False)
        prueba.update({
            "diferencia_medias": round(float(np.mean(vals_corredor) - np.mean(vals_control)), 5),
            "mannwhitney_p": round(float(p_mw), 5),
            "welch_t_p": round(float(p_t), 5),
            "aviso": "p-valores OPTIMISTAS -- pixeles vecinos no son observaciones independientes "
                     "(pseudo-replicacion), ver docstring del modulo. No usar como prueba concluyente sola.",
        })
    else:
        prueba["aviso"] = "insuficientes pixeles validos en corredor y/o control para prueba estadistica"

    perfil_distancia = []
    if distancia_corredor is not None and np.any(mascara_corredor):
        d_validos = distancia_corredor[mascara_corredor]
        bordes = np.unique(np.quantile(d_validos, np.linspace(0, 1, n_bins_distancia + 1)))
        bins_pix = np.digitize(d_validos, bordes[1:-1]) if len(bordes) > 2 else np.zeros_like(d_validos, dtype=int)
        ys, xs = np.where(mascara_corredor)
        for b in sorted(set(bins_pix.tolist())):
            sel = bins_pix == b
            idx_y, idx_x = ys[sel], xs[sel]
            vals = delta[idx_y, idx_x]
            vals = vals[~np.isnan(vals)]
            if len(vals) == 0:
                continue
            perfil_distancia.append({
                "bin_distancia_pasos": int(b),
                "distancia_min_pasos": int(d_validos[sel].min()),
                "distancia_max_pasos": int(d_validos[sel].max()),
                "n_pixeles": int(len(vals)),
                "delta_ndmi_media_corredor": round(float(np.mean(vals)), 5),
            })

    return {
        "resumen_global": resumen_global,
        "prueba_corredor_vs_control": prueba,
        "perfil_por_distancia": perfil_distancia,
    }


# ==============================================================================
# --- EARTH ENGINE: NDMI de UNA ventana, alineado a la malla del terreno ------
# ==============================================================================
def _descargar_ndmi_ventana_alineado(ee, geom_wgs84_visual, fecha_ini, fecha_fin, transform_ref, shape_ref,
                                      utm_crs, nubosidad_max_pct=None, carpeta_tmp=None, etiqueta="ventana"):
    """Compuesto Sentinel-2 (mediana) de UNA ventana de fechas, NDMI =
    normalizedDifference(['B11','B8']), realineado a la misma malla que ya
    usa el terreno 3D de geomatica.py -- mismo patrón exacto que
    core/validacion_incendios._descargar_dnbr_alineado(), pero para una sola
    ventana (no una resta pre/post): aquí la resta post-pre se hace afuera,
    en comparar_ndmi_corredor_control(), sobre dos llamadas independientes a
    esta función. scale=30, igual que el resto del proyecto usa para
    alinear con la malla derivada del SRTM/Hansen (30m) -- Earth Engine
    remuestrea B11 (nativo 20m) internamente, mismo criterio ya documentado
    en validacion_incendios.py para B12."""
    import tempfile

    nubosidad_max_pct = nubosidad_max_pct if nubosidad_max_pct is not None else DEFORESTACION_NUBOSIDAD_MAX_PCT
    carpeta_tmp = carpeta_tmp or tempfile.gettempdir()
    os.makedirs(carpeta_tmp, exist_ok=True)

    aoi = ee.Geometry(geom_wgs84_visual)
    col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
           .filterBounds(aoi).filterDate(fecha_ini, fecha_fin)
           .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", nubosidad_max_pct))
           .select(["B8", "B11"]))
    n = col.size().getInfo()
    log(f"  NDMI {etiqueta} ({fecha_ini} a {fecha_fin}): {n} imágenes Sentinel-2 con <{nubosidad_max_pct}% nubes")
    if n == 0:
        return {"ndmi_alineado": None, "n_imagenes": 0}

    ndmi = col.median().normalizedDifference(["B11", "B8"]).rename("ndmi").clip(aoi)
    url = ndmi.getDownloadURL({"region": aoi, "scale": 30, "crs": "EPSG:4326", "format": "GEO_TIFF"})
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    tif_crudo = os.path.join(carpeta_tmp, f"temp_ndmi_corredor_{etiqueta}.tif")
    with open(tif_crudo, "wb") as f:
        f.write(r.content)

    rows, cols_ = shape_ref
    ndmi_alineado = np.full((rows, cols_), np.nan, dtype=np.float32)  # NaN, no cero -- 0.0 NDMI es un
    # valor real (sin humedad relativa neta), no "sin dato".
    with rasterio.open(tif_crudo) as src:
        reproject(
            source=rasterio.band(src, 1), destination=ndmi_alineado,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=transform_ref, dst_crs=utm_crs,
            resampling=Resampling.bilinear, src_nodata=src.nodata,
        )
    os.remove(tif_crudo)
    return {"ndmi_alineado": ndmi_alineado, "n_imagenes": n}


# ==============================================================================
# --- MAPA 3D: cicatriz / corredor / control sobre el mismo terreno D8 --------
# ==============================================================================
def construir_capas_corredor_3d(mascara_quemado, mascara_corredor, mascara_control, delta_ndmi, hidrologia,
                                 utm_crs=None):
    """Tres capas go.Scatter3d (cicatriz=rojo, corredor=azul, control=gris),
    listas para pasarse como `capas_extra` a geomatica.generar_mapa_3d() --
    mismo patrón que construir_capa_validacion_incendio_3d() de
    core/validacion_incendios.py. El hover de cada punto trae su ΔNDMI real,
    para poder inspeccionar visualmente si el corredor se ve distinto del
    control sin tener que salir del mapa."""
    import plotly.graph_objects as go
    from core.geomatica import calcular_grid_latlon

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
        ("Cicatriz quemada (referencia)", "red", mascara_quemado),
        ("Corredor descendente (aguas abajo, nunca quemado)", "dodgerblue", mascara_corredor),
        ("Control (terreno similar, sin flujo de la cicatriz)", "lightgray", mascara_control),
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

        delta_pts = delta_ndmi[fy, fx]
        customdata_cols = [Z_raw[fy, fx], delta_pts]
        hovertemplate = f"{etiqueta}<br>Altitud: %{{customdata[0]:.0f}} msnm<br>ΔNDMI: %{{customdata[1]:.3f}}"
        if lat_grid is not None:
            customdata_cols += [lat_grid[fy, fx], lon_grid[fy, fx]]
            hovertemplate += "<br>Lat: %{customdata[2]:.5f}<br>Lon: %{customdata[3]:.5f}"
        hovertemplate += "<extra></extra>"
        customdata = np.column_stack(customdata_cols)

        capas.append(go.Scatter3d(
            x=x_km, y=y_km, z=z_km, mode="markers",
            marker=dict(size=2.6, color=color, opacity=0.85),
            customdata=customdata, hovertemplate=hovertemplate, name=etiqueta,
        ))
        log(f"Capa '{etiqueta}': {len(fx)} píxeles.")
    return capas


# ==============================================================================
# --- ORQUESTADOR CON DATOS REALES ---------------------------------------------
# ==============================================================================
def procesar_corredor_descendente(geojson_path, id_proyecto, fecha_evento, anio_hansen, zonas_m=None,
                                   ventana_dias=60, umbral_dnbr_quemado=None, max_pasos_corredor=300,
                                   n_bins_elev=5, n_bins_pendiente=3, factor_control=1,
                                   nubosidad_max_pct=None, percentil_cauce=None, carpeta_salida=None,
                                   carpeta_srtm=None, proyecto_gee=None, mapa_3d=False):
    """Pipeline completo. Reusa terreno/D8/SRTM en caché (igual que
    validacion_incendios.py) y gasta cupo GEE en solo TRES descargas nuevas:
    lossyear (barato, banda categórica), dNBR del evento (para ubicar la
    cicatriz -- ya lo pagas si corriste validacion_incendios.py para el
    mismo evento, aquí se vuelve a pagar porque este módulo no depende de
    que ya lo hayas corrido) y dos ventanas NDMI (pre/post). Nada de esto
    reprocesa años de más -- a diferencia de --historial en
    validacion_incendios.py, aquí SIEMPRE es un evento puntual.

    Las ventanas NDMI pre/post: `ventana_dias` días ANTES de `fecha_evento`
    (vegetación previa al incendio) contra la MISMA ventana calendario un
    AÑO después (no inmediatamente después del incendio) -- comparar la
    misma temporada un año aparte cancela el ciclo estacional normal
    (seco/lluvias) que de otro modo se confundiría con el efecto que se
    busca medir. Además le da tiempo a la vegetación de mostrar cualquier
    efecto de humedad más allá del shock inmediato de la quema."""
    from core import geomatica
    from core.deforestacion import _descargar_lossyear_alineado
    from core.validacion_incendios import _descargar_dnbr_alineado, validar_perdida_contra_incendio
    import pyproj
    from shapely.ops import transform as shp_transform
    from shapely.geometry import mapping
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
    umbral_dnbr_quemado = (umbral_dnbr_quemado if umbral_dnbr_quemado is not None
                            else VALIDACION_INCENDIO_UMBRAL_DNBR_QUEMADO)
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    carpeta_srtm = carpeta_srtm or CARPETA_SRTM
    os.makedirs(carpeta_salida, exist_ok=True)

    log("Cargando terreno y D8 (reusa SRTM en caché)...")
    geom_utm_nucleo, dst_array, meta_utm, utm_crs = geomatica.cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        percentil_cauce, carpeta_srtm, id_proyecto,
    )
    if "flow_dir" not in hidrologia:
        raise RuntimeError(
            "calcular_hidrologia_d8() no devolvió 'flow_dir' -- este módulo requiere una versión de "
            "core/geomatica.py que exponga la dirección de flujo D8 cruda, no solo cauces/acumulación."
        )
    flow_dir = hidrologia["flow_dir"]
    Z_raw = hidrologia["Z_raw"]
    zona_de_pixel = hidrologia["zona_de_pixel"]
    mascara_valida = (zona_de_pixel != "fuera") & ~np.isnan(Z_raw)

    pw_v, ph_v = hidrologia["pw_v"], hidrologia["ph_v"]
    z_grad_base = np.where(mascara_valida, Z_raw, np.nanmean(Z_raw[mascara_valida]))
    dy, gx = np.gradient(z_grad_base, ph_v, pw_v)
    pendiente_deg = np.degrees(np.arctan(np.sqrt(gx ** 2 + dy ** 2)))

    a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True).transform
    geom_visual_geojson = mapping(shp_transform(a_wgs84, geom_utm_nucleo.buffer(max(zonas_m))))

    log("Descargando Hansen lossyear alineado...")
    lossyear_alineado = _descargar_lossyear_alineado(
        ee, geom_visual_geojson, hidrologia["transform"], Z_raw.shape, utm_crs, carpeta_tmp=carpeta_srtm,
    )
    lossyear_alineado = np.where(mascara_valida, lossyear_alineado, 255)

    log(f"Descargando dNBR (evento confirmado {fecha_evento}) para ubicar la cicatriz real...")
    dnbr_resultado = _descargar_dnbr_alineado(
        ee, geom_visual_geojson, fecha_evento, hidrologia["transform"], Z_raw.shape, utm_crs,
        nubosidad_max_pct=nubosidad_max_pct, carpeta_tmp=carpeta_srtm,
    )
    if dnbr_resultado["dnbr_alineado"] is None:
        raise RuntimeError(
            "No se pudo calcular dNBR para ubicar la cicatriz -- sin suficientes imágenes Sentinel-2 "
            "limpias en la ventana pre/post del evento. Prueba con --nubosidad-max más alto."
        )
    resultado_val = validar_perdida_contra_incendio(
        lossyear_alineado, dnbr_resultado["dnbr_alineado"], anio_hansen, umbral_dnbr_quemado,
    )
    mascara_quemado = resultado_val["mascara_confirmados_dnbr"] & mascara_valida
    n_quemado = int(mascara_quemado.sum())
    log(f"Cicatriz confirmada por dNBR: {n_quemado} píxeles dentro de la zona visual.")
    if n_quemado == 0:
        raise RuntimeError(
            "Cero píxeles de cicatriz confirmada -- no hay semilla para trazar el corredor descendente. "
            "Revisa --fecha-evento/--anio-hansen o baja --umbral-dnbr."
        )

    log(f"Trazando corredor descendente (D8, máx {max_pasos_corredor} pasos)...")
    mascara_corredor, distancia = trazar_corredor_descendente(
        flow_dir, mascara_quemado, mascara_valida, max_pasos=max_pasos_corredor,
    )
    n_corredor = int(mascara_corredor.sum())
    log(f"  -> {n_corredor} píxeles en el corredor descendente (nunca quemados, aguas abajo de la cicatriz).")
    if n_corredor == 0:
        log("Corredor descendente vacío -- puede pasar con cicatrices muy pequeñas o si la cicatriz está "
            "en el borde exacto de la zona visual sin nada aguas abajo dentro de ella. No hay nada que "
            "comparar; el resto del pipeline seguirá pero con corredor=0.", nivel="WARN")

    # Excluye del control cualquier píxel con pérdida Hansen propia en los últimos ~10 años (no solo el
    # año del evento) -- terreno recién talado por otra razón tampoco es un buen "control sin disturbio".
    codigo_min_reciente = max(1, anio_hansen - 10 - 2000)
    mascara_perdida_reciente = (lossyear_alineado >= codigo_min_reciente) & (lossyear_alineado <= (anio_hansen - 2000))
    mascara_excluir_control = mascara_quemado | mascara_corredor | mascara_perdida_reciente

    log("Seleccionando grupo control (pareado por elevación y pendiente)...")
    mascara_control, resumen_bins = seleccionar_grupo_control(
        mascara_corredor, mascara_excluir_control, Z_raw, pendiente_deg,
        n_bins_elev=n_bins_elev, n_bins_pendiente=n_bins_pendiente, factor_control=factor_control,
    )
    n_control = int(mascara_control.sum())
    log(f"  -> {n_control} píxeles de control.")
    for b in resumen_bins:
        if b["n_tomado_control"] < b["n_pedido_control"]:
            log(f"  Bin elev={b['bin_elevacion']}/pend={b['bin_pendiente']}: solo {b['n_disponible_control']} "
                f"candidatos disponibles para {b['n_pedido_control']} pedidos ({b['n_corredor']} px de "
                f"corredor) -- control incompleto en ese bin, no se inventan píxeles.", nivel="WARN")

    fecha_dt = datetime.strptime(fecha_evento, "%Y-%m-%d")
    pre_ini = (fecha_dt - timedelta(days=ventana_dias)).strftime("%Y-%m-%d")
    pre_fin = fecha_evento
    post_ini = (fecha_dt + timedelta(days=365 - ventana_dias)).strftime("%Y-%m-%d")
    post_fin = (fecha_dt + timedelta(days=365)).strftime("%Y-%m-%d")
    log(f"Ventanas NDMI -- PRE: {pre_ini} a {pre_fin} (antes del evento) | POST: {post_ini} a {post_fin} "
        f"(misma temporada, 1 año después -- para no confundir estacionalidad con efecto real).")

    ndmi_pre_res = _descargar_ndmi_ventana_alineado(
        ee, geom_visual_geojson, pre_ini, pre_fin, hidrologia["transform"], Z_raw.shape, utm_crs,
        nubosidad_max_pct=nubosidad_max_pct, carpeta_tmp=carpeta_srtm, etiqueta="pre",
    )
    ndmi_post_res = _descargar_ndmi_ventana_alineado(
        ee, geom_visual_geojson, post_ini, post_fin, hidrologia["transform"], Z_raw.shape, utm_crs,
        nubosidad_max_pct=nubosidad_max_pct, carpeta_tmp=carpeta_srtm, etiqueta="post",
    )
    if ndmi_pre_res["ndmi_alineado"] is None or ndmi_post_res["ndmi_alineado"] is None:
        raise RuntimeError(
            "No se pudo calcular NDMI pre y/o post -- sin suficientes imágenes Sentinel-2 limpias en "
            "alguna de las dos ventanas. Prueba con --ventana-dias más grande o --nubosidad-max más alto."
        )

    resultado_comparacion = comparar_ndmi_corredor_control(
        ndmi_pre_res["ndmi_alineado"], ndmi_post_res["ndmi_alineado"],
        mascara_corredor, mascara_control, mascara_quemado, distancia_corredor=distancia,
    )

    pixel_area_ha = pw_v * ph_v / 10000.0
    resumen_ha = {
        "ha_cicatriz_quemada": round(n_quemado * pixel_area_ha, 3),
        "ha_corredor_descendente": round(n_corredor * pixel_area_ha, 3),
        "ha_control": round(n_control * pixel_area_ha, 3),
    }
    log(f"Áreas: cicatriz={resumen_ha['ha_cicatriz_quemada']}ha, "
        f"corredor={resumen_ha['ha_corredor_descendente']}ha, control={resumen_ha['ha_control']}ha")
    for fila in resultado_comparacion["resumen_global"]:
        log(f"  {fila['grupo']}: n={fila['n_pixeles']}, ΔNDMI media={fila['delta_ndmi_media']}, "
            f"mediana={fila['delta_ndmi_mediana']}")
    prueba = resultado_comparacion["prueba_corredor_vs_control"]
    if "mannwhitney_p" in prueba:
        log(f"  Corredor vs control: diferencia de medias={prueba['diferencia_medias']}, "
            f"Mann-Whitney p={prueba['mannwhitney_p']}, Welch t p={prueba['welch_t_p']} -- {prueba['aviso']}")
    else:
        log(f"  Corredor vs control: {prueba.get('aviso')}", nivel="WARN")
    if resultado_comparacion["perfil_por_distancia"]:
        log("  Perfil ΔNDMI del corredor por distancia (pasos D8) desde la cicatriz:")
        for fila in resultado_comparacion["perfil_por_distancia"]:
            log(f"    {fila['distancia_min_pasos']}-{fila['distancia_max_pasos']} pasos: "
                f"n={fila['n_pixeles']}, ΔNDMI media={fila['delta_ndmi_media_corredor']}")

    fila_resumen_csv = {
        "id_proyecto": id_proyecto, "fecha_evento": fecha_evento, "anio_hansen": anio_hansen,
        **resumen_ha,
        "pre_rango": f"{pre_ini} a {pre_fin}", "post_rango": f"{post_ini} a {post_fin}",
        "n_imagenes_s2_pre": ndmi_pre_res["n_imagenes"], "n_imagenes_s2_post": ndmi_post_res["n_imagenes"],
    }
    for fila in resultado_comparacion["resumen_global"]:
        prefijo = fila["grupo"].split("_")[0]
        fila_resumen_csv[f"n_pixeles_{prefijo}"] = fila["n_pixeles"]
        fila_resumen_csv[f"delta_ndmi_media_{prefijo}"] = fila["delta_ndmi_media"]
        fila_resumen_csv[f"delta_ndmi_mediana_{prefijo}"] = fila["delta_ndmi_mediana"]
    fila_resumen_csv.update({f"prueba_{k}": v for k, v in prueba.items()})

    df_resumen = pd.DataFrame([fila_resumen_csv])
    csv_resumen_path = os.path.join(
        carpeta_salida, f"corredor_descendente_resumen_{id_proyecto.lower()}_{anio_hansen}.csv"
    )
    df_resumen.to_csv(csv_resumen_path, index=False)
    log(f"CSV resumen guardado en: {csv_resumen_path}", nivel="OK")

    csv_perfil_path = None
    df_perfil = pd.DataFrame(resultado_comparacion["perfil_por_distancia"])
    if not df_perfil.empty:
        csv_perfil_path = os.path.join(
            carpeta_salida, f"corredor_descendente_perfil_distancia_{id_proyecto.lower()}_{anio_hansen}.csv"
        )
        df_perfil.to_csv(csv_perfil_path, index=False)
        log(f"CSV perfil por distancia guardado en: {csv_perfil_path}", nivel="OK")

    html_path = None
    if mapa_3d:
        delta_ndmi = ndmi_post_res["ndmi_alineado"] - ndmi_pre_res["ndmi_alineado"]
        capas_extra = construir_capas_corredor_3d(
            mascara_quemado, mascara_corredor, mascara_control, delta_ndmi, hidrologia, utm_crs=utm_crs,
        )
        delta_corredor = next((f["delta_ndmi_media"] for f in resultado_comparacion["resumen_global"]
                                if f["grupo"].startswith("corredor")), None)
        delta_control = next((f["delta_ndmi_media"] for f in resultado_comparacion["resumen_global"]
                               if f["grupo"].startswith("control")), None)
        subtitulo = (f"Corredor descendente D8 desde incendio {fecha_evento}: {n_corredor}px "
                     f"({resumen_ha['ha_corredor_descendente']}ha) vs {n_control}px control -- "
                     f"ΔNDMI corredor={delta_corredor}, ΔNDMI control={delta_control} "
                     f"(n=1 caso de estudio -- ver CSV para limitaciones estadísticas)")
        html_path = os.path.join(
            carpeta_salida, f"{id_proyecto.lower()}_3d_corredor_descendente_{anio_hansen}.html"
        )
        titulo_base = f"{id_proyecto} -- Corredor hidrológico descendente desde incendio {fecha_evento}"
        geomatica.generar_mapa_3d(
            hidrologia, id_proyecto, html_path, subtitulo=subtitulo, utm_crs=utm_crs,
            capas_extra=capas_extra, titulo_base=titulo_base,
        )
        log(f"Mapa 3D: {html_path}")

    return resultado_comparacion, csv_resumen_path, csv_perfil_path, html_path


# ==============================================================================
# --- MODO DEMO: terreno/flujo sintéticos, sin Earth Engine, sin red ----------
# ==============================================================================
def demo():
    """Prueba las funciones puras (trazar_corredor_descendente,
    seleccionar_grupo_control, comparar_ndmi_corredor_control) sobre una
    ladera sintética con flujo uniforme hacia el sur y una caída de NDMI
    artificial SOLO en el corredor esperado -- para verificar, sin red ni
    Earth Engine, que el módulo: (1) traza el corredor correcto excluyendo
    la cicatriz, (2) arma un control sin traslape, y (3) detecta la caída
    sintética con un p bajo."""
    log("=== core.corredor_descendente --demo (sintético, sin Earth Engine) ===")
    rng = np.random.default_rng(7)
    size = 40
    y, x = np.mgrid[0:size, 0:size]
    elevacion = (3800 - 3.0 * y + rng.normal(0, 1.0, size=(size, size))).astype(np.float64)
    pendiente = np.full((size, size), 15.0) + rng.normal(0, 1.0, size=(size, size))
    mascara_valida = np.ones((size, size), dtype=bool)

    # flow_dir sintético: todo fluye hacia el sur (código 4) -- coherente con la ladera de arriba
    # (elevación decrece con la fila, o sea "sur" = fila creciente = cuesta abajo).
    flow_dir = np.full((size, size), 4, dtype=np.uint8)

    mascara_quemado = np.zeros((size, size), dtype=bool)
    mascara_quemado[5:9, 15:20] = True  # bloque de "cicatriz" cerca de arriba de la ladera

    mascara_corredor, distancia = trazar_corredor_descendente(
        flow_dir, mascara_quemado, mascara_valida, max_pasos=100,
    )
    n_esperado = 5 * (size - 9)  # 5 columnas, desde la fila justo debajo de la cicatriz hasta el borde
    log(f"Corredor sintético: {mascara_corredor.sum()} píxeles (esperado {n_esperado})")
    assert mascara_corredor.sum() == n_esperado, "el corredor sintético no tiene el tamaño esperado"
    assert not np.any(mascara_corredor & mascara_quemado), "el corredor NO debe incluir píxeles quemados"

    mascara_excluir = mascara_quemado | mascara_corredor
    mascara_control, resumen_bins = seleccionar_grupo_control(
        mascara_corredor, mascara_excluir, elevacion, pendiente, n_bins_elev=3, n_bins_pendiente=2,
    )
    log(f"Control sintético: {mascara_control.sum()} píxeles en {len(resumen_bins)} bins")
    assert not np.any(mascara_control & mascara_corredor), "control y corredor no deben traslaparse"
    assert not np.any(mascara_control & mascara_quemado), "control no debe incluir píxeles quemados"

    ndmi_pre = np.full((size, size), 0.30) + rng.normal(0, 0.01, size=(size, size))
    ndmi_post = ndmi_pre.copy()
    ndmi_post[mascara_corredor] -= 0.05  # caída artificial SOLO en el corredor -- lo que se busca detectar

    resultado = comparar_ndmi_corredor_control(
        ndmi_pre, ndmi_post, mascara_corredor, mascara_control, mascara_quemado, distancia_corredor=distancia,
    )
    print("\n--- Resumen por grupo (demo, sintético) ---")
    for fila in resultado["resumen_global"]:
        print(f"  {fila['grupo']}: n={fila['n_pixeles']}, ΔNDMI media={fila['delta_ndmi_media']}")
    prueba = resultado["prueba_corredor_vs_control"]
    print(f"  Prueba corredor vs control: {prueba}")
    assert prueba["mannwhitney_p"] < 0.05, (
        "con una caída sintética clara y separada del control, se espera un p bajo -- si no sale así, "
        "hay un bug en la comparación."
    )
    log("Demo OK -- el módulo distingue corredor de control y detecta una caída sintética clara.", nivel="OK")
    return resultado


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Corredor hidrológico descendente (D8) desde una cicatriz de incendio/pérdida -- "
                    "prueba si hay un efecto de humedad detectable laderas abajo -- Motor Nacional"
    )
    ap.add_argument("--demo", action="store_true", help="Corre con terreno/flujo sintéticos, sin Earth Engine ni red")
    ap.add_argument("--geojson", type=str, help="Ruta al GeoJSON del polígono núcleo")
    ap.add_argument("--id-proyecto", type=str, help="Nombre identificador del sitio (para nombres de archivo)")
    ap.add_argument("--fecha-evento", type=str, help="Fecha real del incendio/evento, formato YYYY-MM-DD")
    ap.add_argument("--anio-hansen", type=int, help="Año de pérdida Hansen del evento (ej. 2025)")
    ap.add_argument("--zonas", type=str, default=None, help="Buffers en metros separados por coma, ej. '0,500,1000'")
    ap.add_argument("--ventana-dias", type=int, default=60,
                     help="Días de la ventana NDMI antes del evento (y de la ventana espejo, un año "
                          "después) -- default 60")
    ap.add_argument("--umbral-dnbr", type=float, default=None,
                     help="dNBR mínimo para considerar 'quemado' (default: config.VALIDACION_INCENDIO_UMBRAL_DNBR_QUEMADO)")
    ap.add_argument("--max-pasos-corredor", type=int, default=300,
                     help="Tope de pasos D8 al trazar el corredor descendente (protección, no debería activarse)")
    ap.add_argument("--bins-elevacion", type=int, default=5, help="Número de bins de elevación para el control")
    ap.add_argument("--bins-pendiente", type=int, default=3, help="Número de bins de pendiente para el control")
    ap.add_argument("--factor-control", type=int, default=1,
                     help="Múltiplo del tamaño del corredor a buscar por bin en el control (default 1x)")
    ap.add_argument("--nubosidad-max", type=float, default=None)
    ap.add_argument("--percentil-cauce", type=float, default=None)
    ap.add_argument("--carpeta-srtm", type=str, default=None)
    ap.add_argument("--carpeta-salida", type=str, default=None)
    ap.add_argument("--proyecto-gee", type=str, default=None,
                     help="ID de proyecto de Google Cloud para ee.Initialize(project=...)")
    ap.add_argument("--mapa-3d", action="store_true", help="Genera el mapa 3D (cicatriz/corredor/control)")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not args.geojson or not args.id_proyecto or not args.fecha_evento or not args.anio_hansen:
        ap.error("--geojson, --id-proyecto, --fecha-evento y --anio-hansen son obligatorios fuera de --demo")

    zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
    procesar_corredor_descendente(
        geojson_path=args.geojson, id_proyecto=args.id_proyecto, fecha_evento=args.fecha_evento,
        anio_hansen=args.anio_hansen, zonas_m=zonas_m, ventana_dias=args.ventana_dias,
        umbral_dnbr_quemado=args.umbral_dnbr, max_pasos_corredor=args.max_pasos_corredor,
        n_bins_elev=args.bins_elevacion, n_bins_pendiente=args.bins_pendiente, factor_control=args.factor_control,
        nubosidad_max_pct=args.nubosidad_max, percentil_cauce=args.percentil_cauce,
        carpeta_srtm=args.carpeta_srtm, carpeta_salida=args.carpeta_salida, proyecto_gee=args.proyecto_gee,
        mapa_3d=args.mapa_3d,
    )


if __name__ == "__main__":
    main()

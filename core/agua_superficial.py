#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extensión de agua superficial VISIBLE por año (JRC Global Surface Water),
por zona -- para responder con evidencia satelital directa "¿hay más o menos
agua año con año?", en vez de inferirlo del D8 (core/geomatica.py +
core/validacion_hidrologica.py), que es una foto fija del relieve (un solo
SRTM) sin ningún eje temporal.

POR QUÉ ESTE MÓDULO ESTÁ SEPARADO (mismo principio que carbono.py,
deforestacion.py, validacion_hidrologica.py): una fuente de datos externa
por módulo -- aquí, JRC Global Surface Water vía Earth Engine.

LIMITACIÓN REAL, NO ESCONDIDA -- léela ANTES de interpretar cualquier
resultado de este módulo:
    JRC Global Surface Water usa Landsat a 30m de resolución y está
    calibrado para AGUA ABIERTA -- lagos, presas, ríos anchos. Un cauce de
    montaña de unos metros de ancho (como las barrancas de Cofre de Perote
    que SÍ detecta el D8, ver core/validacion_hidrologica.py) puede ser
    demasiado angosto para que Landsat lo clasifique como agua. Un
    resultado en 0 o casi 0 hectáreas NO significa "no hay agua" -- puede
    significar solo que este sensor no la puede ver a esa escala. Si el
    sitio tiene un cuerpo de agua más ancho (laguna, bordo, presa, tramo
    ancho de río), ahí sí es la herramienta correcta.

    Rango de años: JRC_ANIO_MIN-JRC_ANIO_MAX (ver constantes abajo) --
    versión JRC_GSW1_4, la más reciente disponible en el catálogo público
    de Earth Engine al escribir esto (verificado 2026-08-31, fuente:
    https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_4_YearlyHistory).
    JRC ya publicó una versión más nueva (hasta 2024) fuera de Earth Engine
    (https://global-surface-water.appspot.com/download) -- si en el futuro
    también aparece como colección de Earth Engine, actualizar COLECCION_JRC
    y JRC_ANIO_MAX abajo.

QUÉ SÍ CALCULA:
    - Hectáreas de agua ESTACIONAL (waterClass=2) y PERMANENTE
      (waterClass=3) por año, por zona (nucleo + buffers acumulativos,
      misma convención que core.deforestacion) -- nunca sumadas de entrada,
      para poder ver si un año se volvió más "estacional" en vez de
      desaparecer del todo.
    - Un resumen por ANILLO EXCLUSIVO (sin traslape, sumable) de
      hectáreas de agua TOTAL (estacional+permanente), por año -- ver
      generar_resumen_agua_sin_traslape().
    - Una gráfica de LÍNEA 2D (no un mapa 3D -- ver generar_grafica_agua_anual)
      con la tendencia año a año.

QUÉ NO CALCULA:
    - Volumen de agua, caudal, ni nada en litros/m3/s -- esto es SUPERFICIE
      clasificada como agua (hectáreas), no un aforo hidrológico.
    - Ningún ajuste ni corrección del modelo D8 -- son dos fuentes
      independientes, cada una con su propia limitación documentada.
"""

import argparse
import os

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.warp import reproject, Resampling
import requests
from shapely.geometry import mapping

from config import ZONAS_ANALISIS_M, PERCENTIL_CAUCE_HIDROLOGIA, CARPETA_SRTM, log
from core.deforestacion import _reproyectores_utm, geom_zona_precisa

COLECCION_JRC = "JRC/GSW1_4/YearlyHistory"
JRC_ANIO_MIN = 1984
JRC_ANIO_MAX = 2021


# ==============================================================================
# --- JRC GLOBAL SURFACE WATER: HISTORIAL ANUAL PARA UN POLÍGONO ---
# ==============================================================================
def calcular_agua_anual_gee(ee, geom_wgs84_geojson, anio_inicio=None, anio_fin=None):
    """Hectáreas de agua ESTACIONAL y PERMANENTE por año, para UN polígono
    (dict geojson, WGS84), vía JRC Global Surface Water. Un solo
    reduceRegion() para TODOS los años del rango a la vez -- se apila cada
    año como una banda separada de una sola imagen y se pide un
    frequencyHistogram por banda en una sola llamada (mismo truco que
    Hansen en core.deforestacion.calcular_deforestacion_anual_gee), en vez
    de un reduceRegion por año -- mucho más barato en EECU (créditos de
    Earth Engine) que pedir un compuesto por año como hace Sentinel-2 en
    calcular_indices_anuales_gee.

    anio_inicio/anio_fin se acotan SIEMPRE al rango real de JRC -- nunca se
    inventa un año fuera de lo que el dataset realmente cubre. Devuelve
    area_ha, el historial {año: {ha_estacional, ha_permanente}} y el rango
    de años realmente usado (que puede ser más angosto que el pedido)."""
    anio_inicio_pedido = anio_inicio or JRC_ANIO_MIN
    anio_fin_pedido = anio_fin or JRC_ANIO_MAX
    anio_inicio = max(anio_inicio_pedido, JRC_ANIO_MIN)
    anio_fin = min(anio_fin_pedido, JRC_ANIO_MAX)
    if anio_inicio_pedido < JRC_ANIO_MIN or anio_fin_pedido > JRC_ANIO_MAX:
        log(f"JRC Global Surface Water solo cubre {JRC_ANIO_MIN}-{JRC_ANIO_MAX} -- el rango pedido "
            f"({anio_inicio_pedido}-{anio_fin_pedido}) se acota a {anio_inicio}-{anio_fin}.", nivel="WARN")
    if anio_fin < anio_inicio:
        raise ValueError(f"Rango de años inválido tras acotar a JRC ({JRC_ANIO_MIN}-{JRC_ANIO_MAX}): "
                          f"{anio_inicio}-{anio_fin}")

    aoi = ee.Geometry(geom_wgs84_geojson)
    area_ha = ee.Number(aoi.area(1)).divide(10000).getInfo()

    coleccion = ee.ImageCollection(COLECCION_JRC)
    bandas_por_anio = []
    for anio in range(anio_inicio, anio_fin + 1):
        img_anio = coleccion.filterDate(f"{anio}-01-01", f"{anio + 1}-01-01").first()
        bandas_por_anio.append(img_anio.select("waterClass").rename(f"y{anio}"))

    imagen_apilada = ee.Image.cat(bandas_por_anio)
    stats = imagen_apilada.reduceRegion(
        reducer=ee.Reducer.frequencyHistogram(), geometry=aoi, scale=30, maxPixels=1e9, bestEffort=True,
    ).getInfo()

    historial = {}
    for anio in range(anio_inicio, anio_fin + 1):
        hist = stats.get(f"y{anio}") or {}
        n_estacional = hist.get("2", 0)
        n_permanente = hist.get("3", 0)
        # 900 m2/pixel (30x30m de JRC) -> hectareas
        historial[anio] = {
            "ha_estacional": round(n_estacional * 900 / 10000, 4),
            "ha_permanente": round(n_permanente * 900 / 10000, 4),
        }

    return {"area_ha": area_ha, "historial": historial, "anio_inicio": anio_inicio, "anio_fin": anio_fin}


# ==============================================================================
# --- ORQUESTADOR TABULAR (por zona -- núcleo + buffers) ---
# ==============================================================================
def calcular_agua_por_zona_real(ee, geom_wgs84_nucleo, zonas_m, anio_inicio=None, anio_fin=None):
    """JRC Global Surface Water por año para cada zona (núcleo + buffers
    ACUMULATIVOS -- mismo patrón que core.deforestacion.calcular_deforestacion_por_zona_real:
    buffer_500m SIEMPRE incluye completo al núcleo, no es traslape/muñeca
    rusa, es geometría -- ver ese docstring). Devuelve (df_resumen, df_historial):
      - df_resumen: una fila por zona (área, rango de años, promedio/mínimo/
        máximo de hectáreas de agua total en el periodo -- el resumen
        rápido para "¿hay más o menos agua?" sin tener que leer el CSV
        largo entero).
      - df_historial: formato LARGO, una fila por (zona, año)."""
    a_utm, a_wgs84, _ = _reproyectores_utm(geom_wgs84_nucleo)
    filas_resumen, filas_historial = [], []

    for buf_m in zonas_m:
        etiqueta = "nucleo" if buf_m == 0 else f"buffer_{buf_m}m"
        geom_zona_wgs84, _ = geom_zona_precisa(geom_wgs84_nucleo, buf_m, a_utm, a_wgs84)

        log(f"Zona {etiqueta} -- consultando JRC Global Surface Water por año...")
        agua = calcular_agua_anual_gee(ee, mapping(geom_zona_wgs84), anio_inicio, anio_fin)

        anios_orden = sorted(agua["historial"])
        ha_totales = [agua["historial"][a]["ha_estacional"] + agua["historial"][a]["ha_permanente"]
                      for a in anios_orden]
        anio_min_agua = anios_orden[int(np.argmin(ha_totales))]
        anio_max_agua = anios_orden[int(np.argmax(ha_totales))]
        log(f"  -> {agua['anio_inicio']}-{agua['anio_fin']}: promedio {np.mean(ha_totales):.2f} ha de agua "
            f"(mínimo {min(ha_totales):.2f} ha en {anio_min_agua}, máximo {max(ha_totales):.2f} ha en {anio_max_agua})")

        filas_resumen.append({
            "zona": etiqueta, "buffer_m": buf_m, "area_ha": round(agua["area_ha"], 3),
            "anio_inicio": agua["anio_inicio"], "anio_fin": agua["anio_fin"],
            "ha_agua_promedio": round(float(np.mean(ha_totales)), 4),
            "ha_agua_minimo": round(float(min(ha_totales)), 4), "anio_minimo": anio_min_agua,
            "ha_agua_maximo": round(float(max(ha_totales)), 4), "anio_maximo": anio_max_agua,
        })

        for anio in anios_orden:
            vals = agua["historial"][anio]
            filas_historial.append({
                "zona": etiqueta, "buffer_m": buf_m, "anio": anio,
                "ha_agua_estacional": vals["ha_estacional"], "ha_agua_permanente": vals["ha_permanente"],
                "ha_agua_total": round(vals["ha_estacional"] + vals["ha_permanente"], 4),
            })

    return pd.DataFrame(filas_resumen), pd.DataFrame(filas_historial)


def generar_resumen_agua_sin_traslape(df_historial, id_proyecto, carpeta_salida):
    """Anexo de lectura rápida a agua_superficial_historial_anual_*.csv --
    NO lo reemplaza. 'ha_agua_total' ahí es ACUMULATIVA POR ZONA
    (buffer_500m SIEMPRE incluye completo al núcleo, misma convención que
    core.deforestacion.generar_resumen_no_traslapado -- ver ese docstring
    para la prueba geométrica de por qué la resta entre buffers nunca da
    negativo).

    DIFERENCIA A PROPÓSITO con ese mismo resumen de deforestación: pérdida
    de bosque es un EVENTO (ocurre una vez, sí se puede sumar a través de
    los años para un total del periodo). Agua es un ESTADO: sumar
    hectáreas de agua del año 1990 más las del año 1991 no da un "total"
    que signifique algo -- por eso este resumen, a diferencia del de
    deforestación, NO trae una fila "TOTAL {periodo}". Solo trae, para
    CADA año por separado, el anillo exclusivo (restando cada buffer del
    siguiente más grande) y una fila TOTAL de ESE año (= agua en todo el
    sitio ese año, sin traslape).

    Solo opera sobre 'ha_agua_total' (estacional+permanente ya sumados) --
    si en el futuro hace falta el desglose estacional/permanente también
    por anillo exclusivo, se puede correr esta misma lógica dos veces más
    (una por columna), no está incluido aquí para no complicar la primera
    versión de este módulo."""
    faltantes = [c for c in ["zona", "buffer_m", "anio", "ha_agua_total"] if c not in df_historial.columns]
    if faltantes:
        log(f"generar_resumen_agua_sin_traslape: faltan columnas {faltantes} -- no se genera el resumen.", nivel="WARN")
        return None

    buffers_ordenados = sorted(df_historial["buffer_m"].unique())
    nombre_anillo = {}
    for i, buf_m in enumerate(buffers_ordenados):
        if i == 0:
            nombre_anillo[buf_m] = "nucleo (0m)" if buf_m == 0 else f"anillo_0-{buf_m}m"
        else:
            nombre_anillo[buf_m] = f"anillo_{buffers_ordenados[i - 1]}-{buf_m}m"

    pivot = df_historial.pivot_table(index="anio", columns="buffer_m", values="ha_agua_total", aggfunc="first")
    filas = []
    for anio in sorted(pivot.index):
        acumulado_prev, total_anio = 0.0, 0.0
        for buf_m in buffers_ordenados:
            valor = pivot.loc[anio, buf_m]
            valor = 0.0 if pd.isna(valor) else float(valor)
            valor_anillo = round(valor - acumulado_prev, 4)
            if valor_anillo < -0.01:
                # Margen de tolerancia por ruido de borde entre geometrías (cada buffer se calcula
                # con su propio reduceRegion, no compartiendo máscara de píxel) -- un negativo aquí
                # sería un artefacto de redondeo, no una zona con agua "negativa". Se deja en 0.0 y
                # se avisa, en vez de propagar un número sin sentido físico.
                log(f"Anillo {nombre_anillo[buf_m]}, año {anio}: valor negativo ({valor_anillo} ha) -- "
                    f"probable ruido de borde entre geometrías, se deja en 0.0.", nivel="WARN")
                valor_anillo = 0.0
            acumulado_prev = valor
            total_anio += valor_anillo
            filas.append({"anillo": nombre_anillo[buf_m], "anio": anio, "ha_agua": valor_anillo})
        filas.append({"anillo": "TOTAL (sitio completo, sin traslape)", "anio": anio, "ha_agua": round(total_anio, 4)})

    df_out = pd.DataFrame(filas).sort_values(["anio", "anillo"]).reset_index(drop=True)
    nombre_out = f"agua_superficial_resumen_sin_traslape_{id_proyecto.lower()}.csv"
    csv_path = os.path.join(carpeta_salida, nombre_out)
    df_out.to_csv(csv_path, index=False)
    log(f"Resumen SIN traslape (ha_agua por anillo exclusivo, sumable, por año) guardado en: {csv_path}", nivel="OK")
    return csv_path


def procesar_sitio_real(geojson_path, id_proyecto, zonas_m=None, anio_inicio=None, anio_fin=None,
                         carpeta_salida=None, proyecto_gee=None):
    """Pipeline tabular completo: inicializa Earth Engine, consulta JRC
    Global Surface Water por zona y por año, guarda los dos CSV (historial
    largo + resumen sin traslape). Devuelve (df_resumen, df_historial,
    csv_historial, csv_sin_traslape)."""
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

    df_resumen, df_historial = calcular_agua_por_zona_real(ee, geom_wgs84, zonas_m, anio_inicio, anio_fin)

    csv_historial = os.path.join(carpeta_salida, f"agua_superficial_historial_anual_{id_proyecto.lower()}.csv")
    df_historial.to_csv(csv_historial, index=False)
    log(f"CSV historial anual (largo, zona+año, estacional/permanente por separado): {csv_historial}")

    csv_sin_traslape = generar_resumen_agua_sin_traslape(df_historial, id_proyecto, carpeta_salida)

    return df_resumen, df_historial, csv_historial, csv_sin_traslape


# ==============================================================================
# --- GRÁFICA 2D: HECTÁREAS DE AGUA POR AÑO (línea, NO mapa 3D) ---
# ==============================================================================
def generar_grafica_agua_anual(id_proyecto, carpeta_salida=None):
    """Gráfica de línea 2D -- a propósito NO un mapa 3D. Una tendencia en
    el tiempo ("change over time") se lee mejor con un eje temporal simple;
    un terreno 3D obligaría a fijar un solo año para dibujar el relieve y
    escondería justo la tendencia que es el punto de esta gráfica (ver
    skill de dataviz: la forma la elige el trabajo del dato, no al revés).

    Lee el CSV sin traslape YA generado por procesar_sitio_real() -- no
    vuelve a consultar Earth Engine, así que se puede regenerar la gráfica
    (ej. para ajustar el estilo) sin gastar cuota de GEE otra vez."""
    import plotly.graph_objects as go
    from core import geomatica

    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    csv_path = os.path.join(carpeta_salida, f"agua_superficial_resumen_sin_traslape_{id_proyecto.lower()}.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"No se encontró {csv_path} -- corre antes core.agua_superficial (sin --grafica) "
            "para generar el CSV sin traslape."
        )
    df = pd.read_csv(csv_path)

    anillos = [a for a in df["anillo"].unique() if not str(a).startswith("TOTAL")]
    # Colores categóricos en orden fijo (paleta validada del skill de dataviz, slots 1-3) -- nunca
    # ciclados. TOTAL se dibuja aparte, en gris punteado: es un agregado derivado, no una zona más
    # (mismo criterio que "color sigue a la entidad" del skill -- TOTAL no es una entidad geográfica).
    paleta = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]

    fig = go.Figure()
    for i, anillo in enumerate(anillos):
        sub = df[df["anillo"] == anillo].sort_values("anio")
        fig.add_trace(go.Scatter(
            x=sub["anio"], y=sub["ha_agua"], mode="lines+markers", name=anillo,
            line=dict(color=paleta[i % len(paleta)], width=2), marker=dict(size=6),
            hovertemplate=f"{anillo}<br>Año: %{{x}}<br>%{{y:.2f}} ha de agua<extra></extra>",
        ))

    sub_total = df[df["anillo"].astype(str).str.startswith("TOTAL")].sort_values("anio")
    fig.add_trace(go.Scatter(
        x=sub_total["anio"], y=sub_total["ha_agua"], mode="lines", name="TOTAL sitio (sin traslape)",
        line=dict(color="#4d4d4d", width=2, dash="dash"),
        hovertemplate="TOTAL sitio<br>Año: %{x}<br>%{y:.2f} ha de agua<extra></extra>",
    ))

    anio_min, anio_max = int(df["anio"].min()), int(df["anio"].max())
    titulo = (f"{id_proyecto} -- Agua superficial visible por año ({anio_min}-{anio_max})<br>"
              "<sub>JRC Global Surface Water, Landsat 30m -- agua ESTACIONAL + PERMANENTE."
              "<br>Cerca de 0 no significa que no haya agua: cauces angostos de montaña pueden ser "
              "invisibles a esta resolución (ver docstring del módulo).</sub>")
    margin_t = geomatica.calcular_margen_top_titulo(titulo)

    fig.update_layout(
        title=titulo, xaxis_title="Año", yaxis_title="Hectáreas de agua superficial",
        hovermode="x unified",
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.75)"),
        plot_bgcolor="white", margin=dict(l=60, r=20, t=margin_t, b=50),
    )
    fig.update_xaxes(showgrid=True, gridcolor="#e5e5e5")
    fig.update_yaxes(showgrid=True, gridcolor="#e5e5e5", rangemode="tozero")

    html_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_agua_superficial_anual.html")
    fig.write_html(html_path)
    log(f"Gráfica de agua superficial anual: {html_path}")
    return html_path


# ==============================================================================
# --- DIAGNÓSTICO: ¿AGUA REAL, SOMBRA TOPOGRÁFICA O ROCA OSCURA? ---
# ==============================================================================
def _descargar_waterclass_alineado(ee, geom_wgs84_visual, anio, transform_ref, shape_ref, utm_crs,
                                    carpeta_tmp=None):
    """Descarga la banda 'waterClass' de JRC para UN año, recortada al
    polígono visual, y la realinea exacto a la misma malla (transform_ref/
    shape_ref) que ya usa el terreno 3D de geomatica.py -- mismo patrón que
    core.deforestacion._descargar_lossyear_alineado (mismo criterio:
    Resampling.nearest, NUNCA bilinear, porque waterClass es un CÓDIGO
    categórico -- 0=sin dato, 1=no agua, 2=estacional, 3=permanente --
    interpolar entre el código 1 y el código 3 daría un 'agua 2.5' sin
    significado)."""
    import tempfile

    carpeta_tmp = carpeta_tmp or tempfile.gettempdir()
    os.makedirs(carpeta_tmp, exist_ok=True)

    aoi = ee.Geometry(geom_wgs84_visual)
    img_anio = ee.ImageCollection(COLECCION_JRC).filterDate(f"{anio}-01-01", f"{anio + 1}-01-01").first()
    img = ee.Image(img_anio).select("waterClass").clip(aoi)
    url = img.getDownloadURL({"region": aoi, "scale": 30, "crs": "EPSG:4326", "format": "GEO_TIFF"})
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    tif_crudo = os.path.join(carpeta_tmp, f"temp_jrc_waterclass_{anio}_crudo.tif")
    with open(tif_crudo, "wb") as f:
        f.write(r.content)

    rows, cols = shape_ref
    waterclass_alineado = np.zeros((rows, cols), dtype=np.uint8)
    with rasterio.open(tif_crudo) as src:
        reproject(
            source=rasterio.band(src, 1), destination=waterclass_alineado,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=transform_ref, dst_crs=utm_crs,
            resampling=Resampling.nearest,
        )
    os.remove(tif_crudo)
    return waterclass_alineado


def _generar_mapa_diagnostico_2d(Z_raw, pendiente_deg, hillshade, mascara_nucleo, mascara_agua_nucleo,
                                  lats, lons, id_proyecto, anio_diagnostico, carpeta_salida):
    """Mapa 2D -- a propósito NO un mapa 3D: el punto aquí es cruzar tres
    capas alineadas pixel a pixel (hillshade, pendiente, clasificación JRC)
    sobre la MISMA vista, algo que un terreno 3D solo estorbaría (ver
    generar_grafica_agua_anual para el mismo razonamiento aplicado a la
    tendencia anual). Recortado al bounding box de la zona núcleo. Reusa
    geomatica.calcular_margen_top_titulo() para el margen del título."""
    import plotly.graph_objects as go
    from core import geomatica

    filas_nucleo, cols_nucleo = np.where(mascara_nucleo)
    if len(filas_nucleo) == 0:
        raise ValueError("La máscara de la zona núcleo está vacía -- no se puede recortar el mapa de diagnóstico.")
    r0, r1 = int(filas_nucleo.min()), int(filas_nucleo.max()) + 1
    c0, c1 = int(cols_nucleo.min()), int(cols_nucleo.max()) + 1

    hillshade_recorte = hillshade[r0:r1, c0:c1]
    mascara_nucleo_recorte = mascara_nucleo[r0:r1, c0:c1]
    hillshade_vis = np.where(mascara_nucleo_recorte, hillshade_recorte, np.nan)

    fig = go.Figure()
    # reversescale=True: la escala "Greys" de Plotly va blanco(0)->negro(1) por defecto: sin invertir,
    # un pixel MAS iluminado (hillshade alto) saldria MAS OSCURO en el mapa, al reves de la intuicion
    # de una carta de sombras real (iluminado=claro, en sombra=oscuro). Verificado con captura antes
    # de entregar -- sin este flag el mapa engañaba visualmente aunque el veredicto numerico ya era correcto.
    fig.add_trace(go.Heatmap(z=hillshade_vis, colorscale="Greys", reversescale=True,
                              showscale=False, zmin=0, zmax=1, hoverinfo="skip"))

    fa, ca = np.where(mascara_agua_nucleo)
    if len(fa):
        sel = (fa >= r0) & (fa < r1) & (ca >= c0) & (ca < c1)
        fa, ca = fa[sel], ca[sel]
        fig.add_trace(go.Scatter(
            x=ca - c0, y=fa - r0, mode="markers",
            marker=dict(symbol="x", size=7, color="#d62728", line=dict(width=1, color="#7f0000")),
            name=f"JRC agua {anio_diagnostico} (núcleo)",
            customdata=np.column_stack([
                Z_raw[fa, ca], pendiente_deg[fa, ca], hillshade[fa, ca], lats[fa, ca], lons[fa, ca],
            ]),
            hovertemplate=("Altitud: %{customdata[0]:.0f} msnm<br>Pendiente: %{customdata[1]:.1f}°"
                            "<br>Hillshade: %{customdata[2]:.2f}<br>Lat: %{customdata[3]:.5f}, "
                            "Lon: %{customdata[4]:.5f}<extra></extra>"),
        ))

    titulo = (f"{id_proyecto} -- Diagnóstico agua JRC {anio_diagnostico} en zona núcleo<br>"
              "<sub>Fondo=hillshade ilustrativo (NO la geometría solar real de Landsat) | "
              "X roja=píxel clasificado como agua por JRC ese año."
              "<br>Cruce de evidencia indirecta (pendiente + sombra) -- no reemplaza verificación de campo, "
              "ver docstring de generar_diagnostico_nucleo.</sub>")
    margin_t = geomatica.calcular_margen_top_titulo(titulo)

    fig.update_layout(
        title=titulo, xaxis_title="columna (píxel, 30m)", yaxis_title="fila (píxel, 30m)",
        yaxis=dict(autorange="reversed", scaleanchor="x", scaleratio=1),
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.75)"),
        margin=dict(l=60, r=20, t=margin_t, b=50), plot_bgcolor="white",
    )

    html_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_diagnostico_nucleo_{anio_diagnostico}.html")
    fig.write_html(html_path)
    log(f"Mapa de diagnóstico núcleo ({anio_diagnostico}): {html_path}")
    return html_path


def generar_diagnostico_nucleo(geojson_path, id_proyecto, anio_diagnostico=None, zonas_m=None,
                                carpeta_salida=None, carpeta_srtm=None, proyecto_gee=None,
                                sol_azimut_deg=150.0, sol_elevacion_deg=45.0):
    """Diagnóstico puntual para UN año: cruza los píxeles de la zona núcleo
    que JRC clasifica como agua (estacional o permanente) ese año contra
    (a) la pendiente del terreno (del mismo SRTM que ya usa geomatica.py,
    sin descargar nada nuevo) y (b) un hillshade ILUSTRATIVO -- para
    distinguir entre tres hipótesis frente a una señal de agua
    sospechosamente plana año con año en el núcleo:
        1. Agua real (pendiente baja, sin ser sistemáticamente más oscura
           que el resto del núcleo).
        2. Sombra topográfica (pendiente alta, sistemáticamente más
           oscura -- una sombra recurre en la misma geometría cada año y
           Landsat puede confundirla con agua).
        3. Roca volcánica oscura (andesita/basalto) mal clasificada -- OJO:
           esto NO se puede distinguir de sombra solo con pendiente+hillshade
           (dan la misma firma); cuando el veredicto descarta agua real,
           dice 'sombra o roca oscura', nunca afirma sombra con certeza.

    LIMITACIÓN HONESTA DEL HILLSHADE: se calcula con la fórmula estándar
    (Lambertian) a partir del MISMO SRTM del terreno y un azimut/elevación
    solar FIJOS (ajustables via sol_azimut_deg/sol_elevacion_deg) -- esto
    es una ilustración de "¿este píxel normalmente cae en sombra con luz
    de esa dirección?", NO la geometría solar real de la imagen Landsat
    que JRC usó para clasificar agua ese año específico (esa geometría no
    la expone el dataset). El veredicto de esta función es evidencia
    indirecta para decidir si vale la pena verificar en campo o con
    imagen de alta resolución -- no reemplaza esa verificación.

    Reusa el SRTM ya cacheado por core/geomatica.py (no se descarga dos
    veces si ya existe en carpeta_srtm). Descarga de JRC SOLO la banda
    waterClass del año pedido (no el historial completo), ya alineada
    píxel a píxel a la malla del terreno.

    Devuelve un dict: veredicto, n_pixeles_agua_nucleo, pct_pendiente_baja,
    pct_pendiente_alta, hillshade_agua_promedio, hillshade_nucleo_promedio,
    csv_path, html_path, anio_diagnostico."""
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
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    carpeta_srtm = carpeta_srtm or CARPETA_SRTM
    os.makedirs(carpeta_salida, exist_ok=True)
    anio_diagnostico = anio_diagnostico or JRC_ANIO_MAX
    if anio_diagnostico < JRC_ANIO_MIN or anio_diagnostico > JRC_ANIO_MAX:
        raise ValueError(f"anio_diagnostico={anio_diagnostico} fuera del rango de JRC "
                          f"({JRC_ANIO_MIN}-{JRC_ANIO_MAX}).")

    geom_utm_nucleo, dst_array, meta_utm, utm_crs = geomatica.cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        PERCENTIL_CAUCE_HIDROLOGIA, carpeta_srtm, id_proyecto,
    )

    Z_raw = hidrologia["Z_raw"]
    transform_ref = hidrologia["transform"]
    pw_v, ph_v = hidrologia["pw_v"], hidrologia["ph_v"]
    mascara_nucleo = hidrologia["zona_de_pixel"] == "nucleo"
    rows, cols = Z_raw.shape

    a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True).transform
    geom_visual_utm = geom_utm_nucleo.buffer(max(zonas_m))
    geom_visual_wgs84 = shp_transform(a_wgs84, geom_visual_utm)

    log(f"Descargando JRC waterClass {anio_diagnostico} alineado a la malla del terreno...")
    waterclass_alineado = _descargar_waterclass_alineado(
        ee, mapping(geom_visual_wgs84), anio_diagnostico, transform_ref, Z_raw.shape, utm_crs,
        carpeta_tmp=carpeta_srtm,
    )

    # Pendiente (grados) y hillshade ILUSTRATIVO -- ver limitación en el docstring de esta función.
    Z_para_gradiente = np.where(np.isnan(Z_raw), np.nanmean(Z_raw), Z_raw)
    gy, gx = np.gradient(Z_para_gradiente, ph_v, pw_v)
    pendiente_rad = np.arctan(np.sqrt(gx ** 2 + gy ** 2))
    pendiente_deg = np.degrees(pendiente_rad)
    aspecto_rad = np.arctan2(gy, -gx)
    zenith_rad = np.radians(90.0 - sol_elevacion_deg)
    azimut_rad = np.radians(sol_azimut_deg)
    hillshade = (np.cos(zenith_rad) * np.cos(pendiente_rad)
                 + np.sin(zenith_rad) * np.sin(pendiente_rad) * np.cos(azimut_rad - aspecto_rad))
    hillshade = np.clip(hillshade, 0.0, 1.0)

    mascara_agua_nucleo = mascara_nucleo & np.isin(waterclass_alineado, [2, 3])
    n_pixeles_agua_nucleo = int(mascara_agua_nucleo.sum())
    mascara_nucleo_valida = mascara_nucleo & ~np.isnan(Z_raw)

    lats, lons = geomatica.calcular_grid_latlon(transform_ref, utm_crs, rows, cols)

    pct_pendiente_baja = pct_pendiente_alta = hillshade_agua_promedio = None
    hillshade_nucleo_promedio = float(np.nanmean(hillshade[mascara_nucleo_valida])) \
        if mascara_nucleo_valida.any() else None

    if n_pixeles_agua_nucleo == 0:
        veredicto = (f"El año {anio_diagnostico}, JRC no clasifica NINGÚN píxel de la zona núcleo como agua "
                     "(ni estacional ni permanente) -- no hay nada que diagnosticar para este año en concreto; "
                     "prueba con otro año de la serie donde sí aparezca la señal plana.")
    else:
        pendiente_agua = pendiente_deg[mascara_agua_nucleo]
        hillshade_agua = hillshade[mascara_agua_nucleo]
        pct_pendiente_baja = float((pendiente_agua < 10.0).mean())
        pct_pendiente_alta = float((pendiente_agua > 25.0).mean())
        hillshade_agua_promedio = float(hillshade_agua.mean())
        diferencia_hillshade = (hillshade_nucleo_promedio - hillshade_agua_promedio
                                 if hillshade_nucleo_promedio is not None else 0.0)  # >0 = agua mas oscura

        log(f"  Núcleo, {anio_diagnostico}: {n_pixeles_agua_nucleo} píxeles JRC=agua -- "
            f"{pct_pendiente_baja*100:.0f}% en pendiente <10°, {pct_pendiente_alta*100:.0f}% en pendiente >25°, "
            f"hillshade agua={hillshade_agua_promedio:.2f} vs núcleo={hillshade_nucleo_promedio:.2f} "
            f"(diferencia {diferencia_hillshade:+.2f}).")

        if pct_pendiente_alta > 0.5 and diferencia_hillshade > 0.15:
            veredicto = (
                f"CONSISTENTE CON SOMBRA TOPOGRÁFICA (o roca oscura en pendiente pronunciada): "
                f"{pct_pendiente_alta*100:.0f}% de los píxeles 'agua' están en pendiente >25° y son en promedio "
                f"{diferencia_hillshade:.2f} más oscuros (hillshade) que el resto del núcleo. Esto NO distingue "
                "sombra de roca volcánica oscura -- ambas dan esta misma firma. Recomendación: NO presentar "
                "esta cifra al público sin verificación adicional (imagen de alta resolución o visita a la "
                "coordenada del CSV)."
            )
        elif pct_pendiente_baja > 0.5 and diferencia_hillshade < 0.05:
            veredicto = (
                f"CONSISTENTE CON AGUA REAL: {pct_pendiente_baja*100:.0f}% de los píxeles 'agua' están en "
                f"pendiente <10° (plausible para un cuerpo de agua) y no son sistemáticamente más oscuros que "
                f"el resto del núcleo (diferencia de hillshade {diferencia_hillshade:+.2f}). No es una "
                "confirmación de campo, pero no hay señal de sombra ni roca oscura en pendiente."
            )
        else:
            veredicto = (
                f"RESULTADO MIXTO/AMBIGUO: {pct_pendiente_baja*100:.0f}% en pendiente baja, "
                f"{pct_pendiente_alta*100:.0f}% en pendiente alta, diferencia de hillshade "
                f"{diferencia_hillshade:+.2f} -- no permite descartar ni confirmar sombra/roca oscura vs. agua "
                "real con este cruce. Revisar el CSV punto por punto (columnas lat/lon) en campo o en imagen "
                "de alta resolución antes de usar esta cifra."
            )
    log(veredicto, nivel="WARN" if ("SOMBRA" in veredicto or "MIXTO" in veredicto) else "OK")

    filas_csv = []
    fa, ca = np.where(mascara_agua_nucleo)
    for f, c in zip(fa, ca):
        filas_csv.append({
            "fila": int(f), "columna": int(c), "altitud_msnm": round(float(Z_raw[f, c]), 1),
            "pendiente_deg": round(float(pendiente_deg[f, c]), 1),
            "hillshade": round(float(hillshade[f, c]), 3),
            "jrc_waterclass": int(waterclass_alineado[f, c]),
            "lat": round(float(lats[f, c]), 6), "lon": round(float(lons[f, c]), 6),
        })
    df_csv = pd.DataFrame(filas_csv)
    csv_path = os.path.join(carpeta_salida, f"diagnostico_nucleo_{anio_diagnostico}_{id_proyecto.lower()}.csv")
    df_csv.to_csv(csv_path, index=False)
    log(f"CSV punto por punto (para revisar en Google Earth): {csv_path}")

    html_path = _generar_mapa_diagnostico_2d(
        Z_raw, pendiente_deg, hillshade, mascara_nucleo, mascara_agua_nucleo, lats, lons,
        id_proyecto, anio_diagnostico, carpeta_salida,
    )

    return {
        "veredicto": veredicto, "n_pixeles_agua_nucleo": n_pixeles_agua_nucleo,
        "pct_pendiente_baja": pct_pendiente_baja, "pct_pendiente_alta": pct_pendiente_alta,
        "hillshade_agua_promedio": hillshade_agua_promedio, "hillshade_nucleo_promedio": hillshade_nucleo_promedio,
        "csv_path": csv_path, "html_path": html_path, "anio_diagnostico": anio_diagnostico,
    }


# ==============================================================================
# --- MODO DEMO: sin Earth Engine, valores sintéticos deterministas ---
# ==============================================================================
def demo():
    """Prueba TODA la lógica que no depende de Earth Engine: la resta a
    anillo exclusivo (generar_resumen_agua_sin_traslape) y la gráfica 2D
    completa (generar_grafica_agua_anual), con un df_historial sintético
    (3 zonas, 6 años, valores acumulativos monótonos por diseño -- igual
    que las zonas reales de geomatica.py) -- así se prueba sin tocar Earth
    Engine ni descargar nada."""
    log("=== core.agua_superficial --demo (sin Earth Engine, valores sintéticos) ===")

    rng = np.random.default_rng(7)
    anios = list(range(2018, 2024))
    zonas = [(0, "nucleo"), (500, "buffer_500m"), (1000, "buffer_1000m")]

    filas = []
    for anio in anios:
        base = max(0.0, 2.0 + 0.3 * rng.normal())
        acumulado = 0.0
        for buf_m, etiqueta in zonas:
            acumulado += abs(rng.normal(1.0, 0.4)) if buf_m > 0 else base
            filas.append({
                "zona": etiqueta, "buffer_m": buf_m, "anio": anio,
                "ha_agua_estacional": round(acumulado * 0.4, 4),
                "ha_agua_permanente": round(acumulado * 0.6, 4),
                "ha_agua_total": round(acumulado, 4),
            })
    df_historial = pd.DataFrame(filas)
    print("\n--- Historial sintético (acumulativo por zona, como en el pipeline real) ---")
    print(df_historial.to_string(index=False))

    id_proyecto = "DEMO_AGUA"
    carpeta_tmp = os.path.expanduser("~/resultados_demo_agua_superficial")
    os.makedirs(carpeta_tmp, exist_ok=True)
    csv_sin_traslape = generar_resumen_agua_sin_traslape(df_historial, id_proyecto, carpeta_tmp)

    df_sin_traslape = pd.read_csv(csv_sin_traslape)
    print("\n--- Resumen SIN traslape (anillo exclusivo, sumable, por año) ---")
    print(df_sin_traslape.to_string(index=False))

    # Verificación de que el anillo exclusivo SÍ suma al total real de cada año (sin muñeca rusa):
    # el TOTAL de cada año debe ser exactamente el valor bruto de buffer_1000m (el buffer más
    # grande) de ese mismo año -- porque la resta en cadena es una suma telescópica.
    for anio in anios:
        total_calculado = df_sin_traslape[
            (df_sin_traslape["anio"] == anio) & (df_sin_traslape["anillo"].str.startswith("TOTAL"))
        ]["ha_agua"].iloc[0]
        total_esperado = df_historial[(df_historial["anio"] == anio) & (df_historial["buffer_m"] == 1000)][
            "ha_agua_total"
        ].iloc[0]
        assert abs(total_calculado - total_esperado) < 1e-6, (
            f"Año {anio}: TOTAL sin traslape ({total_calculado}) no coincide con buffer_1000m bruto "
            f"({total_esperado}) -- la resta en cadena debería ser exacta."
        )
    print(f"\nOK: el TOTAL sin traslape de cada año coincide exacto con buffer_1000m bruto de ese año "
          f"({len(anios)} años verificados) -- sin muñeca rusa.")

    try:
        html_path = generar_grafica_agua_anual(id_proyecto, carpeta_tmp)
        log(f"Gráfica demo generada en: {html_path}")
    except ImportError as e:
        log(f"plotly no instalado, se omite la gráfica del demo: {e}", nivel="WARN")

    return df_sin_traslape


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Agua superficial visible por año (JRC Global Surface Water) por zona -- Motor Nacional"
    )
    ap.add_argument("--demo", action="store_true", help="Corre con valores sintéticos, sin Earth Engine ni red")
    ap.add_argument("--geojson", type=str, help="Ruta al GeoJSON del polígono núcleo")
    ap.add_argument("--id-proyecto", type=str, help="Nombre identificador (para nombres de archivo)")
    ap.add_argument("--zonas", type=str, default=None, help="Buffers en metros separados por coma, ej. '0,500,1000'")
    ap.add_argument("--anio-inicio", type=int, default=None,
                     help=f"default: config del sitio o {JRC_ANIO_MIN} (JRC no cubre antes)")
    ap.add_argument("--anio-fin", type=int, default=None, help=f"default: {JRC_ANIO_MAX} (última cobertura de JRC v1.4)")
    ap.add_argument("--proyecto-gee", type=str, default=None, help="ID de proyecto de Google Cloud para ee.Initialize(project=...)")
    ap.add_argument("--carpeta-salida", type=str, default=None)
    ap.add_argument("--grafica", action="store_true",
                     help="Solo regenera la gráfica 2D desde el CSV sin traslape YA existente -- no consulta "
                          "Earth Engine otra vez (útil para ajustar el estilo sin gastar cuota).")
    ap.add_argument("--diagnostico-nucleo", action="store_true",
                     help="Corre el diagnóstico sombra/agua/roca oscura para UN año en la zona núcleo "
                          "(ver docstring de generar_diagnostico_nucleo) en vez del historial completo.")
    ap.add_argument("--anio-diagnostico", type=int, default=None,
                     help=f"Año a diagnosticar con --diagnostico-nucleo (default: {JRC_ANIO_MAX}, el más "
                          "reciente disponible en JRC).")
    ap.add_argument("--sol-azimut", type=float, default=150.0,
                     help="Azimut solar en grados (0=norte, 90=este) para el hillshade ILUSTRATIVO del "
                          "diagnóstico (default: 150 -- ver limitación en el docstring, no es la geometría "
                          "solar real de Landsat).")
    ap.add_argument("--sol-elevacion", type=float, default=45.0,
                     help="Elevación solar en grados sobre el horizonte para el hillshade ILUSTRATIVO "
                          "(default: 45).")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not args.id_proyecto:
        ap.error("--id-proyecto es obligatorio fuera de --demo")

    carpeta_salida = args.carpeta_salida or os.path.expanduser(f"~/resultados_{args.id_proyecto.lower()}")

    if args.grafica:
        generar_grafica_agua_anual(args.id_proyecto, carpeta_salida=carpeta_salida)
        return

    if not args.geojson:
        ap.error("--geojson es obligatorio fuera de --demo/--grafica")

    zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None

    if args.diagnostico_nucleo:
        resultado = generar_diagnostico_nucleo(
            geojson_path=args.geojson, id_proyecto=args.id_proyecto, anio_diagnostico=args.anio_diagnostico,
            zonas_m=zonas_m, carpeta_salida=carpeta_salida, proyecto_gee=args.proyecto_gee,
            sol_azimut_deg=args.sol_azimut, sol_elevacion_deg=args.sol_elevacion,
        )
        print(f"\n=== VEREDICTO ({resultado['anio_diagnostico']}) ===\n{resultado['veredicto']}\n")
        return

    procesar_sitio_real(
        geojson_path=args.geojson, id_proyecto=args.id_proyecto, zonas_m=zonas_m,
        anio_inicio=args.anio_inicio, anio_fin=args.anio_fin, carpeta_salida=carpeta_salida,
        proyecto_gee=args.proyecto_gee,
    )
    generar_grafica_agua_anual(args.id_proyecto, carpeta_salida=carpeta_salida)


if __name__ == "__main__":
    main()

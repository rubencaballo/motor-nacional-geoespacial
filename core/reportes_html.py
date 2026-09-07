#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tarjetas HTML + página envolvente para los mapas 3D de esta plataforma.

POR QUÉ EXISTE ESTE MÓDULO:
    core/carbono.py (CO2e ALMACENADO por zona) y core/carbono_perdida.py
    (CO2e LIBERADO por zona) necesitan mostrar exactamente el mismo tipo
    de tarjeta -- número grande y legible por zona, con el color del
    anillo como acento, fuera de la escena 3D -- porque meter esos
    números como texto flotando en el modelo 3D se probó y falló: se
    cortan con los ejes, se tapan con el tooltip nativo de Plotly, o
    quedan ilegibles según el ángulo de cámara ("aquí se pierden los
    números", feedback real del usuario). En vez de que cada módulo
    reimplemente su propio HTML/CSS de tarjetas (y que un ajuste de
    diseño haya que repetirlo en dos archivos, arriesgando que se
    desincronicen), este módulo es el ÚNICO lugar donde vive ese diseño.

    Pensado para reusarse tal cual en cualquier ANP futura -- no depende
    de Cofre de Perote ni de ningún dato específico, solo RENDERIZA
    valores que ya le pasan calculados. Si un número sale mal, el bug
    está en quien llamó a esta función con el dato equivocado, no aquí.

CONVENCIÓN DE COLOR (igual en todo el proyecto, desde geomatica.py):
    núcleo=rojo, buffer_500m=naranja, buffer_1000m=dorado, total=verde."""

import pandas as pd

NOMBRES_LEGIBLES_ZONA = {"nucleo": "Núcleo", "buffer_500m": "Buffer 500 m", "buffer_1000m": "Buffer 1000 m"}
COLORES_ZONA_HEX = {"nucleo": "#d62728", "buffer_500m": "#e8892b", "buffer_1000m": "#c9a227"}
COLOR_TOTAL_HEX = "#2b6b3f"


def tabla_desglose_anual_deforestacion_html(df_desglose, anio_inicio, anio_fin):
    """Tablita año por año (2026-09-06, pedida por el usuario viendo el
    mapa 3D de deforestación: "que se coloque a un costado, los años...
    y ponerle al ladito la cantidad de hectáreas deforestadas por año").

    Pensada para ir AL COSTADO del mapa 3D (ver generar_mapa_3d_deforestacion
    en core/deforestacion.py) -- a propósito NO se intentó meter esto en
    las etiquetas de la barra de color de Plotly: esa barra ya venía con
    ticks cada 2 años (nunca los 25 años completos) precisamente para no
    saturarla de texto ilegible -- 25 etiquetas ahí se hubieran encimado.
    Una tabla aparte, con scroll, es el mismo patrón que ya se usó para la
    tabla de especies en riesgo (ver tabla_especies_riesgo_html) cuando esa
    misma clase de problema -- demasiadas filas para caber legible sobre
    la imagen -- ya se resolvió una vez en este proyecto.

    Dos columnas de datos, porque el usuario pidió las dos cosas a la vez:
    "para que se vea en qué año se deforestó más" (barra + ha de ESE año)
    y "que se sume las hectáreas... por año" (acumulado corriendo, desde
    anio_inicio hasta ese año). La barra horizontal de cada fila es
    relativa al año con MÁS pérdida del rango completo (no a un máximo
    fijo), así que el año peor SIEMPRE se ve con la barra al 100%.

    `df_desglose`: el DataFrame que devuelve
    core.deforestacion.generar_desglose_anual_visual() (columnas
    anio/anillo/perdida_ha_visual) -- se usa tal cual, solo se filtran las
    filas anillo == 'TOTAL (area visual, sin traslape)' para no sumar cada
    anillo por separado aquí (ya vienen sumados en esa fila)."""
    valores_por_anio = {}
    for anio in range(anio_inicio, anio_fin + 1):
        sub = df_desglose[(df_desglose["anio"] == anio) & (df_desglose["anillo"] == "TOTAL (area visual, sin traslape)")]
        valores_por_anio[anio] = float(sub["perdida_ha_visual"].iloc[0]) if len(sub) else 0.0
    max_ha = max(valores_por_anio.values()) if valores_por_anio else 0.0

    filas_html = []
    acumulado = 0.0
    for anio in range(anio_inicio, anio_fin + 1):
        ha = valores_por_anio[anio]
        acumulado += ha
        pct_barra = (ha / max_ha * 100) if max_ha > 0 else 0.0
        filas_html.append(
            f'<tr><td>{anio}</td>'
            f'<td>{ha:,.2f} ha'
            f'<div style="background:#d62728;height:5px;border-radius:3px;'
            f'width:{pct_barra:.0f}%;margin-top:3px"></div></td>'
            f'<td>{acumulado:,.2f} ha</td></tr>'
        )
    # 2026-09-06 (pedido por Ruben): un renglón de GRAN TOTAL separado y
    # resaltado al fondo de la tabla -- antes el total solo se podía leer
    # como el último valor de la columna "Acumulado", fácil de perder de
    # vista en una tabla con scroll de hasta 25 años. Es el mismo número
    # (el acumulado ya llegó al último año del ciclo de arriba), solo que
    # ahora con su propia fila fija, en negritas y con fondo distinto, para
    # que no dependa de hacer scroll hasta el final para verlo.
    filas_html.append(
        f'<tr style="font-weight:700;background:#f4f4f4;border-top:2px solid #999">'
        f'<td colspan="2">GRAN TOTAL {anio_inicio}-{anio_fin}</td>'
        f'<td>{acumulado:,.2f} ha</td></tr>'
    )
    return f'''<div class="tabla-riesgo-wrap" style="max-width:380px">
      <table class="tabla-riesgo">
        <thead><tr><th>Año</th><th>Perdido ese año</th><th>Acumulado desde {anio_inicio}</th></tr></thead>
        <tbody>{"".join(filas_html)}</tbody>
      </table>
    </div>'''


def script_click_deforestacion_timelapse_html(div_id, prefijo_nombre_capa=None):
    """<script> compartido: al picar un punto dentro del mapa Plotly con id
    `div_id`, abre Google Earth Web centrado en ese lat/lon, en pestaña
    nueva.

    CORREGIDO 2026-09-06 (auditoría de git log): esto YA existía en el
    proyecto desde antes -- ver `CLICK_ABRE_GOOGLE_EARTH_JS` en
    core/geomatica.py, usado ahí vía `fig.write_html(..., post_script=...)`
    para los mapas 3D que se escriben directo a disco (hidrología, carbono,
    incendios). Esa versión pasó por varias iteraciones reales:
      - 636e364: primera versión, formato `/web/@lat,lon,...`
      - 0427aa6: ajustó la distancia de cámara (600d -> 150d)
      - 6ed8aad: descubrió que la pose de cámara (`/web/@lat,lon,alt,dist,
        tilt,heading,tilt,roll`) hacía que Google Earth "reacomodara a su
        gusto" la vista (pantalla negra / no aterriza exacto) y la cambió
        a `/web/search/lat,lon` -- EL MISMO mecanismo que teclear la
        coordenada a mano en el buscador de Google Earth, que es el que sí
        aterriza confiable.
    La primera versión de ESTA función (escrita hoy, antes de encontrar lo
    anterior) reinventó el formato viejo ya descartado (`/web/@...`) y
    luego pasó por Earth Engine Timelapse -- ninguno de los dos era
    necesario. Esta versión ya usa el formato correcto y ya validado, y
    generaliza el mismo criterio de geomatica.py: reacciona a CUALQUIER
    trazo cuyo customdata traiga lat/lon en sus ÚLTIMAS dos posiciones (no
    solo el de deforestación) -- así el mapa de biodiversidad se comporta
    igual que el resto de la plataforma (hidrología/carbono/incendios ya
    hacen esto en TODOS sus puntos, no solo en uno). Si se quiere
    restringir a un solo tipo de trazo, pasar `prefijo_nombre_capa` (ej.
    "Deforestación") -- por default (None) es genérico, igual que
    geomatica.CLICK_ABRE_GOOGLE_EARTH_JS.

    LIMITACIÓN CONOCIDA (heredada, no resuelta): esta URL abre Google
    Earth Web centrado en el punto, pero no activa el modo "Imágenes
    históricas" (línea de tiempo) -- eso requeriría Earth Engine Timelapse,
    que el usuario ya probó y prefirió no usar."""
    filtro_nombre = (
        f'if (pt.data.name.indexOf("{prefijo_nombre_capa}") !== 0) return;'
        if prefijo_nombre_capa else ""
    )
    return f'''<script>
    (function() {{
      var _div = document.getElementById("{div_id}");
      if (!_div) return;
      _div.on("plotly_click", function(datos) {{
        var pt = datos.points && datos.points[0];
        if (!pt || !pt.data) return;
        {filtro_nombre}
        var cd = pt.customdata;
        if (!cd || cd.length < 2) return;
        var lat = cd[cd.length - 2], lon = cd[cd.length - 1];
        if (typeof lat !== "number" || typeof lon !== "number") return;
        // MISMO formato validado que geomatica.CLICK_ABRE_GOOGLE_EARTH_JS y
        // _liga_coordenadas_html() en este archivo -- si se cambia uno, cambiar los tres.
        var url = "https://earth.google.com/web/search/" + lat.toFixed(6) + "," + lon.toFixed(6);
        window.open(url, "_blank");
      }});
    }})();
    </script>'''


def linea_secundaria_html(texto_negritas, texto_normal=""):
    """UNA línea secundaria dentro de una tarjeta (ej. el dato de GEDI, o
    un desglose por causa) -- mismo estilo de tipografía en toda la
    plataforma: el dato en negritas, la fuente/nota en gris al lado."""
    extra = f' <span class="fuente">{texto_normal}</span>' if texto_normal else ""
    return f'<div class="secundaria"><b>{texto_negritas}</b>{extra}</div>'


def tabla_especies_riesgo_html(df_riesgo, caption=None):
    """Tabla de coordenadas de especies en riesgo, numerada IGUAL que los
    marcadores del mapa satelital 2D (columna `numero`) -- pensada para ir
    DEBAJO de esa imagen, con scroll interno, en vez de encimada sobre ella.

    Por qué existe: hasta el 2026-09-04, esta tabla se dibujaba directo
    sobre la imagen satelital (matplotlib fig.text). Funcionaba con un
    puñado de especies en riesgo, pero en cuanto el cruce NOM-059 empezó a
    encontrar de verdad (core/gbif_biodiversidad.py) el conteo subió a
    70-80 filas en un sitio real y la tabla tapaba buena parte del mapa
    -- feedback directo del usuario. Sacarla a HTML normal (con scroll)
    resuelve esto para cualquier cantidad de filas, no solo para el caso
    que se probó ese día.

    `df_riesgo`: DataFrame con columnas numero/NOMBRE_CIENTIFICO/
    zona_legible/NIVEL_RIESGO/COLOR_RIESGO/LAT/LON (mismo df que ya arma
    core.biodiversidad._mapa_satelital_riesgo() para dibujar los puntos,
    reusado tal cual para que los números de mapa y tabla NUNCA se
    desincronicen). Columnas opcionales, cada una agrega su propia
    columna en la tabla SOLO si está presente (mismo criterio ya usado
    para ENDEMICO_MX -- una tabla vieja sin esa columna sigue funcionando
    igual, nunca truena por columna faltante):

    - ENDEMICO_MX: pastilla de color -- verde "SI" (confirmado a mano),
      azul "PROBABLE_CALCULADO" (indicio calculado vía GBIF, NO es lo
      mismo que confirmado -- ver NOTA_ENDEMISMO en el CSV completo),
      gris "NO_EVALUADO" (sin dato, nunca "confirmado que no es
      endémica").
    - GBIF_ID (2026-09-06): liga "Ver" al registro público de ESE
      avistamiento exacto en gbif.org (a veces trae foto) -- para que
      alguien en la conferencia pueda ver la especie con un clic, no solo
      leer su nombre en la tabla. Si el ID no es el numérico real de GBIF
      (ej. "demo-..." en --demo), no se dibuja liga rota -- se deja "s/d".
    - CERCA_DEFORESTACION/ANIO_DEFORESTACION_CERCANA (2026-09-06, ver
      core.biodiversidad._marcar_deforestacion_en_puntos): pastilla roja
      si el punto cae en un píxel Hansen (30m) marcado con pérdida de
      bosque -- NUNCA se afirma que ESA fue la causa de que la especie
      esté en riesgo, solo que geográficamente coincide con una mancha de
      deforestación detectada; ver nota al pie del reporte para el
      detalle de esta limitación.

    Lat/Lon siempre se vuelve una liga a Google Maps (vista satelital,
    universal, carga rápido) con una segunda liga chica a Google Earth
    Web al lado (más vistoso pero a veces tarda/pide permisos del
    navegador) -- clic en cualquiera de las dos lleva al punto exacto."""
    tiene_endemismo = "ENDEMICO_MX" in df_riesgo.columns
    tiene_gbif_id = "GBIF_ID" in df_riesgo.columns
    tiene_deforestacion = "CERCA_DEFORESTACION" in df_riesgo.columns
    filas_html = []
    for _, r in df_riesgo.sort_values("numero").iterrows():
        col_endemismo = f'<td>{_pastilla_endemismo_html(r.get("ENDEMICO_MX", ""))}</td>' if tiene_endemismo else ""
        col_deforestacion = (
            f'<td>{_pastilla_deforestacion_html(r.get("CERCA_DEFORESTACION"), r.get("ANIO_DEFORESTACION_CERCANA"))}</td>'
            if tiene_deforestacion else ""
        )
        col_ver = f'<td>{_liga_ver_especie_html(r.get("GBIF_ID"))}</td>' if tiene_gbif_id else ""
        filas_html.append(
            f'<tr><td>{r["numero"]}</td>'
            f'<td>{r["NOMBRE_CIENTIFICO"]}</td>'
            f'<td>{r.get("zona_legible", "")}</td>'
            f'<td><span class="pill-riesgo" style="background:{r["COLOR_RIESGO"]}">{r["NIVEL_RIESGO"]}</span></td>'
            f'{col_endemismo}'
            f'{col_deforestacion}'
            f'<td>{_liga_coordenadas_html(r["LAT"], r["LON"])}</td>'
            f'{col_ver}</tr>'
        )
    caption_html = f'<p class="tabla-riesgo-caption">{caption}</p>' if caption else ""
    th_endemismo = "<th>Endémica</th>" if tiene_endemismo else ""
    th_deforestacion = "<th>Deforestación</th>" if tiene_deforestacion else ""
    th_ver = "<th>Ver</th>" if tiene_gbif_id else ""
    return f'''{caption_html}<div class="tabla-riesgo-wrap">
      <table class="tabla-riesgo">
        <thead><tr><th>#</th><th>Especie</th><th>Zona</th><th>Nivel</th>{th_endemismo}{th_deforestacion}
        <th>Lat, Lon</th>{th_ver}</tr></thead>
        <tbody>{"".join(filas_html)}</tbody>
      </table>
    </div>'''


def _liga_coordenadas_html(lat, lon):
    """Lat/Lon como texto (5 decimales, igual que antes -- para quien las
    quiera copiar a mano, ej. a un GPS de campo) MÁS un par de ligas
    cortas: 'Mapa' abre Google Maps en vista satelital centrada en el
    punto exacto (universal, carga rápido, no depende de WebGL); 'Earth'
    abre Google Earth Web en el mismo punto. Se dan las dos para no
    depender de que una falle.

    CORREGIDO 2026-09-06: la liga 'Earth' traía el formato de POSE DE
    CÁMARA (`/web/@lat,lon,alt,dist,tilt,heading,tilt,roll`) -- ya se
    había descubierto y corregido este mismo problema antes en este
    proyecto (ver commit 6ed8aad, `CLICK_ABRE_GOOGLE_EARTH_JS` en
    core/geomatica.py): esa pose hace que Google Earth "reacomode a su
    gusto" la vista, a veces dejando la pantalla en negro. El formato
    `/web/search/lat,lon` (mismo mecanismo que teclear la coordenada a
    mano en el buscador de Earth) es el que en la práctica aterriza
    confiable. Si se ajusta aquí, ajustar también
    script_click_deforestacion_timelapse_html() en este mismo archivo y
    CLICK_ABRE_GOOGLE_EARTH_JS en core/geomatica.py -- los tres deben
    quedar iguales."""
    try:
        lat_f, lon_f = float(lat), float(lon)
    except (TypeError, ValueError):
        return ""
    texto = f"{lat_f:.5f}, {lon_f:.5f}"
    liga_maps = f"https://www.google.com/maps/place/{lat_f:.6f},{lon_f:.6f}/@{lat_f:.6f},{lon_f:.6f},700m/data=!3m1!1e3"
    liga_earth = f"https://earth.google.com/web/search/{lat_f:.6f},{lon_f:.6f}"
    return (f'{texto}<br>'
            f'<a href="{liga_maps}" target="_blank" rel="noopener" class="liga-coord">Mapa</a> · '
            f'<a href="{liga_earth}" target="_blank" rel="noopener" class="liga-coord">Earth</a>')


def _liga_ver_especie_html(gbif_id):
    """Liga al registro público de ESE avistamiento en gbif.org (a veces
    trae foto subida por quien lo observó) -- para que alguien viendo el
    reporte pueda ver la especie con un clic. `gbif_id` de un --demo (ej.
    'demo-Ambystoma velasci-12-8', no numérico) NO genera una liga rota:
    GBIF_ID real de la API siempre es un entero, así que se valida antes
    de armar la liga."""
    gbif_id_str = str(gbif_id or "").strip()
    if not gbif_id_str or gbif_id_str.lower() == "nan" or not gbif_id_str.replace(".", "", 1).isdigit():
        return '<span style="color:#999">s/d</span>'
    gbif_id_str = gbif_id_str.split(".")[0]  # por si pandas lo trajo como float (123.0)
    return (f'<a href="https://www.gbif.org/occurrence/{gbif_id_str}" target="_blank" rel="noopener" '
            f'class="liga-coord">GBIF ↗</a>')


_COLORES_PASTILLA_ENDEMISMO = {
    "SI": ("#1a7a3c", "Endémica (confirmada)"),
    "PROBABLE_CALCULADO": ("#2a5db0", "Posible (calculado GBIF, sin confirmar)"),
    "NO_EVALUADO": ("#8a8a8a", "Sin dato"),
}


def _pastilla_endemismo_html(estatus):
    """Pastilla de color para la columna 'Endémica' de la tabla de
    coordenadas -- ver _COLORES_PASTILLA_ENDEMISMO arriba para el
    criterio de color por estatus. A PROPÓSITO usa colores distintos
    para "SI" (verde, confirmado) y "PROBABLE_CALCULADO" (azul, indicio
    calculado) -- nunca deben verse igual en la tabla."""
    color, etiqueta = _COLORES_PASTILLA_ENDEMISMO.get(str(estatus or "").strip().upper(), ("#8a8a8a", "Sin dato"))
    return f'<span class="pill-riesgo" style="background:{color}" title="{etiqueta}">{etiqueta}</span>'


def _pastilla_deforestacion_html(cerca_deforestacion, anio):
    """Pastilla roja si el punto de este registro cae en un píxel Hansen
    (30m de resolución) marcado con pérdida de bosque -- ver docstring de
    tabla_especies_riesgo_html() arriba para la limitación importante: es
    una COINCIDENCIA GEOGRÁFICA a nivel de píxel, nunca una afirmación de
    causa ("esta especie está en riesgo POR la deforestación"). Verde
    "Sin pérdida" si se evaluó y no cayó en ningún píxel con pérdida.
    Gris "No evaluado" si esta corrida no cruzó contra deforestación
    (--con-deforestacion no se usó, o el punto cayó fuera del área visual
    del terreno 3D)."""
    if cerca_deforestacion is True:
        anio_txt = f" ({int(anio)})" if anio is not None and pd.notna(anio) else ""
        return (f'<span class="pill-riesgo" style="background:#8b0000" '
                f'title="Píxel Hansen con pérdida de bosque">Pérdida de bosque{anio_txt}</span>')
    if cerca_deforestacion is False:
        return '<span class="pill-riesgo" style="background:#2ca02c" title="Sin pérdida Hansen en este píxel">Sin pérdida</span>'
    return '<span class="pill-riesgo" style="background:#8a8a8a" title="No evaluado contra deforestación">No evaluado</span>'


def tarjeta_html(zona, color, valor_principal_fmt, unidad_principal, nota_principal=None, lineas_secundarias=None,
                  nombre_mostrado=None, area_ha=None, es_total=False):
    """UNA tarjeta -- valor grande en negritas (nunca coloreado: el color
    solo vive en el borde/acento, así el texto siempre lleva tinta neutra
    y el color lleva la identidad de la zona). Hectáreas (si se dan) van
    chico arriba del valor principal -- responde "cuánto mide" antes de
    "cuánto guarda/libera". `nota_principal` va chica y gris justo bajo el
    valor (incertidumbre, o cualquier nota corta); `lineas_secundarias` es
    una lista de HTML ya armado (ver linea_secundaria_html) para datos
    adicionales (GEDI, desglose por causa, etc.) -- puede ir vacía.

    `valor_principal_fmt`: el número YA formateado como string (ej.
    f"{x:,.0f}") -- esta función no decide redondeo ni separador de miles,
    eso es decisión de quien arma los datos de la tarjeta."""
    clase = "tarjeta tarjeta-total" if es_total else "tarjeta"
    linea_area = f'<div class="area">{area_ha:,.1f} ha</div>' if area_ha is not None and pd.notna(area_ha) else ""
    nota_html = f'<div class="nota-principal">{nota_principal}</div>' if nota_principal else ""
    secundarias_html = "".join(lineas_secundarias or [])
    return f'''<div class="{clase}" style="--color:{color}">
      <div class="zona">{nombre_mostrado or NOMBRES_LEGIBLES_ZONA.get(zona, zona)}</div>
      {linea_area}
      <div class="valor">{valor_principal_fmt}<span class="unidad">{unidad_principal}</span></div>
      {nota_html}
      {secundarias_html}
    </div>'''


CSS_TARJETAS = '''
  * { box-sizing: border-box; }
  body { margin:0; padding:22px 28px 36px; background:#f6f6f4; color:#1c1c1c;
         font-family: -apple-system, "Segoe UI", Roboto, Arial, sans-serif; }
  h1 { font-size:19px; font-weight:600; margin:0 0 4px; }
  .subtitulo { font-size:13px; color:#666; margin:0 0 18px; }
  .tarjetas { display:flex; gap:12px; flex-wrap:wrap; margin-bottom:20px; }
  .tarjeta { flex:1 1 210px; background:#fff; border-radius:10px; border-left:5px solid var(--color);
             box-shadow:0 1px 3px rgba(0,0,0,.08); padding:14px 16px; }
  .tarjeta-total { background:#eef6f0; }
  .tarjeta .zona { font-size:11.5px; font-weight:600; text-transform:uppercase; letter-spacing:.03em;
                    color:#666; margin-bottom:4px; }
  .tarjeta .area { font-size:14px; font-weight:500; color:#444; margin-bottom:8px; }
  .tarjeta .valor { font-size:26px; font-weight:600; line-height:1.15; }
  .tarjeta .valor .unidad { font-size:13px; font-weight:500; color:#666; margin-left:4px; }
  .tarjeta .nota-principal { font-size:11.5px; color:#888; margin:2px 0 8px; }
  .tarjeta .secundaria { font-size:12px; color:#444; border-top:1px solid #eee; padding-top:7px; margin-top:2px; }
  .tarjeta .secundaria + .secundaria { border-top:none; padding-top:2px; margin-top:0; }
  .tarjeta .secundaria .fuente { color:#888; font-weight:400; }
  .mapa { background:#fff; border-radius:10px; box-shadow:0 1px 3px rgba(0,0,0,.08); padding:4px; }
  .nota { font-size:11px; color:#999; margin-top:8px; max-width:900px; }
  .tabla-riesgo-caption { font-size:11.5px; color:#888; margin:10px 0 4px; }
  .tabla-riesgo-wrap { max-height:420px; overflow-y:auto; border:1px solid #e6e6e3; border-radius:8px;
                        background:#fff; box-shadow:0 1px 3px rgba(0,0,0,.06); }
  table.tabla-riesgo { width:100%; border-collapse:collapse; font-size:12.5px; }
  table.tabla-riesgo thead th { position:sticky; top:0; background:#f6f6f4; text-align:left;
                                 padding:8px 10px; border-bottom:1px solid #ddd; font-weight:600;
                                 color:#444; z-index:1; }
  table.tabla-riesgo td { padding:6px 10px; border-bottom:1px solid #f2f2f0; color:#333; }
  table.tabla-riesgo tbody tr:nth-child(even) { background:#fafafa; }
  table.tabla-riesgo tbody tr:hover { background:#f0f3ff; }
  .pill-riesgo { display:inline-block; padding:2px 9px; border-radius:10px; color:#fff; font-size:11px;
                 font-weight:600; letter-spacing:.02em; }
  a.liga-coord { color:#2a5db0; text-decoration:none; font-size:11.5px; font-weight:600; }
  a.liga-coord:hover { text-decoration:underline; }
'''


def pagina_html_con_tarjetas(titulo_pagina, h1, subtitulo, tarjetas_html, div_mapa, nota_pie):
    """Página completa: encabezado + fila de tarjetas + mapa 3D (ya
    convertido a HTML embebible con fig.to_html(full_html=False)) + nota
    al pie. Si algún día se quiere cambiar el look de las tarjetas para
    TODOS los mapas 3D de la plataforma (los de hoy y los de próximas
    ANP), este es el único archivo donde hay que tocarlo."""
    return f'''<!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<title>{titulo_pagina}</title>
<style>{CSS_TARJETAS}</style>
</head>
<body>
  <h1>{h1}</h1>
  <p class="subtitulo">{subtitulo}</p>
  <div class="tarjetas">
    {tarjetas_html}
  </div>
  <div class="mapa">{div_mapa}</div>
  <p class="nota">{nota_pie}</p>
</body>
</html>'''

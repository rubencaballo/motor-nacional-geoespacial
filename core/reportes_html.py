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


def linea_secundaria_html(texto_negritas, texto_normal=""):
    """UNA línea secundaria dentro de una tarjeta (ej. el dato de GEDI, o
    un desglose por causa) -- mismo estilo de tipografía en toda la
    plataforma: el dato en negritas, la fuente/nota en gris al lado."""
    extra = f' <span class="fuente">{texto_normal}</span>' if texto_normal else ""
    return f'<div class="secundaria"><b>{texto_negritas}</b>{extra}</div>'


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

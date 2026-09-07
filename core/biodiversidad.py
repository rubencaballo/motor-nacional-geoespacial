#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
core/biodiversidad.py
=======================

Módulo público de biodiversidad de la plataforma Salamandra -- distinto
de `core/V19_4_HIBRIDO_CORREGIDO.py` (que arma un dictamen técnico para
CONANP/SEDEMA). Este es para el módulo que va a ver la comunidad/la
conferencia: qué especies hay por zona, cuáles están en algún nivel de
riesgo (NOM-059/IUCN), un mapa 3D de "zonas a vigilar" sobre el mismo
terreno que ya usa el módulo de deforestación, hectáreas reales por zona
(2D plano vs 3D corregido por relieve), y un mapa satelital 2D (Esri World
Imagery vía contextily) con esos mismos puntos de riesgo.

DE DÓNDE VIENE CADA PIEZA (para no reinventar nada, ver RUNBOOK.md):
    - Descarga GBIF por zona (núcleo/buffer_500m/buffer_1000m, anillo
      EXCLUSIVO) y validación taxonómica CLASE vs FAMILIA: ya existían y
      funcionan de verdad en V19_4_HIBRIDO_CORREGIDO.py -- se sacaron a
      core/gbif_biodiversidad.py el 2026-09-04 para que ese script y
      este las compartan sin duplicar código.
    - Cruce NOM-059/IUCN/endemismo: nuevo en core/gbif_biodiversidad.py
      (2026-09-04), usando la tabla nacional de 2,678 especies
      digitalizada de NOM-059-SEMARNAT-2010 (DOF 14/11/2019 + Fe de
      erratas 04/03/2020) -- antes de esto, NINGÚN script del repo tenía
      este cruce (verificado con grep exhaustivo).
    - Terreno 3D + malla (DEM, hidrología D8): core/geomatica.py, EXACTO
      mismo mecanismo que ya usa core/deforestacion.py para dibujar las
      manchas de pérdida Hansen -- este módulo solo agrega una capa
      Scatter3d nueva sobre esa misma malla, no descarga/calcula terreno
      por su cuenta.
    - Tarjetas HTML por zona: core/reportes_html.py, mismo diseño que ya
      usan core/carbono.py y core/carbono_perdida.py.

ZONAS, NO MUNICIPIOS: el script viejo de Texolo (irdcloudbiodiversidadgbf.py,
nunca integrado a este repo) dividía por municipio (Xico/Teocelo), lo que
requería un script de pre-proceso que no existe aquí
(irdcludramsar_v2_split_municipal.py). Este módulo usa la MISMA unidad
que ya usa deforestación/carbono/hidrología: núcleo/buffer_500m/
buffer_1000m -- comparable entre sí y sin depender de nada externo.

MAPA 3D: SOLO especies en algún nivel de riesgo (CRITICO/ALTO/MEDIO/
VIGILAR) -- ésa es la capa "zonas a vigilar" pensada para community/
conferencia, no un mapa de "todas las especies" (ese es el CSV completo).

MAPA SATELITAL 2D (`_mapa_satelital_riesgo()`, 2026-09-04): adaptado del
patrón de `generar_mapa_satelital_riesgo()` del script viejo de Texolo
(irdcloudbiodiversidadgbf.py, nunca integrado a este repo), pero con
anillo núcleo/buffer_500m/buffer_1000m dibujado en vez de contorno por
municipio (mismo criterio "ZONAS, NO MUNICIPIOS" de arriba) y reusando los
colores de zona de core/reportes_html.py (mismo `#d62728`/`#e8892b`/
`#c9a227` que ya ve el usuario en las tarjetas). Solo especies en riesgo,
numeradas, con inset de zoom sobre el cluster más denso y tabla de
coordenadas para ir a campo -- igual que el mapa 3D, pero sobre imagen
satelital real (Esri World Imagery) en vez de terreno DEM. `contextily`
es dependencia OPCIONAL: si no está instalado, o si no hay salida a
internet a las teselas en el momento de correr, esta función AVISA (WARN)
y el resto del módulo sigue sin el mapa 2D -- no truena (mismo patrón que
`verificador_biodiversidad.py`/`rasterio` para el área 3D del dictamen
CONANP). NO SE PUDO VERIFICAR el render real (con teselas descargadas)
desde el entorno donde se escribió esta función, por no tener contextily
instalado ni salida a internet ahí -- sí se verificó que degrada
correctamente (WARN, sin truenar) cuando falta la dependencia.

LO QUE ESTE MÓDULO **NO** HACE TODAVÍA, A PROPÓSITO, NO ESCONDIDO:
    - Informe PDF técnico -- la entrega CSV + página HTML ya cubre lo
      pedido para este primer corte; un PDF tipo FPDF se puede agregar
      después si hace falta para un público distinto (biólogos/CONANP).
    - Lista nacional de endemismo (sigue siendo una semilla de 4
      especies, ver core/gbif_biodiversidad.py) -- construirla es un
      trabajo de digitalización aparte, del tamaño del de NOM-059.

No está conectado a run_pipeline.sh todavía -- se corre aparte
(`python3 -m core.biodiversidad ...`) hasta validarlo con un sitio real.
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd

from config import ZONAS_ANALISIS_M, CARPETA_SRTM, PERCENTIL_CAUCE_HIDROLOGIA, DEFORESTACION_ANIO_INICIO_DEFAULT, log

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gbif_biodiversidad import (  # noqa: E402
    zonas_anillo_exclusivo, descargar_biodiversidad_gbif,
    cargar_tabla_nom059, cargar_notas_nom059, cargar_tabla_endemismo, enriquecer_con_riesgo,
)
from salamandra_biodiversidad import validate_taxonomic_class, FAMILY_CLASS_MAP  # noqa: E402
# Rarefacción + Chao1 (2026-09-06, pedido explícito del usuario: "que no nos
# puedan decir que está sesgado" por diferencia de esfuerzo de muestreo entre
# zonas). Se eligió verificador_biodiversidad.py sobre
# salamandra_biodiversidad.rarefy_richness() -- ésta última solo compara un
# punto fijo por zona; verificador_biodiversidad da la curva completa
# (riqueza esperada por tamaño de muestra, con IC 95%) MÁS el estimador
# Chao1 (especies no detectadas todavía) con los mismos datos -- más
# defendible frente a una audiencia académica. Ninguna de las dos estaba
# conectada a este módulo antes de hoy (verificado con grep).
from verificador_biodiversidad import rarefaccion_monte_carlo, chao1_riqueza  # noqa: E402

NOMBRES_LEGIBLES_NIVEL = {
    "NUCLEO_DESGLOSE_29": "nucleo",
    "BUFFER_500m": "buffer_500m",
    "BUFFER_1000m": "buffer_1000m",
}
NIVELES_EN_RIESGO = ("CRITICO", "ALTO", "MEDIO", "VIGILAR")
MARCADORES_CLASE = {
    "Aves": "o", "Mammalia": "^", "Amphibia": "s", "Reptilia": "D",
    "Insecta": "P", "Arachnida": "X", "Chondrichthyes": "v",
}


def _separar_coordenadas_para_dibujo(lons, lats, radio_grados=0.0009):
    """Cuando 2+ registros caen en la misma coordenada (o casi), matplotlib
    los dibuja apilados uno sobre otro (colores mezclados por transparencia,
    números encimados). Esta función separa las posiciones de DIBUJO en
    círculo alrededor del punto original -- no toca lat/lon reales (esas
    son las que van al CSV/tabla de campo, sin modificar). radio_grados
    default (~90-100m a latitudes templadas) pensado para un mapa que cubre
    varios km (sitio completo + buffers), no para el mapa 3D de detalle."""
    import math
    from collections import defaultdict

    grupos = defaultdict(list)
    for i, (lon, lat) in enumerate(zip(lons, lats)):
        grupos[(round(lon, 5), round(lat, 5))].append(i)

    lons_dibujo, lats_dibujo = list(lons), list(lats)
    for indices in grupos.values():
        n = len(indices)
        if n <= 1:
            continue
        for k, idx in enumerate(indices):
            angulo = 2 * math.pi * k / n
            lons_dibujo[idx] = lons[idx] + radio_grados * math.cos(angulo)
            lats_dibujo[idx] = lats[idx] + radio_grados * math.sin(angulo)
    return lons_dibujo, lats_dibujo


def _envolver_texto_hover(texto, ancho=58):
    """Inserta <br> cada ~`ancho` caracteres (en el espacio más cercano, sin
    cortar palabras) -- pensado para notas largas (NOTA_ENDEMISMO/NOTA_RIESGO)
    que se meten en un hovertemplate de Plotly.

    OJO 2026-09-06 (pedido real del usuario viendo el HTML de Texolo): sin
    esto, una nota de una sola línea sin saltos (ej. la cita académica de
    Charadrahyla taeniopus, ~280 caracteres) hace que Plotly dibuje el cuadro
    de hover TAN ANCHO como esa línea -- en el screenshot que mandó se veía
    como un bloque rojo que atravesaba toda la figura, tapando la leyenda y
    el mapa. `textwrap` ya resuelve esto bien (no corta palabras a la mitad),
    solo hace falta unir con '<br>' en vez de '\\n' porque esto va dentro de
    HTML de Plotly, no a una terminal."""
    if not texto:
        return texto
    import textwrap
    return "<br>".join(textwrap.wrap(str(texto), width=ancho))


def _separar_para_inset(xs, ys, distancia_min_m=45.0, iteraciones=60):
    """Complemento de _separar_coordenadas_para_dibujo(), pensado
    específicamente para las coordenadas que se dibujan DENTRO de un
    recuadro de zoom (inset). Esa función de arriba solo separa
    duplicados EXACTOS (mismo lon/lat redondeado a ~1m) -- suficiente
    para el mapa completo (varios km de ancho), pero un inset de zoom
    muestra apenas unos cientos de metros: dos puntos a, digamos, 20-30m
    de distancia real se ven perfectamente separados en el mapa
    principal, pero dentro del recuadro (círculos de ~15-20px de
    diámetro) pueden terminar pegados o tapándose uno a otro -- feedback
    real del usuario ("si notas en el circulo hay unos ocultitos").

    MÉTODO: relajación iterativa tipo "resolver solapes de etiquetas" --
    en cada pasada, cualquier par de puntos a menos de `distancia_min_m`
    se empuja en direcciones opuestas (mitad de la distancia que falta
    cada uno) a lo largo de la línea que los une; se repite hasta que
    ningún par quede más cerca del mínimo, o se agoten las iteraciones.
    Solo cambia la posición de DIBUJO dentro del inset -- xs/ys reales
    (CSV/tabla de coordenadas) nunca se tocan, igual que la función
    hermana de arriba.

    Por qué no un reparto en círculo por grupo (primer intento, ya
    descartado): agrupar por vecindad y repartir cada grupo en un
    círculo separa bien a los miembros DE UN MISMO grupo entre sí, pero
    dos grupos VECINOS pueden terminar con miembros que se acercan de
    nuevo entre grupos -- verificado con el cluster de 25 puntos del
    sitio real (Ramsar 1601 Texolo): quedaban 2 pares aún encimados sin
    importar el radio elegido. La relajación iterativa, al comparar
    TODOS los pares en cada pasada (no solo dentro de un grupo), no
    tiene ese punto ciego -- con el mismo cluster de 25 puntos deja
    CERO pares por debajo del mínimo."""
    import math

    n = len(xs)
    xs_out, ys_out = list(xs), list(ys)
    for _ in range(iteraciones):
        movido = False
        for a in range(n):
            for b in range(a + 1, n):
                dx, dy = xs_out[b] - xs_out[a], ys_out[b] - ys_out[a]
                d = (dx * dx + dy * dy) ** 0.5
                if d >= distancia_min_m:
                    continue
                movido = True
                if d < 1e-6:
                    # Coordenada idéntica -- separar en una dirección
                    # arbitraria pero DISTINTA por par (a,b), para que dos
                    # pares distintos de duplicados no se muevan igual.
                    ang = 2 * math.pi * ((a * 7 + b * 13) % 360) / 360.0
                    ux, uy = math.cos(ang), math.sin(ang)
                else:
                    ux, uy = dx / d, dy / d
                empuje = (distancia_min_m - d) / 2.0
                xs_out[a] -= ux * empuje
                ys_out[a] -= uy * empuje
                xs_out[b] += ux * empuje
                ys_out[b] += uy * empuje
        if not movido:
            break
    return xs_out, ys_out


# ==============================================================================
# --- MAPA SATELITAL 2D -- ver docstring del módulo ("MAPA SATELITAL 2D")
#     para de dónde sale el patrón y por qué contextily es opcional.
# ==============================================================================
def _mapa_satelital_riesgo(df_bio, geom_zonas_wgs84, id_proyecto, carpeta_salida, subtitulo_extra=""):
    try:
        import contextily as ctx
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import geopandas as gpd
        from shapely.geometry import Point
    except ImportError as exc:
        log(f"Mapa satelital 2D OMITIDO -- falta una dependencia ({exc}). Instala con: "
            f"pip install contextily --break-system-packages", nivel="WARN")
        return None, None

    col_clase = "CLASE_VALIDADA" if "CLASE_VALIDADA" in df_bio.columns else ("CLASE" if "CLASE" in df_bio.columns else None)
    df_riesgo = df_bio[df_bio["NIVEL_RIESGO"].isin(NIVELES_EN_RIESGO)].copy()
    df_riesgo = df_riesgo.dropna(subset=["LAT", "LON"]) if {"LAT", "LON"}.issubset(df_riesgo.columns) else df_riesgo.iloc[0:0]
    df_riesgo = df_riesgo.drop_duplicates(subset=["NOMBRE_CIENTIFICO", "LAT", "LON"]).reset_index(drop=True)

    if len(df_riesgo) == 0:
        log("Sin especies en riesgo con coordenadas -- se omite el mapa satelital 2D.", nivel="WARN")
        return None, None

    df_riesgo["numero"] = df_riesgo.index + 1
    df_riesgo["zona_legible"] = (
        df_riesgo["NIVEL"].map(lambda n: NOMBRES_LEGIBLES_NIVEL.get(n, str(n)))
        if "NIVEL" in df_riesgo.columns else ""
    )
    lon_dibujo, lat_dibujo = _separar_coordenadas_para_dibujo(df_riesgo["LON"].tolist(), df_riesgo["LAT"].tolist())
    gdf_puntos = gpd.GeoDataFrame(
        df_riesgo, geometry=[Point(xy) for xy in zip(lon_dibujo, lat_dibujo)], crs="EPSG:4326",
    ).to_crs(epsg=3857)

    fig, ax = plt.subplots(figsize=(13, 13))

    from core import reportes_html
    for zona, geom in (geom_zonas_wgs84 or {}).items():
        if geom is None:
            continue
        color = reportes_html.COLORES_ZONA_HEX.get(zona, "#ffff00")
        gpd.GeoDataFrame(geometry=[geom], crs="EPSG:4326").to_crs(epsg=3857).boundary.plot(
            ax=ax, color=color, linewidth=2.4, zorder=4, label=zona,
        )

    for _, row in gdf_puntos.iterrows():
        clase = row.get(col_clase, "") if col_clase else ""
        marker = MARCADORES_CLASE.get(clase, "o")
        ax.scatter(row.geometry.x, row.geometry.y, c=row["COLOR_RIESGO"], s=220, marker=marker,
                   edgecolors="black", linewidths=1.6, zorder=6)
        ax.annotate(str(row["numero"]), (row.geometry.x, row.geometry.y),
                    color="white", fontsize=10, fontweight="bold", ha="center", va="center", zorder=7)

    try:
        # OJO 2026-09-06 (pedido por el usuario -- "en algunos el letrerito se pierde"):
        # el texto de atribución que contextily dibuja solo (largo: "Tiles (C) Esri --
        # Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye...") caía en la MISMA esquina
        # inferior izquierda que nuestro propio aviso de "N especies numeradas..." (ver
        # más abajo) -- los dos textos se encimaban y se volvían ilegibles. Se apaga la
        # atribución automática aquí y se agrega UNA sola, corta, combinada con nuestro
        # aviso -- sigue dando el crédito a Esri (obligatorio por sus términos de uso),
        # solo que en un solo texto en vez de dos compitiendo por el mismo rincón.
        ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom="auto", attribution=False)
    except Exception as e:
        log(f"No se pudo descargar el fondo satelital ({str(e)[:120]}). Verifica conexión a internet -- "
            f"el mapa se genera solo con los puntos/anillos, sin imagen de fondo.", nivel="WARN")

    # --- insets de zoom sobre TODOS los clusters densos (no solo el más denso) ---
    # NOTA 2026-09-05: antes solo se buscaba el cluster #1 (el de mas puntos) y se
    # le ponia un solo recuadro de zoom -- feedback real del usuario: en el mismo
    # mapa suele haber un segundo (o tercer) amontonamiento de puntos en otra parte
    # del sitio (ej. el cluster de 5 puntos junto al cauce central) que se quedaba
    # sin su propio zoom y se veia encimado en el mapa principal. Ahora se buscan
    # clusters de forma golosa (greedy): se toma el punto con mas vecinos dentro
    # del radio, se arma su cluster, se quitan esos puntos de la bolsa, y se repite
    # hasta agotar los puntos densos o llegar al maximo de recuadros que caben sin
    # taparse entre si (ni tapar la leyenda de zonas, que vive en "lower right").
    # NOTA 2026-09-06: el arreglo de arriba (multiples insets) trajo un problema
    # nuevo que el usuario detecto viendo el resultado -- los recuadros son
    # imagenes OPACAS dibujadas encima del mapa principal, y por pura coincidencia
    # de posicion algunos puntos sueltos (que no son parte de ningun cluster)
    # quedaban tapados detras de un recuadro que vive en esa misma esquina. Antes
    # de fijar la posicion de cada inset, ahora se revisa cuantos puntos ajenos a
    # su propio cluster taparia en cada esquina candidata, y se elige la que tape
    # menos (idealmente cero) -- ver _rect_datos_de_inset()/_n_puntos_tapados()
    # mas abajo.
    xs = gdf_puntos.geometry.x.values
    ys = gdf_puntos.geometry.y.values
    # NOTA 2026-09-06: max_insets subió de 3 a 4 -- un sitio real (Ramsar
    # 1601 Texolo) demostró tener 4 clusters densos de verdad, no 3: con
    # el límite viejo, el 4o cluster se quedaba sin su propio recuadro. Si
    # algún día hubiera un 5o cluster denso, el aviso WARN de abajo ("Se
    # detectaron más de N clusters densos...") lo reporta en vez de fallar
    # en silencio.
    radio_vecindad_m, min_puntos_para_inset, max_insets = 350, 4, 4
    # 6 posiciones candidatas (más que max_insets=4 A PROPÓSITO) -- ver
    # razonamiento completo en la búsqueda de asignación más abajo: entre
    # más posiciones libres para elegir, más fácil encontrar una
    # combinación sin colisiones para TODOS los clusters a la vez.
    posiciones_inset = [
        dict(loc="upper right", bbox_to_anchor=(0, 0, 0.98, 0.97)),
        dict(loc="upper left", bbox_to_anchor=(0.02, 0, 1, 0.97)),
        dict(loc="lower left", bbox_to_anchor=(0.02, 0.03, 1, 1)),
        dict(loc="lower right", bbox_to_anchor=(0, 0.03, 0.98, 1)),
        dict(loc="center left", bbox_to_anchor=(0.02, 0, 1, 1)),
        dict(loc="center right", bbox_to_anchor=(0, 0, 0.98, 1)),
    ]

    indices_disponibles = set(range(len(xs)))
    clusters = []
    while indices_disponibles:
        idxs_rest = np.array(sorted(indices_disponibles))
        mejor_idx_local, mejor_n_vecinos = None, 0
        for i in idxs_rest:
            dist = ((xs[idxs_rest] - xs[i]) ** 2 + (ys[idxs_rest] - ys[i]) ** 2) ** 0.5
            n_vecinos = int((dist <= radio_vecindad_m).sum())
            if n_vecinos > mejor_n_vecinos:
                mejor_n_vecinos, mejor_idx_local = n_vecinos, i
        if mejor_idx_local is None or mejor_n_vecinos < min_puntos_para_inset:
            break
        dist_centro = ((xs[idxs_rest] - xs[mejor_idx_local]) ** 2 + (ys[idxs_rest] - ys[mejor_idx_local]) ** 2) ** 0.5
        en_cluster_local = idxs_rest[dist_centro <= radio_vecindad_m]
        clusters.append((mejor_n_vecinos, en_cluster_local))
        indices_disponibles -= set(en_cluster_local.tolist())
        if len(clusters) >= max_insets:
            break
    clusters.sort(key=lambda c: c[0], reverse=True)

    if clusters:
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
        import itertools

        def _rect_datos_de_inset(pos, ancho_pct):
            """Crea el inset en la posición candidata `pos`, con tamaño de
            pantalla `ancho_pct`% (sin dibujar nada todavía) y devuelve
            (axins, rect) donde rect=(x0,x1,y0,y1) es el área que ese
            recuadro ocupa, expresada en las mismas coordenadas de datos
            que xs/ys (Web Mercator) -- así se puede comparar directo contra
            las coordenadas de los puntos del mapa principal. `ancho_pct`
            solo cambia qué tan grande se ve el recuadro EN LA PÁGINA, no
            el rango de datos que muestra (eso lo fija axins.set_xlim/
            set_ylim más abajo, aparte) -- por eso encoger el tamaño es una
            forma válida de reducir el área que un recuadro tapa en el
            mapa principal, sin perder ningún punto del cluster."""
            axins_ = inset_axes(ax, width=f"{ancho_pct}%", height=f"{ancho_pct}%", loc=pos["loc"],
                                 bbox_to_anchor=pos["bbox_to_anchor"], bbox_transform=ax.transAxes, borderpad=0)
            fig.canvas.draw()
            bbox_disp = axins_.get_window_extent()
            esquinas = ax.transData.inverted().transform(
                [(bbox_disp.x0, bbox_disp.y0), (bbox_disp.x1, bbox_disp.y1)])
            (dx0, dy0), (dx1, dy1) = esquinas
            return axins_, (min(dx0, dx1), max(dx0, dx1), min(dy0, dy1), max(dy0, dy1))

        def _n_puntos_tapados(rect, excluir_idx):
            x0, x1, y0, y1 = rect
            dentro = (xs >= x0) & (xs <= x1) & (ys >= y0) & (ys <= y1)
            return set(np.where(dentro)[0].tolist()) - set(excluir_idx.tolist())

        # NOTA 2026-09-06 (segunda vuelta): la versión anterior de este
        # bloque asignaba posiciones de forma VORAZ -- el cluster más
        # grande elegía primero la posición que menos le tapara A ÉL, sin
        # ver el efecto en los clusters que faltaban. En un caso real
        # (Ramsar 1601 Texolo) eso obligaba al último cluster a quedarse
        # con la única posición libre restante, aunque esa posición
        # tapara justo al punto 74 (Aquiloeurycea cafetalera, la especie
        # insignia) -- verificado reproduciendo la corrida completa contra
        # el CSV real del usuario. Dos cambios lo resuelven de raíz:
        #   1. TAMAÑO variable por recuadro (ver _rect_datos_de_inset):
        #      antes del 30% fijo, un recuadro en una esquina podía tapar
        #      un punto que a un tamaño más chico (25%/20%/15%) ya no
        #      tapa -- se prueba de más grande a más chico y se usa el
        #      primero sin colisión.
        #   2. BÚSQUEDA de la asignación óptima, no voraz: se prueban
        #      TODAS las formas de repartir las posiciones candidatas
        #      entre los clusters detectados (permutaciones -- con
        #      max_insets<=4 clusters y 6 posiciones candidatas son como
        #      mucho unos cientos de combinaciones, trivial en tiempo) y
        #      se elige la que tapa MENOS puntos en total, con corte
        #      temprano en cuanto se encuentra una combinación con cero.
        #      Verificado con el CSV real: esto baja el total de puntos
        #      tapados de 1 (cafetalera, con el método voraz de antes) a
        #      0 -- ver bitácora de esta corrida.
        TAMANOS_INSET_PCT = (30, 25, 20, 15)
        n_clusters, n_pos = len(clusters), len(posiciones_inset)

        # El rectángulo de datos que ocupa (posición, tamaño) NO depende de
        # qué cluster lo use -- solo de los límites actuales del mapa. Se
        # calcula UNA sola vez por combinación (antes se recalculaba por
        # cada cluster, repitiendo el mismo fig.canvas.draw() -- el paso
        # caro -- hasta n_clusters veces de más por combinación).
        rects_por_posicion_tamano = {}
        for j in range(n_pos):
            for pct in TAMANOS_INSET_PCT:
                axins_probar, rect = _rect_datos_de_inset(posiciones_inset[j], pct)
                axins_probar.remove()
                rects_por_posicion_tamano[(j, pct)] = rect

        matriz_tapados = [[None] * n_pos for _ in range(n_clusters)]
        matriz_pct = [[None] * n_pos for _ in range(n_clusters)]
        for i in range(n_clusters):
            for j in range(n_pos):
                mejor_tapa, mejor_pct = None, None
                for pct in TAMANOS_INSET_PCT:
                    tapa = _n_puntos_tapados(rects_por_posicion_tamano[(j, pct)], clusters[i][1])
                    if mejor_tapa is None or len(tapa) < len(mejor_tapa):
                        mejor_tapa, mejor_pct = tapa, pct
                        if not mejor_tapa:
                            break
                matriz_tapados[i][j] = mejor_tapa
                matriz_pct[i][j] = mejor_pct

        mejor_asignacion, mejor_total = None, None
        for combo in itertools.permutations(range(n_pos), n_clusters):
            total = sum(len(matriz_tapados[i][j]) for i, j in enumerate(combo))
            if mejor_total is None or total < mejor_total:
                mejor_asignacion, mejor_total = combo, total
                if mejor_total == 0:
                    break

        for i, (n_vecinos, en_cluster) in enumerate(clusters):
            j = mejor_asignacion[i]
            pos_elegida, pct_elegido, tapados = posiciones_inset[j], matriz_pct[i][j], matriz_tapados[i][j]
            axins, _rect_final = _rect_datos_de_inset(pos_elegida, pct_elegido)
            if tapados:
                log(f"El recuadro de zoom del cluster {i + 1} tapa {len(tapados)} punto(s) sueltos del "
                    f"mapa principal en todas las combinaciones de posición/tamaño probadas -- se dejó en la "
                    f"que menos tapa. Considera revisar esos puntos en la tabla de coordenadas.", nivel="WARN")

            # NOTA 2026-09-06: dentro del recuadro de zoom, dos puntos a
            # unos 20-30m reales uno del otro (perfectamente separados en
            # el mapa principal, que cubre kilómetros) pueden verse
            # pegados/tapados uno sobre otro -- feedback real del usuario
            # ("si notas en el circulo hay unos ocultitos"). Se calcula una
            # separación de DIBUJO adicional SOLO para este recuadro (no
            # afecta la posición en el mapa principal ni las coordenadas
            # reales de la tabla) -- ver _separar_para_inset() arriba.
            xs_inset, ys_inset = _separar_para_inset(xs[en_cluster].tolist(), ys[en_cluster].tolist())

            # El límite del recuadro se calcula sobre la UNIÓN de las
            # coordenadas originales del cluster y las ya separadas -- en
            # clusters grandes _separar_para_inset() puede necesitar un
            # radio de separación mayor al margen fijo de antes, y si el
            # límite del recuadro no lo cubre, el punto separado quedaría
            # CORTADO fuera de la vista (peor que el traslape que se
            # intentaba arreglar). Unir ambos rangos garantiza que todo lo
            # que se dibuja cabe dentro del recuadro, sin importar cuánto
            # haya separado ese cluster en particular.
            margen_m = 110
            x0 = min(xs[en_cluster].min(), min(xs_inset)) - margen_m
            x1 = max(xs[en_cluster].max(), max(xs_inset)) + margen_m
            y0 = min(ys[en_cluster].min(), min(ys_inset)) - margen_m
            y1 = max(ys[en_cluster].max(), max(ys_inset)) + margen_m
            for pos_dibujo, idx in enumerate(en_cluster):
                row = gdf_puntos.iloc[idx]
                clase = row.get(col_clase, "") if col_clase else ""
                marker = MARCADORES_CLASE.get(clase, "o")
                x_dibujo, y_dibujo = xs_inset[pos_dibujo], ys_inset[pos_dibujo]
                axins.scatter(x_dibujo, y_dibujo, c=row["COLOR_RIESGO"], s=260, marker=marker,
                              edgecolors="black", linewidths=1.6, zorder=6)
                axins.annotate(str(row["numero"]), (x_dibujo, y_dibujo),
                               color="white", fontsize=10, fontweight="bold", ha="center", va="center", zorder=7)
            try:
                ctx.add_basemap(axins, source=ctx.providers.Esri.WorldImagery, zoom="auto")
            except Exception:
                pass
            axins.set_xlim(x0, x1)
            axins.set_ylim(y0, y1)
            axins.set_xticks([])
            axins.set_yticks([])
            for spine in axins.spines.values():
                spine.set_edgecolor("yellow")
                spine.set_linewidth(2.5)
            axins.set_title(f"Zoom -- cluster {i + 1} ({n_vecinos} puntos)", fontsize=8.5,
                            fontweight="bold", color="white", bbox=dict(facecolor="black", alpha=0.7, pad=3))
            mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec="yellow", linewidth=1.5, zorder=5)

        # si se agotaron los recuadros disponibles pero todavia queda un cluster
        # denso sin mostrar, se avisa en el log en vez de fallar en silencio
        if indices_disponibles:
            idxs_rest = np.array(sorted(indices_disponibles))
            queda_denso = False
            for i in idxs_rest:
                dist = ((xs[idxs_rest] - xs[i]) ** 2 + (ys[idxs_rest] - ys[i]) ** 2) ** 0.5
                if int((dist <= radio_vecindad_m).sum()) >= min_puntos_para_inset:
                    queda_denso = True
                    break
            if queda_denso:
                log(f"Se detectaron más de {max_insets} clusters densos -- solo se dibujan los "
                    f"{max_insets} más grandes para no saturar el mapa con recuadros.", nivel="WARN")

    ax.set_axis_off()
    ax.set_title(f"{id_proyecto} -- ZONAS A VIGILAR (satelital)\nEspecies en riesgo NOM-059/IUCN"
                 f"{subtitulo_extra}", fontsize=13, fontweight="bold", color="white",
                 bbox=dict(facecolor="black", alpha=0.6, pad=6))
    if geom_zonas_wgs84:
        # NOTA 2026-09-06: la leyenda de zonas vivía en loc="lower right",
        # DENTRO de la imagen -- esa es justo una de las 4 esquinas que
        # ahora puede ocupar un recuadro de zoom (ver max_insets=4 arriba).
        # Con 4 clusters densos reales las 4 esquinas quedan ocupadas por
        # insets, así que dejar la leyenda ahí la taparía. Se saca
        # COMPLETAMENTE fuera del área de la imagen (a la derecha) para
        # que nunca compita por espacio con ningún recuadro sin importar
        # cuántos haya -- bbox_inches="tight" en el savefig de abajo
        # expande el lienzo para que quepa completa, no se corta.
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=9, framealpha=0.85, borderaxespad=0)

    # NOTA 2026-09-04: antes, aquí se dibujaba la tabla de coordenadas
    # completa ENCIMA de la imagen satelital (fig.text en la esquina
    # inferior izquierda). Funcionaba con ~10 especies en riesgo, pero en
    # cuanto el cruce NOM-059 empezó a detectar de verdad (35+ especies,
    # ver core/gbif_biodiversidad.py) la tabla creció a 70-80 filas y tapaba
    # buena parte del mapa -- feedback real del usuario ("sale amontonada
    # la tabla... se ve un poquito amontonadito"). La tabla completa ahora
    # se arma aparte, como HTML debajo de esta imagen (ver
    # reportes_html.tabla_especies_riesgo_html(), llamada desde
    # _procesar_biodiversidad() con el mismo df_riesgo/numero que aquí, así
    # los números del mapa y los de la tabla siempre coinciden). Sobre la
    # imagen solo queda un aviso corto de dónde encontrar la tabla.
    ax.text(0.015, 0.015, f"{len(df_riesgo)} especies numeradas -- tabla de coordenadas completa debajo de esta "
            f"imagen | Imagery (C) Esri, Maxar, Earthstar Geographics, and the GIS User Community",
            transform=ax.transAxes, fontsize=9.5, color="white", va="bottom", ha="left",
            bbox=dict(facecolor="black", alpha=0.6, pad=5))

    png_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_satelital_riesgo.png")
    plt.tight_layout()
    plt.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log(f"Mapa satelital 2D (zonas a vigilar): {png_path} ({len(df_riesgo)} puntos numerados)", nivel="OK")
    return png_path, df_riesgo


def _marcar_deforestacion_en_puntos(df, hidrologia, utm_crs, lossyear_alineado):
    """Agrega CERCA_DEFORESTACION (True/False/None) y
    ANIO_DEFORESTACION_CERCANA (año real o NaN) a `df`, cruzando cada
    registro (por LAT/LON) contra `lossyear_alineado` -- la banda Hansen
    YA alineada a la MISMA malla de `hidrologia` (ver
    core.deforestacion._descargar_lossyear_alineado(), reusada tal cual
    por _descargar_lossyear_para_biodiversidad() más abajo, sin
    reimplementar la descarga/realineación).

    Mismo mecanismo de transformación lat/lon -> fila/columna de malla
    que ya usa el bucle del mapa 3D de riesgo (ver
    _procesar_biodiversidad) -- un solo lugar donde vive esa cuenta.

    LIMITACIÓN A PROPÓSITO, NO ESCONDIDA: esto es una COINCIDENCIA
    GEOGRÁFICA a nivel de píxel Hansen (30m de resolución), nunca una
    prueba de que la deforestación causó el riesgo de esa especie --
    dos cosas cercanas en el mapa no implican relación causal, eso
    requeriría un estudio biológico real. Se documenta también en la
    nota al pie del reporte y en el título de la pastilla de la tabla.

    Siempre agrega las dos columnas, incluso si `lossyear_alineado` es
    None (--con-deforestacion no se usó) -- en ese caso quedan en
    None/NaN para TODOS los registros, que la pastilla de
    reportes_html.py dibuja como "No evaluado" (gris), nunca como "Sin
    pérdida" (verde) -- no evaluar no es lo mismo que evaluar y no
    encontrar nada. Un punto sin LAT/LON, o que cae fuera de la malla
    visual del terreno 3D (buffer máximo), también queda en None/NaN por
    la misma razón."""
    import pyproj

    out = df.copy()
    out["CERCA_DEFORESTACION"] = None
    out["ANIO_DEFORESTACION_CERCANA"] = np.nan
    if lossyear_alineado is None or not {"LAT", "LON"}.issubset(out.columns):
        return out

    Z_smooth = hidrologia["Z_smooth"]
    transform = hidrologia["transform"]
    rows, cols = Z_smooth.shape
    inv = ~transform
    a_utm = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True).transform

    cerca_col, anio_col = [], []
    for lat, lon in zip(out["LAT"], out["LON"]):
        if pd.isna(lat) or pd.isna(lon):
            cerca_col.append(None)
            anio_col.append(np.nan)
            continue
        x_utm, y_utm = a_utm(float(lon), float(lat))
        col_f, row_f = inv * (x_utm, y_utm)
        row_i, col_i = int(round(row_f)), int(round(col_f))
        if not (0 <= row_i < rows and 0 <= col_i < cols) or np.isnan(Z_smooth[row_i, col_i]):
            cerca_col.append(None)  # fuera de la malla visual -- no evaluado, nunca "sin pérdida"
            anio_col.append(np.nan)
            continue
        codigo = int(lossyear_alineado[row_i, col_i])
        if codigo > 0:
            cerca_col.append(True)
            anio_col.append(2000 + codigo)
        else:
            cerca_col.append(False)
            anio_col.append(np.nan)
    out["CERCA_DEFORESTACION"] = cerca_col
    out["ANIO_DEFORESTACION_CERCANA"] = anio_col
    return out


def _descargar_lossyear_para_biodiversidad(geom_utm_nucleo, zonas_m, hidrologia, utm_crs, carpeta_srtm, proyecto_gee):
    """Descarga la banda Hansen 'lossyear' alineada a la MISMA malla que ya
    usa el terreno 3D de este módulo -- reusa core.deforestacion (misma
    función/patrón que ya corre en producción para el mapa 3D dedicado de
    deforestación, `generar_mapa_3d_deforestacion()`), sin reimplementar
    nada. Devuelve None (con un WARN, NUNCA truena el resto del reporte)
    si Earth Engine no está instalado, no se pudo autenticar, o la
    descarga falla por cualquier razón -- mismo espíritu "degrada, no
    truena" que ya usa este módulo para contextily/el mapa satelital 2D
    (ver docstring de _mapa_satelital_riesgo)."""
    try:
        import ee
        from shapely.geometry import mapping
        from shapely.ops import transform as shp_transform
        import pyproj
        from core import deforestacion
    except ImportError as exc:
        log(f"Cruce con deforestación OMITIDO -- falta una dependencia ({exc}). Instala earthengine-api "
            f"con: pip install earthengine-api --break-system-packages", nivel="WARN")
        return None

    try:
        if proyecto_gee:
            ee.Initialize(project=proyecto_gee)
        else:
            ee.Initialize()
    except Exception as e:
        log(f"Cruce con deforestación OMITIDO -- no se pudo inicializar Earth Engine ({e}). Corre "
            f"'earthengine authenticate' una vez en esta máquina, o revisa --proyecto-gee.", nivel="WARN")
        return None

    try:
        a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True).transform
        geom_visual_utm = geom_utm_nucleo.buffer(max(zonas_m))
        geom_visual_wgs84 = shp_transform(a_wgs84, geom_visual_utm)
        log("Descargando Hansen lossyear alineado al terreno (para cruzar con los puntos de biodiversidad)...")
        lossyear_alineado = deforestacion._descargar_lossyear_alineado(
            ee, mapping(geom_visual_wgs84), hidrologia["transform"], hidrologia["Z_raw"].shape, utm_crs,
            carpeta_tmp=carpeta_srtm,
        )
        log("Hansen lossyear descargado y alineado -- listo para cruzar con los puntos de biodiversidad y "
            "para dibujarse en el mismo mapa 3D.", nivel="OK")
        return lossyear_alineado
    except Exception as e:
        log(f"Cruce con deforestación OMITIDO -- la descarga de Hansen lossyear falló ({e}). El resto del "
            f"reporte de biodiversidad se genera igual, sin esta capa.", nivel="WARN")
        return None


# ==============================================================================
# --- RAREFACCIÓN + CHAO1 POR ZONA (2026-09-06) -- ver nota de import arriba
#     para por qué se eligió verificador_biodiversidad.py sobre
#     salamandra_biodiversidad.rarefy_richness(). Corre sobre TODOS los
#     registros de cada zona (no solo especies en riesgo): la pregunta que
#     responde es "¿alcanzó el muestreo?", no "¿cuántas en riesgo hay?".
# ==============================================================================
MIN_REGISTROS_RAREFACCION = 5  # por debajo de esto, la curva/Chao1 no dice nada útil -- se omite (WARN), no se fuerza


def _rarefaccion_por_zona(df_bio, carpeta_salida, id_proyecto, n_iter=200, seed=42):
    """Para cada zona (NIVEL -> nombre legible), calcula la curva de
    rarefacción individual-based (riqueza esperada por tamaño de muestra,
    con IC 95%) y el estimador Chao1 (especies no detectadas estimadas),
    ambos de verificador_biodiversidad.py, sobre la lista completa de
    NOMBRE_CIENTIFICO de esa zona (todos los registros, no solo riesgo).

    Devuelve (df_curvas, resumen_por_zona: dict zona->{...}, png_path o None).
    Zonas con menos de MIN_REGISTROS_RAREFACCION registros se excluyen de
    la curva/Chao1 (no hay suficiente para que el número signifique algo)
    pero SÍ aparecen en el resumen con sus valores en None -- para que se
    note la omisión en vez de faltar en silencio."""
    if "NIVEL" not in df_bio.columns:
        log("Sin columna NIVEL -- se omite rarefacción/Chao1 por zona.", nivel="WARN")
        return pd.DataFrame(), {}, None

    curvas = []
    resumen = {}
    for nivel_gbif, nombre_zona in NOMBRES_LEGIBLES_NIVEL.items():
        especies = df_bio.loc[df_bio["NIVEL"] == nivel_gbif, "NOMBRE_CIENTIFICO"].dropna().tolist()
        n_regs = len(especies)
        if n_regs < MIN_REGISTROS_RAREFACCION:
            log(f"[RAREFACCION] Zona '{nombre_zona}' con solo {n_regs} registro(s) -- se omite curva/Chao1 "
                f"(mínimo {MIN_REGISTROS_RAREFACCION}, el número no sería confiable).", nivel="WARN")
            resumen[nombre_zona] = {
                "rarefaccion_n_registros": n_regs, "rarefaccion_s_observada": None,
                "rarefaccion_chao1_estimado": None, "rarefaccion_especies_no_detectadas_est": None,
                "rarefaccion_singletons_f1": None,
            }
            continue
        curva = rarefaccion_monte_carlo(especies, n_iter=n_iter, seed=seed)
        curva.insert(0, "zona", nombre_zona)
        curvas.append(curva)
        chao1 = chao1_riqueza(especies)
        resumen[nombre_zona] = {
            "rarefaccion_n_registros": n_regs, "rarefaccion_s_observada": chao1["s_observada"],
            "rarefaccion_chao1_estimado": chao1["chao1_estimado"],
            "rarefaccion_especies_no_detectadas_est": chao1["especies_no_detectadas_estimadas"],
            "rarefaccion_singletons_f1": chao1["singletons_f1"],
        }
        log(f"[RAREFACCION] Zona '{nombre_zona}': {n_regs} registros, {chao1['s_observada']} spp observadas, "
            f"Chao1={chao1['chao1_estimado']} (≈{chao1['especies_no_detectadas_estimadas']} spp no detectadas "
            f"estimadas).", nivel="OK")

    df_curvas = pd.concat(curvas, ignore_index=True) if curvas else pd.DataFrame()
    if not df_curvas.empty:
        curvas_csv_path = os.path.join(carpeta_salida, f"biodiversidad_rarefaccion_{id_proyecto.lower()}.csv")
        df_curvas.to_csv(curvas_csv_path, index=False, encoding="utf-8-sig")
        log(f"CSV de curvas de rarefacción guardado en: {curvas_csv_path}", nivel="OK")

    png_path = None
    if not df_curvas.empty:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError as exc:
            log(f"Gráfica de rarefacción OMITIDA -- falta matplotlib ({exc}). El CSV/resumen de rarefacción "
                f"se generan igual, sin la imagen.", nivel="WARN")
            return df_curvas, resumen, None

        from core import reportes_html
        fig, ax = plt.subplots(figsize=(8, 5.5))
        for nombre_zona, grupo in df_curvas.groupby("zona"):
            color = reportes_html.COLORES_ZONA_HEX.get(nombre_zona, "#666")
            ax.plot(grupo["n_muestra"], grupo["riqueza_media"], color=color, linewidth=2, label=nombre_zona)
            ax.fill_between(grupo["n_muestra"], grupo["riqueza_ci_low"], grupo["riqueza_ci_high"],
                             color=color, alpha=0.18, linewidth=0)
        ax.set_xlabel("Registros muestreados (n)")
        ax.set_ylabel("Especies esperadas (riqueza rarefectada)")
        ax.set_title(f"{id_proyecto} -- Curva de rarefacción por zona (sombreado = IC 95%)")
        ax.legend(frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        png_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_rarefaccion.png")
        fig.savefig(png_path, dpi=150)
        plt.close(fig)
        log(f"Gráfica de rarefacción guardada en: {png_path}", nivel="OK")

    return df_curvas, resumen, png_path


# ==============================================================================
# --- NÚCLEO COMPARTIDO: recibe un df_bio (real, de GBIF, o sintético de
#     --demo) + una hidrología ya calculada, y produce CSV + tarjetas +
#     mapa 3D + página HTML. Ambos caminos (real y demo) pasan por AQUÍ,
#     así el --demo prueba de verdad el mismo código que corre en un
#     sitio real -- lo único sintético es el origen de las especies y
#     del terreno, no la lógica que los procesa (mismo criterio que
#     core.deforestacion.demo(), que solo simula el array 'lossyear',
#     no la alineación/el mapa). ---
# ==============================================================================
def _procesar_biodiversidad(df_bio, hidrologia, utm_crs, id_proyecto, carpeta_salida,
                             tabla_nom059, tabla_endemismo, consultar_iucn, max_consultas_iucn,
                             mapa_3d, subtitulo_extra="", metricas_zona=None,
                             mapa_2d=True, geom_zonas_wgs84=None,
                             calcular_endemismo_heuristico=True, max_consultas_endemismo=300,
                             usar_cache_disco=True, ruta_cache_disco=None,
                             lossyear_alineado=None, anio_inicio_deforestacion=None, anio_fin_deforestacion=None,
                             calcular_rarefaccion=True, n_iter_rarefaccion=200):
    """`metricas_zona`: DataFrame opcional de geomatica.calcular_metricas_por_zona()
    (columnas zona/area_2d_ha/area_3d_ha/factor_relieve/...) -- mismo cálculo
    YA usado por core.analizar_sitio (terreno) para el resto de la plataforma,
    aquí reusado sin volver a cargar DEM ni tocar disco/red de más: el
    llamador (generar_biodiversidad()/demo()) ya tiene geom_utm_nucleo/
    dst_array/meta_utm en memoria por el propio cálculo de hidrología D8.
    Si no se da (None), las tarjetas/CSV simplemente no muestran hectáreas."""
    os.makedirs(carpeta_salida, exist_ok=True)

    # OJO 2026-09-06 (bug real, detectado por el usuario viendo el reporte de
    # Texolo: el TOTAL de hectáreas sumaba ~5,000+ ha cuando el sitio real es
    # mucho más chico -- "hay muñeca rusa en las hectáreas"): geomatica.
    # calcular_metricas_por_zona() devuelve, A PROPÓSITO, el área ACUMULADA
    # de cada disco (geom_nucleo.buffer(0m), .buffer(500m), .buffer(1000m) --
    # cada disco CONTIENE por completo al anterior), no el área del ANILLO
    # exclusivo. Esto es el mismo "zona" que ya usan los conteos de especies
    # (anillo EXCLUSIVO, cada registro cuenta en una sola zona), así que
    # sumar las 3 filas de metricas_zona tal cual TRIPLE-cuenta el núcleo y
    # DOBLE-cuenta el anillo de 500m -- exactamente el patrón "muñeca rusa"
    # que describió el usuario. core.deforestacion YA resuelve esto restando
    # el acumulado anterior (ver 'valor_anillo = valor - acumulado_prev' en
    # ese archivo) -- aquí faltaba la misma resta. Se corrige aquí, ANTES de
    # construir metricas_por_zona_map, para que el resto del código
    # (tarjetas, TOTAL, CSV) reciba área ya-exclusiva sin tener que tocarse.
    if metricas_zona is not None:
        metricas_zona = metricas_zona.sort_values("buffer_m").reset_index(drop=True)
        suma_discos_2d_antes = round(metricas_zona["area_2d_ha"].sum(), 1)  # lo que se hubiera sumado mal
        filas_anillo = []
        prev_a2d, prev_a3d = 0.0, 0.0
        for _, fila in metricas_zona.iterrows():
            a2d_anillo = fila["area_2d_ha"] - prev_a2d
            a3d_anillo = fila["area_3d_ha"] - prev_a3d
            fila_nueva = fila.to_dict()
            fila_nueva["area_2d_ha"] = round(a2d_anillo, 3)
            fila_nueva["area_3d_ha"] = round(a3d_anillo, 3)
            fila_nueva["factor_relieve"] = round(a3d_anillo / a2d_anillo, 4) if a2d_anillo > 0 else float("nan")
            filas_anillo.append(fila_nueva)
            prev_a2d, prev_a3d = fila["area_2d_ha"], fila["area_3d_ha"]
        metricas_zona = pd.DataFrame(filas_anillo)
        suma_anillos_2d_despues = round(metricas_zona["area_2d_ha"].sum(), 1)
        log(f"Métricas de zona convertidas de disco acumulado a anillo exclusivo -- área 2D total real del "
            f"sitio: {suma_anillos_2d_despues:,.1f} ha (antes de esta corrección se hubiera sumado "
            f"{suma_discos_2d_antes:,.1f} ha por error, triple/doble-contando núcleo y el anillo de 500m).",
            nivel="OK")

    metricas_por_zona_map = (
        metricas_zona.set_index("zona").to_dict("index") if metricas_zona is not None else {}
    )

    if "NOMBRE_CIENTIFICO" not in df_bio.columns:
        sys.exit(f"ERROR: el dataframe de biodiversidad no tiene columna NOMBRE_CIENTIFICO. "
                  f"Columnas: {list(df_bio.columns)}")

    # Validación taxonómica CLASE vs FAMILIA -- mismo mecanismo ya
    # obligatorio en V19_4_HIBRIDO_CORREGIDO.py (ver ese script para el
    # origen del bug histórico que esto corrige).
    if "FAMILIA" in df_bio.columns:
        df_bio, reporte_clase = validate_taxonomic_class(df_bio, class_col="CLASE", family_col="FAMILIA")
        col_clase = "CLASE_VALIDADA"
        log(f"Validación CLASE vs FAMILIA: {reporte_clase.registros_corregidos}/{reporte_clase.total_registros} "
            f"registros ({reporte_clase.porcentaje_corregido:.1f}%) corregidos "
            f"(tabla de {len(FAMILY_CLASS_MAP)} familias).", nivel="OK")
    else:
        col_clase = "CLASE" if "CLASE" in df_bio.columns else None
        log("Sin columna FAMILIA -- no se puede validar CLASE. Se usa la columna CLASE cruda, sin garantía.",
            nivel="WARN")

    # Cruce NOM-059/IUCN/endemismo -- lo nuevo de este módulo.
    # notas_nom059: SIEMPRE de la tabla nacional default (no del override
    # local por sitio) -- es solo para el aviso de "posible población
    # restringida" (ver PATRON_POBLACION_RESTRINGIDA_REGEX en
    # gbif_biodiversidad.py), independiente de qué categoría haya quedado
    # vigente para este sitio en particular.
    df_bio = enriquecer_con_riesgo(
        df_bio, col_especie="NOMBRE_CIENTIFICO", tabla_nom059=tabla_nom059, tabla_endemismo=tabla_endemismo,
        consultar_iucn=consultar_iucn, max_consultas_iucn=max_consultas_iucn,
        calcular_endemismo_heuristico=calcular_endemismo_heuristico, max_consultas_endemismo=max_consultas_endemismo,
        usar_cache_disco=usar_cache_disco, ruta_cache_disco=ruta_cache_disco,
        notas_nom059=cargar_notas_nom059(),
    )

    # Cruce con deforestación (Hansen) -- opcional, ver
    # _descargar_lossyear_para_biodiversidad()/generar_biodiversidad()
    # (--con-deforestacion). Siempre se llama, incluso si
    # lossyear_alineado es None, para que las columnas CERCA_DEFORESTACION/
    # ANIO_DEFORESTACION_CERCANA existan siempre (en None/NaN cuando no se
    # evaluó) -- así el resto del código (CSV, tabla, mapas) no necesita
    # ramas "si existe la columna" para esto en particular.
    df_bio = _marcar_deforestacion_en_puntos(df_bio, hidrologia, utm_crs, lossyear_alineado)
    hay_deforestacion = lossyear_alineado is not None
    if hay_deforestacion:
        n_cerca = int((df_bio["CERCA_DEFORESTACION"] == True).sum())  # noqa: E712
        n_evaluados = int(df_bio["CERCA_DEFORESTACION"].notna().sum())
        log(f"[DEFORESTACION] {n_cerca}/{n_evaluados} registro(s) de biodiversidad caen en un píxel Hansen "
            f"con pérdida de bosque detectada ({len(df_bio) - n_evaluados} de {len(df_bio)} registros totales "
            f"quedaron fuera de la malla visual del terreno 3D -- no evaluados).", nivel="OK")

    csv_path = os.path.join(carpeta_salida, f"biodiversidad_{id_proyecto.lower()}.csv")
    df_bio.to_csv(csv_path, index=False, encoding="utf-8-sig")
    log(f"CSV de biodiversidad ({len(df_bio)} registros, {df_bio['NOMBRE_CIENTIFICO'].nunique()} especies "
        f"únicas) guardado en: {csv_path}", nivel="OK")

    # ------------------------------------------------------------------
    # Rarefacción + Chao1 por zona (2026-09-06) -- ANTES del resumen por
    # zona para poder fusionar sus columnas ahí mismo. Sobre TODOS los
    # registros de cada zona (no solo especies en riesgo) -- responde
    # "¿alcanzó el muestreo?", pregunta distinta a "¿cuántas en riesgo hay?".
    # ------------------------------------------------------------------
    rarefaccion_resumen_por_zona = {}
    png_rarefaccion_path = None
    if calcular_rarefaccion:
        _df_curvas_rarefaccion, rarefaccion_resumen_por_zona, png_rarefaccion_path = _rarefaccion_por_zona(
            df_bio, carpeta_salida, id_proyecto, n_iter=n_iter_rarefaccion,
        )
    else:
        log("Rarefacción/Chao1 OMITIDA (--sin-rarefaccion).", nivel="WARN")

    # ------------------------------------------------------------------
    # Resumen por zona -- para las tarjetas HTML. Mismo criterio de anillo
    # EXCLUSIVO que ya usa todo el resto de la plataforma: cada registro
    # cuenta en UNA sola zona (viene así desde zonas_anillo_exclusivo()).
    # ------------------------------------------------------------------
    resumen_por_zona = []
    for nivel_gbif, nombre_zona in NOMBRES_LEGIBLES_NIVEL.items():
        sub = df_bio[df_bio.get("NIVEL", "") == nivel_gbif] if "NIVEL" in df_bio.columns else df_bio.iloc[0:0]
        n_regs = len(sub)
        n_spp = int(sub["NOMBRE_CIENTIFICO"].nunique()) if n_regs else 0
        n_riesgo = int(sub.loc[sub["NIVEL_RIESGO"].isin(NIVELES_EN_RIESGO), "NOMBRE_CIENTIFICO"].nunique()) if n_regs else 0
        n_endemicas = int(sub.loc[sub["ENDEMICO_MX"] == "SI", "NOMBRE_CIENTIFICO"].nunique()) if n_regs else 0
        # "PROBABLE_CALCULADO" = indicio de la heurística GBIF (ver
        # core/gbif_biodiversidad.py::obtener_endemismo_heuristico_gbif) --
        # A PROPÓSITO se cuenta y se muestra APARTE de las confirmadas
        # ("SI"), nunca sumado con ellas: es un indicio calculado, no una
        # confirmación oficial de catálogo.
        n_endemicas_probables = int(
            sub.loc[sub["ENDEMICO_MX"] == "PROBABLE_CALCULADO", "NOMBRE_CIENTIFICO"].nunique()
        ) if n_regs else 0
        n_deforestacion = int(
            sub.loc[sub["CERCA_DEFORESTACION"] == True, "NOMBRE_CIENTIFICO"].nunique()  # noqa: E712
        ) if n_regs and hay_deforestacion else 0
        m = metricas_por_zona_map.get(nombre_zona, {})
        fila_zona = {
            "zona": nombre_zona, "registros": n_regs, "especies": n_spp,
            "especies_en_riesgo": n_riesgo, "especies_endemicas": n_endemicas,
            "especies_endemicas_probables": n_endemicas_probables,
            "especies_en_zona_deforestada": n_deforestacion,
            "area_2d_ha": m.get("area_2d_ha"), "area_3d_ha": m.get("area_3d_ha"),
            "factor_relieve": m.get("factor_relieve"),
        }
        fila_zona.update(rarefaccion_resumen_por_zona.get(nombre_zona, {}))
        resumen_por_zona.append(fila_zona)
    df_resumen = pd.DataFrame(resumen_por_zona)
    resumen_csv_path = os.path.join(carpeta_salida, f"biodiversidad_resumen_zona_{id_proyecto.lower()}.csv")
    df_resumen.to_csv(resumen_csv_path, index=False, encoding="utf-8-sig")

    salidas = {"csv": csv_path, "resumen_csv": resumen_csv_path, "html": None}

    # ------------------------------------------------------------------
    # Mapa 3D "zonas a vigilar" -- SOLO especies en riesgo, sobre la
    # MISMA malla de terreno que ya recibe este módulo (hidrologia),
    # reusando geomatica.generar_mapa_3d() -- ningún DEM/malla propio.
    # ------------------------------------------------------------------
    fig = None
    n_puntos_riesgo = 0
    tabla_anual_deforestacion_html = ""
    if mapa_3d:
        from core import geomatica
        import pyproj

        df_riesgo = df_bio[df_bio["NIVEL_RIESGO"].isin(NIVELES_EN_RIESGO)].copy()
        df_riesgo = df_riesgo.dropna(subset=["LAT", "LON"]) if {"LAT", "LON"}.issubset(df_riesgo.columns) else df_riesgo.iloc[0:0]
        # OJO 2026-09-06 (pregunta real del usuario: "porque en el 3d hay pocos puntos y
        # abajo en la lista tenemos 74"): el mapa satelital 2D y la tabla de especies en
        # riesgo (_mapa_satelital_riesgo() más arriba, tabla_especies_riesgo_html() en
        # reportes_html.py) SIEMPRE deduplican por (NOMBRE_CIENTIFICO, LAT, LON) antes de
        # numerar/dibujar -- varios registros GBIF de la misma especie en el mismo punto
        # exacto (mismo avistamiento, o el mismo sitio visitado varias veces) cuentan una
        # sola vez. Este bloque del mapa 3D NO lo hacía -- iteraba sobre cada REGISTRO
        # crudo. Con los datos reales de Texolo (2026-09-06) eso eran 242 registros en
        # riesgo vs. 74 combinaciones únicas especie+coordenada en la tabla -- ya la mitad
        # de la diferencia era esto. Se agrega el mismo dedup aquí para que el conteo
        # coincida con lo que ya ve el usuario en la tabla/mapa satelital.
        n_registros_riesgo_crudos = len(df_riesgo)
        df_riesgo = df_riesgo.drop_duplicates(subset=["NOMBRE_CIENTIFICO", "LAT", "LON"])

        if len(df_riesgo) > 0:
            Z_smooth = hidrologia["Z_smooth"]
            pw_v, ph_v = hidrologia["pw_v"], hidrologia["ph_v"]
            transform = hidrologia["transform"]
            rows, cols = Z_smooth.shape
            inv = ~transform
            a_utm = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True).transform

            xs_km, ys_km, zs_km, colores, tamanos, textos = [], [], [], [], [], []
            lats_reales, lons_reales = [], []
            for _, r in df_riesgo.iterrows():
                x_utm, y_utm = a_utm(float(r["LON"]), float(r["LAT"]))
                col_f, row_f = inv * (x_utm, y_utm)
                row_i, col_i = int(round(row_f)), int(round(col_f))
                if not (0 <= row_i < rows and 0 <= col_i < cols) or np.isnan(Z_smooth[row_i, col_i]):
                    continue  # punto fuera de la malla visual (buffer máximo) -- se omite, no se fuerza
                xs_km.append(col_f * pw_v / 1000.0)
                ys_km.append((rows - 1 - row_f) * ph_v / 1000.0)
                # OJO 2026-09-07 (pedido por el usuario viendo el mapa real de Texolo): estos puntos
                # se "perdían" visualmente entre las manchas de deforestación -- dos causas, no una:
                # (1) misma altura exacta que construir_capa_deforestacion_3d() (ambas usaban
                # max(pw_v, ph_v)*0.6), así que a ciertos ángulos de cámara quedaban en el mismo plano
                # y se encimaban; (2) COLORES_RIESGO (rojo/naranja/amarillo por nivel) es casi la MISMA
                # paleta que el colorscale de año de deforestación (amarillo->rojo oscuro) -- a
                # propósito no se cambia esa paleta aquí (es la convención semáforo de riesgo, ya
                # usada en toda la plataforma: tarjetas, tabla, mapa 2D), así que la corrección es por
                # ALTURA + FORMA, no por color: estos puntos ahora flotan bastante más arriba del
                # terreno que la capa de deforestación (que sigue en *0.6), y usan diamante en vez de
                # círculo -- mismo símbolo que ya usa "Punto más alto real" en geomatica.py para
                # marcadores que deben distinguirse del resto a cualquier ángulo.
                zs_km.append(float(Z_smooth[row_i, col_i]) + max(pw_v, ph_v) * 3.0)
                colores.append(r["COLOR_RIESGO"])
                tamanos.append(9.5 if r["NIVEL_RIESGO"] == "CRITICO" else 7.5)
                # OJO 2026-09-06: esta capa no traía customdata -- a diferencia de la
                # capa de deforestación (construir_capa_deforestacion_3d en
                # deforestacion.py) y de geomatica.CLICK_ABRE_GOOGLE_EARTH_JS, el click
                # en un diamante de especie en riesgo no hacía NADA (sin lat/lon que
                # leer). Se guarda aquí el LAT/LON REAL (no el de dibujo, que abajo se
                # separa un poco para que los puntos no se tapen entre sí) -- mismo
                # convenio de "últimas dos posiciones = lat, lon" que ya usan las demás
                # capas clicables de la plataforma.
                lats_reales.append(float(r["LAT"]))
                lons_reales.append(float(r["LON"]))
                deforestacion_txt = (
                    f"<br>⚠ Píxel Hansen con pérdida de bosque ({int(r['ANIO_DEFORESTACION_CERCANA'])})"
                    if r.get("CERCA_DEFORESTACION") is True and pd.notna(r.get("ANIO_DEFORESTACION_CERCANA"))
                    else ""
                )
                textos.append(f"{r['NOMBRE_CIENTIFICO']}<br>NOM-059: {r.get('CATEGORIA_NOM059') or 's/d'} | "
                               f"IUCN: {r.get('CATEGORIA_IUCN') or 's/d'} | Nivel: {r['NIVEL_RIESGO']} | "
                               f"Endémica: {r.get('ENDEMICO_MX', 'NO_EVALUADO')}"
                               + (f"<br><i>{_envolver_texto_hover(r['NOTA_ENDEMISMO'])}</i>" if r.get("NOTA_ENDEMISMO") else "")
                               + deforestacion_txt)
            n_puntos_riesgo = len(xs_km)

            # OJO 2026-09-06 (segunda mitad de la misma pregunta del usuario): aun ya
            # deduplicado por especie+coordenada, la malla D8 tiene resolución de
            # ~pw_v/ph_v metros por píxel (30m típico con SRTM) -- dos puntos reales a,
            # digamos, 10-15m de distancia (normal en avistamientos de un mismo grupo/
            # sendero) redondean al MISMO row_i/col_i y quedan con exactamente el mismo
            # x,y,z: un diamante tapa al otro por completo, aunque haya varios "apilados"
            # en ese lugar. Verificado con los datos reales de Texolo: 74 puntos
            # deduplicados caían en solo ~35 píxeles distintos de la malla (un solo
            # píxel absorbía hasta 6-7 especies distintas del mismo punto de avistamiento
            # popular). Se reutiliza la misma relajación iterativa que ya separa puntos
            # encimados en los insets de zoom del mapa 2D (_separar_para_inset) -- aquí
            # en metros reales sobre el plano x,y de la escena 3D, con una distancia
            # mínima más chica que un píxel (para no "mover" un punto a otro lugar del
            # terreno, solo despegarlo visualmente de sus vecinos exactos).
            if n_puntos_riesgo > 1:
                xs_m = [v * 1000.0 for v in xs_km]
                ys_m = [v * 1000.0 for v in ys_km]
                distancia_min_m = max(pw_v, ph_v) * 0.4
                xs_m, ys_m = _separar_para_inset(xs_m, ys_m, distancia_min_m=distancia_min_m, iteraciones=60)
                xs_km = [v / 1000.0 for v in xs_m]
                ys_km = [v / 1000.0 for v in ys_m]

            capas_extra = []
            if n_puntos_riesgo > 0:
                import plotly.graph_objects as go
                capas_extra.append(go.Scatter3d(
                    x=xs_km, y=ys_km, z=zs_km, mode="markers",
                    marker=dict(size=tamanos, color=colores, opacity=0.95, symbol="diamond",
                                line=dict(width=1.4, color="black")),
                    text=textos, hovertemplate="%{text}<extra></extra>",
                    customdata=np.column_stack([lats_reales, lons_reales]),
                    name="Especies en riesgo (NOM-059/IUCN)",
                ))
                log(f"{n_puntos_riesgo} puntos de especies en riesgo colocados sobre la malla 3D "
                    f"({n_registros_riesgo_crudos} registros crudos -> {n_puntos_riesgo} especie+coordenada "
                    f"únicas; {len(df_riesgo) - n_puntos_riesgo} quedaron fuera del área visual del buffer "
                    f"máximo).", nivel="OK")
            else:
                log("Ninguna especie en riesgo cayó dentro del área visual de la malla 3D -- se genera el "
                    "mapa solo con el terreno, sin capa de riesgo.", nivel="WARN")

            # Capa de deforestación (Hansen) EN EL MISMO mapa 3D que las especies en
            # riesgo -- justo lo que se pidió: ver de un vistazo dónde caen los
            # avistamientos y dónde se perdió bosque, sobre el mismo terreno, sin
            # tener que abrir dos mapas separados y adivinar si coinciden.
            tabla_anual_deforestacion_html = ""
            if lossyear_alineado is not None:
                from core import deforestacion
                from core import reportes_html as _reportes_html_tabla_anual
                # Mismos defaults que usa construir_capa_deforestacion_3d() internamente
                # (config.DEFORESTACION_ANIO_INICIO_DEFAULT / año actual - 1) -- se resuelven
                # aquí también porque generar_desglose_anual_visual() los necesita explícitos
                # para la tabla, y así los dos (mapa + tabla) quedan seguro con el MISMO rango.
                anio_ini_resuelto = anio_inicio_deforestacion or DEFORESTACION_ANIO_INICIO_DEFAULT
                anio_fin_resuelto = anio_fin_deforestacion or (datetime.now().year - 1)
                capa_deforestacion = deforestacion.construir_capa_deforestacion_3d(
                    lossyear_alineado, hidrologia, anio_inicio_deforestacion, anio_fin_deforestacion, utm_crs=utm_crs,
                )
                if capa_deforestacion is not None:
                    capas_extra.append(capa_deforestacion)
                # 2026-09-06 (pedido por el usuario): la misma tablita año-por-año +
                # GRAN TOTAL que ya tiene el mapa 3D standalone de core.deforestacion --
                # antes solo vivía ahí, y para verla había que correr ese módulo aparte.
                # Se agrega aquí también para que UNA sola corrida de core.biodiversidad
                # (con --con-deforestacion) ya la incluya, sin depender de un segundo comando.
                try:
                    _csv_desglose, df_desglose = deforestacion.generar_desglose_anual_visual(
                        lossyear_alineado, hidrologia, anio_ini_resuelto, anio_fin_resuelto,
                        id_proyecto, carpeta_salida,
                    )
                    tabla_anual_deforestacion_html = _reportes_html_tabla_anual.tabla_desglose_anual_deforestacion_html(
                        df_desglose, anio_ini_resuelto, anio_fin_resuelto
                    )
                except Exception as e:
                    log(f"No se pudo generar la tabla anual de deforestación ({e}) -- el mapa 3D con la "
                        f"capa Hansen se genera igual, solo se omite esta tabla.", nivel="WARN")

            html_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_3d_biodiversidad.html")
            subtitulo = (f"Zonas a vigilar: {n_puntos_riesgo} puntos de especies en NOM-059/IUCN "
                         f"(crítico/alto/medio/vigilar)"
                         + (", con capa de deforestación Hansen encima -- mismo terreno." if lossyear_alineado is not None
                            else ", sobre el mismo terreno que el mapa de deforestación.")
                         + f"{subtitulo_extra}")
            titulo_base = f"{id_proyecto} -- Biodiversidad: especies en riesgo (NOM-059/IUCN)"
            fig = geomatica.generar_mapa_3d(
                hidrologia, id_proyecto, html_path, subtitulo=subtitulo, utm_crs=utm_crs,
                capas_extra=capas_extra, titulo_base=titulo_base, devolver_fig=True,
            )
        else:
            log("Sin coordenadas LAT/LON en el dataframe de biodiversidad (o ninguna especie en riesgo) -- "
                "se omite el mapa 3D.", nivel="WARN")

    # ------------------------------------------------------------------
    # Página HTML final: tarjetas por zona (core/reportes_html.py, mismo
    # estilo que carbono/carbono_perdida) + mapa 3D embebido si se generó.
    # ------------------------------------------------------------------
    from core import reportes_html

    tarjetas_html = []
    for fila in resumen_por_zona:
        zona = fila["zona"]
        color = reportes_html.COLORES_ZONA_HEX.get(zona, "#666")
        lineas_sec = [
            reportes_html.linea_secundaria_html(f"{fila['especies_en_riesgo']}", "en riesgo (NOM-059/IUCN)"),
            reportes_html.linea_secundaria_html(f"{fila['especies_endemicas']}", "endémicas confirmadas"),
        ]
        if fila.get("especies_endemicas_probables"):
            lineas_sec.append(reportes_html.linea_secundaria_html(
                f"{fila['especies_endemicas_probables']}", "posibles endémicas (calculado GBIF, sin confirmar)"
            ))
        if hay_deforestacion and fila.get("especies_en_zona_deforestada"):
            lineas_sec.append(reportes_html.linea_secundaria_html(
                f"{fila['especies_en_zona_deforestada']}", "en un píxel con pérdida de bosque (Hansen)"
            ))
        # Área 3D real (corregida por relieve, célula por célula vía DEM) al
        # lado del área 2D (la que ya va arriba de la tarjeta, `area_ha=`) --
        # la comparación es el punto: en terreno con pendiente, la superficie
        # real siempre es mayor a la proyección plana (factor_relieve >= 1.0).
        if fila.get("area_3d_ha") is not None:
            lineas_sec.append(reportes_html.linea_secundaria_html(
                f"{fila['area_3d_ha']:,.1f} ha", f"área 3D real (relieve) · factor {fila['factor_relieve']:.2f}×"
            ))
        # Chao1 (2026-09-06): estimador de riqueza mínima esperada, contando
        # especies aún no detectadas -- ver _rarefaccion_por_zona(). None si
        # la zona no tuvo suficientes registros (MIN_REGISTROS_RAREFACCION).
        if fila.get("rarefaccion_chao1_estimado") is not None:
            lineas_sec.append(reportes_html.linea_secundaria_html(
                f"≈{fila['rarefaccion_chao1_estimado']:.0f} spp",
                f"estimado Chao1 (+{fila['rarefaccion_especies_no_detectadas_est']:.0f} sin detectar aún)"
            ))
        tarjetas_html.append(reportes_html.tarjeta_html(
            zona, color, f"{fila['especies']}", "especies",
            nota_principal=f"{fila['registros']:,} registros GBIF",
            lineas_secundarias=lineas_sec, area_ha=fila.get("area_2d_ha"),
        ))

    total_spp = int(df_bio["NOMBRE_CIENTIFICO"].nunique())
    total_riesgo = int(df_bio.loc[df_bio["NIVEL_RIESGO"].isin(NIVELES_EN_RIESGO), "NOMBRE_CIENTIFICO"].nunique())
    total_endemicas = int(df_bio.loc[df_bio["ENDEMICO_MX"] == "SI", "NOMBRE_CIENTIFICO"].nunique())
    total_endemicas_probables = int(
        df_bio.loc[df_bio["ENDEMICO_MX"] == "PROBABLE_CALCULADO", "NOMBRE_CIENTIFICO"].nunique()
    )
    total_deforestacion = int(
        df_bio.loc[df_bio["CERCA_DEFORESTACION"] == True, "NOMBRE_CIENTIFICO"].nunique()  # noqa: E712
    ) if hay_deforestacion else 0
    areas_2d = [f["area_2d_ha"] for f in resumen_por_zona if f.get("area_2d_ha") is not None]
    areas_3d = [f["area_3d_ha"] for f in resumen_por_zona if f.get("area_3d_ha") is not None]
    lineas_sec_total = [reportes_html.linea_secundaria_html(f"{total_riesgo}", "especies en riesgo") ,
                         reportes_html.linea_secundaria_html(f"{total_endemicas}", "especies endémicas confirmadas")]
    if total_endemicas_probables:
        lineas_sec_total.append(reportes_html.linea_secundaria_html(
            f"{total_endemicas_probables}", "posibles endémicas (calculado GBIF, sin confirmar)"
        ))
    if hay_deforestacion and total_deforestacion:
        lineas_sec_total.append(reportes_html.linea_secundaria_html(
            f"{total_deforestacion}", "especies en algún píxel con pérdida de bosque (Hansen)"
        ))
    total_area_2d = round(sum(areas_2d), 1) if areas_2d else None
    if areas_3d and total_area_2d:
        total_area_3d = round(sum(areas_3d), 1)
        factor_total = round(total_area_3d / total_area_2d, 2) if total_area_2d else None
        lineas_sec_total.append(reportes_html.linea_secundaria_html(
            f"{total_area_3d:,.1f} ha", f"área 3D real (relieve) · factor {factor_total:.2f}×"
        ))
    nota_principal_total = f"{total_riesgo} en riesgo | {total_endemicas} endémicas confirmadas"
    if total_endemicas_probables:
        nota_principal_total += f" | {total_endemicas_probables} posibles (calculado)"
    if hay_deforestacion and total_deforestacion:
        nota_principal_total += f" | {total_deforestacion} en zona deforestada (Hansen)"
    tarjeta_total = reportes_html.tarjeta_html(
        "TOTAL", reportes_html.COLOR_TOTAL_HEX, f"{total_spp}", "especies (sitio completo)",
        nota_principal=nota_principal_total,
        lineas_secundarias=lineas_sec_total, area_ha=total_area_2d,
        nombre_mostrado="TOTAL (zonas exclusivas, sí sumable)", es_total=True,
    )

    div_mapa_3d = (
        (fig.to_html(full_html=False, include_plotlyjs=True, config={"displaylogo": False}, div_id="mapa3d")
         # 2026-09-06 (pedido por el usuario): mismo click-a-Timelapse que ya tiene
         # el mapa 3D de core/deforestacion.py -- solo dispara si el punto picado
         # es de la capa de deforestación (ver reportes_html para la nota honesta
         # sobre el formato de la liga). Si esta corrida no trajo --con-deforestacion,
         # el script queda inofensivo (no hay trazo "Deforestación" que lo dispare).
         + reportes_html.script_click_deforestacion_timelapse_html("mapa3d"))
        if fig is not None else
        "<p style='color:#666'>Mapa 3D no generado en esta corrida "
        "(sin especies en riesgo dentro del área visual, o --sin-mapa-3d).</p>"
    )

    png_2d_path, df_riesgo_2d = None, None
    if mapa_2d:
        png_2d_path, df_riesgo_2d = _mapa_satelital_riesgo(
            df_bio, geom_zonas_wgs84, id_proyecto, carpeta_salida, subtitulo_extra
        )
        salidas["png_2d"] = png_2d_path
    tabla_2d_html = (
        reportes_html.tabla_especies_riesgo_html(
            df_riesgo_2d, caption=(f"Tabla de coordenadas -- mismos números que el mapa de arriba, "
                                    f"{len(df_riesgo_2d)} especies en riesgo, para ir a campo.")
        ) if df_riesgo_2d is not None and len(df_riesgo_2d) else ""
    )
    bloque_2d = (
        f'<h2 style="margin:26px 0 8px;font-size:15px;">Mapa satelital 2D (Esri World Imagery)</h2>'
        f'<img src="{os.path.basename(png_2d_path)}" alt="Mapa satelital 2D -- zonas a vigilar" '
        f'style="max-width:100%;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.15);display:block">'
        f'{tabla_2d_html}'
        if png_2d_path else
        '<h2 style="margin:26px 0 8px;font-size:15px;">Mapa satelital 2D (Esri World Imagery)</h2>'
        '<p style="color:#666">No generado en esta corrida (falta contextily, sin internet a las teselas '
        'en el momento de correr, sin especies en riesgo con coordenadas, o --sin-mapa-2d).</p>'
    )
    tabla_rarefaccion_filas = "".join(
        f"<tr><td>{f['zona']}</td><td>{f['registros']:,}</td>"
        f"<td>{f['rarefaccion_s_observada']}</td>"
        f"<td>{f['rarefaccion_chao1_estimado']:.0f}</td>"
        f"<td>+{f['rarefaccion_especies_no_detectadas_est']:.0f}</td></tr>"
        if f.get("rarefaccion_chao1_estimado") is not None else
        f"<tr><td>{f['zona']}</td><td>{f['registros']:,}</td>"
        f"<td colspan='3' style='color:#666'>menos de {MIN_REGISTROS_RAREFACCION} registros -- omitido</td></tr>"
        for f in resumen_por_zona
    )
    bloque_rarefaccion = (
        f'<h2 style="margin:26px 0 8px;font-size:15px;">Rarefacción y Chao1 (¿alcanzó el muestreo?)</h2>'
        f'<img src="{os.path.basename(png_rarefaccion_path)}" alt="Curva de rarefacción por zona" '
        f'style="max-width:100%;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.15);display:block">'
        f'<table style="border-collapse:collapse;margin-top:10px;font-size:13px">'
        f'<thead><tr style="text-align:left;border-bottom:1px solid #ccc">'
        f'<th style="padding:4px 10px">Zona</th><th style="padding:4px 10px">Registros</th>'
        f'<th style="padding:4px 10px">Spp observadas</th><th style="padding:4px 10px">Chao1 estimado</th>'
        f'<th style="padding:4px 10px">Sin detectar (est.)</th></tr></thead>'
        f'<tbody>{tabla_rarefaccion_filas}</tbody></table>'
        if png_rarefaccion_path else
        '<h2 style="margin:26px 0 8px;font-size:15px;">Rarefacción y Chao1 (¿alcanzó el muestreo?)</h2>'
        '<p style="color:#666">No generado en esta corrida (--sin-rarefaccion, falta matplotlib, o ninguna '
        f'zona alcanzó el mínimo de {MIN_REGISTROS_RAREFACCION} registros).</p>'
    )
    # 2026-09-06 (pedido por el usuario): la tablita año-por-año + GRAN TOTAL de
    # deforestación, generada arriba (junto con el mapa 3D) cuando --con-deforestacion
    # sí trajo datos -- antes esto solo existía si se corría core.deforestacion --mapa-3d
    # como un comando aparte; ahora una sola corrida de core.biodiversidad ya la incluye.
    _titulo_bloque_defo_anual = '<h2 style="margin:26px 0 8px;font-size:15px;">Deforestación por año (Hansen)</h2>'
    if tabla_anual_deforestacion_html:
        bloque_deforestacion_anual = _titulo_bloque_defo_anual + tabla_anual_deforestacion_html
    else:
        bloque_deforestacion_anual = (
            _titulo_bloque_defo_anual
            + '<p style="color:#666">No generado en esta corrida (falta --con-deforestacion, o la descarga de '
              'Hansen falló -- ver el WARN correspondiente arriba en el log).</p>'
        )
    div_mapa = div_mapa_3d + bloque_2d + bloque_rarefaccion + bloque_deforestacion_anual

    nota_pie = (
        "Área 2D (arriba de cada tarjeta) = proyección plana del polígono; área 3D real = suma célula por "
        "célula del DEM corregida por pendiente local (área_píxel / cos(pendiente)) -- mismo cálculo que ya "
        "usa el terreno de core.analizar_sitio. En terreno con relieve, el área 3D SIEMPRE es mayor o igual "
        "a la 2D (factor de relieve >= 1.0); la diferencia es superficie real de hábitat que un mapa plano "
        "no muestra. "
        "Especies por zona en anillo EXCLUSIVO (cada registro cuenta en una sola zona, sin traslape -- "
        "mismo criterio que deforestación/carbono). NOM-059 tiene prioridad sobre IUCN por ser la norma "
        "mexicana aplicable (ver core/gbif_biodiversidad.py). La tabla NOM-059 es la lista nacional completa "
        "(2,678 especies, DOF 14/11/2019 + Fe de erratas 04/03/2020); la tabla de endemismo sigue siendo una "
        "lista corta confirmada a mano (semilla + catálogo local). 'endémicas confirmadas' cuenta SOLO esa "
        "lista corta. 'posibles endémicas (calculado GBIF)' es un indicio adicional, calculado automáticamente "
        "para las especies que no estaban en esa lista corta: si ~95% o más de los registros MUNDIALES de la "
        "especie en GBIF caen dentro de México, se marca como PROBABLE_CALCULADO -- un indicio, NO una "
        "confirmación oficial de un catálogo experto (CONABIO/Enciclovida); requiere verificación antes de "
        "usarse como dato duro para un dictamen o publicación (ver NOTA_ENDEMISMO en el CSV para el detalle "
        "del cálculo por especie, y --sin-endemismo-heuristico para desactivarlo). 'NO_EVALUADO' en endemismo "
        "significa 'sin indicio ni dato', nunca 'confirmado que no es endémica'. La ausencia de registros "
        "recientes de una especie NO significa que ya no esté presente: "
        "es sesgo de muestreo de GBIF (esfuerzo de observación desigual entre zonas/años), no evidencia de "
        "extirpación. El mapa satelital 2D usa teselas públicas Esri World Imagery (vía contextily) -- "
        "misma fecha/resolución que Esri publique al momento de generarse, no necesariamente reciente. "
        "Variedades domesticadas/cultivos conocidos (columna NOTA_RIESGO en el CSV, ej. café cultivado o "
        "aves de traspatio) NO heredan el nivel de riesgo de la población silvestre de referencia -- el "
        "registro se conserva, solo no se cuenta como 'especie en riesgo' en zonas/mapas. La tabla NOM-059 "
        "nacional a veces aplica a una subespecie/población geográfica específica (ver notas de "
        "data/nom059/tabla_nom059_nacional.csv); un CSV local con categoria_nom059='NINGUNA' para una "
        "especie puntual anula esa categoría para este sitio (ver --tabla-nom059-local). "
        + (
            "La columna/pastilla 'Deforestación' de la tabla de coordenadas (y la capa amarilla/roja del mapa "
            "3D) cruza cada registro contra la banda Hansen 'lossyear' a 30 m de resolución: 'Pérdida de "
            "bosque (año)' significa que ESE PÍXEL fue detectado con pérdida de cobertura arbórea ese año -- "
            "es una COINCIDENCIA GEOGRÁFICA, NUNCA una prueba de que la deforestación causó el riesgo de esa "
            "especie ni de que el individuo observado murió o se fue; 'No evaluado' significa que el punto "
            "cayó fuera del área visual del terreno 3D, nunca 'sin pérdida'. "
            if hay_deforestacion else
            "Esta corrida no cruzó contra deforestación (falta --con-deforestacion) -- la columna "
            "'Deforestación' de la tabla, si aparece, sale 'No evaluado' en todas las filas. "
        )
        + "Las ligas 'Mapa'/'Earth' junto a cada coordenada abren Google Maps/Earth en ese punto exacto "
        "(vista satelital de terceros, no verificada por este proyecto); la liga 'GBIF' abre el registro "
        "público de ESE avistamiento en gbif.org (puede o no traer foto, según quien lo subió). "
        "Rarefacción/Chao1 (sección de abajo): la curva compara la riqueza esperada de cada zona AL MISMO "
        "tamaño de muestra, para que una zona con más registros no parezca 'más biodiversa' solo por tener "
        "más esfuerzo de observación; Chao1 estima cuántas especies existen pero aún no se detectaron, a "
        "partir de cuántas especies solo aparecen 1 o 2 veces en los datos (mientras más 'singletons', mayor "
        "la estimación de faltantes). Esto reduce el sesgo de muestreo desigual, NO lo elimina: Chao1 asume "
        "un muestreo razonablemente aleatorio, y las observaciones de GBIF están sesgadas hacia sitios "
        "accesibles/populares (senderos, miradores) -- sigue siendo la mejor defensa disponible con estos "
        "datos, no una garantía estadística absoluta. Zonas con menos de "
        f"{MIN_REGISTROS_RAREFACCION} registros se omiten de esta sección por instrucción explícita (ver tabla). "
        "Este mapa y CSV no sustituyen una verificación de campo por biólogo certificado."
    )
    html_final = reportes_html.pagina_html_con_tarjetas(
        f"{id_proyecto} -- Biodiversidad por zona", f"{id_proyecto} -- Biodiversidad",
        f"{total_spp} especies detectadas (GBIF) | {total_riesgo} en algún nivel de riesgo NOM-059/IUCN",
        "".join(tarjetas_html) + tarjeta_total, div_mapa, nota_pie,
    )
    html_final_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_biodiversidad.html")
    with open(html_final_path, "w", encoding="utf-8") as f:
        f.write(html_final)
    log(f"Página de biodiversidad (tarjetas + mapa 3D + mapa 2D): {html_final_path}", nivel="OK")
    salidas["html"] = html_final_path
    salidas["mapa_3d_puntos_riesgo"] = n_puntos_riesgo
    salidas["resumen_por_zona"] = resumen_por_zona
    salidas["png_rarefaccion"] = png_rarefaccion_path
    return salidas


# ==============================================================================
# --- CAMINO REAL: GBIF + terreno real (Earth Engine no hace falta aquí --
#     el DEM sale de SRTM local, igual que geomatica.py/deforestacion.py) ---
# ==============================================================================
def generar_biodiversidad(geojson_path, id_proyecto, zonas_m=None, carpeta_salida=None, carpeta_srtm=None,
                           percentil_cauce=None, gbif_max_regs_por_zona=2000, consultar_iucn=True,
                           max_consultas_iucn=500, nom059_tabla_local=None, endemismo_tabla_local=None,
                           mapa_3d=True, mapa_2d=True,
                           calcular_endemismo_heuristico=True, max_consultas_endemismo=300,
                           usar_cache_disco=True, ruta_cache_disco=None,
                           con_deforestacion=False, anio_inicio_deforestacion=None,
                           anio_fin_deforestacion=None, proyecto_gee=None,
                           calcular_rarefaccion=True, n_iter_rarefaccion=200):
    """Pipeline completo para un sitio real: descarga GBIF por zona
    (anillo exclusivo), cruza NOM-059/IUCN/endemismo, y arma el mapa 3D
    "zonas a vigilar" sobre el mismo terreno que usa deforestación.

    NO SE PUDO PROBAR CONTRA GBIF/RED REAL desde donde se escribió este
    módulo (sin salida a internet, verificado). La lógica de descarga es
    LITERAL la misma que ya corre en producción en
    V19_4_HIBRIDO_CORREGIDO.py (se movió a core/gbif_biodiversidad.py sin
    cambiarla) -- lo nuevo (cruce de riesgo + capa 3D) está probado con
    --demo (ver demo() más abajo), pero la validación end-to-end contra
    un sitio real queda pendiente de correr en tu máquina."""
    from core import geomatica
    import pyproj
    from shapely.ops import transform as shp_transform
    import geopandas as gpd

    zonas_m = sorted(zonas_m if zonas_m is not None else ZONAS_ANALISIS_M)
    if len(zonas_m) > 3 or zonas_m[0] != 0:
        sys.exit("ERROR: core.biodiversidad solo soporta el patrón estándar de zonas [0, buffer1, buffer2] "
                  "(hasta 3 niveles: nucleo/buffer_500m/buffer_1000m) -- mismo límite que "
                  "zonas_anillo_exclusivo() en core/gbif_biodiversidad.py.")
    percentil_cauce = percentil_cauce if percentil_cauce is not None else PERCENTIL_CAUCE_HIDROLOGIA
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    carpeta_srtm = carpeta_srtm or CARPETA_SRTM
    os.makedirs(carpeta_salida, exist_ok=True)

    geom_utm_nucleo, dst_array, meta_utm, utm_crs = geomatica.cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        percentil_cauce, carpeta_srtm, id_proyecto,
    )

    # Hectáreas reales por zona (2D plano vs 3D corregido por relieve) --
    # mismo cálculo de core.analizar_sitio (geomatica.calcular_metricas_por_zona),
    # reusando el DEM/geom ya cargados arriba, sin descarga/cálculo aparte.
    df_metricas_zona = geomatica.calcular_metricas_por_zona(geom_utm_nucleo, dst_array, meta_utm, zonas_m)

    a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True).transform

    def _gdf_zona_wgs84(buf_m):
        geom_utm = geom_utm_nucleo.buffer(buf_m) if buf_m > 0 else geom_utm_nucleo
        geom_wgs84 = shp_transform(a_wgs84, geom_utm)
        return gpd.GeoDataFrame({"zona": ["z"]}, geometry=[geom_wgs84], crs="EPSG:4326")

    gdf_nucleo = _gdf_zona_wgs84(0)
    gdf_500 = _gdf_zona_wgs84(zonas_m[1]) if len(zonas_m) > 1 else gpd.GeoDataFrame()
    gdf_1000 = _gdf_zona_wgs84(zonas_m[2]) if len(zonas_m) > 2 else gpd.GeoDataFrame()

    # Mismos polígonos de zona (WGS84) para el mapa satelital 2D -- ya se
    # calculan arriba para la descarga GBIF, se reusan aquí sin recalcular.
    geom_zonas_wgs84 = {"nucleo": gdf_nucleo.geometry.iloc[0]}
    if len(gdf_500) > 0:
        geom_zonas_wgs84["buffer_500m"] = gdf_500.geometry.iloc[0]
    if len(gdf_1000) > 0:
        geom_zonas_wgs84["buffer_1000m"] = gdf_1000.geometry.iloc[0]

    log(f"Descargando biodiversidad GBIF por zona (núcleo + anillos exclusivos, "
        f"máx {gbif_max_regs_por_zona:,} regs/zona, todas las clases taxonómicas)...")
    df_bio = descargar_biodiversidad_gbif(
        zonas_anillo_exclusivo(gdf_nucleo, gdf_500, gdf_1000), gbif_max_regs_por_zona,
    )
    gbif_csv_path = os.path.join(carpeta_salida, f"biodiversidad_gbif_bruto_{id_proyecto.lower()}_"
                                  f"{datetime.now().strftime('%Y-%m-%d')}.csv")
    df_bio.to_csv(gbif_csv_path, index=False, encoding="utf-8-sig")
    log(f"Descarga GBIF cruda guardada (fuente trazable): {gbif_csv_path}", nivel="OK")

    tabla_nom059 = cargar_tabla_nom059(ruta_local_override=nom059_tabla_local)
    tabla_endemismo = cargar_tabla_endemismo(ruta_local=endemismo_tabla_local)

    # Cruce con deforestación (Hansen) -- opcional (--con-deforestacion),
    # reusa el MISMO terreno/hidrología ya calculados arriba para el mapa
    # 3D de riesgo, así que no descarga ni recalcula el DEM de nuevo. Si
    # falla por cualquier razón (sin Earth Engine, sin autenticar, sin
    # internet a Hansen), degrada con un WARN y sigue sin esta capa -- ver
    # _descargar_lossyear_para_biodiversidad().
    lossyear_alineado = None
    if con_deforestacion:
        lossyear_alineado = _descargar_lossyear_para_biodiversidad(
            geom_utm_nucleo, zonas_m, hidrologia, utm_crs, carpeta_srtm, proyecto_gee,
        )

    return _procesar_biodiversidad(
        df_bio, hidrologia, utm_crs, id_proyecto, carpeta_salida,
        tabla_nom059, tabla_endemismo, consultar_iucn, max_consultas_iucn, mapa_3d,
        metricas_zona=df_metricas_zona, mapa_2d=mapa_2d, geom_zonas_wgs84=geom_zonas_wgs84,
        calcular_endemismo_heuristico=calcular_endemismo_heuristico, max_consultas_endemismo=max_consultas_endemismo,
        usar_cache_disco=usar_cache_disco, ruta_cache_disco=ruta_cache_disco,
        lossyear_alineado=lossyear_alineado, anio_inicio_deforestacion=anio_inicio_deforestacion,
        anio_fin_deforestacion=anio_fin_deforestacion,
        calcular_rarefaccion=calcular_rarefaccion, n_iter_rarefaccion=n_iter_rarefaccion,
    )


# ==============================================================================
# --- DEMO: sin GBIF, sin Earth Engine, sin red -- mismo espíritu que
#     --demo en el resto del proyecto (deforestacion.py, geomatica.py,
#     carbono_perdida.py). Corre el mismo _procesar_biodiversidad() de
#     arriba, solo con un df_bio sintético y el DEM sintético ya
#     existente en geomatica._dem_sintetico(). Las especies sintéticas
#     son NOMBRES REALES de la tabla NOM-059 que ya digitalizamos y
#     verificamos (no inventados), así el cruce de riesgo se prueba con
#     categorías reales, no con datos de juguete. ---
# ==============================================================================
def demo():
    log("=== core.biodiversidad --demo (sin GBIF/Earth Engine, valores sintéticos) ===")
    from core import geomatica
    import pyproj

    # 0/500/1000 -- antes era [0, 300, 600], que no coincidía con las
    # etiquetas NUCLEO_DESGLOSE_29/BUFFER_500m/BUFFER_1000m que sí usan las
    # especies sintéticas de abajo (y todo el resto de la plataforma). Con
    # 300/600 el resumen por zona igual salía "correcto" porque el nombre de
    # zona no depende de zonas_m aquí, pero el corte real del DEM ocurría en
    # radios distintos a los que el resto del demo asume -- ahora coincide.
    zonas_m = [0, 500, 1000]
    id_proyecto = "DEMO_BIODIVERSIDAD"
    carpeta_salida = os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    os.makedirs(carpeta_salida, exist_ok=True)

    dst_array, meta_utm, geom_utm_nucleo, utm_crs = geomatica._dem_sintetico()
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        22, carpeta_salida, id_proyecto,
    )
    df_metricas_zona = geomatica.calcular_metricas_por_zona(geom_utm_nucleo, dst_array, meta_utm, zonas_m)
    rows, cols = hidrologia["Z_raw"].shape
    a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True).transform

    # Mismos anillos que ve el mapa 3D, en WGS84, para que el mapa satelital
    # 2D dibuje los mismos contornos núcleo/buffer_500m/buffer_1000m (aquí
    # sí es geometría sintética, como el resto del --demo).
    from shapely.ops import transform as shp_transform
    geom_zonas_wgs84 = {
        "nucleo": shp_transform(a_wgs84, geom_utm_nucleo),
        "buffer_500m": shp_transform(a_wgs84, geom_utm_nucleo.buffer(zonas_m[1])),
        "buffer_1000m": shp_transform(a_wgs84, geom_utm_nucleo.buffer(zonas_m[2])),
    }

    # Especies reales (nombre + categoría ya verificados hoy contra la
    # tabla nacional NOM-059, ver nom059_digitalizacion/) repartidas en
    # puntos sintéticos DENTRO de la malla del DEM sintético, una por
    # nivel/zona, para que el resumen por zona salga con datos en las 3.
    especies_demo = [
        # (nombre_cientifico, clase, familia, nivel_gbif, row, col)
        ("Cupressus guadalupensis", "Pinopsida", "Cupressaceae", "NUCLEO_DESGLOSE_29", rows // 2, cols // 2),
        ("Pristis pectinata", "Chondrichthyes", "Pristidae", "NUCLEO_DESGLOSE_29", rows // 2 - 5, cols // 2 + 4),
        ("Brachypelma hamorii", "Arachnida", "Theraphosidae", "BUFFER_500m", rows // 3, cols // 3),
        ("Halodule wrightii", "Magnoliopsida", "Cymodoceaceae", "BUFFER_500m", rows // 3 + 6, cols // 3 - 3),
        ("Geonoma interrupta", "Liliopsida", "Arecaceae", "BUFFER_1000m", 8, 8),
        ("Especie Sin Categoria Nom059", "Insecta", "Formicidae", "BUFFER_1000m", 12, 15),
    ]
    filas = []
    for nombre, clase, familia, nivel, row_i, col_i in especies_demo:
        row_i = min(max(row_i, 0), rows - 1)
        col_i = min(max(col_i, 0), cols - 1)
        x_utm, y_utm = hidrologia["transform"] * (col_i, row_i)
        lon, lat = a_wgs84(x_utm, y_utm)
        for _ in range(3):  # varios "registros" por especie, como en una descarga real
            filas.append({
                "NIVEL": nivel, "CLASE": clase, "ORDEN": "", "FAMILIA": familia,
                "NOMBRE_CIENTIFICO": nombre, "LAT": lat, "LON": lon,
                "ANIO": 2024, "BASE": "HUMAN_OBSERVATION", "GBIF_ID": f"demo-{nombre}-{row_i}-{col_i}",
            })
    df_bio_sintetico = pd.DataFrame(filas)
    print(f"\n--- Especies sintéticas del demo ({len(especies_demo)} spp, {len(df_bio_sintetico)} 'regs') ---")
    print(df_bio_sintetico[["NIVEL", "NOMBRE_CIENTIFICO", "FAMILIA"]].drop_duplicates().to_string(index=False))

    tabla_nom059 = cargar_tabla_nom059()  # tabla nacional real, sin override local
    tabla_endemismo = cargar_tabla_endemismo()

    salidas = _procesar_biodiversidad(
        df_bio_sintetico, hidrologia, utm_crs, id_proyecto, carpeta_salida,
        tabla_nom059, tabla_endemismo, consultar_iucn=False, max_consultas_iucn=0, mapa_3d=True,
        subtitulo_extra=" (DEMO -- especies reales de NOM-059, coordenadas sintéticas, sin descarga GBIF real)",
        metricas_zona=df_metricas_zona, mapa_2d=True, geom_zonas_wgs84=geom_zonas_wgs84,
        # sin red en --demo (mismo espíritu que consultar_iucn=False arriba) -- la
        # heurística de endemismo también necesita GBIF, así que se apaga aquí.
        calcular_endemismo_heuristico=False, max_consultas_endemismo=0,
        # sin caché de disco ni deforestación en --demo -- no hay nada que cachear
        # sin red, y Earth Engine tampoco aplica aquí (mismo espíritu que arriba).
        usar_cache_disco=False,
    )
    log(f"Demo completo. CSV: {salidas['csv']} | HTML: {salidas['html']} | "
        f"puntos de riesgo en el mapa 3D: {salidas['mapa_3d_puntos_riesgo']}", nivel="OK")
    return salidas


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Módulo de Biodiversidad Salamandra -- GBIF por zona + cruce NOM-059/IUCN/endemismo "
                    "+ mapa 3D de especies en riesgo, para cualquier polígono (ANP, predio, sitio Ramsar)."
    )
    ap.add_argument("--demo", action="store_true", help="Corre con valores sintéticos, sin GBIF ni Earth Engine")
    ap.add_argument("--geojson", type=str, help="Ruta al GeoJSON del polígono (núcleo -- los buffers se calculan)")
    ap.add_argument("--id-proyecto", type=str, help="Nombre identificador (para nombres de archivo)")
    ap.add_argument("--zonas", type=str, default=None, help="Buffers en metros separados por coma, ej. '0,500,1000'")
    ap.add_argument("--carpeta-salida", type=str, default=None)
    ap.add_argument("--carpeta-srtm", type=str, default=None)
    ap.add_argument("--percentil-cauce", type=float, default=None)
    ap.add_argument("--gbif-max-regs-por-zona", type=int, default=2000)
    ap.add_argument("--sin-iucn", action="store_true", help="No consulta IUCN vía GBIF (más rápido, solo cruce NOM-059)")
    ap.add_argument("--max-consultas-iucn", type=int, default=500)
    ap.add_argument("--tabla-nom059-local", type=str, default=None,
                     help="CSV local opcional (nombre_cientifico,categoria_nom059,notas) que sobreescribe/agrega "
                          "casos particulares sobre la tabla nacional")
    ap.add_argument("--tabla-endemismo-local", type=str, default=None,
                     help="CSV local opcional (nombre_cientifico,endemico_mx,nota)")
    ap.add_argument("--sin-endemismo-heuristico", action="store_true",
                     help="No calcula el indicio de endemismo vía GBIF (distribución mundial por país) para "
                          "especies fuera de la lista confirmada -- deja esas especies en NO_EVALUADO, sin "
                          "gastar consultas extra a GBIF. Ver core/gbif_biodiversidad.py::"
                          "obtener_endemismo_heuristico_gbif().")
    ap.add_argument("--max-consultas-endemismo", type=int, default=300,
                     help="Presupuesto de consultas a GBIF para el indicio de endemismo (mismo espíritu que "
                          "--max-consultas-iucn) -- default 300.")
    ap.add_argument("--sin-mapa-3d", action="store_true", help="Solo CSV + resumen, sin generar el mapa 3D")
    ap.add_argument("--sin-mapa-2d", action="store_true",
                     help="No genera el mapa satelital 2D (Esri World Imagery vía contextily). Sin esta "
                          "bandera, el mapa 2D se intenta igual y se omite solo con un aviso si falta "
                          "contextily o no hay internet a las teselas en ese momento -- nunca truena.")
    ap.add_argument("--cache-gbif", type=str, default=None,
                     help="Ruta al archivo de caché de disco compartido de IUCN/endemismo (default: "
                          "cache/gbif_cache.json dentro del repo). Especies ya resueltas en corridas "
                          "anteriores (de este sitio o de cualquier otro) no se vuelven a consultar a GBIF.")
    ap.add_argument("--sin-cache-gbif", action="store_true",
                     help="No usa ni actualiza el caché de disco de IUCN/endemismo -- cada corrida vuelve a "
                          "consultar todo desde cero (más lento, pero útil si sospechas que el caché quedó "
                          "con un dato viejo/erróneo).")
    ap.add_argument("--con-deforestacion", action="store_true",
                     help="Cruza cada registro de biodiversidad contra la banda Hansen 'lossyear' (misma "
                          "fuente que core.deforestacion) y agrega esa capa al mismo mapa 3D -- requiere "
                          "Earth Engine autenticado en esta máquina (earthengine-api instalado + "
                          "'earthengine authenticate' corrido al menos una vez). Si falla por cualquier "
                          "razón, se omite con un WARN -- el resto del reporte se genera igual.")
    ap.add_argument("--anio-inicio-deforestacion", type=int, default=None,
                     help="Primer año Hansen a considerar para el cruce de deforestación (default: mismo "
                          "default que core.deforestacion.DEFORESTACION_ANIO_INICIO_DEFAULT).")
    ap.add_argument("--anio-fin-deforestacion", type=int, default=None,
                     help="Último año Hansen a considerar (default: el año pasado).")
    ap.add_argument("--proyecto-gee", type=str, default=None,
                     help="ID del proyecto de Google Earth Engine a usar con --con-deforestacion (ee.Initialize"
                          "(project=...)) -- solo hace falta si tu cuenta de GEE lo requiere.")
    ap.add_argument("--sin-rarefaccion", action="store_true",
                     help="No calcula la curva de rarefacción ni el estimador Chao1 por zona (más rápido, "
                          "menos defendible frente a una crítica de sesgo de muestreo).")
    ap.add_argument("--n-iter-rarefaccion", type=int, default=200,
                     help="Remuestreos Monte Carlo por punto de la curva de rarefacción (default 200; sube a "
                          "999+ para una corrida final que vaya a presentación/publicación).")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not args.geojson or not args.id_proyecto:
        ap.error("--geojson y --id-proyecto son obligatorios fuera de --demo")

    zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
    generar_biodiversidad(
        geojson_path=args.geojson, id_proyecto=args.id_proyecto, zonas_m=zonas_m,
        carpeta_salida=args.carpeta_salida, carpeta_srtm=args.carpeta_srtm, percentil_cauce=args.percentil_cauce,
        gbif_max_regs_por_zona=args.gbif_max_regs_por_zona, consultar_iucn=not args.sin_iucn,
        max_consultas_iucn=args.max_consultas_iucn, nom059_tabla_local=args.tabla_nom059_local,
        endemismo_tabla_local=args.tabla_endemismo_local, mapa_3d=not args.sin_mapa_3d,
        mapa_2d=not args.sin_mapa_2d,
        calcular_endemismo_heuristico=not args.sin_endemismo_heuristico,
        max_consultas_endemismo=args.max_consultas_endemismo,
        usar_cache_disco=not args.sin_cache_gbif, ruta_cache_disco=args.cache_gbif,
        con_deforestacion=args.con_deforestacion, anio_inicio_deforestacion=args.anio_inicio_deforestacion,
        anio_fin_deforestacion=args.anio_fin_deforestacion, proyecto_gee=args.proyecto_gee,
        calcular_rarefaccion=not args.sin_rarefaccion, n_iter_rarefaccion=args.n_iter_rarefaccion,
    )


if __name__ == "__main__":
    main()

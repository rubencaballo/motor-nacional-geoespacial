#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
core/gbif_biodiversidad.py
============================

Lugar ÚNICO donde vive la descarga GBIF por zona y el cruce contra
NOM-059/IUCN/endemismo -- mismo principio que ya usa core/reportes_html.py
para las tarjetas HTML (ver su docstring): si esto se reimplementara
adentro de cada script que necesita biodiversidad, un ajuste (ej. el
reintento ante un 503 de GBIF, o una corrección a la tabla NOM-059)
habría que repetirlo en cada copia, arriesgando que se desincronicen.

POR QUÉ EXISTE ESTE MÓDULO:
    `core/V19_4_HIBRIDO_CORREGIDO.py` (el dictamen CO2+Biodiversidad para
    CONANP/SEDEMA) ya traía su propia `zonas_anillo_exclusivo()` y
    `descargar_biodiversidad_gbif()` funcionando de verdad contra GBIF.
    Cuando se diseñó el nuevo `core/biodiversidad.py` (el módulo público
    de la plataforma Salamandra, para mapas 2D/3D de especies en riesgo)
    hacía falta exactamente la misma descarga por zona -- así que en vez
    de escribirla otra vez, se saca aquí y los dos la importan.

    De regalo, este módulo agrega lo que NINGÚN script de este repo tenía
    todavía (verificado con grep exhaustivo el 2026-09-04): cruce contra
    NOM-059-SEMARNAT-2010 (tabla nacional completa, 2,678 especies,
    digitalizada del DOF 14/11/2019 + Fe de erratas 04/03/2020) e IUCN vía
    GBIF, con la misma regla de prioridad ya usada en el script viejo de
    Texolo (irdcloudbiodiversidadgbf.py, nunca integrado a este repo):
    NOM-059 manda sobre IUCN, por ser la norma mexicana aplicable.

LIMITACIÓN CONOCIDA, A PROPÓSITO NO ESCONDIDA:
    La tabla de endemismo sigue siendo una lista local mínima (4 especies
    confirmadas a mano para Texolo, heredadas del script viejo) -- NO hay
    todavía una lista nacional de endemismo equivalente a la de NOM-059.
    Construir esa lista es un trabajo de digitalización tan grande como
    el de NOM-059 y queda fuera de este módulo por ahora. `ENDEMICO_MX`
    para cualquier especie fuera de esa lista corta sale "NO_EVALUADO",
    nunca "NO" (nunca se afirma la ausencia de endemismo sin dato real).

MÉTODO NO CARO DE INDICIO DE ENDEMISMO (agregado 2026-09-06,
`obtener_endemismo_heuristico_gbif()`): mientras no exista una lista
nacional digitalizada, se aprovecha que GBIF ya sabe -- por cada
ocurrencia MUNDIAL de una especie, no solo las de este sitio -- en qué
país cayó. Si casi el 100% de los registros mundiales de una especie
caen dentro de México, es un indicio fuerte (calculado, gratis, con la
misma infraestructura que ya usa `obtener_categoria_iucn()`) de que esa
especie podría ser endémica. ESTO NO ES UNA CONFIRMACIÓN OFICIAL -- por
eso el resultado usa un estatus DISTINTO ("PROBABLE_CALCULADO"), nunca
se mezcla con "SI" (que solo sale de una fuente confirmada a mano:
CONABIO, un catálogo experto, o la semilla incrustada). Un catálogo
oficial (CONABIO/Enciclovida) sigue siendo la única fuente para pasar de
"PROBABLE_CALCULADO" a "SI" de verdad.

Dependencias: requests, pandas, geopandas (esta última solo para el
fallback bbox+recorte de descargar_biodiversidad_gbif).
"""

from __future__ import annotations

import json
import os
import re
import sys
import time

import pandas as pd
import requests

GBIF_OCCURRENCE_URL = "https://api.gbif.org/v1/occurrence/search"
GBIF_SPECIES_MATCH_URL = "https://api.gbif.org/v1/species/match"

# Ruta default a la tabla nacional -- relativa a este archivo, no al cwd
# desde donde se invoque el script, para que funcione sin importar desde
# dónde se corra `python3 -m core.biodiversidad` o
# `python3 core/V19_4_HIBRIDO_CORREGIDO.py`.
_AQUI = os.path.dirname(os.path.abspath(__file__))
TABLA_NOM059_NACIONAL_DEFAULT = os.path.join(_AQUI, "..", "data", "nom059", "tabla_nom059_nacional.csv")

# Ruta default del caché de disco de IUCN/endemismo -- ver
# cargar_cache_disco()/guardar_cache_disco() más abajo. Es UN caché
# compartido por especie, no por sitio: una especie que ya se resolvió
# corriendo Texolo sirve gratis para Cofre de Perote o cualquier otro
# sitio futuro que la vuelva a encontrar.
CACHE_GBIF_DEFAULT = os.path.join(_AQUI, "..", "cache", "gbif_cache.json")

# Semilla de endemismo heredada del script viejo de Texolo
# (irdcloudbiodiversidadgbf.py) -- ver limitación documentada arriba.
SEMILLA_ENDEMISMO_INCRUSTADA = [
    ("Microtus quasiater", "SI", "Endémico de México (SEMARNAT/CONABIO)."),
    ("Handleyomys chapmani", "SI", "Endémico de México, bosque mesófilo de montaña."),
    ("Oryzomys chapmani", "SI", "Sinónimo de Handleyomys chapmani, mismo estatus."),
    ("Lineatriton lineola", "SI",
     "Endémico, con registro histórico específico en Barranca de Teocelo (espécimen MVZ, "
     "ver Parra-Olea et al.) -- ligado a este micro-sitio, no asumir presencia en otro sin verificar."),
]

COLORES_RIESGO = {
    "CRITICO": "#7f0000",
    "ALTO": "#d62728",
    "MEDIO": "#ff7f0e",
    "VIGILAR": "#ffcc00",
    "BAJO": "#2ca02c",
    "SIN_DATO": "#7f7f7f",
}

# Semilla de especies cuya categoría NOM-059/IUCN de referencia en realidad
# describe a OTRA población -- cultivada/doméstica, o silvestre pero de una
# subespecie/región geográficamente distinta -- que NO es la que normalmente
# aparece en GBIF cerca del sitio analizado. Mismo espíritu que
# SEMILLA_ENDEMISMO_INCRUSTADA (lista corta, a mano, ampliable). Dos casos
# reales encontrados auditando corridas de Texolo:
#   - 2026-09-04: 'Coffea arabica' salía "ALTO" porque IUCN global clasifica
#     la población SILVESTRE ancestral de café (Etiopía/Sudán del Sur) como
#     "EN" por cambio climático -- pero el café detectado vía GBIF cerca de
#     Texolo es cultivo, no esa población silvestre.
#   - 2026-09-07: 'Procyon lotor' (mapache) salía "CRÍTICO" (NOM-059 "P")
#     porque la tabla nacional lista a la especie completa como en peligro
#     -- pero el propio texto de esa entrada aclara que es por el "mapache
#     de Islas Marías / Tres Marías" (subespecie insular, Nayarit), no por
#     el mapache continental común que aparece en GBIF en Veracruz (y en
#     casi cualquier sitio del país). Un mapache real avistado en Texolo no
#     es la subespecie protegida.
#   - 2026-09-06 (auditando la corrida real de Texolo con el mapache ya
#     corregido, buscando el MISMO patrón en el resto de la tabla nacional):
#     tres casos más, mismo patrón exacto -- especie común/extendida en
#     GBIF que hereda la categoría de una población geográficamente
#     restringida mencionada en el propio texto de la norma:
#       'Troglodytes aedon' (chivirín/house wren): NOM-059 "Pr" es por el
#       "chivirín saltapared DE COZUMEL" (subespecie/población endémica de
#       la isla) -- el chivirín continental es de los pájaros más comunes
#       y extendidos de América.
#       'Leptotila verreauxi' (paloma arroyera): NOM-059 "Pr" es por la
#       "paloma arroyera DE TRES MARÍAS" (población insular, Nayarit) -- la
#       paloma continental es común en casi todo México/Centroamérica.
#       'Setophaga coronata' (chipe coronado): NOM-059 "A" es por el "chipe
#       coronado GUATEMALTECO" (población residente de tierras altas de
#       Chiapas/Guatemala) -- las aves migratorias que invernan en Veracruz
#       (y en gran parte de México) son la forma común, no esa población.
#     Un grep sistemático de la tabla completa contra palabras clave
#     geográficas (Cozumel/insular/Tres Marías/guatemalteco/peninsular)
#     encontró 64 entradas con este patrón en total -- ver
#     ADVERTIR_PATRON_POBLACION_RESTRINGIDA_REGEX más abajo para un aviso
#     automático (no exclusión automática) cuando aparezca una especie que
#     matchee ese patrón y todavía no esté en esta lista -- así un sitio
#     nuevo (Cofre de Perote, Pico de Orizaba) no vuelve a agarrar esto por
#     sorpresa solo hasta la siguiente auditoría manual.
# En todos los casos esto NO oculta el registro (sigue intacto en el CSV
# completo), solo evita que el mapa/tarjetas de "especies en riesgo" trate
# a la población equivocada como si estuviera en peligro.
SEMILLA_POBLACION_DISTINTA_EXCLUIR_RIESGO = {
    "coffea arabica": "IUCN 'EN' es por la población SILVESTRE ancestral (Etiopía/Sudán del Sur, "
                       "amenazada por cambio climático) -- el café detectado vía GBIF cerca de zonas "
                       "habitadas casi siempre es cultivo, no esa población silvestre.",
    "procyon lotor": "NOM-059 'P' (peligro de extinción) es por la subespecie insular 'mapache de Islas "
                      "Marías / Tres Marías' (Nayarit) -- el propio texto de la norma lo aclara. Un mapache "
                      "avistado en tierra firme (Veracruz u otro estado continental) es la población común "
                      "y extendida, no la subespecie protegida.",
    "troglodytes aedon": "NOM-059 'Pr' es por el 'chivirín saltapared de Cozumel' (subespecie/población "
                          "endémica de la isla de Cozumel) -- el chivirín continental (Veracruz u otro "
                          "estado continental) es de los pájaros más comunes y extendidos del continente, "
                          "no la población insular protegida.",
    "leptotila verreauxi": "NOM-059 'Pr' es por la 'paloma arroyera de Tres Marías' (población insular, "
                            "Nayarit) -- la paloma continental (Veracruz u otro estado continental) es "
                            "común y extendida en casi todo México y Centroamérica, no la población insular.",
    "setophaga coronata": "NOM-059 'A' es por el 'chipe coronado guatemalteco' (población residente de "
                           "tierras altas de Chiapas/Guatemala) -- las aves migratorias que invernan en "
                           "Veracruz (y en gran parte de México) son la forma común y extendida, no esa "
                           "población residente restringida.",
}

# 2026-09-06: patrón para AVISAR (no excluir automático) cuando una especie
# en el mapa/tarjetas de riesgo tiene, en su nota de la tabla NOM-059, una
# palabra clave de población geográficamente restringida y TODAVÍA no está
# en SEMILLA_POBLACION_DISTINTA_EXCLUIR_RIESGO -- ver enriquecer_con_riesgo()
# más abajo. A propósito NO excluye solo: casos como 'Setophaga pitiayumi'
# ("MULTIPLE ENTRIES: subsp. graysoni=Pr, subsp. insularis=P") son
# ambiguos -- puede que la subespecie mexicana continental (graysoni)
# también esté protegida, no solo la insular -- eso requiere criterio
# humano, no debe auto-excluirse con la misma confianza que los 4 casos
# confirmados arriba.
PATRON_POBLACION_RESTRINGIDA_REGEX = re.compile(
    r"cozumel|\bisla\b|insular|guatemalteco|tres mar[ií]as|peninsular", re.IGNORECASE
)

_PATRON_VARIEDAD_DOMESTICADA = re.compile(r"\bvar\.?\s*domestic|\bf\.?\s*domestic|domesticus\b", re.IGNORECASE)


def _nota_exclusion_riesgo_cultivo(nombre_cientifico, nombre_canon):
    """Devuelve una nota (string) si `nombre_cientifico` es una variedad
    domesticada detectada por patrón (ej. 'Cairina moschata var. domestica')
    o está en SEMILLA_POBLACION_DISTINTA_EXCLUIR_RIESGO -- o None si no
    aplica. Ver docstring de SEMILLA_POBLACION_DISTINTA_EXCLUIR_RIESGO
    arriba para los casos reales que motivaron esto (Cairina moschata var.
    domestica heredaba el "P" de la pata real SILVESTRE nacional; Procyon
    lotor heredaba el "P" de la subespecie insular de Islas Marías)."""
    if _PATRON_VARIEDAD_DOMESTICADA.search(str(nombre_cientifico or "")):
        return ("variedad/raza domesticada (nombre incluye 'var. domestica' o equivalente) -- no se le "
                "aplica el estatus NOM-059/IUCN de la población silvestre de la especie.")
    if nombre_canon in SEMILLA_POBLACION_DISTINTA_EXCLUIR_RIESGO:
        return SEMILLA_POBLACION_DISTINTA_EXCLUIR_RIESGO[nombre_canon]
    return None


def nombre_canonico(nombre_cientifico):
    """
    GBIF suele traer el nombre con autor y año pegados, ej.
    'Amazona albifrons (Sparrman, 1788)'. Para cruzar contra la tabla
    NOM-059 (y para /species/match) nos quedamos solo con el binomio
    'genus species' en minúsculas -- MISMA función y mismo criterio que
    ya usaba irdcloudbiodiversidadgbf.py, para que una tabla NOM-059
    hecha para ese script siga funcionando igual aquí.
    """
    if not nombre_cientifico:
        return ""
    match = re.match(r"^([A-ZÁÉÍÓÚÑ][a-záéíóúñ\-]+(?:\s+[a-záéíóúñ\-]+)?)", str(nombre_cientifico).strip())
    if match:
        return match.group(1).strip().lower()
    return str(nombre_cientifico).strip().lower()


# ------------------------------------------------------------------
# Zonas en anillo exclusivo (mismo patrón que ya usa deforestacion.py/
# geomatica.py vía zona_de_pixel: cada registro cuenta en UNA sola zona,
# nunca se duplica entre núcleo/buffer_500m/buffer_1000m).
# ------------------------------------------------------------------
def zonas_anillo_exclusivo(gdf_nucleo, gdf_500=None, gdf_1000=None):
    """
    Devuelve [(NIVEL, geometria), ...] en ANILLO EXCLUSIVO (no acumulado),
    para consultar GBIF sin duplicar registros entre zonas. Los polígonos
    de entrada son discos ACUMULADOS (0-500m incluye el núcleo, 0-1000m
    incluye 0-500m) -- convención estándar de este repo (ver
    revisar_area() en V19_4_HIBRIDO_CORREGIDO.py, o zona_de_pixel en
    geomatica.py para el equivalente raster). Si se consultara GBIF con
    esos polígonos acumulados tal cual, un registro dentro del núcleo
    saldría también en las descargas de "500m" y "1000m".

    `gdf_nucleo`/`gdf_500`/`gdf_1000`: GeoDataFrames (o None/vacíos si esa
    zona no aplica). Idéntica firma lógica a la función original de
    V19_4_HIBRIDO_CORREGIDO.py, solo que ahora recibe los GeoDataFrames
    como parámetros explícitos en vez de leerlos de variables globales
    del script -- así la puede llamar cualquier módulo, no solo ese.
    """
    geom_nucleo = gdf_nucleo.geometry.union_all()
    geom_500_acum = gdf_500.geometry.union_all() if (gdf_500 is not None and len(gdf_500) > 0) else None
    geom_1000_acum = gdf_1000.geometry.union_all() if (gdf_1000 is not None and len(gdf_1000) > 0) else None
    anillo_500 = geom_500_acum.difference(geom_nucleo) if geom_500_acum is not None else None
    if geom_1000_acum is not None:
        base_1000 = geom_500_acum if geom_500_acum is not None else geom_nucleo
        anillo_1000 = geom_1000_acum.difference(base_1000)
    else:
        anillo_1000 = None
    return [
        ("NUCLEO_DESGLOSE_29", geom_nucleo),
        ("BUFFER_500m", anillo_500),
        ("BUFFER_1000m", anillo_1000),
    ]


def descargar_biodiversidad_gbif(zonas_geom, max_regs_por_zona):
    """
    Descarga ocurrencias de GBIF para cada zona, SIN filtrar por clase
    taxonómica (todas las clases -- aves, insectos, plantas, hongos,
    arácnidos, etc.). Usa geometría WKT primero; si GBIF la rechaza
    (HTTP 400, típico con polígonos complejos), cae a bbox + recorte
    local con GeoPandas.

    CLASE se toma de o.get("class") -- lo que GBIF REALMENTE reporta por
    registro -- nunca de una variable de clase que se le pidió a la API
    (ese fue el bug histórico de los scripts V16/V16.5 de Texolo, que
    corrompió CLASE durante 4 versiones sucesivas: escribían la clase
    PEDIDA, no la devuelta). Aguas abajo, de todas formas pasa por
    validate_taxonomic_class() (salamandra_biodiversidad.py) como
    segunda red de seguridad -- nunca hay que confiar en una sola capa.

    Reintenta con backoff creciente (hasta 6 veces, 10s/20s/.../60s) ante
    503/502/500/504 -- caídas cortas de servidor de GBIF, no errores del
    polígono ni de la consulta.

    Idéntica a la función que vivía dentro de V19_4_HIBRIDO_CORREGIDO.py
    (incluido el endurecimiento de reintentos del 2026-09-02) -- se movió
    aquí tal cual, sin cambiar comportamiento, para que ese script y
    core/biodiversidad.py compartan una sola copia.
    """
    todos = []
    for nombre_zona, geom in zonas_geom:
        if geom is None or geom.is_empty:
            print(f"[GBIF] {nombre_zona}: geometría vacía, se omite.")
            continue
        wkt = geom.simplify(0.002).wkt
        minx, miny, maxx, maxy = geom.bounds
        print(f"[GBIF] {nombre_zona}: consultando (bbox {minx:.4f},{miny:.4f},{maxx:.4f},{maxy:.4f})...")

        def fetch(use_geometry):
            res = []
            off = 0
            reintentos_5xx = 0
            MAX_REINTENTOS_5XX = 6
            while len(res) < max_regs_por_zona:
                params = {"hasCoordinate": "true", "hasGeospatialIssue": "false", "limit": 300, "offset": off}
                if use_geometry:
                    params["geometry"] = wkt
                else:
                    params["decimalLongitude"] = f"{minx:.4f},{maxx:.4f}"
                    params["decimalLatitude"] = f"{miny:.4f},{maxy:.4f}"
                try:
                    r = requests.get(GBIF_OCCURRENCE_URL, params=params, timeout=90)
                except Exception as e:
                    print(f"  [GBIF] {nombre_zona} EXCEPCIÓN off={off}: {e}")
                    break
                if r.status_code == 400 and use_geometry:
                    return None
                if r.status_code == 429:
                    print(f"  [GBIF] {nombre_zona} rate limit (429), esperando 20s...")
                    time.sleep(20)
                    continue
                if r.status_code in (500, 502, 503, 504):
                    reintentos_5xx += 1
                    if reintentos_5xx > MAX_REINTENTOS_5XX:
                        print(f"  [GBIF] {nombre_zona} HTTP {r.status_code} off={off}: sigue caído tras "
                              f"{MAX_REINTENTOS_5XX} reintentos -- sigo con los {len(res)} regs que ya junté "
                              f"para esta zona (si es 0, GBIF está caído de verdad, vuelve a correr en unos minutos).")
                        break
                    espera = min(10 * reintentos_5xx, 60)
                    print(f"  [GBIF] {nombre_zona} HTTP {r.status_code} off={off} "
                          f"(GBIF caído momentáneamente, reintento {reintentos_5xx}/{MAX_REINTENTOS_5XX}) -- "
                          f"esperando {espera}s...")
                    time.sleep(espera)
                    continue
                if r.status_code != 200:
                    print(f"  [GBIF] {nombre_zona} HTTP {r.status_code} off={off}: {r.text[:200]}")
                    break
                reintentos_5xx = 0
                data = r.json()
                rr = data.get("results", [])
                if not rr:
                    break
                for o in rr:
                    if o.get("decimalLatitude") is None:
                        continue
                    res.append({
                        "NIVEL": nombre_zona,
                        "CLASE": o.get("class", ""),
                        "ORDEN": o.get("order", ""),
                        "FAMILIA": o.get("family", ""),
                        "NOMBRE_CIENTIFICO": o.get("scientificName"),
                        "LAT": o.get("decimalLatitude"),
                        "LON": o.get("decimalLongitude"),
                        "ANIO": o.get("year"),
                        "BASE": o.get("basisOfRecord"),
                        "GBIF_ID": o.get("gbifID"),
                    })
                print(f"  [GBIF] {nombre_zona} off={off} +{len(rr)} -> {len(res)}/{max_regs_por_zona} "
                      f"(modo {'WKT' if use_geometry else 'BBOX'})")
                if data.get("endOfRecords") or len(rr) < 300:
                    break
                off += 300
                time.sleep(1.0)
            return res[:max_regs_por_zona]

        res = fetch(use_geometry=True)
        if res is None:
            print(f"  [GBIF] {nombre_zona}: geometría WKT rechazada por GBIF (HTTP 400) -- cayendo a bbox + recorte local.")
            res_bbox = fetch(use_geometry=False)
            if res_bbox:
                import geopandas as gpd
                df_tmp = pd.DataFrame(res_bbox)
                gdf_pts = gpd.GeoDataFrame(df_tmp, geometry=gpd.points_from_xy(df_tmp.LON, df_tmp.LAT), crs="EPSG:4326")
                dentro = gdf_pts[gdf_pts.geometry.within(geom)]
                print(f"  [GBIF] {nombre_zona}: bbox trajo {len(df_tmp)}, dentro del polígono real {len(dentro)}.")
                res = dentro.drop(columns="geometry").to_dict("records")
            else:
                res = []
        print(f"[GBIF] {nombre_zona}: {len(res)} registros finales.")
        todos.extend(res)

    if not todos:
        sys.exit("ERROR: la descarga de GBIF no trajo ningún registro para ninguna zona -- revisa tu conexión a "
                  "internet o que el polígono tenga coordenadas válidas dentro de México/el área esperada. Este "
                  "script no genera un dictamen/reporte con biodiversidad vacía.")
    df = pd.DataFrame(todos)
    antes = len(df)
    if "GBIF_ID" in df.columns:
        df = df.drop_duplicates("GBIF_ID")
    print(f"[GBIF] TOTAL descargado: {antes} regs brutos -> {len(df)} tras quitar duplicados por GBIF_ID.")
    return df


# ------------------------------------------------------------------
# NOM-059 / IUCN / endemismo -- lo nuevo de este módulo (2026-09-04).
# ------------------------------------------------------------------
def cargar_tabla_nom059(ruta_nacional=None, ruta_local_override=None):
    """
    Combina dos fuentes, en este orden de prioridad:
    1. La tabla NACIONAL (`ruta_nacional`, default: data/nom059/tabla_nom059_nacional.csv
       -- las 2,678 especies digitalizadas de NOM-059-SEMARNAT-2010, DOF
       14/11/2019 + Fe de erratas 04/03/2020) -- SIEMPRE se aplica si el
       archivo existe.
    2. Un CSV local opcional por sitio (`ruta_local_override`, columnas
       nombre_cientifico/categoria_nom059/notas -- mismo formato que la
       tabla nacional) -- para casos particulares de un sitio específico;
       si una especie está en ambos, el CSV local sobreescribe la
       categoría de la tabla nacional.

    Devuelve dict[nombre_canonico -> categoria_str] (valores "P"/"A"/"Pr"/"E").
    Si ninguna de las dos rutas existe, devuelve {} (vacío) e imprime un
    WARN -- nunca inventa categorías.
    """
    ruta_nacional = ruta_nacional or TABLA_NOM059_NACIONAL_DEFAULT
    tabla = {}
    if os.path.isfile(ruta_nacional):
        df_nac = pd.read_csv(ruta_nacional, encoding="utf-8-sig")
        for _, r in df_nac.iterrows():
            nc = nombre_canonico(r.get("nombre_cientifico", ""))
            cat = str(r.get("categoria_nom059", "")).strip()
            if nc and cat and cat.lower() != "nan":
                tabla[nc] = cat
        print(f"[NOM-059] tabla nacional cargada: {len(tabla)} especies desde {ruta_nacional}")
    else:
        print(f"[NOM-059][WARN] no encuentro la tabla nacional en {ruta_nacional} -- el cruce NOM-059 "
              f"saldrá vacío (todas las especies sin categoría) hasta que exista ese archivo.")

    if ruta_local_override and os.path.isfile(ruta_local_override):
        df_local = pd.read_csv(ruta_local_override, encoding="utf-8-sig")
        n_nuevas = 0
        n_sobreescritas = 0
        n_anuladas = 0
        for _, r in df_local.iterrows():
            nc = nombre_canonico(r.get("nombre_cientifico", ""))
            cat_raw = str(r.get("categoria_nom059", "")).strip()
            if not nc or not cat_raw or cat_raw.lower() == "nan":
                continue
            # Sentinela especial: "NINGUNA"/"NO_APLICA" en el CSV local ANULA
            # explícitamente la categoría nacional para esta especie en este
            # sitio (deja tabla[nc] = "" en vez de saltarse la fila). Caso de
            # uso real: la tabla nacional a veces lista una categoría que en
            # realidad aplica solo a una subespecie/población geográfica
            # distinta (ej. NOM-059 lista 'Procyon lotor' como 'P' -- peligro
            # de extinción -- pero la nota de esa fila aclara que es por la
            # población insular de Islas Marías/Tres Marías, no por el
            # mapache continental común que aparece en la mayoría de sitios).
            # Sin este sentinela no había forma de decirle al sistema "sí ya
            # vi la categoría nacional, y para ESTE sitio no aplica" -- una
            # fila con categoria_nom059 vacía simplemente se ignoraba (ver
            # 'if not cat_raw' arriba).
            if cat_raw.upper() in ("NINGUNA", "NO_APLICA", "N/A", "NA"):
                if nc in tabla:
                    n_anuladas += 1
                tabla[nc] = ""
                continue
            cat = cat_raw
            if nc in tabla and tabla[nc] != cat:
                n_sobreescritas += 1
            elif nc not in tabla:
                n_nuevas += 1
            tabla[nc] = cat
        print(f"[NOM-059] override local aplicado desde {ruta_local_override}: "
              f"{n_nuevas} especies nuevas, {n_sobreescritas} categorías sobreescritas, "
              f"{n_anuladas} categorías nacionales anuladas para este sitio.")
    return tabla


def cargar_notas_nom059(ruta_nacional=None):
    """Devuelve dict[nombre_canonico -> notas] con el texto de la columna
    'notas' de la tabla nacional (familia + nombres comunes) -- separado
    de cargar_tabla_nom059() porque esa función solo devuelve la categoría
    (P/A/Pr/E), no el texto, y aquí SÍ hace falta el texto: ver
    enriquecer_con_riesgo()/PATRON_POBLACION_RESTRINGIDA_REGEX arriba, el
    aviso de "posible población geográficamente restringida" necesita leer
    la nota para buscar palabras como 'Cozumel'/'Tres Marías'/'insular'."""
    ruta_nacional = ruta_nacional or TABLA_NOM059_NACIONAL_DEFAULT
    notas = {}
    if os.path.isfile(ruta_nacional):
        df_nac = pd.read_csv(ruta_nacional, encoding="utf-8-sig")
        for _, r in df_nac.iterrows():
            nc = nombre_canonico(r.get("nombre_cientifico", ""))
            nota = str(r.get("notas", "")).strip()
            if nc and nota and nota.lower() != "nan":
                notas[nc] = nota
    return notas


def cargar_tabla_endemismo(ruta_local=None):
    """
    Semilla incrustada (heredada de irdcloudbiodiversidadgbf.py, ver
    limitación documentada en el docstring del módulo) + CSV local
    opcional (columnas nombre_cientifico/endemico_mx/nota) que puede
    agregar o sobreescribir. Devuelve dict[nombre_canonico -> (estatus, nota)].
    """
    tabla = {nombre_canonico(n): (estatus, nota) for n, estatus, nota in SEMILLA_ENDEMISMO_INCRUSTADA}
    if ruta_local and os.path.isfile(ruta_local):
        df_local = pd.read_csv(ruta_local, encoding="utf-8-sig")
        for _, r in df_local.iterrows():
            nc = nombre_canonico(r.get("nombre_cientifico", ""))
            if nc:
                tabla[nc] = (str(r.get("endemico_mx", "")).strip(), str(r.get("nota", "")).strip())
        print(f"[ENDEMISMO] override local aplicado desde {ruta_local} -- tabla total: {len(tabla)} especies.")
    return tabla


_cache_iucn = {}
_contador_consultas_iucn = 0


def obtener_categoria_iucn(scientific_name, max_consultas=500):
    """
    GBIF /species/match -> usageKey -> /species/{usageKey}/iucnRedListCategory,
    lee el campo `code` corto (CR/EN/VU/NT/LC/DD/...), no `category` (el
    nombre largo). Cachea en memoria por nombre_canonico. Presupuesto
    global `max_consultas` para no saturar GBIF en sitios con muchas
    especies -- una vez agotado, devuelve None para el resto sin avisar
    especie por especie (mismo comportamiento que el script viejo de
    Texolo).
    """
    global _contador_consultas_iucn
    nc = nombre_canonico(scientific_name)
    if not nc:
        return None
    if nc in _cache_iucn:
        return _cache_iucn[nc]
    if _contador_consultas_iucn >= max_consultas:
        return None
    _contador_consultas_iucn += 1
    try:
        r_match = requests.get(GBIF_SPECIES_MATCH_URL, params={"name": scientific_name}, timeout=20)
        usage_key = r_match.json().get("usageKey") if r_match.status_code == 200 else None
        if not usage_key:
            _cache_iucn[nc] = None
            return None
        r_iucn = requests.get(f"https://api.gbif.org/v1/species/{usage_key}/iucnRedListCategory", timeout=20)
        codigo = r_iucn.json().get("code") if r_iucn.status_code == 200 else None
        _cache_iucn[nc] = codigo
        return codigo
    except Exception:
        _cache_iucn[nc] = None
        return None


_cache_endemismo_gbif = {}
_contador_consultas_endemismo_gbif = 0

# Umbral y muestra mínima -- ver docstring de obtener_endemismo_heuristico_gbif()
# para el razonamiento completo de por qué existen estos dos números.
UMBRAL_PCT_MEXICO_ENDEMISMO_DEFAULT = 95.0
MINIMO_REGISTROS_ENDEMISMO_DEFAULT = 8


def obtener_endemismo_heuristico_gbif(scientific_name, max_consultas=300,
                                       umbral_pct_mexico=UMBRAL_PCT_MEXICO_ENDEMISMO_DEFAULT,
                                       minimo_registros=MINIMO_REGISTROS_ENDEMISMO_DEFAULT):
    """
    "Método no caro" de indicio de endemismo (discutido en el chat el
    2026-09-06): en vez de digitalizar una lista nacional de endemismo
    desde cero (trabajo del tamaño del de NOM-059), se pregunta a GBIF en
    qué países cayó CADA registro MUNDIAL de esta especie (no solo los de
    este sitio) vía `occurrence/search?facet=country&limit=0` -- una sola
    llamada de conteo agrupado por país, no se descarga ninguna ocurrencia
    individual. Mismo patrón de caché/presupuesto que
    obtener_categoria_iucn(): species/match -> usageKey, cachea en
    memoria por nombre_canonico, presupuesto global `max_consultas`.

    Devuelve SIEMPRE una tupla (estatus, nota, pct_mexico, total_mundial):

    - Sin red / sin match / presupuesto agotado / GBIF no contestó 200:
      (None, "", None, None) -- el llamador debe tratar esto exactamente
      igual que "no se pudo evaluar" (NO_EVALUADO), nunca como negativo.
    - Menos de `minimo_registros` registros mundiales en GBIF: sigue
      ("NO_EVALUADO", nota de muestra insuficiente, pct, total) -- con
      pocos registros mundiales un "100% en México" no es confiable (GBIF
      simplemente no tiene datos de otros países, no significa que la
      especie no viva ahí).
    - >= `umbral_pct_mexico`% de los registros mundiales caen en México:
      ("PROBABLE_CALCULADO", nota con el % y el total, pct, total) --
      INDICIO calculado, a propósito un estatus distinto de "SI"
      confirmado (ver LIMITACIÓN/MÉTODO NO CARO en el docstring del
      módulo). Nunca se debe mostrar igual que una especie ya confirmada
      a mano por un catálogo oficial.
    - Cualquier otro caso (especie con presencia real fuera de México):
      ("NO_EVALUADO", nota con el % real para trazabilidad, pct, total)
      -- sigue sin afirmar "NO es endémica", solo "esta heurística no
      encontró indicio con este método".

    GBIF usa códigos ISO 3166-1 alfa-2 para el campo country (ej. "MX")
    tanto para filtrar como para las cuentas de un facet -- se revisan
    varias formas ("MX"/"MEXICO"/"Mexico") por si la API cambia de
    convención entre versiones, para no fallar en silencio por un
    detalle de formato.
    """
    global _contador_consultas_endemismo_gbif
    nc = nombre_canonico(scientific_name)
    if not nc:
        return None, "", None, None
    if nc in _cache_endemismo_gbif:
        return _cache_endemismo_gbif[nc]
    if _contador_consultas_endemismo_gbif >= max_consultas:
        return None, "", None, None
    _contador_consultas_endemismo_gbif += 1
    try:
        r_match = requests.get(GBIF_SPECIES_MATCH_URL, params={"name": scientific_name}, timeout=20)
        usage_key = r_match.json().get("usageKey") if r_match.status_code == 200 else None
        if not usage_key:
            resultado = (None, "", None, None)
            _cache_endemismo_gbif[nc] = resultado
            return resultado
        r_fac = requests.get(GBIF_OCCURRENCE_URL, params={
            "taxonKey": usage_key, "facet": "country", "facetLimit": 300, "limit": 0,
        }, timeout=30)
        if r_fac.status_code != 200:
            resultado = (None, "", None, None)
            _cache_endemismo_gbif[nc] = resultado
            return resultado
        data = r_fac.json()
        conteos_pais = {}
        for f in data.get("facets", []):
            if f.get("field", "").upper() == "COUNTRY":
                for c in f.get("counts", []):
                    conteos_pais[str(c.get("name", "")).strip().upper()] = int(c.get("count", 0))
        total_mundial = sum(conteos_pais.values())
        n_mexico = conteos_pais.get("MX", 0) or conteos_pais.get("MEXICO", 0)
        pct_mexico = (100.0 * n_mexico / total_mundial) if total_mundial else None

        if not total_mundial:
            resultado = ("NO_EVALUADO", "GBIF no devolvió registros mundiales para calcular distribución por país.",
                         pct_mexico, total_mundial)
        elif total_mundial < minimo_registros:
            resultado = ("NO_EVALUADO",
                         f"Solo {total_mundial} registro(s) mundiales en GBIF -- muestra insuficiente para un "
                         f"indicio de endemismo confiable (mínimo {minimo_registros}).",
                         pct_mexico, total_mundial)
        elif pct_mexico >= umbral_pct_mexico:
            resultado = ("PROBABLE_CALCULADO",
                         f"Heurística GBIF: {pct_mexico:.0f}% de {total_mundial} registros mundiales caen dentro "
                         f"de México -- indicio de endemismo calculado, NO es una confirmación oficial "
                         f"(CONABIO/catálogo experto). Verificar antes de tratarlo como dato duro.",
                         pct_mexico, total_mundial)
        else:
            resultado = ("NO_EVALUADO",
                         f"Heurística GBIF: {pct_mexico:.0f}% de {total_mundial} registros mundiales caen dentro "
                         f"de México -- por debajo del umbral ({umbral_pct_mexico:.0f}%), sin indicio de "
                         f"endemismo por este método.",
                         pct_mexico, total_mundial)
        _cache_endemismo_gbif[nc] = resultado
        return resultado
    except Exception:
        resultado = (None, "", None, None)
        _cache_endemismo_gbif[nc] = resultado
        return resultado


def cargar_cache_disco(ruta=None):
    """Carga a memoria (_cache_iucn/_cache_endemismo_gbif) los resultados de
    IUCN/endemismo ya consultados en corridas ANTERIORES -- de este sitio o
    de cualquier otro, el caché es por especie, no por sitio. Sin esto,
    CADA corrida completa le vuelve a preguntar a GBIF absolutamente todo
    desde cero, aunque ya se le haya preguntado ayer por las mismas
    especies y la respuesta no cambie (la categoría IUCN o la distribución
    mundial de una especie no cambia de un día para otro) -- motivado por
    quejas reales del usuario de corridas repetidas lentas mientras
    depurábamos el cruce de endemismo (2026-09-06).

    Solo importa (y solo se guarda, ver guardar_cache_disco()) resultados
    ÚTILES -- nunca un None (fallo de red / sin match) -- así un hipo de
    conexión de una corrida no se vuelve un "no se pudo evaluar"
    permanente para esa especie; la siguiente corrida la vuelve a
    intentar sola. No falla si el archivo no existe o está corrupto --
    en ese caso arranca vacío, exactamente como si no hubiera caché."""
    global _cache_iucn, _cache_endemismo_gbif
    ruta = ruta or CACHE_GBIF_DEFAULT
    if not os.path.isfile(ruta):
        return
    try:
        with open(ruta, "r", encoding="utf-8") as f:
            data = json.load(f)
        n_iucn = len(data.get("iucn", {}))
        _cache_iucn.update(data.get("iucn", {}))
        n_endemismo = 0
        for nc, valor in data.get("endemismo", {}).items():
            if isinstance(valor, list) and len(valor) == 4:
                _cache_endemismo_gbif[nc] = tuple(valor)
                n_endemismo += 1
        print(f"[CACHE GBIF] cargado desde {ruta}: {n_iucn} especie(s) con IUCN ya conocida, "
              f"{n_endemismo} especie(s) con indicio de endemismo ya conocido -- esas NO se vuelven "
              f"a consultar a GBIF en esta corrida.")
    except Exception as e:
        print(f"[CACHE GBIF][WARN] no se pudo leer {ruta} ({e}) -- se arranca sin caché de disco, "
              f"sin afectar el resultado de esta corrida, solo tardará más de lo normal.")


def guardar_cache_disco(ruta=None):
    """Guarda a disco lo que se sabe en memoria (_cache_iucn/
    _cache_endemismo_gbif) al terminar una corrida, para que la SIGUIENTE
    (de este sitio o de cualquier otro) no tenga que volver a preguntarle
    a GBIF por las mismas especies -- ver cargar_cache_disco() para el
    razonamiento completo. Solo persiste resultados útiles, nunca un None
    (fallo de red o especie sin match en GBIF)."""
    ruta = ruta or CACHE_GBIF_DEFAULT
    try:
        os.makedirs(os.path.dirname(ruta), exist_ok=True)
        data = {
            "iucn": {nc: v for nc, v in _cache_iucn.items() if v},
            "endemismo": {nc: list(v) for nc, v in _cache_endemismo_gbif.items() if v and v[0] is not None},
        }
        with open(ruta, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2, sort_keys=True)
        print(f"[CACHE GBIF] guardado en {ruta}: {len(data['iucn'])} especie(s) IUCN, "
              f"{len(data['endemismo'])} especie(s) con indicio de endemismo -- disponibles sin "
              f"volver a consultar GBIF en la próxima corrida (de este sitio o de cualquier otro).")
    except Exception as e:
        print(f"[CACHE GBIF][WARN] no se pudo guardar el caché en {ruta} ({e}) -- no afecta el "
              f"resultado de esta corrida, solo que la siguiente no se beneficiará del ahorro.")


def asignar_nivel_riesgo(categoria_nom059, categoria_iucn):
    """
    Combina NOM-059 (prioridad, por ser la norma mexicana aplicable) con
    IUCN como respaldo. Devuelve (nivel, color_hex). Mismo criterio ya
    validado en el script viejo de Texolo (irdcloudbiodiversidadgbf.py).
    """
    nom = (categoria_nom059 or "").strip().upper()
    iucn = (categoria_iucn or "").strip().upper()
    if nom == "P" or iucn in ("CR", "EW", "EX"):
        return "CRITICO", COLORES_RIESGO["CRITICO"]
    if nom == "A" or iucn == "EN":
        return "ALTO", COLORES_RIESGO["ALTO"]
    if nom == "PR" or iucn == "VU":
        return "MEDIO", COLORES_RIESGO["MEDIO"]
    if iucn == "NT":
        return "VIGILAR", COLORES_RIESGO["VIGILAR"]
    if iucn == "LC":
        return "BAJO", COLORES_RIESGO["BAJO"]
    return "SIN_DATO", COLORES_RIESGO["SIN_DATO"]


def enriquecer_con_riesgo(df, col_especie="NOMBRE_CIENTIFICO", tabla_nom059=None, tabla_endemismo=None,
                           consultar_iucn=True, max_consultas_iucn=500,
                           calcular_endemismo_heuristico=True, max_consultas_endemismo=300,
                           usar_cache_disco=True, ruta_cache_disco=None, notas_nom059=None):
    """
    Agrega columnas CATEGORIA_NOM059, CATEGORIA_IUCN, NIVEL_RIESGO,
    COLOR_RIESGO, ENDEMICO_MX, NOTA_ENDEMISMO a `df` (una fila por
    registro/ocurrencia o por especie única, funciona igual en ambos
    casos), cruzando `col_especie` contra las tablas cargadas con
    cargar_tabla_nom059()/cargar_tabla_endemismo().

    `consultar_iucn=False` salta las consultas de red a GBIF por
    completo (útil para --demo o para no gastar presupuesto de
    consultas si solo interesa el cruce NOM-059). Las consultas IUCN se
    hacen UNA VEZ por especie única (no por registro/ocurrencia), para
    no repetir la misma consulta cientos de veces en un polígono con
    muchos registros de la misma especie.

    `calcular_endemismo_heuristico=True` hace lo mismo pero para el
    "método no caro" de endemismo (ver obtener_endemismo_heuristico_gbif()
    arriba): SOLO para especies que no están ya en `tabla_endemismo`
    (semilla confirmada + CSV local por sitio) -- una especie ya
    confirmada a mano ("SI") o ya evaluada explícitamente en el CSV local
    nunca se toca ni se sobreescribe con el resultado de la heurística.
    `max_consultas_endemismo` es el mismo tipo de presupuesto global que
    `max_consultas_iucn`, para no saturar GBIF en sitios con muchas
    especies.

    `usar_cache_disco=True` (default) carga, ANTES de consultar nada, el
    caché de disco compartido (ver cargar_cache_disco()) con especies ya
    resueltas en corridas anteriores -- de este sitio o de cualquier
    otro -- y lo actualiza en disco al terminar (ver guardar_cache_disco()),
    para que la siguiente corrida no vuelva a preguntarle a GBIF por las
    mismas especies. `ruta_cache_disco=None` usa CACHE_GBIF_DEFAULT.

    Devuelve una copia de `df` (no modifica in-place).
    """
    if usar_cache_disco:
        cargar_cache_disco(ruta_cache_disco)

    tabla_nom059 = tabla_nom059 if tabla_nom059 is not None else {}
    tabla_endemismo = tabla_endemismo if tabla_endemismo is not None else {}
    notas_nom059 = notas_nom059 if notas_nom059 is not None else {}
    out = df.copy()

    especies_unicas = out[col_especie].dropna().unique().tolist()
    cat_nom059_por_especie = {}
    cat_iucn_por_especie = {}
    endemico_por_especie = {}
    nota_endemismo_por_especie = {}
    nota_riesgo_por_especie = {}
    n_excluidas_cultivo = 0
    n_endemismo_confirmado = 0
    n_endemismo_calculado = 0
    especies_a_revisar_poblacion_restringida = []  # ver PATRON_POBLACION_RESTRINGIDA_REGEX arriba
    for nombre in especies_unicas:
        nc = nombre_canonico(nombre)
        nota_cultivo = _nota_exclusion_riesgo_cultivo(nombre, nc)
        if nota_cultivo:
            # Variedad domesticada/cultivo conocido -- NO se le aplica la
            # categoría NOM-059/IUCN de la población silvestre de referencia
            # (ver SEMILLA_POBLACION_DISTINTA_EXCLUIR_RIESGO arriba). El registro sigue
            # intacto en el CSV, solo sale SIN_DATO en vez de heredar un
            # nivel de riesgo que describe a otra población.
            cat_nom059_por_especie[nombre] = ""
            cat_iucn_por_especie[nombre] = ""
            nota_riesgo_por_especie[nombre] = nota_cultivo
            n_excluidas_cultivo += 1
        else:
            cat_nom059_por_especie[nombre] = tabla_nom059.get(nc, "")
            cat_iucn_por_especie[nombre] = (
                obtener_categoria_iucn(nombre, max_consultas=max_consultas_iucn) if consultar_iucn else None
            ) or ""
            nota_riesgo_por_especie[nombre] = ""
            # 2026-09-06: la especie SÍ tiene una categoría NOM-059 asignada
            # (no está en la lista de exclusión confirmada) -- pero si la
            # nota de esa fila de la tabla nacional menciona una población
            # geográficamente restringida (Cozumel/Tres Marías/insular/
            # guatemalteco/peninsular), es candidata al MISMO patrón que ya
            # se confirmó 4 veces (mapache, chivirín, paloma, chipe) -- se
            # avisa (no se excluye solo) para que un humano decida, mismo
            # espíritu que el WARN de endemismo heurístico "no es
            # confirmación oficial" más abajo.
            nota_tabla = notas_nom059.get(nc, "")
            if (cat_nom059_por_especie[nombre] and nota_tabla
                    and PATRON_POBLACION_RESTRINGIDA_REGEX.search(nota_tabla)):
                especies_a_revisar_poblacion_restringida.append((nombre, cat_nom059_por_especie[nombre], nota_tabla))

        entrada_endemismo = tabla_endemismo.get(nc)  # None = no está en semilla ni en CSV local
        if entrada_endemismo is not None:
            estatus_endemismo, nota_endemismo = entrada_endemismo
            if estatus_endemismo == "SI":
                n_endemismo_confirmado += 1
        elif calcular_endemismo_heuristico:
            estatus_calc, nota_calc, _pct, _total = obtener_endemismo_heuristico_gbif(
                nombre, max_consultas=max_consultas_endemismo,
            )
            if estatus_calc is None:
                estatus_endemismo, nota_endemismo = "NO_EVALUADO", ""
            else:
                estatus_endemismo, nota_endemismo = estatus_calc, nota_calc
                if estatus_calc == "PROBABLE_CALCULADO":
                    n_endemismo_calculado += 1
        else:
            estatus_endemismo, nota_endemismo = "NO_EVALUADO", ""
        endemico_por_especie[nombre] = estatus_endemismo
        nota_endemismo_por_especie[nombre] = nota_endemismo

    out["CATEGORIA_NOM059"] = out[col_especie].map(cat_nom059_por_especie)
    out["CATEGORIA_IUCN"] = out[col_especie].map(cat_iucn_por_especie)
    out["ENDEMICO_MX"] = out[col_especie].map(endemico_por_especie)
    out["NOTA_ENDEMISMO"] = out[col_especie].map(nota_endemismo_por_especie)
    out["NOTA_RIESGO"] = out[col_especie].map(nota_riesgo_por_especie)
    niveles_colores = out.apply(
        lambda r: asignar_nivel_riesgo(r["CATEGORIA_NOM059"], r["CATEGORIA_IUCN"]), axis=1
    )
    out["NIVEL_RIESGO"] = niveles_colores.map(lambda t: t[0])
    out["COLOR_RIESGO"] = niveles_colores.map(lambda t: t[1])

    n_con_nom059 = sum(1 for v in cat_nom059_por_especie.values() if v)
    n_en_riesgo = sum(1 for n in out["NIVEL_RIESGO"] if n in ("CRITICO", "ALTO", "MEDIO"))
    print(f"[RIESGO] {len(especies_unicas)} especies únicas evaluadas: {n_con_nom059} con categoría NOM-059, "
          f"{sum(1 for v in cat_iucn_por_especie.values() if v)} con categoría IUCN, "
          f"{n_excluidas_cultivo} excluidas por describir otra población (cultivo/doméstico, o subespecie "
          f"geográfica distinta -- ver NOTA_RIESGO). {n_en_riesgo} registros en algún nivel de riesgo "
          f"(CRITICO/ALTO/MEDIO).")
    print(f"[ENDEMISMO] {n_endemismo_confirmado} especie(s) endémica(s) CONFIRMADA(S) (SI, catálogo/semilla), "
          f"{n_endemismo_calculado} especie(s) con indicio PROBABLE_CALCULADO (heurística GBIF -- distribución "
          f"mundial, no es confirmación oficial{' -- desactivada con calcular_endemismo_heuristico=False' if not calcular_endemismo_heuristico else ''}).")

    if especies_a_revisar_poblacion_restringida:
        print(f"[RIESGO][WARN] {len(especies_a_revisar_poblacion_restringida)} especie(s) con categoría NOM-059 "
              f"cuya nota en la tabla nacional menciona una población geográficamente restringida (Cozumel/Tres "
              f"Marías/insular/guatemalteco/peninsular) -- MISMO patrón ya confirmado 4 veces en este proyecto "
              f"(mapache, chivirín, paloma arroyera, chipe coronado: la categoría describe una subespecie/"
              f"población distinta, no la que normalmente aparece en GBIF en tierra firme). Revisar a mano si "
              f"aplica agregarlas a SEMILLA_POBLACION_DISTINTA_EXCLUIR_RIESGO -- NO se excluyen solas, esto es "
              f"solo un aviso:")
        for nombre, cat, nota in especies_a_revisar_poblacion_restringida:
            print(f"    - {nombre}: NOM-059='{cat}' | nota: {nota}")

    if usar_cache_disco:
        guardar_cache_disco(ruta_cache_disco)

    return out

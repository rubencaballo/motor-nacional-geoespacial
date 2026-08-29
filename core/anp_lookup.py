#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Busca un Área Natural Protegida por nombre dentro del shapefile nacional
de CONANP (232 ANPs, descargado una sola vez) y exporta su polígono como
GeoJSON individual -- para no tener que repetir manualmente el proceso de
"entrar al portal, descargar, filtrar con geopandas" cada vez que se
quiere analizar un sitio nuevo (ej. Pico de Orizaba, después de Cofre de
Perote).

FUENTE DE DATOS:
    Shapefile nacional de CONANP, descargado del portal:
        https://sig.conanp.gob.mx/Shape
    ("Áreas Naturales Protegidas", el paquete completo -- no hace falta
    volver a descargarlo por cada sitio nuevo, ya trae las 232 ANPs).

QUÉ HACE:
    - Busca por nombre (insensible a mayúsculas y acentos) en la columna
      NOMBRE del shapefile.
    - Si hay una sola coincidencia, exporta su geometría a GeoJSON
      (reproyectada a WGS84) y devuelve la ruta.
    - Si hay varias coincidencias, las lista (con NUM_ANP para desambiguar)
      en vez de adivinar cuál quiso decir el usuario -- mejor pedir que
      sea más específico que exportar el ANP equivocado.
    - Si no hay ninguna, sugiere los nombres más parecidos (por si hay un
      error de tecleo o el nombre oficial es distinto al esperado).

QUÉ NO HACE:
    - No descarga el shapefile nacional -- eso sigue siendo manual, una
      sola vez (ver FUENTE DE DATOS arriba). Si el archivo no está en la
      ruta esperada, este módulo avisa dónde conseguirlo, no intenta
      descargarlo solo (el portal de CONANP no tiene una URL de descarga
      directa estable para automatizar).
"""

import argparse
import difflib
import os
import re
import unicodedata

import geopandas as gpd
from shapely.geometry import Polygon

from config import ANP_SHAPEFILE_NACIONAL, log


def _normalizar(texto):
    """Quita acentos y pasa a minúsculas, para comparar nombres sin que
    un acento de más o de menos cause un falso 'no encontrado'."""
    return unicodedata.normalize("NFKD", str(texto)).encode("ascii", "ignore").decode().lower()


def _nombre_archivo_seguro(nombre):
    """Convierte el nombre del ANP en un nombre de archivo válido:
    sin acentos, sin espacios ni símbolos raros."""
    limpio = _normalizar(nombre)
    limpio = re.sub(r"[^a-z0-9]+", "_", limpio).strip("_")
    return limpio.title().replace("_", "_")  # se deja con guiones bajos, estilo Cofre_De_Perote


def buscar_anp(nombre_busqueda=None, num_anp=None, shapefile_nacional=None, carpeta_salida=None):
    """Busca un ANP por nombre (o por NUM_ANP exacto, para desambiguar) en
    el shapefile nacional de CONANP. Devuelve:
        - (geojson_path, fila) si hay una sola coincidencia -- ya exportado.
        - (None, DataFrame_de_candidatos) si hay varias coincidencias.
        - (None, None) si no hay ninguna (revisa log de sugerencias).
    """
    shapefile_nacional = shapefile_nacional or ANP_SHAPEFILE_NACIONAL
    if not os.path.exists(shapefile_nacional):
        raise FileNotFoundError(
            f"No se encontró el shapefile nacional de CONANP en: {shapefile_nacional}\n"
            f"Descárgalo una sola vez desde https://sig.conanp.gob.mx/Shape "
            f"(paquete 'Áreas Naturales Protegidas') y ajusta la ruta con --shapefile-nacional "
            f"o la variable de entorno IRDCLOUD_ANP_SHAPEFILE."
        )

    gdf = gpd.read_file(shapefile_nacional)
    log(f"Shapefile nacional cargado: {len(gdf)} ANPs")

    if num_anp is not None:
        match = gdf[gdf["NUM_ANP"] == num_anp]
    else:
        if not nombre_busqueda:
            raise ValueError("Se requiere --buscar (nombre) o --num-anp (identificador exacto).")
        objetivo = _normalizar(nombre_busqueda)
        match = gdf[gdf["NOMBRE"].apply(_normalizar).str.contains(objetivo, na=False)]

    if len(match) == 0:
        nombres = gdf["NOMBRE"].tolist()
        sugerencias = difflib.get_close_matches(nombre_busqueda or "", nombres, n=5, cutoff=0.4)
        log(f"No se encontró ningún ANP que coincida con '{nombre_busqueda}'.", nivel="WARN")
        if sugerencias:
            log("¿Quisiste decir alguno de estos?", nivel="WARN")
            for s in sugerencias:
                log(f"  - {s}", nivel="WARN")
        return None, None

    if len(match) > 1:
        log(f"Se encontraron {len(match)} ANPs que coinciden con '{nombre_busqueda}' -- "
            f"sé más específico o usa --num-anp con el NUM_ANP exacto:", nivel="WARN")
        for _, fila in match.iterrows():
            log(f"  NUM_ANP={fila['NUM_ANP']}  |  {fila['NOMBRE']}  |  {fila['ESTADOS']}  |  {fila['SUPERFICIE']:.1f} ha", nivel="WARN")
        return None, match[["NUM_ANP", "NOMBRE", "ESTADOS", "SUPERFICIE"]]

    fila = match.iloc[0]
    log(f"Encontrado: {fila['NOMBRE']} (NUM_ANP={fila['NUM_ANP']}, {fila['ESTADOS']}, {fila['SUPERFICIE']:.1f} ha)")

    carpeta_salida = carpeta_salida or os.getcwd()
    os.makedirs(carpeta_salida, exist_ok=True)
    nombre_archivo = _nombre_archivo_seguro(fila["NOMBRE"])
    geojson_path = os.path.join(carpeta_salida, f"{nombre_archivo}.geojson")

    gdf_uno = gpd.GeoDataFrame(geometry=[fila.geometry], crs=gdf.crs).to_crs("EPSG:4326")
    gdf_uno.to_file(geojson_path, driver="GeoJSON")
    log(f"GeoJSON exportado: {geojson_path}")

    return geojson_path, fila


# ==============================================================================
# --- MODO DEMO: shapefile sintético en memoria, sin depender del archivo real ---
# ==============================================================================
def demo():
    """Prueba la lógica de búsqueda (coincidencia única, ambigua, y sin
    resultados) sobre un GeoDataFrame sintético de 3 ANPs inventadas, sin
    depender del shapefile nacional real."""
    log("=== core.anp_lookup --demo (ANPs sintéticas, sin archivo real) ===")

    poligono = Polygon([(-97.0, 19.0), (-97.0, 19.1), (-96.9, 19.1), (-96.9, 19.0)])
    gdf_demo = gpd.GeoDataFrame({
        "NUM_ANP": [901, 902, 903],
        "NOMBRE": ["Parque Nacional Demo Uno", "Parque Nacional Demo Dos", "Reserva Sintética Ejemplo"],
        "ESTADOS": ["Veracruz", "Veracruz", "Oaxaca"],
        "SUPERFICIE": [123.4, 567.8, 999.0],
        "geometry": [poligono, poligono, poligono],
    }, crs="EPSG:4326")

    import tempfile
    with tempfile.TemporaryDirectory() as tmp_shp_dir, tempfile.TemporaryDirectory() as tmp_salida:
        shp_demo_path = os.path.join(tmp_shp_dir, "anp_demo.shp")
        gdf_demo.to_file(shp_demo_path)

        print("\n--- Caso 1: coincidencia única ('Demo Uno') ---")
        path1, fila1 = buscar_anp("Demo Uno", shapefile_nacional=shp_demo_path, carpeta_salida=tmp_salida)
        print(f"  -> geojson generado: {path1 is not None} | nombre encontrado: {fila1['NOMBRE'] if fila1 is not None else None}")

        print("\n--- Caso 2: coincidencia ambigua ('Demo') ---")
        path2, candidatos2 = buscar_anp("Demo", shapefile_nacional=shp_demo_path, carpeta_salida=tmp_salida)
        print(f"  -> geojson generado: {path2 is not None} (debe ser False/None) | candidatos: {len(candidatos2) if candidatos2 is not None else 0}")

        print("\n--- Caso 3: sin coincidencias ('Inexistente XYZ') ---")
        path3, res3 = buscar_anp("Inexistente XYZ", shapefile_nacional=shp_demo_path, carpeta_salida=tmp_salida)
        print(f"  -> geojson generado: {path3 is not None} (debe ser False/None)")


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(description="Busca un ANP por nombre en el shapefile nacional de CONANP -- Motor Nacional")
    ap.add_argument("--demo", action="store_true", help="Corre con ANPs sintéticas, sin depender del shapefile real")
    ap.add_argument("--buscar", type=str, help="Nombre (o parte del nombre) del ANP a buscar, ej. 'Pico de Orizaba'")
    ap.add_argument("--num-anp", type=int, default=None, help="NUM_ANP exacto, para desambiguar si --buscar da varios resultados")
    ap.add_argument("--shapefile-nacional", type=str, default=None, help="Ruta al shapefile nacional de CONANP")
    ap.add_argument("--salida", type=str, default=None, help="Carpeta donde guardar el .geojson (default: carpeta actual)")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not args.buscar and args.num_anp is None:
        ap.error("Se requiere --buscar 'nombre' o --num-anp (fuera de --demo)")

    buscar_anp(
        nombre_busqueda=args.buscar, num_anp=args.num_anp,
        shapefile_nacional=args.shapefile_nacional, carpeta_salida=args.salida,
    )


if __name__ == "__main__":
    main()

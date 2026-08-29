#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Orquestador: corre el análisis completo de un sitio con un solo comando,
en vez de tener que ejecutar core.geomatica, core.carbono y
core.validacion_hidrologica por separado, copiando rutas de archivo entre
uno y otro a mano.

Acepta CUALQUIER polígono, no solo Áreas Naturales Protegidas:
    --anp "Pico de Orizaba"    -> lo busca en el shapefile nacional de
                                   CONANP (core.anp_lookup) y arma el geojson.
    --geojson mi_finca.geojson -> usa directamente el geojson que le des
                                   (para predios privados, fincas de
                                   pequeños o medianos productores, etc.
                                   que no están en el catálogo de ANPs).

Qué genera, siempre que se le den los datos necesarios para cada parte
(cada pieza es opcional según qué argumentos se pasen -- ver más abajo):
    1. Terreno (siempre): area_2d/3d, factor de relieve, pendiente, mapa 3D
       con coordenadas reales (lat/lon) al pasar el mouse.
    2. Carbono/CO2e (si se da --proyecto-gee): biomasa vía Earth Engine,
       CO2e por zona (acumulado Y por anillo, sumable sin doble conteo),
       mapa 3D con el CO2e en el título, y el CSV combinado
       "resumen_terreno_y_carbono" (terreno + carbono en un solo archivo).
    3. Validación hidrológica (si se da --shapefile-inegi): % de cauces D8
       que coinciden con la Red Hidrográfica oficial de INEGI, CSV y mapa
       3D verde/rojo/azul.

Si no se dan --proyecto-gee ni --shapefile-inegi, solo corre el paso 1
(terreno) -- no truena, simplemente no genera lo que no tiene con qué.
"""

import argparse
import os

from config import ZONAS_ANALISIS_M, PERCENTIL_CAUCE_HIDROLOGIA, log


def analizar_sitio(id_proyecto=None, anp=None, geojson_path=None, num_anp=None,
                    shapefile_nacional=None, shapefile_inegi=None, proyecto_gee=None,
                    zonas_m=None, percentil_cauce=None, carpeta_salida=None, layer="corriente_ag_l"):
    """Corre el pipeline completo sobre un sitio. Devuelve un dict con las
    rutas de todo lo que se generó (algunas claves pueden faltar si no se
    dieron los datos necesarios para esa parte)."""
    from core import geomatica

    zonas_m = zonas_m if zonas_m is not None else ZONAS_ANALISIS_M
    percentil_cauce = percentil_cauce if percentil_cauce is not None else PERCENTIL_CAUCE_HIDROLOGIA
    resultados = {}

    # --- Paso 0: conseguir el geojson (búsqueda de ANP, o el que ya se dio) ---
    if geojson_path:
        log(f"Usando geojson dado directamente: {geojson_path}")
    elif anp:
        from core.anp_lookup import buscar_anp
        log(f"Buscando ANP: '{anp}'...")
        geojson_path, fila = buscar_anp(
            nombre_busqueda=anp, num_anp=num_anp,
            shapefile_nacional=shapefile_nacional, carpeta_salida=carpeta_salida or os.getcwd(),
        )
        if geojson_path is None:
            raise ValueError(
                "No se pudo resolver el ANP a un único geojson -- revisa el log de arriba "
                "(o no hubo coincidencias, o hubo varias y hace falta --num-anp para desambiguar)."
            )
        if not id_proyecto:
            id_proyecto = fila["NOMBRE"].replace(" ", "_")
    else:
        raise ValueError("Se requiere --anp (buscar por nombre) o --geojson (polígono propio).")

    if not id_proyecto:
        id_proyecto = os.path.splitext(os.path.basename(geojson_path))[0]

    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    os.makedirs(carpeta_salida, exist_ok=True)
    log(f"=== Analizando: {id_proyecto} (salida en {carpeta_salida}) ===")

    # --- Paso 1: terreno (siempre) ---
    log("\n--- [1/3] Terreno (geomatica.py) ---")
    df_zonas, csv_geo, html_geo = geomatica.procesar_sitio_real(
        geojson_path=geojson_path, id_proyecto=id_proyecto, zonas_m=zonas_m,
        percentil_cauce=percentil_cauce, carpeta_salida=carpeta_salida,
    )
    resultados["terreno_csv"] = csv_geo
    resultados["terreno_mapa_3d"] = html_geo

    # --- Paso 2: carbono (solo si hay proyecto de Earth Engine) ---
    if proyecto_gee:
        log("\n--- [2/3] Carbono/CO2e (carbono.py, vía Earth Engine) ---")
        from core import carbono
        df_carbono, csv_carbono = carbono.procesar_sitio_real(
            geojson_path=geojson_path, id_proyecto=id_proyecto, zonas_m=zonas_m,
            carpeta_salida=carpeta_salida, proyecto_gee=proyecto_gee,
        )
        resultados["carbono_csv"] = csv_carbono

        html_carbono = carbono.generar_mapa_3d_con_carbono(
            geojson_path=geojson_path, id_proyecto=id_proyecto, zonas_m=zonas_m,
            percentil_cauce=percentil_cauce, carpeta_salida=carpeta_salida,
        )
        resultados["carbono_mapa_3d"] = html_carbono
        resultados["resumen_terreno_y_carbono_csv"] = os.path.join(
            carpeta_salida, f"resumen_terreno_y_carbono_{id_proyecto.lower()}.csv"
        )
    else:
        log("\n--- [2/3] Carbono: SALTADO (no se dio --proyecto-gee) ---")

    # --- Paso 3: validación hidrológica (solo si hay shapefile de INEGI) ---
    if shapefile_inegi:
        log("\n--- [3/3] Validación hidrológica D8 vs INEGI (validacion_hidrologica.py) ---")
        from core import validacion_hidrologica as vh
        df_val, csv_val = vh.validar_sitio_real(
            geojson_path=geojson_path, shapefile_inegi_path=shapefile_inegi, id_proyecto=id_proyecto,
            zonas_m=zonas_m, percentil_cauce=percentil_cauce, carpeta_salida=carpeta_salida, layer=layer,
        )
        resultados["validacion_csv"] = csv_val

        html_val = vh.generar_mapa_3d_validacion(
            geojson_path=geojson_path, shapefile_inegi_path=shapefile_inegi, id_proyecto=id_proyecto,
            zonas_m=zonas_m, percentil_cauce=percentil_cauce, carpeta_salida=carpeta_salida, layer=layer,
        )
        resultados["validacion_mapa_3d"] = html_val
    else:
        log("\n--- [3/3] Validación hidrológica: SALTADA (no se dio --shapefile-inegi) ---")

    log("\n=== ANÁLISIS COMPLETO ===")
    for k, v in resultados.items():
        log(f"  {k}: {v}")

    return resultados


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(description="Analiza un sitio completo (terreno + carbono + validación) con un solo comando -- Motor Nacional")
    ap.add_argument("--anp", type=str, help="Nombre de un ANP a buscar automáticamente, ej. 'Pico de Orizaba'")
    ap.add_argument("--geojson", type=str, help="Ruta a un geojson propio (predio, finca -- no necesita ser un ANP)")
    ap.add_argument("--num-anp", type=int, default=None, help="NUM_ANP exacto, si --anp da varias coincidencias")
    ap.add_argument("--id-proyecto", type=str, default=None, help="Nombre identificador (default: se deriva del ANP o del geojson)")
    ap.add_argument("--shapefile-nacional", type=str, default=None, help="Shapefile nacional de CONANP (solo si se usa --anp)")
    ap.add_argument("--shapefile-inegi", type=str, default=None, help="Si se da, corre también la validación D8 vs INEGI")
    ap.add_argument("--proyecto-gee", type=str, default=None, help="Si se da, corre también carbono/CO2e vía Earth Engine")
    ap.add_argument("--zonas", type=str, default=None, help="Buffers en metros separados por coma, ej. '0,500,1000'")
    ap.add_argument("--percentil-cauce", type=float, default=None)
    ap.add_argument("--carpeta-salida", type=str, default=None)
    ap.add_argument("--layer", type=str, default="corriente_ag_l")
    args = ap.parse_args()

    if not args.anp and not args.geojson:
        ap.error("Se requiere --anp 'nombre' o --geojson ruta.geojson")

    zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
    analizar_sitio(
        id_proyecto=args.id_proyecto, anp=args.anp, geojson_path=args.geojson, num_anp=args.num_anp,
        shapefile_nacional=args.shapefile_nacional, shapefile_inegi=args.shapefile_inegi,
        proyecto_gee=args.proyecto_gee, zonas_m=zonas_m, percentil_cauce=args.percentil_cauce,
        carpeta_salida=args.carpeta_salida, layer=args.layer,
    )


if __name__ == "__main__":
    main()

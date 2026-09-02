#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V19.4 HÍBRIDO CORREGIDO - V17 + V19, con datos reales (sin constantes fijas)
Basado en V19_3_FINAL_HIBRIDO_V17_V19_MIX_DEFINITIVO.py, con 3 correcciones de integridad de datos:

  1) YA NO busca el geojson recorriendo todo $HOME (riesgo de cargar un polígono viejo/incorrecto
     si hay duplicados en el disco). Ahora la ruta es un argumento explícito obligatorio.

  2) YA NO usa constantes de AGB/CO2 escritas a mano en el código (AGB_NUCLEO=76.53,
     CO2_NUCLEO=65983.13, etc. -- que nunca cambiaban sin importar el polígono real).
     Ahora lee los valores reales de anillo exclusivo ("*_incremental_t", ya corregidos
     para no ser "muñeca rusa") desde el CSV resumen_terreno_y_carbono_*.csv que genera
     core/carbono.py. Además cruza el área del CSV contra el área del geojson y avisa si
     no coinciden (eso es justamente lo que puede pasar cuando el script "arrastra" un
     polígono que no es el vigente).

  3) YA NO cae en silencio a valores por defecto si falta una columna o un dato. Si algo
     no cuadra, el script se detiene con un error explícito -- este script alimenta un
     dictamen técnico-científico que se manda a CONANP/SEDEMA, así que es preferible que
     truene a que publique un número inventado.

  BONUS (encontrado al revisar el original): el texto del dictamen (Sección III, IV, V, VI
  y la tabla del PDF) tenía los conteos de especies/registros/densidades por zona escritos
  como texto literal (no eran ni siquiera variables), por lo que el "dictamen" siempre iba a
  decir "359 spp / 12,974 regs / 222 spp núcleo" sin importar los datos reales del CSV de
  biodiversidad que se cargó. Aquí esos números se calculan de verdad desde el CSV.

  También se notó una inconsistencia: el script original dice "Ramsar 1791" en todo el
  dictamen, pero el resto del proyecto (los dashboards 3D) usa "Ramsar 1601". Aquí NO se
  asume cuál es el correcto -- se pasa como argumento (--sitio-ramsar) y se imprime una
  advertencia hasta que se confirme.

  ACTUALIZACIÓN 02/09/2026 -- dos correcciones más, consolidadas aquí mismo (no repartidas
  en scripts separados, por instrucción explícita):

  4) VALIDACIÓN TAXONÓMICA DE CLASE. Se rastreó el origen de que Mammalia/Amphibia/Reptilia
     salieran con conteos absurdos en el CSV de biodiversidad: los scripts de descarga GBIF
     (V16_BBOX_FALLBACK_FINAL.py, V16_5_BBOX.py, V16_5_BBOX_FIX_SSL.py, V16_5_CORTADO.py)
     escribían "CLASE": clase (la clase que el loop le PIDIÓ a GBIF), nunca la que GBIF
     realmente devolvió por registro (o.get("class")) -- mientras que FAMILIA sí se leía bien.
     Ese bug quedó heredado sin corregir hasta esta versión. Ahora V19.4 importa
     validate_taxonomic_class() de salamandra_biodiversidad.py y cruza FAMILIA (confiable)
     contra una tabla de autoridad familia->clase para producir CLASE_VALIDADA -- esa es la
     columna que se usa para TODO conteo, mapa y texto del dictamen por clase taxonómica, no
     la CLASE cruda del CSV. Requiere que salamandra_biodiversidad.py esté en la misma carpeta.

  5) ID Y CRITERIOS RAMSAR. Se confirmó, cruzando el polígono oficial WDPA (campo INT_CRIT)
     contra RSIS y CONANP, que el Sitio Ramsar correcto es el 1601 (1791 era un error
     heredado desde V16) y que sus criterios oficiales son (ii) y (iv). Ya no se deja como
     "confirmar cuál es correcto" cuando --sitio-ramsar es 1601 (default); si se pasa otro
     ID, el dictamen sigue marcando la advertencia porque esos criterios no fueron
     verificados para un sitio distinto.

  ACTUALIZACIÓN 02/09/2026 (2ª parte, mismo día) -- --csv-biodiversidad ahora es OPCIONAL.

  6) DESCARGA GBIF INTEGRADA (para ANP/fincas nuevas, ej. Cofre de Perote, que todavía no
     tienen un CSV de biodiversidad armado a mano como Texolo). Si NO se pasa
     --csv-biodiversidad, el script descarga directo de GBIF usando el polígono de
     --geojson: TODAS las clases taxonómicas (sin filtrar por Aves/Mammalia/Amphibia/
     Reptilia como hacían los V16 -- así sale un dataset tan rico como el de Texolo,
     con insectos/plantas/hongos/arácnidos incluidos), una consulta por zona
     (núcleo/500m/1000m) usando geometría en ANILLO EXCLUSIVO (no los polígonos
     acumulados del geojson) para no duplicar un mismo registro del núcleo en las 3
     zonas -- el mismo principio de "anillo exclusivo, no muñeca rusa" que ya se usa
     para área y CO2 en este script. CLASE se toma de o.get("class") (lo que GBIF
     REALMENTE devuelve por registro), nunca de la clase que se le pidió a la API --
     eso es precisamente el bug que corrompió el CLASE de Texolo en V16/V16.5, y aquí
     se evita desde el origen (aunque de todas formas pasa después por
     validate_taxonomic_class() como segunda red de seguridad). La descarga se GUARDA
     como CSV en --out-dir (no se queda solo "en memoria") para que el dictamen cite
     una fuente trazable y reproducible, igual que si hubiera sido un CSV manual.

Uso (con CSV de biodiversidad ya existente, como Texolo):
  conda activate geo
  python V19_4_HIBRIDO_CORREGIDO.py \
      --geojson /ruta/exacta/Buffer_500_1000m.geojson \
      --csv-carbono /ruta/exacta/resumen_terreno_y_carbono_2025.csv \
      --csv-biodiversidad /ruta/exacta/Biodiversidad_DESGLOSE_29_COMPLETO.csv \
      --sitio-ramsar 1601

Uso (SIN CSV de biodiversidad -- descarga automática de GBIF, para un ANP nuevo):
  conda activate geo
  python V19_4_HIBRIDO_CORREGIDO.py \
      --geojson /ruta/exacta/Buffer_500_1000m_Perote.geojson \
      --csv-carbono /ruta/exacta/resumen_terreno_y_carbono_cofre_de_perote.csv \
      --sitio-ramsar [ID del sitio, si aplica -- Cofre de Perote no es Ramsar, puede omitirse]

Si no sabes las rutas exactas, en WSL corre:
  find ~ -iname "Buffer_500_1000m.geojson" 2>/dev/null
  find ~ -iname "resumen_terreno_y_carbono*.csv" 2>/dev/null
  find ~ -iname "*DESGLOSE_29_COMPLETO*.csv" 2>/dev/null
y copia/pega la ruta que corresponda a la corrida vigente (la más reciente que tú reconozcas,
no la que el script "adivine").
"""
import os, sys, time, argparse
import requests
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime

# ------------------------------------------------------------------
# 0) VALIDACIÓN TAXONÓMICA -- integrada desde salamandra_biodiversidad.py
#    (consolidado aquí en V19.4 en vez de repartirse entre scripts, por
#    instrucción explícita: todo entra por un solo lugar, el pipeline
#    de producción, no cada script de descarga GBIF por separado.)
# ------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from salamandra_biodiversidad import validate_taxonomic_class, FAMILY_CLASS_MAP
except ImportError:
    sys.exit(
        "ERROR: no encuentro salamandra_biodiversidad.py en la misma carpeta que este script.\n"
        "Desde el 02/09/2026 V19.4 depende de él para validar la columna CLASE contra FAMILIA -- "
        "se encontró que el bug histórico de los scripts de descarga GBIF (V16/V16.5: escribían "
        "'CLASE' como la clase PEDIDA a GBIF, no la que GBIF realmente devolvía por registro) "
        "sigue arrastrado en el CSV de biodiversidad que usa esta corrida, aun cuando FAMILIA es "
        "confiable. Copia salamandra_biodiversidad.py junto a este script y vuelve a correr."
    )

try:
    import contextily as ctx
    CTX = True
except Exception:
    CTX = False
    print("[WARN] sin contextily, sin Esri")

try:
    from reportlab.lib.pagesizes import letter
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib import colors
    from reportlab.lib.units import inch
    REPORTLAB = True
except Exception:
    REPORTLAB = False

print("=" * 90)
print("V19.4 HÍBRIDO CORREGIDO - CO2 (anillo exclusivo, real) + BIODIVERSIDAD (real)")
print("=" * 90)

# ------------------------------------------------------------------
# 1) RUTAS EXPLÍCITAS -- se acabó el rglob() por todo $HOME
# ------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--geojson", required=True, help="Ruta EXACTA al geojson de polígonos (núcleo + buffers 500/1000m)")
parser.add_argument("--csv-carbono", required=True, help="Ruta EXACTA al resumen_terreno_y_carbono_*.csv (salida de core/carbono.py, con columnas *_incremental_t)")
parser.add_argument("--csv-biodiversidad", default=None, help="Ruta EXACTA al CSV de biodiversidad (Biodiversidad_DESGLOSE_29_COMPLETO o equivalente). OPCIONAL: si no se pasa, el script descarga automáticamente de GBIF usando el polígono de --geojson (todas las clases taxonómicas, anillo exclusivo por zona) y guarda el CSV descargado en --out-dir.")
parser.add_argument("--gbif-max-regs-por-zona", type=int, default=5000, help="Límite de registros a descargar de GBIF por zona (núcleo/500m/1000m) cuando NO se pasa --csv-biodiversidad. Default 5000.")
parser.add_argument("--gbif-guardar-csv", default=None, help="Ruta donde guardar el CSV descargado de GBIF cuando NO se pasa --csv-biodiversidad (default: dentro de --out-dir, con la fecha de hoy).")
parser.add_argument("--sitio-ramsar", default="1601", help="ID del sitio Ramsar a imprimir en mapas/dictamen. El script original traía '1791' escrito a mano; el resto del proyecto usa '1601'. CONFIRMAR cuál es el correcto -- default aquí: 1601")
parser.add_argument("--ndwi", type=float, default=None, help="Valor NDWI a citar en el dictamen (el original traía -0.681 fijo, sin fuente verificable en este script). Si no se da, se marca como NO VERIFICADO.")
parser.add_argument("--spp-alto-riesgo", type=int, default=None, help="Nº de especies en alto riesgo dentro del núcleo, para el dictamen (el original traía '24' fijo, sin fuente en este script). Si no se da, se marca como NO VERIFICADO.")
parser.add_argument("--out-dir", default=None, help="Carpeta de salida (default: junto al geojson)")
parser.add_argument("--tolerancia-area-pct", type=float, default=2.0, help="Tolerancia %% para avisar si el área del geojson no coincide con la del CSV de carbono (default 2%%)")
args = parser.parse_args()

_archivos_a_checar = [("geojson", args.geojson), ("csv-carbono", args.csv_carbono)]
if args.csv_biodiversidad:
    _archivos_a_checar.append(("csv-biodiversidad", args.csv_biodiversidad))
for etiqueta, ruta in _archivos_a_checar:
    if not os.path.isfile(ruta):
        sys.exit(f"ERROR: no existe el archivo pasado en --{etiqueta}: {ruta}")

print(f"[GEOJSON] {args.geojson}")
print(f"[CSV CARBONO] {args.csv_carbono}")
if args.csv_biodiversidad:
    print(f"[CSV BIODIVERSIDAD] {args.csv_biodiversidad}")
else:
    print(f"[CSV BIODIVERSIDAD] NO se pasó -- se va a descargar automáticamente de GBIF "
          f"(máx {args.gbif_max_regs_por_zona:,} regs/zona, todas las clases taxonómicas).")

# ID y criterios Ramsar -- resuelto el 02/09/2026 cruzando WDPA (campo INT_CRIT del
# polígono oficial "Cascadas de Texolo y su entorno") contra RSIS y el sitio de CONANP:
# el sitio SÍ es el 1601 (1791 era un error que venía arrastrado desde V16), y sus
# criterios oficiales son (ii) y (iv). Ya no se deja como "confirmar cuál es correcto".
if args.sitio_ramsar == "1601":
    RAMSAR_ID_SUFFIX = " (confirmado: RSIS/CONANP, 02/09/2026)"
    RAMSAR_CRITERIOS_TXT = (
        "Criterios Ramsar del sitio: (ii) y (iv) -- confirmado cruzando WDPA (campo INT_CRIT del "
        "polígono oficial 'Cascadas de Texolo y su entorno', SITE_ID 902863) contra RSIS y CONANP."
    )
    print("[OK] Sitio Ramsar 1601 -- confirmado como ID correcto (el script original decía 1791, error "
          "histórico). Criterios oficiales: (ii);(iv).")
else:
    RAMSAR_ID_SUFFIX = " (VERIFICAR ID -- no coincide con el 1601 confirmado para 'Cascadas de Texolo')"
    RAMSAR_CRITERIOS_TXT = (
        f"[Criterios Ramsar NO verificados para --sitio-ramsar={args.sitio_ramsar} -- solo se confirmó "
        f"(ii);(iv) para el Sitio Ramsar 1601 'Cascadas de Texolo y su entorno'. Si esta corrida es para "
        f"ese sitio, usa --sitio-ramsar 1601; si es otro sitio, verifica sus criterios oficiales aparte.]"
    )
    print(f"[ALERTA] --sitio-ramsar={args.sitio_ramsar} no es el 1601 confirmado por RSIS/CONANP para "
          f"'Cascadas de Texolo y su entorno'. Verifica el ID antes de mandar el dictamen.")

# ------------------------------------------------------------------
# GEOJSON: mismo detector de columnas que el original (tipo/NIVEL, area_ha/SUP_HA)
# ------------------------------------------------------------------
gdf = gpd.read_file(args.geojson)
print(f"[GDF] {len(gdf)} features columns={list(gdf.columns)[:10]}")

col_tipo = "tipo" if "tipo" in gdf.columns else "NIVEL" if "NIVEL" in gdf.columns else gdf.columns[0]
col_area_geojson = "area_ha" if "area_ha" in gdf.columns else "SUP_HA" if "SUP_HA" in gdf.columns else None

def get_gdf(tipo_keyword):
    mask = gdf[col_tipo].astype(str).str.contains(tipo_keyword, case=False, na=False)
    return gdf[mask] if mask.any() else gpd.GeoDataFrame()

gdf_nucleo = get_gdf("Nucleo")
if len(gdf_nucleo) == 0:
    gdf_nucleo = get_gdf("DESGLOSE_29")
gdf_500 = get_gdf("500")
gdf_1000 = get_gdf("1000")

if len(gdf_nucleo) == 0:
    print("[WARN] no encontré el núcleo por nombre ('Nucleo'/'DESGLOSE_29'), uso la primera feature del geojson -- verifica que sea la correcta.")
    gdf_nucleo = gdf.iloc[[0]]

print(f"Nucleo: {len(gdf_nucleo)} | 500m: {len(gdf_500)} | 1000m: {len(gdf_1000)}")

def area_geojson(gdf_zona, nombre):
    if col_area_geojson is None or len(gdf_zona) == 0:
        return None
    try:
        return float(gdf_zona[col_area_geojson].iloc[0])
    except Exception:
        print(f"[WARN] no pude leer área de {nombre} desde el geojson")
        return None

area_geo_nucleo = area_geojson(gdf_nucleo, "núcleo")
area_geo_500 = area_geojson(gdf_500, "buffer 500m")
area_geo_1000 = area_geojson(gdf_1000, "buffer 1000m")

# ------------------------------------------------------------------
# 2) DATOS REALES DE CARBONO -- ya NO hay AGB_NUCLEO/CO2_NUCLEO escritos a mano.
#    Se leen del CSV corregido (anillo exclusivo, "*_incremental_t").
# ------------------------------------------------------------------
df_carbono = pd.read_csv(args.csv_carbono)
print(f"[CSV CARBONO] {len(df_carbono)} filas, columnas={list(df_carbono.columns)}")

def buscar_columna(df, candidatas, descripcion):
    for c in candidatas:
        if c in df.columns:
            return c
    sys.exit(
        f"ERROR: no encuentro una columna de '{descripcion}' en {args.csv_carbono}.\n"
        f"Columnas disponibles: {list(df.columns)}\n"
        f"Nombres que busqué: {candidatas}\n"
        f"Este script NO va a inventar un valor por defecto -- corrige el nombre de columna "
        f"arriba en el script (lista 'candidatas') o renombra la columna en el CSV, y vuelve a correr."
    )

col_buffer_m = buscar_columna(df_carbono, ["buffer_m"], "buffer_m (0=núcleo, 500, 1000)")
col_area_c = buscar_columna(df_carbono, ["area_ha"], "área en hectáreas")
col_agb = buscar_columna(df_carbono, ["agb_mgha", "agb_mg_ha", "agb"], "AGB (Mg/ha)")
col_co2_incr = buscar_columna(
    df_carbono,
    ["co2e_t_incremental", "co2e_incremental_t", "co2e_incremental", "co2e_t_incr", "co2e_incr_t"],
    "CO2e incremental de anillo exclusivo (t) -- la columna que genera _agregar_co2e_incremental en core/carbono.py",
)

# IMPORTANTE (encontrado comparando contra el dashboard 3D el 01/09/2026):
# la columna 'area_ha' de resumen_terreno_y_carbono_*.csv SIGUE SIENDO ACUMULATIVA
# (muñeca rusa) aunque 'co2e_incremental_t' ya esté corregida a anillo exclusivo.
# No existe una columna 'area_ha_incremental' en el CSV -- hay que derivarla aquí,
# igual que _agregar_co2e_incremental lo hace para el carbono (inc = valor - anterior,
# ordenado por buffer_m ascendente). Además la fila 'TOTAL (...)' del CSV comparte el
# mismo buffer_m que el último anillo real y siempre trae 'area_ha' vacío (NaN) --
# se descarta filtrando por área no nula, en vez de confiar en el nombre de la fila.
df_zonas = df_carbono[df_carbono[col_area_c].notna()].copy()
df_zonas = df_zonas.sort_values(col_buffer_m).reset_index(drop=True)
if len(df_zonas) == 0:
    sys.exit(f"ERROR: ninguna fila de {args.csv_carbono} tiene un valor no nulo en '{col_area_c}'.")

df_zonas["area_ha_acumulada"] = df_zonas[col_area_c]
df_zonas["area_ha_anillo"] = df_zonas["area_ha_acumulada"].diff()
df_zonas.loc[df_zonas.index[0], "area_ha_anillo"] = df_zonas["area_ha_acumulada"].iloc[0]

def fila_zona(buffer_m_val, nombre):
    sub = df_zonas[df_zonas[col_buffer_m] == buffer_m_val]
    if len(sub) == 0:
        sys.exit(
            f"ERROR: no hay fila con {col_buffer_m}={buffer_m_val} ({nombre}) en {args.csv_carbono} (con '{col_area_c}' no nulo).\n"
            f"Valores encontrados: {sorted(df_zonas[col_buffer_m].dropna().unique().tolist())}"
        )
    if len(sub) > 1:
        print(f"[WARN] {len(sub)} filas con {col_buffer_m}={buffer_m_val} ({nombre}), uso la primera")
    return sub.iloc[0]

fila_nucleo = fila_zona(0, "núcleo")
fila_500 = fila_zona(500, "buffer 500m")
fila_1000 = fila_zona(1000, "buffer 1000m")

# área de ANILLO (exclusiva) -- esta es la que se reporta y se usa para densidades
SUP_NUCLEO = float(fila_nucleo["area_ha_anillo"])
SUP_500 = float(fila_500["area_ha_anillo"])
SUP_1000 = float(fila_1000["area_ha_anillo"])

# área ACUMULADA (tal como viene en el CSV) -- se usa solo para cruzar contra el geojson,
# porque el geojson también guarda los polígonos de buffer como disco acumulado, no como anillo
SUP_NUCLEO_ACUM = float(fila_nucleo["area_ha_acumulada"])
SUP_500_ACUM = float(fila_500["area_ha_acumulada"])
SUP_1000_ACUM = float(fila_1000["area_ha_acumulada"])

AGB_NUCLEO = float(fila_nucleo[col_agb])
AGB_500 = float(fila_500[col_agb])
AGB_1000 = float(fila_1000[col_agb])

CO2_NUCLEO = float(fila_nucleo[col_co2_incr])
CO2_500 = float(fila_500[col_co2_incr])
CO2_1000 = float(fila_1000[col_co2_incr])

# anillo exclusivo -> sí es sumable directamente (a diferencia del script original,
# que sumaba columnas acumulativas disfrazadas de "variable")
CO2_TOTAL = CO2_NUCLEO + CO2_500 + CO2_1000

print(f"[ÁREA] anillo exclusivo (lo que se reporta): núcleo {SUP_NUCLEO:.2f} ha | "
      f"buf500 {SUP_500:.2f} ha | buf1000 {SUP_1000:.2f} ha | total {SUP_NUCLEO+SUP_500+SUP_1000:.2f} ha")
print(f"[ÁREA] acumulada (tal cual el CSV, solo para cruzar contra el geojson): núcleo {SUP_NUCLEO_ACUM:.2f} ha | "
      f"buf500 {SUP_500_ACUM:.2f} ha | buf1000 {SUP_1000_ACUM:.2f} ha")

# ------------------------------------------------------------------
# CORRECCIÓN DE RELIEVE (opcional): si el CSV trae área 3D (superficie real de
# terreno, no la proyección plana), se deriva el mismo anillo exclusivo que para
# el área 2D y se calcula el factor por ZONA (no acumulado -- el factor acumulado
# mezcla el núcleo empinado con el buffer más plano y exagera el efecto en los anillos
# exteriores). La hectárea DECLARADA/oficial sigue siendo la de planta (SUP_*); esto
# es solo para justificar en el dictamen por qué el CO2 por hectárea-planta sale más
# alto que un cálculo 2D ingenuo.
col_area_3d = next((c for c in ["area_3d_ha", "area_ha_3d", "area_relieve_ha"] if c in df_zonas.columns), None)
if col_area_3d:
    df_zonas["area_3d_acumulada"] = df_zonas[col_area_3d]
    df_zonas["area_3d_anillo"] = df_zonas["area_3d_acumulada"].diff()
    df_zonas.loc[df_zonas.index[0], "area_3d_anillo"] = df_zonas["area_3d_acumulada"].iloc[0]
    fila_nucleo_3d = df_zonas[df_zonas[col_buffer_m] == 0].iloc[0]
    fila_500_3d = df_zonas[df_zonas[col_buffer_m] == 500].iloc[0]
    fila_1000_3d = df_zonas[df_zonas[col_buffer_m] == 1000].iloc[0]
    AREA3D_NUCLEO = float(fila_nucleo_3d["area_3d_anillo"])
    AREA3D_500 = float(fila_500_3d["area_3d_anillo"])
    AREA3D_1000 = float(fila_1000_3d["area_3d_anillo"])
    FACTOR_RELIEVE_NUCLEO = AREA3D_NUCLEO / SUP_NUCLEO if SUP_NUCLEO else None
    FACTOR_RELIEVE_500 = AREA3D_500 / SUP_500 if SUP_500 else None
    FACTOR_RELIEVE_1000 = AREA3D_1000 / SUP_1000 if SUP_1000 else None
    print(f"[RELIEVE] superficie real de terreno por anillo (no acumulada): núcleo {AREA3D_NUCLEO:.2f} ha "
          f"(x{FACTOR_RELIEVE_NUCLEO:.4f}) | buf500 {AREA3D_500:.2f} ha (x{FACTOR_RELIEVE_500:.4f}) | "
          f"buf1000 {AREA3D_1000:.2f} ha (x{FACTOR_RELIEVE_1000:.4f})")
else:
    AREA3D_NUCLEO = AREA3D_500 = AREA3D_1000 = None
    FACTOR_RELIEVE_NUCLEO = FACTOR_RELIEVE_500 = FACTOR_RELIEVE_1000 = None
    print(f"[RELIEVE] {args.csv_carbono} no trae columna de área 3D (busqué area_3d_ha/area_ha_3d/area_relieve_ha) "
          f"-- el dictamen se genera sin el factor de corrección de relieve por anillo.")

def texto_relieve(factor, area3d, sup2d):
    if factor is None:
        return "[FACTOR DE RELIEVE NO DISPONIBLE EN EL CSV]"
    return f"x{factor:.2f} ({area3d:.1f} ha de superficie real vs {sup2d:.1f} ha en planta)"

AREA_PLANA_TOTAL = SUP_NUCLEO + SUP_500 + SUP_1000
if AREA3D_NUCLEO is not None:
    AREA3D_TOTAL = AREA3D_NUCLEO + AREA3D_500 + AREA3D_1000
    TEXTO_EXTENSION_REAL = (
        f"Extensión real total (3D): {AREA3D_TOTAL:.2f} ha, frente a {AREA_PLANA_TOTAL:.2f} ha en planta "
        f"({(AREA3D_TOTAL/AREA_PLANA_TOTAL - 1)*100:.1f}% más superficie real) -- la superficie catastral/declarada "
        f"ante CONANP se mantiene en planta por convención registral; esta cifra es la extensión física real del "
        f"terreno y es la que se usa para justificar la captura de carbono por hectárea."
    )
else:
    TEXTO_EXTENSION_REAL = "Extensión real total (3D): no disponible (el CSV de carbono no trae columna de área 3D)."

print(f"CO2 (anillo exclusivo, real): Núcleo {SUP_NUCLEO:.2f} ha = {CO2_NUCLEO:,.0f} t | "
      f"Buf500 {SUP_500:.1f} ha = {CO2_500:,.0f} t | Buf1000 {SUP_1000:.1f} ha = {CO2_1000:,.0f} t | "
      f"TOTAL {CO2_TOTAL:,.0f} t")

# Cruce geojson vs CSV de carbono: si no coinciden, es justo la señal de que se está
# mezclando un polígono con datos de otra corrida.
def revisar_area(nombre, area_geo, area_csv):
    if area_geo is None:
        return
    if area_csv == 0:
        return
    diff_pct = abs(area_geo - area_csv) / area_csv * 100
    if diff_pct > args.tolerancia_area_pct:
        print(f"[ALERTA] Área {nombre}: geojson={area_geo:.2f} ha vs CSV carbono={area_csv:.2f} ha "
              f"(difieren {diff_pct:.1f}%, tolerancia {args.tolerancia_area_pct}%). "
              f"¿Son del geojson y el CSV de la MISMA corrida? Verifica antes de usar esto en un dictamen oficial.")
    else:
        print(f"[OK] Área {nombre}: geojson y CSV carbono coinciden ({diff_pct:.1f}% de diferencia)")

revisar_area("núcleo", area_geo_nucleo, SUP_NUCLEO_ACUM)
revisar_area("buffer 500m (acumulado, disco 0-500m)", area_geo_500, SUP_500_ACUM)
revisar_area("buffer 1000m (acumulado, disco 0-1000m)", area_geo_1000, SUP_1000_ACUM)

base_dir = args.out_dir or os.path.dirname(args.geojson)
out_dir = os.path.join(base_dir, "ENTREGA_V19_4_HIBRIDO_CORREGIDO")
os.makedirs(out_dir, exist_ok=True)

# ------------------------------------------------------------------
# DESCARGA GBIF (solo si NO se pasó --csv-biodiversidad) -- ver punto 6) del
# docstring al inicio del archivo.
# ------------------------------------------------------------------
def zonas_anillo_exclusivo():
    """
    Devuelve [(NIVEL, geometria), ...] en ANILLO EXCLUSIVO (no acumulado), para
    consultar GBIF sin duplicar registros entre zonas. Los polígonos del geojson
    (gdf_nucleo/gdf_500/gdf_1000) son discos ACUMULADOS (0-500m incluye el núcleo,
    0-1000m incluye 0-500m) -- ver el comentario en revisar_area() más arriba, es
    el mismo hecho ya documentado para el área/carbono. Si se consultara GBIF con
    esos polígonos acumulados tal cual, un registro dentro del núcleo saldría
    también en las descargas de "500m" y "1000m". Aquí se resta geométricamente
    cada anillo del disco más chico, igual que ya se hace para SUP_NUCLEO/SUP_500/
    SUP_1000.
    """
    geom_nucleo = gdf_nucleo.geometry.union_all()
    geom_500_acum = gdf_500.geometry.union_all() if len(gdf_500) > 0 else None
    geom_1000_acum = gdf_1000.geometry.union_all() if len(gdf_1000) > 0 else None
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
    Descarga ocurrencias de GBIF para cada zona, SIN filtrar por clase taxonómica
    (todas las clases -- aves, insectos, plantas, hongos, arácnidos, etc., igual
    que el dataset real de Texolo, 12 clases). Usa geometría WKT primero; si GBIF
    la rechaza (HTTP 400, típico con polígonos complejos), cae a bbox + recorte
    local con GeoPandas (mismo patrón que V16_BBOX_FALLBACK_FINAL.py).

    CLASE se toma de o.get("class") -- lo que GBIF REALMENTE reporta por
    registro -- nunca de una variable de clase que se le pidió a la API. Ese fue
    exactamente el bug encontrado en los scripts V16/V16.5 (escribían la clase
    PEDIDA, no la devuelta) que corrompió el CLASE del CSV de Texolo durante 4
    versiones sucesivas. Aquí se evita desde el origen, aunque de todas formas
    aguas abajo pasa por validate_taxonomic_class() como segunda red de
    seguridad -- nunca hay que confiar en una sola capa de validación.
    """
    url = "https://api.gbif.org/v1/occurrence/search"
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
            while len(res) < max_regs_por_zona:
                params = {"hasCoordinate": "true", "hasGeospatialIssue": "false", "limit": 300, "offset": off}
                if use_geometry:
                    params["geometry"] = wkt
                else:
                    params["decimalLongitude"] = f"{minx:.4f},{maxx:.4f}"
                    params["decimalLatitude"] = f"{miny:.4f},{maxy:.4f}"
                try:
                    r = requests.get(url, params=params, timeout=90)
                except Exception as e:
                    print(f"  [GBIF] {nombre_zona} EXCEPCIÓN off={off}: {e}")
                    break
                if r.status_code == 400 and use_geometry:
                    return None
                if r.status_code == 429:
                    print(f"  [GBIF] {nombre_zona} rate limit (429), esperando 20s...")
                    time.sleep(20)
                    continue
                if r.status_code != 200:
                    print(f"  [GBIF] {nombre_zona} HTTP {r.status_code} off={off}: {r.text[:200]}")
                    break
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
                  "script no genera un dictamen con biodiversidad vacía.")

    df = pd.DataFrame(todos)
    antes = len(df)
    if "GBIF_ID" in df.columns:
        df = df.drop_duplicates("GBIF_ID")
    print(f"[GBIF] TOTAL descargado: {antes} regs brutos -> {len(df)} tras quitar duplicados por GBIF_ID.")
    return df


# ------------------------------------------------------------------
# BIODIVERSIDAD: si hay --csv-biodiversidad se lee igual que siempre; si no, se
# descarga de GBIF aquí mismo (ver arriba). En cualquier caso, TODO el texto del
# dictamen se calcula de aquí en vez de estar escrito a mano.
# ------------------------------------------------------------------
if args.csv_biodiversidad:
    print(f"[BIO] Usando CSV ya existente: {args.csv_biodiversidad} (no se descarga nada de GBIF).")
    df_bio = pd.read_csv(args.csv_biodiversidad)
    fuente_biodiversidad = args.csv_biodiversidad
else:
    print(f"[GBIF] Descargando biodiversidad para este polígono (núcleo + anillos exclusivos 500m/1000m, "
          f"todas las clases taxonómicas, límite {args.gbif_max_regs_por_zona:,} regs/zona)...")
    df_bio = descargar_biodiversidad_gbif(zonas_anillo_exclusivo(), args.gbif_max_regs_por_zona)
    fuente_biodiversidad = args.gbif_guardar_csv or os.path.join(
        out_dir, f"Biodiversidad_GBIF_{datetime.now().strftime('%Y-%m-%d')}.csv"
    )
    df_bio.to_csv(fuente_biodiversidad, index=False, encoding="utf-8-sig")
    print(f"[GBIF] {len(df_bio)} regs descargados y guardados en {fuente_biodiversidad} -- esta es la fuente "
          f"citable/trazable del dictamen, no una descarga que se quedó solo 'en memoria'.")

if "NOMBRE_CIENTIFICO" not in df_bio.columns:
    sys.exit(f"ERROR: {fuente_biodiversidad} no tiene columna NOMBRE_CIENTIFICO. Columnas: {list(df_bio.columns)}")
print(f"[BIO] {len(df_bio)} regs, {df_bio.NOMBRE_CIENTIFICO.nunique()} spp, cols={list(df_bio.columns)[:6]}")

# ------------------------------------------------------------------
# VALIDACIÓN TAXONÓMICA DE CLASE (validate_taxonomic_class, de salamandra_biodiversidad.py)
# ------------------------------------------------------------------
# ORIGEN DEL BUG (confirmado 02/09/2026, rastreado hasta los 4 scripts de descarga GBIF
# V16_BBOX_FALLBACK_FINAL.py, V16_5_BBOX.py, V16_5_BBOX_FIX_SSL.py y V16_5_CORTADO.py):
# su función bajar()/fetch_gbif() escribía "CLASE": clase -- la clase que el LOOP le pidió
# a GBIF -- en vez de o.get("class") -- la clase que GBIF realmente devolvió para ese
# registro -- mientras que FAMILIA sí se leía bien de o.get("family"). Ese error quedó
# arrastrado sin corregir en TODAS las versiones posteriores (V19 incluido) porque ninguna
# valida CLASE contra FAMILIA, solo la heredan tal cual del CSV. Antes de este parche,
# V19.4 tampoco lo hacía. A partir de aquí, todo conteo/mapa/texto por clase taxonómica
# usa CLASE_VALIDADA (FAMILIA cruzada contra la tabla de autoridad FAMILY_CLASS_MAP),
# nunca la CLASE cruda del CSV.
if "FAMILIA" in df_bio.columns:
    df_bio, reporte_clase = validate_taxonomic_class(df_bio, class_col="CLASE", family_col="FAMILIA")
    col_clase_usar = "CLASE_VALIDADA"
    print(f"[VALIDACIÓN CLASE] {reporte_clase.registros_corregidos:,}/{reporte_clase.total_registros:,} "
          f"registros ({reporte_clase.porcentaje_corregido:.1f}%) tenían CLASE declarada que NO coincide "
          f"con su FAMILIA real -- corregidos en CLASE_VALIDADA (familia -> clase, tabla de autoridad "
          f"de {len(FAMILY_CLASS_MAP)} familias en salamandra_biodiversidad.py).")
    if reporte_clase.porcentaje_corregido > 10:
        print(f"[ALERTA] >10% de discrepancia CLASE vs FAMILIA -- mismo patrón del bug histórico V16/V16.5 "
              f"(ver arriba). Todas las cifras de esta corrida por clase taxonómica usan CLASE_VALIDADA, "
              f"NO la columna CLASE original del CSV.")
    else:
        print(f"[OK] CLASE vs FAMILIA: {reporte_clase.porcentaje_corregido:.1f}% de discrepancia "
              f"(dentro de lo esperable).")
    if reporte_clase.familias_desconocidas:
        print(f"[WARN] {len(reporte_clase.familias_desconocidas)} familias no están en la tabla de "
              f"referencia de salamandra_biodiversidad.py y se quedaron con su CLASE original SIN "
              f"validar: {reporte_clase.familias_desconocidas[:15]}"
              f"{'...' if len(reporte_clase.familias_desconocidas) > 15 else ''}")
    print("[VALIDACIÓN CLASE] Detalle por clase declarada (total / corregidos / %% corregido):")
    print(reporte_clase.detalle_por_clase_declarada.to_string(index=False))
    TEXTO_VALIDACION_CLASE = (
        f"validada contra FAMILIA (tabla de autoridad familia->clase, salamandra_biodiversidad.py, "
        f"{len(FAMILY_CLASS_MAP)} familias): {reporte_clase.registros_corregidos:,}/"
        f"{reporte_clase.total_registros:,} regs ({reporte_clase.porcentaje_corregido:.1f}%) tenían CLASE "
        f"original incorrecta y se corrigieron a CLASE_VALIDADA."
    )
else:
    print(f"[ALERTA] {fuente_biodiversidad} no tiene columna FAMILIA -- NO se puede validar CLASE contra "
          f"FAMILIA (ver bug histórico V16/V16.5 documentado arriba). Se usará la columna CLASE cruda, "
          f"SIN garantía de que sea correcta. No reportar estas cifras por clase a CONANP/SEDEMA sin "
          f"verificar manualmente contra FAMILIA primero.")
    col_clase_usar = "CLASE"
    TEXTO_VALIDACION_CLASE = (
        "[SIN VALIDAR -- el CSV de biodiversidad no tiene columna FAMILIA; cifras por clase (CLASE cruda) "
        "SIN garantía, ver [ALERTA] impresa arriba]"
    )

gdf_nucleo_3857 = gdf_nucleo.to_crs(epsg=3857)
gdf_500_3857 = gdf_500.to_crs(epsg=3857) if len(gdf_500) > 0 else None
gdf_1000_3857 = gdf_1000.to_crs(epsg=3857) if len(gdf_1000) > 0 else None

def regs_spp(df, mask=None):
    d = df[mask] if mask is not None else df
    return len(d), int(d["NOMBRE_CIENTIFICO"].nunique())

total_regs, total_spp = regs_spp(df_bio)

if "NIVEL" in df_bio.columns:
    mask_nucleo = df_bio["NIVEL"].astype(str).str.contains("NUCLEO|DESGLOSE_29", na=False)
    mask_500 = df_bio["NIVEL"].astype(str).str.contains("500", na=False)
    mask_1000 = df_bio["NIVEL"].astype(str).str.contains("1000", na=False)
else:
    print("[WARN] el CSV de biodiversidad no tiene columna NIVEL -- no puedo desglosar regs/spp por zona, se reportará todo como 'sin zona'.")
    mask_nucleo = pd.Series([False] * len(df_bio))
    mask_500 = pd.Series([False] * len(df_bio))
    mask_1000 = pd.Series([False] * len(df_bio))

regs_nucleo, spp_nucleo = regs_spp(df_bio, mask_nucleo)
regs_500, spp_500 = regs_spp(df_bio, mask_500)
regs_1000, spp_1000 = regs_spp(df_bio, mask_1000)

dens_nucleo = regs_nucleo / SUP_NUCLEO if SUP_NUCLEO else 0.0
dens_500 = regs_500 / SUP_500 if SUP_500 else 0.0
dens_1000 = regs_1000 / SUP_1000 if SUP_1000 else 0.0

if col_clase_usar in df_bio.columns:
    clase_counts = df_bio[col_clase_usar].value_counts()
    clase_pct = (clase_counts / total_regs * 100).round(1) if total_regs else clase_counts
else:
    clase_counts = pd.Series(dtype=int)
    clase_pct = pd.Series(dtype=float)

def texto_clase(nombre):
    if nombre in clase_counts.index:
        return f"{clase_counts[nombre]:,} ({clase_pct[nombre]:.0f}%)"
    return "sin datos"

NDWI_TXT = f"{args.ndwi:.3f}" if args.ndwi is not None else "[NDWI NO VERIFICADO -- pasar --ndwi]"
RIESGO_TXT = f"{args.spp_alto_riesgo}" if args.spp_alto_riesgo is not None else "[Nº SPP ALTO RIESGO NO VERIFICADO -- pasar --spp-alto-riesgo]"

# ------------------------------------------------------------------
# MAPAS (misma lógica visual que el original; solo cambian los datos que citan)
# ------------------------------------------------------------------
def plot_hibrido(gdf_puntos_df, color_col, color_dict, title, out_path, legend_title, show_poligonos=True):
    if gdf_1000_3857 is not None and len(gdf_1000_3857) > 0:
        xmin, ymin, xmax, ymax = gdf_1000_3857.total_bounds
    else:
        xmin, ymin, xmax, ymax = gdf_nucleo_3857.total_bounds
        xmin -= 2000; ymin -= 2000; xmax += 2000; ymax += 2000

    fig, ax = plt.subplots(1, 1, figsize=(14, 14), dpi=300)
    ax.set_xlim(xmin - 500, xmax + 500)
    ax.set_ylim(ymin - 500, ymax + 500)

    if CTX:
        try:
            ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=15, attribution=False, alpha=0.95)
        except Exception as e:
            print(f"[WARN] Esri fail {e}")

    if show_poligonos:
        if gdf_1000_3857 is not None and len(gdf_1000_3857) > 0:
            gdf_1000_3857.boundary.plot(ax=ax, color="red", linewidth=1.8, linestyle=":", alpha=0.8,
                label=f"Buffer 1000m {SUP_1000:.1f} ha | AGB {AGB_1000:.1f} | CO2 {CO2_1000:,.0f} t | "
                      f"{len(df_bio[df_bio.NIVEL=='BUFFER_1000m']) if 'NIVEL' in df_bio.columns and 'BUFFER_1000m' in df_bio.NIVEL.values else 0} regs")
        if gdf_500_3857 is not None and len(gdf_500_3857) > 0:
            gdf_500_3857.boundary.plot(ax=ax, color="orange", linewidth=2.2, linestyle="--", alpha=0.9,
                label=f"Buffer 500m {SUP_500:.1f} ha | AGB {AGB_500:.1f} | CO2 {CO2_500:,.0f} t | "
                      f"{len(df_bio[df_bio.NIVEL=='BUFFER_500m']) if 'NIVEL' in df_bio.columns and 'BUFFER_500m' in df_bio.NIVEL.values else 0} regs")
        gdf_nucleo_3857.plot(ax=ax, facecolor="yellow", edgecolor="black", alpha=0.60, linewidth=2.5)

    gdf_puntos_3857 = gpd.GeoDataFrame(gdf_puntos_df, geometry=gpd.points_from_xy(gdf_puntos_df.LON, gdf_puntos_df.LAT), crs="EPSG:4326").to_crs(epsg=3857)

    for key, color in color_dict.items():
        pts = gdf_puntos_3857[gdf_puntos_3857[color_col] == key]
        if len(pts) == 0:
            continue
        pts.plot(ax=ax, color=color, markersize=18 if color_col != "ANIO" else 14, alpha=0.85,
                  edgecolor="black", linewidth=0.35, label=f"{key} ({len(pts)} regs)")

    ax.text(0.04, 0.92, "▲\nN", transform=ax.transAxes, fontsize=22, fontweight="bold", ha="center", va="top",
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"))

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    sx = x_min + (x_max - x_min) * 0.04
    sy = y_min + (y_max - y_min) * 0.06
    ax.plot([sx, sx + 1000], [sy, sy], color="black", linewidth=7)
    ax.plot([sx, sx + 500], [sy, sy], color="white", linewidth=7)
    ax.text(sx, sy + (y_max - y_min) * 0.018, "1 km | 0.5 km", fontsize=10, fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.95))

    ax.set_title(
        f"{title}\n{SUP_NUCLEO:.2f} ha = {CO2_NUCLEO:,.0f} t | "
        f"{SUP_NUCLEO+SUP_500+SUP_1000:.0f} ha = {CO2_TOTAL:,.0f} t (anillo exclusivo, sí sumable) | "
        f"AGB núcleo {AGB_NUCLEO:.1f} / buf500 {AGB_500:.1f} / buf1000 {AGB_1000:.1f} | NDWI {NDWI_TXT} | GEDI L4B MU + S2",
        fontsize=10, fontweight="bold", pad=12)
    # BUG encontrado y corregido aquí (heredado del script original): geopandas dibuja el
    # relleno del núcleo como PatchCollection, y matplotlib NO registra automáticamente un
    # PatchCollection en la leyenda aunque se le pase label= -- por eso "NÚCLEO ... ha" nunca
    # aparecía en la leyenda de ningún mapa (solo salían los anillos 500m/1000m, que sí son
    # líneas). Se arma un proxy (mpatches.Patch) a mano para que el núcleo sí aparezca.
    handles, labels = ax.get_legend_handles_labels()
    if show_poligonos:
        nucleo_patch = mpatches.Patch(
            facecolor="yellow", edgecolor="black", alpha=0.60,
            label=f"NÚCLEO {SUP_NUCLEO:.2f} ha | AGB {AGB_NUCLEO:.1f} | CO2 {CO2_NUCLEO:,.0f} t | {regs_nucleo} regs")
        handles = [nucleo_patch] + handles
        labels = [nucleo_patch.get_label()] + labels
    ax.legend(handles=handles, labels=labels, title=legend_title, loc="upper right", fontsize=8, framealpha=0.97, facecolor="white")
    ax.set_axis_off()

    cajetin = f"""EPSG:3857 / 4326 - WGS84 UTM 14N
SRC: Esri Sat + GEDI L4B MU + S2 + GBIF iNat
V19.4 HÍBRIDO CORREGIDO - {datetime.now().strftime('%d/%m/%Y')}
Ramsar {args.sitio_ramsar}{RAMSAR_ID_SUFFIX} - CONANP/SEDEMA VER
{SUP_NUCLEO:.2f} ha = {CO2_NUCLEO:,.0f} t CO2e | {total_spp} spp | {total_regs:,} regs
CO2 total paisaje (anillo exclusivo, sí sumable) {CO2_TOTAL:,.0f} t
"""
    fig.text(0.99, 0.02, cajetin, ha="right", va="bottom", fontsize=7,
              bbox=dict(facecolor="white", alpha=0.95, edgecolor="black"), transform=fig.transFigure)
    fig.text(0.5, 0.97, "V19.4 HÍBRIDO CORREGIDO - POLÍGONO CO2 (real) + BIODIVERSIDAD (real) - EUDR / CONANP / SEDEMA",
              ha="center", fontsize=10, fontweight="bold", color="white", bbox=dict(facecolor="black"))

    plt.tight_layout()
    plt.savefig(out_path, dpi=350, bbox_inches="tight")
    plt.close()
    print(f"[MAPA MIX PNG] {out_path}")

col_zona = {"NUCLEO_DESGLOSE_29": "#FFD700", "BUFFER_500m": "#FF8C00", "BUFFER_1000m": "#FF0000"}
if "NIVEL" in df_bio.columns:
    for k in df_bio["NIVEL"].unique():
        if k not in col_zona:
            col_zona[k] = "#FF0000"

col_clase = {"Mammalia": "#8B4513", "Amphibia": "#00CED1", "Aves": "#1E90FF", "Reptilia": "#32CD32"}
if col_clase_usar in df_bio.columns:
    for cl in df_bio[col_clase_usar].unique():
        if cl not in col_clase:
            col_clase[cl] = "#FF00FF"

if "NIVEL" in df_bio.columns:
    plot_hibrido(df_bio, "NIVEL", col_zona,
                 "V19.4 CORREGIDO - CO2 + DISTRIBUCIÓN BIODIVERSIDAD POR ZONA: núcleo / 500m / 1000m",
                 os.path.join(out_dir, "V19_4_MIX_CO2_ZONAS_HD.png"),
                 f"Zona (NIVEL) - {total_regs:,} regs", show_poligonos=True)

if col_clase_usar in df_bio.columns:
    _sufijo_validada = " (CLASE_VALIDADA)" if col_clase_usar == "CLASE_VALIDADA" else " (CLASE cruda, sin validar)"
    plot_hibrido(df_bio, col_clase_usar, col_clase,
                 f"V19.4 CORREGIDO - CO2 + DISTRIBUCIÓN POR CLASE TAXONÓMICA{_sufijo_validada}",
                 os.path.join(out_dir, "V19_4_MIX_CO2_CLASES_HD.png"),
                 f"Clase Taxonómica{_sufijo_validada} - {total_spp} spp", show_poligonos=True)

    df_nucleo_only = df_bio[mask_nucleo] if mask_nucleo.any() else df_bio
    plot_hibrido(df_nucleo_only, col_clase_usar, col_clase,
                 f"V19.4 CORREGIDO - SOLO NÚCLEO + CLASES{_sufijo_validada} - {len(df_nucleo_only)} regs",
                 os.path.join(out_dir, "V19_4_MIX_SOLO_NUCLEO_CLASES_HD.png"),
                 f"Clase Núcleo{_sufijo_validada} - {len(df_nucleo_only)} regs", show_poligonos=True)

col_anio_name = "ANIO" if "ANIO" in df_bio.columns else "YEAR" if "YEAR" in df_bio.columns else None
if col_anio_name:
    anios = sorted(df_bio[col_anio_name].dropna().unique())
    cmap = plt.cm.tab10
    col_anio = {a: cmap(i % 10) for i, a in enumerate(anios)}
    plot_hibrido(df_bio, col_anio_name, col_anio,
                 "V19.4 CORREGIDO - CO2 + DISTRIBUCIÓN POR AÑO",
                 os.path.join(out_dir, "V19_4_MIX_ANIO_HD.png"),
                 "Año", show_poligonos=True)

# ------------------------------------------------------------------
# TABLAS
# ------------------------------------------------------------------
if "NIVEL" in df_bio.columns:
    resumen = df_bio.groupby("NIVEL").agg(REGS=("NOMBRE_CIENTIFICO", "count"), SPP=("NOMBRE_CIENTIFICO", "nunique")).reset_index()
    resumen.to_csv(os.path.join(out_dir, "V19_4_RESUMEN_BIO.csv"), index=False, encoding="utf-8-sig")
    print("[TABLA] Resumen bio guardado")

# ------------------------------------------------------------------
# DICTAMEN -- ahora con datos reales, no texto fijo
# ------------------------------------------------------------------
dictamen_largo = f"""
DICTAMEN TÉCNICO-CIENTÍFICO - V19.4 HÍBRIDO CORREGIDO (V17+V19, datos reales)
RAMSAR {args.sitio_ramsar}{RAMSAR_ID_SUFFIX} - DESGLOSE 29 - {datetime.now().strftime('%d/%m/%Y')}

I. ANTECEDENTES Y OBJETO:
Desglose 29 de {SUP_NUCLEO:.2f} ha dentro de Sitio Ramsar {args.sitio_ramsar}{RAMSAR_ID_SUFFIX}.
{RAMSAR_CRITERIOS_TXT}
Dictamen integrado CO2 + Biodiversidad para CONANP y SEDEMA Veracruz, modelo híbrido V17
(biomasa) + V19 (biodiversidad, {total_regs:,} regs, {total_spp} spp).

II. METODOLOGÍA:
- Base cartográfica: {os.path.basename(args.geojson)} -- Núcleo {SUP_NUCLEO:.2f} ha, Buffer 500m {SUP_500:.1f} ha,
  Buffer 1000m {SUP_1000:.1f} ha. Total paisaje (superficie en planta / catastral) {SUP_NUCLEO+SUP_500+SUP_1000:.2f} ha.
- Superficie real de terreno (corregida por relieve, 3D -- la que efectivamente ocupa el bosque en el
  suelo, siguiendo la pendiente, no su proyección plana): núcleo {texto_relieve(FACTOR_RELIEVE_NUCLEO, AREA3D_NUCLEO, SUP_NUCLEO)},
  buf500 {texto_relieve(FACTOR_RELIEVE_500, AREA3D_500, SUP_500)}, buf1000 {texto_relieve(FACTOR_RELIEVE_1000, AREA3D_1000, SUP_1000)}.
  {TEXTO_EXTENSION_REAL}
- Biomasa/CO2: valores de anillo exclusivo tomados de {os.path.basename(args.csv_carbono)} (columna '{col_co2_incr}'),
  ya corregidos para no duplicar carbono entre zonas anidadas. AGB núcleo {AGB_NUCLEO:.1f} Mg/ha,
  buf500 {AGB_500:.1f} Mg/ha, buf1000 {AGB_1000:.1f} Mg/ha.
- CO2: Núcleo {CO2_NUCLEO:,.0f} t, Buf500 {CO2_500:,.0f} t, Buf1000 {CO2_1000:,.0f} t.
  Total paisaje (anillo exclusivo, sí sumable) {CO2_TOTAL:,.0f} t.
- Biodiversidad: {os.path.basename(fuente_biodiversidad)}{' (descargado de GBIF por este mismo script, todas las clases)' if not args.csv_biodiversidad else ''}. {total_regs:,} registros, {total_spp} especies únicas.
  Núcleo {regs_nucleo:,} regs ({dens_nucleo:.1f} regs/ha, {spp_nucleo} spp), Buf500 {regs_500:,} regs
  ({dens_500:.1f} regs/ha, {spp_500} spp), Buf1000 {regs_1000:,} regs ({dens_1000:.1f} regs/ha, {spp_1000} spp).
- Mapas: Esri World Imagery zoom 15 + contextily + matplotlib 350 DPI.

III. RESULTADOS CO2:
{SUP_NUCLEO:.2f} ha = {CO2_NUCLEO:,.0f} t CO2e. Paisaje total {SUP_NUCLEO+SUP_500+SUP_1000:.2f} ha = {CO2_TOTAL:,.0f} t
(anillo exclusivo -- estas cifras SÍ se pueden sumar entre zonas porque cada una representa
solo su propio anillo, no el área acumulada desde el centro).

IV. RESULTADOS BIODIVERSIDAD:
- Clasificación taxonómica: {TEXTO_VALIDACION_CLASE}
- Mammalia: {texto_clase('Mammalia')} | Amphibia: {texto_clase('Amphibia')} | Aves: {texto_clase('Aves')} | Reptilia: {texto_clase('Reptilia')}
  (clase taxonómica {'VALIDADA contra FAMILIA' if col_clase_usar == 'CLASE_VALIDADA' else 'SIN VALIDAR -- CLASE cruda del CSV'})
- Densidad núcleo: {dens_nucleo:.1f} regs/ha.

V. INTERPRETACIÓN CO2+BIO (cualitativa -- no derivada de los CSV, requiere revisión del autor):
[Esta sección del dictamen original contenía interpretación ecológica narrativa
(gradiente de riqueza vs. biomasa, efecto borde, sesgo de observador urbano, etc.).
Ese análisis no se recalcula automáticamente aquí porque no proviene de una columna
verificable en los CSV de entrada -- Ruben debe revisar si sigue aplicando con los
números reales de esta corrida antes de reinsertarla.]

VI. NOTAS PENDIENTES DE VERIFICACIÓN ANTES DE ENVIAR A CONANP/SEDEMA:
- NDWI citado: {NDWI_TXT}
- Nº especies en alto riesgo (núcleo): {RIESGO_TXT}
- ID de sitio Ramsar: {args.sitio_ramsar}{RAMSAR_ID_SUFFIX}. {RAMSAR_CRITERIOS_TXT}
- Validación taxonómica CLASE vs FAMILIA: {TEXTO_VALIDACION_CLASE}
- Verificar que --geojson y --csv-carbono correspondan a la MISMA corrida (ver [ALERTA]/[OK] de área impresos al inicio de esta ejecución).

VII. PRODUCTOS V19.4:
- V19_4_MIX_CO2_ZONAS_HD.png / V19_4_MIX_CO2_CLASES_HD.png / V19_4_MIX_SOLO_NUCLEO_CLASES_HD.png / V19_4_MIX_ANIO_HD.png
- CSV de resumen y este dictamen.

Fecha: {datetime.now().strftime('%d/%m/%Y')} - {SUP_NUCLEO:.2f} ha = {CO2_NUCLEO:,.0f} t CO2e (fuente: {os.path.basename(args.csv_carbono)})
"""

with open(os.path.join(out_dir, "V19_4_DICTAMEN_HIBRIDO_CORREGIDO.txt"), "w", encoding="utf-8") as f:
    f.write(dictamen_largo)

if REPORTLAB:
    pdf_path = os.path.join(out_dir, "V19_4_DICTAMEN_HIBRIDO_CORREGIDO.pdf")
    doc = SimpleDocTemplate(pdf_path, pagesize=letter, rightMargin=45, leftMargin=45, topMargin=45, bottomMargin=45)
    styles = getSampleStyleSheet()
    story = []
    story.append(Paragraph(
        f"DICTAMEN TÉCNICO-CIENTÍFICO - V19.4 HÍBRIDO CORREGIDO<br/>RAMSAR {args.sitio_ramsar}{RAMSAR_ID_SUFFIX} - DESGLOSE 29 - "
        f"{SUP_NUCLEO:.2f} ha = {CO2_NUCLEO:,.0f} t CO2e - {total_spp} spp", styles["Heading1"]))
    story.append(Spacer(1, 0.2 * inch))
    data = [["ZONA", "HA", "AGB", "CO2e (anillo excl.)", "REGS", "SPP", "DENS regs/ha"]]
    data.append(["NÚCLEO", f"{SUP_NUCLEO:.2f}", f"{AGB_NUCLEO:.1f}", f"{CO2_NUCLEO:,.0f}", f"{regs_nucleo}", f"{spp_nucleo}", f"{dens_nucleo:.1f}"])
    data.append(["BUF 500m", f"{SUP_500:.1f}", f"{AGB_500:.1f}", f"{CO2_500:,.0f}", f"{regs_500}", f"{spp_500}", f"{dens_500:.1f}"])
    data.append(["BUF 1000m", f"{SUP_1000:.1f}", f"{AGB_1000:.1f}", f"{CO2_1000:,.0f}", f"{regs_1000}", f"{spp_1000}", f"{dens_1000:.1f}"])
    data.append(["TOTAL", f"{SUP_NUCLEO+SUP_500+SUP_1000:.1f}", "", f"{CO2_TOTAL:,.0f}", f"{regs_nucleo+regs_500+regs_1000}", f"{total_spp}", ""])
    t = Table(data, colWidths=[65, 45, 35, 75, 40, 35, 55])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#003366")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
    ]))
    story.append(t)
    story.append(Spacer(1, 0.25 * inch))
    for para in dictamen_largo.split("\n\n"):
        if para.strip():
            story.append(Paragraph(para.replace("\n", "<br/>"), styles["Normal"]))
            story.append(Spacer(1, 0.15 * inch))
    doc.build(story)
    print(f"[DICTAMEN PDF] {pdf_path}")

print("\n" + "=" * 90)
print(f"V19.4 HÍBRIDO CORREGIDO COMPLETADO - Carpeta: {out_dir}")
for fn in sorted(os.listdir(out_dir)):
    print(f" - {fn}")
print("=" * 90)

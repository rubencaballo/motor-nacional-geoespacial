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

Uso:
  conda activate geo
  python V19_4_HIBRIDO_CORREGIDO.py \
      --geojson /ruta/exacta/Buffer_500_1000m.geojson \
      --csv-carbono /ruta/exacta/resumen_terreno_y_carbono_2025.csv \
      --csv-biodiversidad /ruta/exacta/Biodiversidad_DESGLOSE_29_COMPLETO.csv \
      --sitio-ramsar 1601

Si no sabes las rutas exactas, en WSL corre:
  find ~ -iname "Buffer_500_1000m.geojson" 2>/dev/null
  find ~ -iname "resumen_terreno_y_carbono*.csv" 2>/dev/null
  find ~ -iname "*DESGLOSE_29_COMPLETO*.csv" 2>/dev/null
y copia/pega la ruta que corresponda a la corrida vigente (la más reciente que tú reconozcas,
no la que el script "adivine").
"""
import os, sys, argparse
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime

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
parser.add_argument("--csv-biodiversidad", required=True, help="Ruta EXACTA al CSV de biodiversidad (Biodiversidad_DESGLOSE_29_COMPLETO o equivalente)")
parser.add_argument("--sitio-ramsar", default="1601", help="ID del sitio Ramsar a imprimir en mapas/dictamen. El script original traía '1791' escrito a mano; el resto del proyecto usa '1601'. CONFIRMAR cuál es el correcto -- default aquí: 1601")
parser.add_argument("--ndwi", type=float, default=None, help="Valor NDWI a citar en el dictamen (el original traía -0.681 fijo, sin fuente verificable en este script). Si no se da, se marca como NO VERIFICADO.")
parser.add_argument("--spp-alto-riesgo", type=int, default=None, help="Nº de especies en alto riesgo dentro del núcleo, para el dictamen (el original traía '24' fijo, sin fuente en este script). Si no se da, se marca como NO VERIFICADO.")
parser.add_argument("--out-dir", default=None, help="Carpeta de salida (default: junto al geojson)")
parser.add_argument("--tolerancia-area-pct", type=float, default=2.0, help="Tolerancia %% para avisar si el área del geojson no coincide con la del CSV de carbono (default 2%%)")
args = parser.parse_args()

for etiqueta, ruta in [("geojson", args.geojson), ("csv-carbono", args.csv_carbono), ("csv-biodiversidad", args.csv_biodiversidad)]:
    if not os.path.isfile(ruta):
        sys.exit(f"ERROR: no existe el archivo pasado en --{etiqueta}: {ruta}")

print(f"[GEOJSON] {args.geojson}")
print(f"[CSV CARBONO] {args.csv_carbono}")
print(f"[CSV BIODIVERSIDAD] {args.csv_biodiversidad}")
if args.sitio_ramsar == "1601":
    print("[AVISO] Usando Ramsar 1601 por default (el script original decía 1791 en el texto del dictamen). Confirma cuál es el ID correcto del sitio antes de mandar el dictamen a CONANP/SEDEMA.")

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
# BIODIVERSIDAD: se sigue leyendo igual, pero ahora TODO el texto del dictamen
# se calcula de aquí en vez de estar escrito a mano.
# ------------------------------------------------------------------
df_bio = pd.read_csv(args.csv_biodiversidad)
if "NOMBRE_CIENTIFICO" not in df_bio.columns:
    sys.exit(f"ERROR: {args.csv_biodiversidad} no tiene columna NOMBRE_CIENTIFICO. Columnas: {list(df_bio.columns)}")
print(f"[BIO] {len(df_bio)} regs, {df_bio.NOMBRE_CIENTIFICO.nunique()} spp, cols={list(df_bio.columns)[:6]}")

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

if "CLASE" in df_bio.columns:
    clase_counts = df_bio["CLASE"].value_counts()
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
Ramsar {args.sitio_ramsar} (verificar) - CONANP/SEDEMA VER
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
if "CLASE" in df_bio.columns:
    for cl in df_bio["CLASE"].unique():
        if cl not in col_clase:
            col_clase[cl] = "#FF00FF"

if "NIVEL" in df_bio.columns:
    plot_hibrido(df_bio, "NIVEL", col_zona,
                 "V19.4 CORREGIDO - CO2 + DISTRIBUCIÓN BIODIVERSIDAD POR ZONA: núcleo / 500m / 1000m",
                 os.path.join(out_dir, "V19_4_MIX_CO2_ZONAS_HD.png"),
                 f"Zona (NIVEL) - {total_regs:,} regs", show_poligonos=True)

if "CLASE" in df_bio.columns:
    plot_hibrido(df_bio, "CLASE", col_clase,
                 "V19.4 CORREGIDO - CO2 + DISTRIBUCIÓN POR CLASE TAXONÓMICA",
                 os.path.join(out_dir, "V19_4_MIX_CO2_CLASES_HD.png"),
                 f"Clase Taxonómica - {total_spp} spp", show_poligonos=True)

    df_nucleo_only = df_bio[mask_nucleo] if mask_nucleo.any() else df_bio
    plot_hibrido(df_nucleo_only, "CLASE", col_clase,
                 f"V19.4 CORREGIDO - SOLO NÚCLEO + CLASES - {len(df_nucleo_only)} regs",
                 os.path.join(out_dir, "V19_4_MIX_SOLO_NUCLEO_CLASES_HD.png"),
                 f"Clase Núcleo - {len(df_nucleo_only)} regs", show_poligonos=True)

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
RAMSAR {args.sitio_ramsar} (VERIFICAR ID) - DESGLOSE 29 - {datetime.now().strftime('%d/%m/%Y')}

I. ANTECEDENTES Y OBJETO:
Desglose 29 de {SUP_NUCLEO:.2f} ha dentro de Sitio Ramsar {args.sitio_ramsar} (confirmar ID exacto).
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
- Biodiversidad: {os.path.basename(args.csv_biodiversidad)}. {total_regs:,} registros, {total_spp} especies únicas.
  Núcleo {regs_nucleo:,} regs ({dens_nucleo:.1f} regs/ha, {spp_nucleo} spp), Buf500 {regs_500:,} regs
  ({dens_500:.1f} regs/ha, {spp_500} spp), Buf1000 {regs_1000:,} regs ({dens_1000:.1f} regs/ha, {spp_1000} spp).
- Mapas: Esri World Imagery zoom 15 + contextily + matplotlib 350 DPI.

III. RESULTADOS CO2:
{SUP_NUCLEO:.2f} ha = {CO2_NUCLEO:,.0f} t CO2e. Paisaje total {SUP_NUCLEO+SUP_500+SUP_1000:.2f} ha = {CO2_TOTAL:,.0f} t
(anillo exclusivo -- estas cifras SÍ se pueden sumar entre zonas porque cada una representa
solo su propio anillo, no el área acumulada desde el centro).

IV. RESULTADOS BIODIVERSIDAD:
- Mammalia: {texto_clase('Mammalia')} | Amphibia: {texto_clase('Amphibia')} | Aves: {texto_clase('Aves')} | Reptilia: {texto_clase('Reptilia')}
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
- ID de sitio Ramsar: {args.sitio_ramsar} (el script original decía 1791; el resto del proyecto usa 1601 -- confirmar cuál es correcto)
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
        f"DICTAMEN TÉCNICO-CIENTÍFICO - V19.4 HÍBRIDO CORREGIDO<br/>RAMSAR {args.sitio_ramsar} (verificar) - DESGLOSE 29 - "
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

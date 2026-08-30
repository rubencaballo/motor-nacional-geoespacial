#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cruza el carbono almacenado por zona (core/carbono.py: biomasa AGB
satelital ESA CCI + muestreo LiDAR real GEDI) con las hectáreas de pérdida
CONFIRMADAS por Hansen (core/deforestacion.py, en anillo exclusivo vía
generar_resumen_no_traslapado) para estimar el CO2e asociado a la
superficie que YA SE PERDIÓ -- no el inventario completo de la zona, solo
la parte que el propio historial de pérdida dice que ya no está.

POR QUÉ UN MÓDULO APARTE (aunque solo cruza dos CSV que ya existen):
    core/carbono.py estima el carbono ALMACENADO en una zona completa hoy
    (inventario). core/deforestacion.py estima las HECTÁREAS perdidas
    desde un año base. Ninguno de los dos, por separado, contesta "¿cuánto
    CO2e representa lo que ya se perdió?" -- eso requiere multiplicar
    densidad de carbono (t CO2e/ha) por hectáreas perdidas, zona por zona,
    y esta es la única función del proyecto que hace esa multiplicación.

MUÑECA RUSA, atendida desde el diseño: tanto el carbono como la pérdida
    existen en dos formas -- acumulativa por disco (zona/buffer_m) y
    exclusiva por anillo (sin traslape). Cruzar el CO2e acumulativo de un
    disco contra la pérdida acumulativa de ese mismo disco daría un número
    "parece que sí pero no" (ambos inflados por lo mismo, se cancelaría
    parcialmente el error de forma impredecible). Por eso esta función
    SIEMPRE trabaja en anillo exclusivo de los dos lados: la densidad de
    carbono del anillo (co2e_incremental_t / área del anillo) multiplicada
    por la pérdida del anillo (ya exclusiva, de
    generar_resumen_no_traslapado) -- así el resultado por anillo SÍ se
    puede sumar directo para un total real.

QUÉ ASUME (documentado, no escondido):
    - Densidad uniforme dentro del anillo: la pérdida confirmada se valora
      con la densidad de carbono PROMEDIO de todo el anillo, porque no
      sabemos si los píxeles talados específicamente eran más o menos
      densos que el resto del anillo (ej. si estaban junto a un arroyo,
      probablemente más densos -- esto podría estar SUBESTIMANDO el CO2e
      real perdido, no sobreestimándolo).
    - Remoción CASI COMPLETA de biomasa aérea en el píxel perdido -- válido
      para tala/desmonte confirmado (lo que este módulo asume por
      defecto), pero NO para pérdida por INCENDIO: un incendio no combuste
      el 100% de la biomasa (ver core/validacion_incendios.py y el factor
      de combustión IPCC ~0.40-0.50 usado ahí). Si la causa confirmada es
      incendio, ese número -- no este módulo de remoción completa -- es el
      correcto.
    - No incluye biomasa subterránea, carbono del suelo, ni necromasa --
      mismo alcance (solo AGB) que core/carbono.py, heredado tal cual.
    - AGB (ESA CCI) y GEDI son dos estimaciones INDEPENDIENTES, nunca se
      promedian entre sí -- se reportan ambas, lado a lado, como ya hace
      core/carbono.py (uno es modelo satelital de cobertura completa, el
      otro es muestreo LiDAR real pero disperso)."""

import argparse
import os

import numpy as np
import pandas as pd

from config import log


# ==============================================================================
# --- ANILLO EXCLUSIVO DE CARBONO (reusa incremental_t si ya viene calculado) --
# ==============================================================================
def _anillo_exclusivo_de_carbono(df_terreno_carbono):
    """Acepta tanto resumen_terreno_y_carbono_*.csv (ya trae
    co2e_incremental_t/gedi_co2e_incremental_t calculados por
    core/carbono.py) como carbono_por_zona_*.csv (formato más viejo, sin
    esas columnas -- en ese caso se calculan aquí por resta consecutiva,
    con la incertidumbre combinada en cuadratura como aproximación
    explícita, NO una propagación rigurosa -- documentado, no una verdad
    exacta).

    Devuelve un dict {buffer_m: {area_ha, co2e_t, co2e_incert_t,
    gedi_co2e_t, gedi_co2e_incert_t}} -- todo YA en anillo exclusivo -- más
    zona_de_buffer y nombre_anillo (mismo esquema de nombres que el resto
    del proyecto)."""
    faltantes = [c for c in ["zona", "buffer_m"] if c not in df_terreno_carbono.columns]
    if faltantes:
        raise ValueError(f"Faltan columnas {faltantes} en el CSV de carbono.")

    col_area = "area_2d_ha" if "area_2d_ha" in df_terreno_carbono.columns else "area_ha"
    tiene_incremental = "co2e_incremental_t" in df_terreno_carbono.columns

    df = df_terreno_carbono.sort_values("buffer_m").reset_index(drop=True)
    buffers_ordenados = df["buffer_m"].tolist()
    zona_de_buffer = {row["buffer_m"]: row["zona"] for _, row in df.iterrows()}
    nombre_anillo = {}
    for i, buf_m in enumerate(buffers_ordenados):
        if i == 0:
            nombre_anillo[buf_m] = "nucleo (0m)" if buf_m == 0 else f"anillo_0-{buf_m}m"
        else:
            nombre_anillo[buf_m] = f"anillo_{buffers_ordenados[i - 1]}-{buf_m}m"

    anillos = {}
    area_prev = co2e_prev = co2e_i_prev = gedi_prev = gedi_i_prev = 0.0
    for _, row in df.iterrows():
        buf_m = row["buffer_m"]
        area = float(row[col_area])
        co2e = float(row["co2e_t"])
        co2e_i = float(row.get("co2e_incertidumbre_t", 0.0) or 0.0)
        gedi = float(row["gedi_co2e_t"]) if "gedi_co2e_t" in row and pd.notna(row["gedi_co2e_t"]) else None
        gedi_i = float(row.get("gedi_co2e_incertidumbre_t", 0.0) or 0.0)

        area_anillo = round(area - area_prev, 4)
        if tiene_incremental:
            co2e_anillo = float(row["co2e_incremental_t"])
            co2e_i_anillo = float(row.get("co2e_incremental_incertidumbre_t", 0.0) or 0.0)
            gedi_anillo = float(row["gedi_co2e_incremental_t"]) if "gedi_co2e_incremental_t" in row else None
            gedi_i_anillo = float(row.get("gedi_co2e_incremental_incertidumbre_t", 0.0) or 0.0)
        else:
            co2e_anillo = round(co2e - co2e_prev, 4)
            co2e_i_anillo = round(float(np.sqrt(co2e_i ** 2 + co2e_i_prev ** 2)), 4)  # aproximacion, ver docstring
            gedi_anillo = round(gedi - gedi_prev, 4) if gedi is not None else None
            gedi_i_anillo = round(float(np.sqrt(gedi_i ** 2 + gedi_i_prev ** 2)), 4)

        anillos[buf_m] = {
            "area_ha": area_anillo, "co2e_t": co2e_anillo, "co2e_incertidumbre_t": co2e_i_anillo,
            "gedi_co2e_t": gedi_anillo, "gedi_co2e_incertidumbre_t": gedi_i_anillo,
        }
        area_prev, co2e_prev, co2e_i_prev = area, co2e, co2e_i
        gedi_prev, gedi_i_prev = (gedi or 0.0), gedi_i

    return anillos, zona_de_buffer, nombre_anillo


# ==============================================================================
# --- CRUCE: densidad de carbono del anillo x hectáreas perdidas del anillo ---
# ==============================================================================
def cruzar_carbono_con_perdida(df_terreno_carbono, df_perdida_sin_traslape, id_proyecto, carpeta_salida=None,
                                causa_default="tala/desmonte u otra causa no confirmada (remocion completa asumida, cota superior)",
                                factor_combustion_por_anio=None, causa_por_anio=None):
    """df_terreno_carbono: resumen_terreno_y_carbono_*.csv o carbono_por_zona_*.csv
    (columnas zona/buffer_m/area_.../co2e_t/gedi_co2e_t...).
    df_perdida_sin_traslape: deforestacion_resumen_sin_traslape_*.csv de
    core.deforestacion.generar_resumen_no_traslapado (columnas
    anillo/anio/perdida_ha, incluye filas 'TOTAL (suma sin traslape)' y
    'TOTAL {anio_min}-{anio_max}').

    factor_combustion_por_anio: dict opcional {anio: factor 0-1} -- para
    AÑOS CONFIRMADOS COMO INCENDIO (ver core/validacion_incendios.py), un
    incendio NO combuste el 100% de la biomasa aérea, a diferencia de tala/
    desmonte (que este módulo asume por defecto como remoción casi
    completa, factor=1.0). Pásale aquí el factor de combustión IPCC
    (~0.40-0.50 para bosque templado, mismo valor ya usado en el análisis
    de Cofre de Perote) para el/los años donde la causa confirmada fue
    incendio -- años no listados usan factor=1.0 (supuesto por defecto,
    cota superior). NUNCA se debe dejar factor=1.0 en un año de incendio
    confirmado -- sobreestimaría el CO2e de ese año.
    causa_por_anio: dict opcional {anio: texto} para anotar la causa
    específica de ese año en la columna 'causa_probable' (ej. "incendio
    confirmado por dNBR") en vez de causa_default.

    Devuelve un DataFrame con, por año y anillo: perdida_ha, densidad de
    carbono del anillo (t CO2e/ha, AGB y GEDI por separado), el factor de
    combustión aplicado, y el CO2e asociado a esa pérdida específica (AGB y
    GEDI) -- todo en anillo exclusivo, así que las filas TOTAL sí se pueden
    sumar directo. NUNCA promedia AGB con GEDI -- se reportan las dos
    estimaciones lado a lado, ver docstring del módulo."""
    factor_combustion_por_anio = factor_combustion_por_anio or {}
    causa_por_anio = causa_por_anio or {}
    anillos, zona_de_buffer, nombre_anillo = _anillo_exclusivo_de_carbono(df_terreno_carbono)

    densidad = {}
    for buf_m, datos in anillos.items():
        area = datos["area_ha"]
        if area <= 0:
            log(f"  anillo '{nombre_anillo[buf_m]}' con área <=0 ha ({area}) -- densidad no calculable, se omite.",
                nivel="WARN")
            continue
        densidad[nombre_anillo[buf_m]] = {
            "densidad_agb_tco2e_ha": round(datos["co2e_t"] / area, 4),
            "densidad_gedi_tco2e_ha": round(datos["gedi_co2e_t"] / area, 4) if datos["gedi_co2e_t"] is not None else None,
        }
        log(f"  {nombre_anillo[buf_m]}: densidad AGB={densidad[nombre_anillo[buf_m]]['densidad_agb_tco2e_ha']:.1f} "
            f"t CO2e/ha | GEDI={densidad[nombre_anillo[buf_m]]['densidad_gedi_tco2e_ha']:.1f} t CO2e/ha "
            f"(área del anillo: {area:.2f} ha)")

    faltantes = [c for c in ["anillo", "anio", "perdida_ha"] if c not in df_perdida_sin_traslape.columns]
    if faltantes:
        raise ValueError(f"Faltan columnas {faltantes} en el CSV de pérdida sin traslape.")

    nombres_anillo_reales = set(nombre_anillo.values())

    # --- PASO 1: solo años reales (enteros) -- las filas 'TOTAL ...' del CSV de entrada se ignoran aquí a
    # propósito. El factor de combustión es por-año, así que la única forma correcta de armar un total es
    # SUMAR los años ya calculados uno por uno -- nunca re-multiplicar hectáreas-ya-sumadas x un solo factor,
    # eso perdería el descuento de combustión de cualquier año de incendio que caiga dentro del periodo
    # (bug real, encontrado y corregido en este mismo módulo antes de entregarlo).
    filas = []
    perdida_ha_total_entrada = {}  # {anillo: perdida_ha de la fila 'TOTAL {min}-{max}' del CSV de entrada, para validar
    for _, row in df_perdida_sin_traslape.iterrows():
        anillo = row["anillo"]
        if anillo not in nombres_anillo_reales:
            continue  # 'TOTAL (suma sin traslape)' del archivo de entrada -- no hace falta, se recalcula abajo
        anio_raw = row["anio"]
        try:
            anio_int = int(anio_raw)
        except (ValueError, TypeError):
            perdida_ha_total_entrada[anillo] = float(row["perdida_ha"])  # fila 'TOTAL {min}-{max}' -- guardar para validar
            continue

        d = densidad.get(anillo)
        if d is None:
            continue
        perdida_ha = float(row["perdida_ha"])
        factor = factor_combustion_por_anio.get(anio_int, 1.0)
        causa = causa_por_anio.get(anio_int, causa_default)
        co2e_agb = round(perdida_ha * d["densidad_agb_tco2e_ha"] * factor, 2)
        co2e_gedi = round(perdida_ha * d["densidad_gedi_tco2e_ha"] * factor, 2) if d["densidad_gedi_tco2e_ha"] is not None else None
        filas.append({
            "anio": anio_int, "anillo": anillo, "perdida_ha": perdida_ha,
            "densidad_agb_tco2e_ha": d["densidad_agb_tco2e_ha"], "factor_combustion_aplicado": factor,
            "co2e_asociado_agb_t": co2e_agb,
            "densidad_gedi_tco2e_ha": d["densidad_gedi_tco2e_ha"], "co2e_asociado_gedi_t": co2e_gedi,
            "causa_probable": causa,
        })

    df_out = pd.DataFrame(filas)
    if df_out.empty:
        raise ValueError("No se encontró ningún año real (entero) en df_perdida_sin_traslape -- nada que cruzar.")
    anio_min, anio_max = int(df_out["anio"].min()), int(df_out["anio"].max())
    etiqueta_periodo = f"TOTAL {anio_min}-{anio_max}"

    # --- PASO 2: fila 'TOTAL (suma sin traslape)' por cada AÑO REAL -- suma directa entre anillos, no se traslapan ---
    filas_total_anio = []
    for anio, grupo in df_out.groupby("anio"):
        causas_del_anio = grupo["causa_probable"].unique()
        causa_total = causas_del_anio[0] if len(causas_del_anio) == 1 else \
            "MIXTO -- ver por anillo (distinta causa/factor entre anillos este año)"
        filas_total_anio.append({
            "anio": anio, "anillo": "TOTAL (suma sin traslape)",
            "perdida_ha": round(grupo["perdida_ha"].sum(), 4),
            "densidad_agb_tco2e_ha": None, "factor_combustion_aplicado": None,
            "co2e_asociado_agb_t": round(grupo["co2e_asociado_agb_t"].sum(), 2),
            "densidad_gedi_tco2e_ha": None,
            "co2e_asociado_gedi_t": round(grupo["co2e_asociado_gedi_t"].sum(), 2) if grupo["co2e_asociado_gedi_t"].notna().all() else None,
            "causa_probable": causa_total,
        })

    # --- PASO 3: GRAN TOTAL por anillo (suma de TODOS los años reales de ese anillo -- respeta el factor de
    # cada año porque cada fila de 'filas' ya lo tiene aplicado) + el gran total general ---
    filas_gran_total = []
    for anillo in nombres_anillo_reales:
        grupo = df_out[df_out["anillo"] == anillo]
        if grupo.empty:
            continue
        causas = grupo["causa_probable"].unique()
        causa_total = causas[0] if len(causas) == 1 else "MIXTO -- ver por año (distinta causa/factor entre años)"
        suma_ha = round(grupo["perdida_ha"].sum(), 4)
        esperado = perdida_ha_total_entrada.get(anillo)
        if esperado is not None and abs(esperado - suma_ha) > 0.01:
            log(f"  AVISO: '{anillo}' -- la suma de los años reales ({suma_ha} ha) no coincide con la fila "
                f"'TOTAL {anio_min}-{anio_max}' del CSV de pérdida de entrada ({esperado} ha) -- revisa si "
                "los dos archivos (carbono y pérdida) cubren exactamente el mismo periodo.", nivel="WARN")
        filas_gran_total.append({
            "anio": etiqueta_periodo, "anillo": anillo, "perdida_ha": suma_ha,
            "densidad_agb_tco2e_ha": None, "factor_combustion_aplicado": None,
            "co2e_asociado_agb_t": round(grupo["co2e_asociado_agb_t"].sum(), 2),
            "densidad_gedi_tco2e_ha": None,
            "co2e_asociado_gedi_t": round(grupo["co2e_asociado_gedi_t"].sum(), 2) if grupo["co2e_asociado_gedi_t"].notna().all() else None,
            "causa_probable": causa_total,
        })
    df_gran_total_anillos = pd.DataFrame(filas_gran_total)
    filas_gran_total.append({
        "anio": etiqueta_periodo, "anillo": "TOTAL (suma sin traslape)",
        "perdida_ha": round(df_gran_total_anillos["perdida_ha"].sum(), 4),
        "densidad_agb_tco2e_ha": None, "factor_combustion_aplicado": None,
        "co2e_asociado_agb_t": round(df_gran_total_anillos["co2e_asociado_agb_t"].sum(), 2),
        "densidad_gedi_tco2e_ha": None,
        "co2e_asociado_gedi_t": round(df_gran_total_anillos["co2e_asociado_gedi_t"].sum(), 2),
        "causa_probable": "MIXTO -- ver por año/anillo" if len(set(df_gran_total_anillos["causa_probable"])) > 1
                           else df_gran_total_anillos["causa_probable"].iloc[0],
    })

    df_out = pd.concat([df_out, pd.DataFrame(filas_total_anio), pd.DataFrame(filas_gran_total)], ignore_index=True)
    df_out["anio"] = df_out["anio"].astype(str)
    df_out = df_out.sort_values(["anio", "anillo"]).reset_index(drop=True)

    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    os.makedirs(carpeta_salida, exist_ok=True)
    csv_out = os.path.join(carpeta_salida, f"co2e_asociado_perdida_{id_proyecto.lower()}.csv")
    df_out.to_csv(csv_out, index=False)
    log(f"CSV cruce carbono x pérdida (anillo exclusivo) guardado en: {csv_out}", nivel="OK")

    # --- resumen en log, solo sobre las filas de años reales (no las de período 'TOTAL 20XX-20XX') ---
    es_anio_real = pd.to_numeric(df_out["anio"], errors="coerce").notna()
    total_agb = df_out.loc[es_anio_real & (df_out["anillo"] == "TOTAL (suma sin traslape)"), "co2e_asociado_agb_t"].sum()
    total_gedi = df_out.loc[es_anio_real & (df_out["anillo"] == "TOTAL (suma sin traslape)"), "co2e_asociado_gedi_t"].sum()
    log(f"CO2e asociado al total de la pérdida confirmada (suma de todos los años, todos los anillos): "
        f"AGB={total_agb:,.0f} t CO2e | GEDI={total_gedi:,.0f} t CO2e", nivel="OK")

    return df_out, csv_out


# ==============================================================================
# --- ORQUESTADOR CON ARCHIVOS REALES ------------------------------------------
# ==============================================================================
def procesar_carbono_perdida_real(carbono_csv, perdida_sin_traslape_csv, id_proyecto, carpeta_salida=None,
                                   causa_default="tala/desmonte u otra causa no confirmada (remocion completa asumida, cota superior)",
                                   factor_combustion_por_anio=None, causa_por_anio=None):
    if not os.path.exists(carbono_csv):
        raise FileNotFoundError(f"No se encontró el CSV de carbono: {carbono_csv}")
    if not os.path.exists(perdida_sin_traslape_csv):
        raise FileNotFoundError(f"No se encontró el CSV de pérdida sin traslape: {perdida_sin_traslape_csv}")

    df_carbono = pd.read_csv(carbono_csv)
    df_perdida = pd.read_csv(perdida_sin_traslape_csv)
    return cruzar_carbono_con_perdida(df_carbono, df_perdida, id_proyecto, carpeta_salida=carpeta_salida,
                                       causa_default=causa_default,
                                       factor_combustion_por_anio=factor_combustion_por_anio,
                                       causa_por_anio=causa_por_anio)


# ==============================================================================
# --- DEMO: sintético, sin archivos ---
# ==============================================================================
def demo():
    log("=== core.carbono_perdida --demo (sintético, sin archivos reales) ===")
    df_carbono = pd.DataFrame([
        {"zona": "nucleo", "buffer_m": 0, "area_2d_ha": 100.0, "co2e_t": 20000.0, "co2e_incertidumbre_t": 5000.0,
         "gedi_co2e_t": 22000.0, "gedi_co2e_incertidumbre_t": 2000.0,
         "co2e_incremental_t": 20000.0, "co2e_incremental_incertidumbre_t": 5000.0,
         "gedi_co2e_incremental_t": 22000.0, "gedi_co2e_incremental_incertidumbre_t": 2000.0},
        {"zona": "buffer_500m", "buffer_m": 500, "area_2d_ha": 300.0, "co2e_t": 50000.0, "co2e_incertidumbre_t": 9000.0,
         "gedi_co2e_t": 55000.0, "gedi_co2e_incertidumbre_t": 3500.0,
         "co2e_incremental_t": 30000.0, "co2e_incremental_incertidumbre_t": 10295.6,
         "gedi_co2e_incremental_t": 33000.0, "gedi_co2e_incremental_incertidumbre_t": 4031.1},
    ])
    df_perdida = pd.DataFrame([
        {"anillo": "nucleo (0m)", "anio": 2020, "perdida_ha": 2.0},
        {"anillo": "anillo_0-500m", "anio": 2020, "perdida_ha": 5.0},
        {"anillo": "TOTAL (suma sin traslape)", "anio": 2020, "perdida_ha": 7.0},
        {"anillo": "nucleo (0m)", "anio": "TOTAL 2020-2020", "perdida_ha": 2.0},
        {"anillo": "anillo_0-500m", "anio": "TOTAL 2020-2020", "perdida_ha": 5.0},
        {"anillo": "TOTAL (suma sin traslape)", "anio": "TOTAL 2020-2020", "perdida_ha": 7.0},
    ])
    df_out, _ = cruzar_carbono_con_perdida(df_carbono, df_perdida, "demo", carpeta_salida="/tmp/mn_demo_carbono_perdida")
    print(df_out.to_string(index=False))
    return df_out


# ==============================================================================
# --- CLI ---
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(description="Cruza carbono por zona con pérdida Hansen confirmada -- Motor Nacional")
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--carbono-csv", type=str,
                     help="resumen_terreno_y_carbono_*.csv o carbono_por_zona_*.csv (core/carbono.py)")
    ap.add_argument("--perdida-sin-traslape-csv", type=str,
                     help="deforestacion_resumen_sin_traslape_*.csv (core/deforestacion.py)")
    ap.add_argument("--id-proyecto", type=str)
    ap.add_argument("--carpeta-salida", type=str, default=None)
    ap.add_argument("--causa", type=str,
                     default="tala/desmonte u otra causa no confirmada (remocion completa asumida, cota superior)",
                     help="Causa por defecto para los años SIN validación de incendio -- se guarda como nota en "
                          "el CSV. Para el/los años CON incendio confirmado, usa --anios-incendio (ver abajo), "
                          "no cambies este default.")
    ap.add_argument("--anios-incendio", type=str, default=None,
                     help="Años (separados por coma) donde la causa confirmada es incendio (ver "
                          "core/validacion_incendios.py) -- ej. '2025'. Esos años usan --factor-combustion en vez "
                          "de remoción completa, porque un incendio no consume el 100% de la biomasa.")
    ap.add_argument("--factor-combustion", type=float, default=0.45,
                     help="Factor de combustión IPCC (0-1) aplicado SOLO a los años en --anios-incendio -- default "
                          "0.45, punto medio del rango 0.40-0.50 para bosque templado ya usado en el análisis de "
                          "incendio de Cofre de Perote.")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if not (args.carbono_csv and args.perdida_sin_traslape_csv and args.id_proyecto):
        ap.error("--carbono-csv, --perdida-sin-traslape-csv e --id-proyecto son obligatorios fuera de --demo")

    anios_incendio = [int(a) for a in args.anios_incendio.split(",")] if args.anios_incendio else []
    factor_combustion_por_anio = {a: args.factor_combustion for a in anios_incendio}
    causa_por_anio = {a: f"incendio confirmado (factor de combustion={args.factor_combustion}, "
                          "ver core/validacion_incendios.py)" for a in anios_incendio}

    procesar_carbono_perdida_real(
        carbono_csv=args.carbono_csv, perdida_sin_traslape_csv=args.perdida_sin_traslape_csv,
        id_proyecto=args.id_proyecto, carpeta_salida=args.carpeta_salida, causa_default=args.causa,
        factor_combustion_por_anio=factor_combustion_por_anio or None,
        causa_por_anio=causa_por_anio or None,
    )


if __name__ == "__main__":
    main()

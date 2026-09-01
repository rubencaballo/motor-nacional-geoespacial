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

from config import log, ZONAS_ANALISIS_M, PERCENTIL_CAUCE_HIDROLOGIA, CARPETA_SRTM


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

    # BUG REAL encontrado y corregido aquí: resumen_terreno_y_carbono_*.csv SIEMPRE trae, al final, una fila
    # 'TOTAL (anillo exclusivo, 0-Xm, sí sumable)' (la agrega core.carbono._agregar_fila_total_carbono) -- esa
    # fila comparte buffer_m con el ÚLTIMO anillo real (ej. ambas buffer_m=1000). Si esta función recibe el CSV
    # tal cual (sin que el llamador la filtre primero -- procesar_carbono_perdida_real() no lo hacía), el
    # buffer_m duplicado rompe zona_de_buffer/nombre_anillo: el nombre del ÚLTIMO anillo real queda sobrescrito
    # por algo como 'anillo_1000-1000m' en vez de 'anillo_500-1000m', y ese anillo real desaparece
    # silenciosamente de todo lo que dependa de este nombre (cruzar_carbono_con_perdida, el balance, el mapa 3D
    # de pérdida) -- exactamente el síntoma real reportado: el anillo de 1000m sin tarjeta ni hover. Se filtra
    # aquí, en la ÚNICA función reusada por todo el módulo, para que ningún llamador (presente o futuro) tenga
    # que acordarse de hacerlo por su cuenta.
    df_terreno_carbono = df_terreno_carbono[~df_terreno_carbono["zona"].astype(str).str.startswith("TOTAL")]

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
# --- BALANCE: stock actual (carbono.py) vs. CO2e ya liberado (este módulo) ---
# ==============================================================================
def generar_balance_stock_vs_perdida(df_terreno_carbono, df_co2e_perdida, id_proyecto, carpeta_salida=None):
    """Compara, por anillo, el carbono ALMACENADO hoy (stock, de
    core/carbono.py, anillo exclusivo) contra el CO2e ya LIBERADO por
    pérdida confirmada en el periodo evaluado (flujo acumulado, salida de
    cruzar_carbono_con_perdida() de este mismo módulo) -- para poder decir
    en una sola tabla "este anillo almacena tanto, y de eso, tanto por
    ciento ya se liberó según la deforestación confirmada".

    Es una lectura de dos resultados ya calculados, NO una fuente nueva de
    datos: reusa _anillo_exclusivo_de_carbono() (mismo stock que ya usa
    cruzar_carbono_con_perdida más arriba, así que no se puede desincronizar
    del número que ya reporta ese cruce) y las filas 'TOTAL {min}-{max}' que
    cruzar_carbono_con_perdida() ya deja en df_co2e_perdida (el gran total
    por anillo del periodo evaluado, ya en anillo exclusivo).

    OJO -- el porcentaje puede superar 100: el stock es una FOTO de un año
    (config.CARBONO_ANIO, ej. 2022) y la pérdida puede cubrir años
    posteriores a esa foto (ej. 2010-2025 incluye 2023-2025, después de la
    foto de stock). Eso no es un error del cálculo, es información real --
    dice que la pérdida acumulada ya superó el inventario más reciente que
    tenemos. Se deja tal cual, sin recortar el porcentaje a 100, a
    propósito (ver docstring del módulo -- no escondemos números feos)."""
    anillos, zona_de_buffer, nombre_anillo = _anillo_exclusivo_de_carbono(df_terreno_carbono)

    df_periodo = df_co2e_perdida[df_co2e_perdida["anio"].astype(str).str.startswith("TOTAL ")].copy()
    if df_periodo.empty:
        raise ValueError("df_co2e_perdida no trae ninguna fila 'TOTAL {min}-{max}' -- ¿es la salida de "
                          "cruzar_carbono_con_perdida() (o el CSV co2e_asociado_perdida_*.csv que ese produce)?")
    etiqueta_periodo = df_periodo["anio"].iloc[0]

    filas = []
    anillos_omitidos = []
    for buf_m, datos in anillos.items():
        anillo = nombre_anillo[buf_m]
        fila_perdida = df_periodo[df_periodo["anillo"] == anillo]
        if fila_perdida.empty:
            anillos_omitidos.append(anillo)
            log(f"  AVISO: '{anillo}' no tiene fila '{etiqueta_periodo}' en el CSV de pérdida -- se OMITE del "
                f"balance (no es 0% liberado, es que no hay dato). Revisa que co2e_asociado_perdida_*.csv "
                f"realmente incluya este anillo (¿corriste cruzar_carbono_con_perdida()/deforestacion.py con las "
                f"mismas zonas_m que carbono.py, ej. --zonas 0,500,1000 en los tres módulos?).", nivel="WARN")
            continue
        fila_perdida = fila_perdida.iloc[0]
        stock_agb, stock_gedi = datos["co2e_t"], datos["gedi_co2e_t"]
        perdido_agb = fila_perdida["co2e_asociado_agb_t"]
        perdido_gedi = fila_perdida["co2e_asociado_gedi_t"]
        filas.append({
            "anillo": anillo, "area_ha": datos["area_ha"],
            "stock_agb_tco2e": stock_agb, "stock_agb_incertidumbre_tco2e": datos["co2e_incertidumbre_t"],
            "co2e_liberado_agb_t": perdido_agb,
            "pct_stock_agb_liberado": round(perdido_agb / stock_agb * 100, 2) if stock_agb else None,
            "stock_gedi_tco2e": stock_gedi, "stock_gedi_incertidumbre_tco2e": datos["gedi_co2e_incertidumbre_t"],
            "co2e_liberado_gedi_t": perdido_gedi,
            "pct_stock_gedi_liberado": round(perdido_gedi / stock_gedi * 100, 2) if (stock_gedi and pd.notna(perdido_gedi)) else None,
            "periodo_perdida_evaluado": etiqueta_periodo,
        })

    if not filas:
        raise ValueError("Ningún anillo de carbono coincidió con los anillos del CSV de pérdida -- revisa que "
                          "ambos archivos sean del mismo id_proyecto/zonas_m.")

    df_balance = pd.DataFrame(filas)

    # El TOTAL de este balance SIEMPRE se recalcula sumando las filas de df_balance que sí se pudieron armar
    # (las de arriba, una por anillo con dato) -- nunca se copia la fila 'TOTAL (suma sin traslape)' que ya
    # trae df_co2e_perdida tal cual, aunque exista. Por qué: esa fila se calculó en cruzar_carbono_con_perdida()
    # en su momento, y si algún anillo queda fuera de ESTE balance (anillos_omitidos, arriba) porque no tiene
    # fila de pérdida, el TOTAL tiene que verse EXACTAMENTE como la suma de las tarjetas que sí se muestran --
    # cualquier otra cosa sería un número que no cuadra con lo que el usuario ve enfrente (la propia duda de
    # "¿no tendrá muñeca rusa?" que motivó este chequeo). La fila de la CSV solo se usa como VALIDACIÓN cruzada
    # (si difiere del recálculo, es una señal real de inconsistencia entre archivos -- se avisa, no se oculta).
    fila_total_perdida = df_periodo[df_periodo["anillo"] == "TOTAL (suma sin traslape)"]
    stock_agb_total = df_balance["stock_agb_tco2e"].sum()
    stock_gedi_total = df_balance["stock_gedi_tco2e"].sum()
    perdido_agb_total = df_balance["co2e_liberado_agb_t"].sum()
    perdido_gedi_total = (df_balance["co2e_liberado_gedi_t"].sum()
                           if df_balance["co2e_liberado_gedi_t"].notna().all() else None)
    if not fila_total_perdida.empty:
        agb_csv = fila_total_perdida["co2e_asociado_agb_t"].iloc[0]
        if anillos_omitidos:
            log(f"  El TOTAL de este balance ({perdido_agb_total:,.2f} t CO2e) es la suma SOLO de los anillos "
                f"con dato -- la fila 'TOTAL (suma sin traslape)' del CSV de pérdida ({agb_csv:,.2f} t CO2e) "
                f"probablemente incluye más anillos y NO se usa aquí a propósito.", nivel="WARN")
        elif abs(agb_csv - perdido_agb_total) > max(1.0, abs(agb_csv) * 0.001):
            log(f"  AVISO: el TOTAL recalculado aquí ({perdido_agb_total:,.2f} t CO2e, suma de los anillos de "
                f"este balance) no coincide con la fila 'TOTAL (suma sin traslape)' del CSV de pérdida "
                f"({agb_csv:,.2f} t CO2e) -- revisa si los anillos de carbono_csv y de perdida_sin_traslape_csv "
                "son exactamente los mismos.", nivel="WARN")
    unc_agb_total = float(np.sqrt((df_balance["stock_agb_incertidumbre_tco2e"].dropna() ** 2).sum()))
    unc_gedi_vals = df_balance["stock_gedi_incertidumbre_tco2e"].dropna()
    unc_gedi_total = float(np.sqrt((unc_gedi_vals ** 2).sum())) if len(unc_gedi_vals) else None

    df_balance = pd.concat([df_balance, pd.DataFrame([{
        "anillo": "TOTAL (suma sin traslape)", "area_ha": round(df_balance["area_ha"].sum(), 4),
        "stock_agb_tco2e": round(stock_agb_total, 1), "stock_agb_incertidumbre_tco2e": round(unc_agb_total, 1),
        "co2e_liberado_agb_t": perdido_agb_total,
        "pct_stock_agb_liberado": round(perdido_agb_total / stock_agb_total * 100, 2) if stock_agb_total else None,
        "stock_gedi_tco2e": round(stock_gedi_total, 1),
        "stock_gedi_incertidumbre_tco2e": round(unc_gedi_total, 1) if unc_gedi_total is not None else None,
        "co2e_liberado_gedi_t": perdido_gedi_total,
        "pct_stock_gedi_liberado": (round(perdido_gedi_total / stock_gedi_total * 100, 2)
                                     if (stock_gedi_total and pd.notna(perdido_gedi_total)) else None),
        "periodo_perdida_evaluado": etiqueta_periodo,
    }])], ignore_index=True)

    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    os.makedirs(carpeta_salida, exist_ok=True)
    csv_out = os.path.join(carpeta_salida, f"balance_stock_vs_perdida_{id_proyecto.lower()}.csv")
    df_balance.to_csv(csv_out, index=False)
    if anillos_omitidos:
        log(f"  Balance INCOMPLETO -- faltan {len(anillos_omitidos)} anillo(s) sin dato de pérdida: "
            f"{', '.join(anillos_omitidos)}. El CSV/mapa que use este balance lo mostrará tal cual (sin ese "
            "anillo, y el TOTAL sin ese anillo) -- no se rellena con 0 para no fingir 'sin pérdida' donde en "
            "realidad es 'sin dato'.", nivel="WARN")
    log(f"Balance stock (carbono.py) vs. CO2e liberado ({etiqueta_periodo}) por anillo guardado en: {csv_out}",
        nivel="OK")
    # anillos_omitidos va en .attrs (metadata de pandas, no una columna) para no romper a quien ya
    # desempaqueta esta función como (df_balance, csv_out) -- generar_mapa_3d_perdida_carbono() lo lee de aquí
    # para poder avisar del hueco directamente en la página, en vez de solo en el log.
    df_balance.attrs["anillos_omitidos"] = anillos_omitidos
    return df_balance, csv_out


# ==============================================================================
# --- DESGLOSE POR CAUSA, POR ANILLO (tala/desmonte vs. incendio confirmado) ---
# ==============================================================================
def _desglose_causa_por_anillo(df_co2e_perdida, nombres_anillo_reales):
    """Agrupa las filas de AÑO REAL (no las 'TOTAL...') de
    co2e_asociado_perdida_*.csv por (anillo, categoría de causa), sumando
    hectáreas y CO2e -- para poder mostrar, en la tarjeta de cada anillo,
    cuánto de lo liberado vino de tala/desmonte vs. de incendio confirmado.

    Por qué no basta con la fila 'TOTAL {min}-{max}' que ya deja
    cruzar_carbono_con_perdida(): en cuanto un anillo tuvo años de las dos
    causas (como Cofre de Perote: 2010-2024 tala, 2025 incendio), esa fila
    queda etiquetada 'MIXTO -- ver por año' a propósito (no inventa un
    promedio de dos factores de combustión distintos) -- aquí sí se separa,
    sumando cada categoría por separado, que es justo lo que hace falta
    para mostrarlo en una tarjeta.

    categoria_causa: colapsa el texto largo de causa_probable (que incluye
    el factor de combustión aplicado, distinto según el año) a categorías
    legibles -- 'Tala/desmonte', 'Incendio confirmado', o 'Otra causa' si
    algún día se agrega una tercera. No pierde información: el texto
    completo, con el factor exacto, sigue disponible en el CSV crudo.

    Devuelve {anillo: [{"anillo", "categoria_causa", "perdida_ha",
    "co2e_agb_t", "co2e_gedi_t"}, ...]}, cada lista ordenada de mayor a
    menor CO2e AGB (para que la causa dominante salga primero en la
    tarjeta)."""
    es_anio_real = pd.to_numeric(df_co2e_perdida["anio"], errors="coerce").notna()
    df_real = df_co2e_perdida[es_anio_real & df_co2e_perdida["anillo"].isin(nombres_anillo_reales)].copy()
    if df_real.empty:
        return {}

    def _categoria(causa):
        c = str(causa).lower()
        if c.startswith("incendio"):
            return "Incendio confirmado"
        if c.startswith("tala"):
            return "Tala/desmonte"
        return "Otra causa"

    df_real["categoria_causa"] = df_real["causa_probable"].apply(_categoria)
    agrupado = df_real.groupby(["anillo", "categoria_causa"], as_index=False).agg(
        perdida_ha=("perdida_ha", "sum"),
        co2e_agb_t=("co2e_asociado_agb_t", "sum"),
        co2e_gedi_t=("co2e_asociado_gedi_t", "sum"),
    )

    desglose = {}
    for anillo, grupo in agrupado.groupby("anillo"):
        desglose[anillo] = grupo.sort_values("co2e_agb_t", ascending=False).to_dict("records")
    return desglose


# ==============================================================================
# --- ORQUESTADOR CON ARCHIVOS REALES ------------------------------------------
# ==============================================================================
def procesar_carbono_perdida_real(carbono_csv, perdida_sin_traslape_csv, id_proyecto, carpeta_salida=None,
                                   causa_default="tala/desmonte u otra causa no confirmada (remocion completa asumida, cota superior)",
                                   factor_combustion_por_anio=None, causa_por_anio=None, generar_balance=False):
    if not os.path.exists(carbono_csv):
        raise FileNotFoundError(f"No se encontró el CSV de carbono: {carbono_csv}")
    if not os.path.exists(perdida_sin_traslape_csv):
        raise FileNotFoundError(f"No se encontró el CSV de pérdida sin traslape: {perdida_sin_traslape_csv}")

    df_carbono = pd.read_csv(carbono_csv)
    df_perdida = pd.read_csv(perdida_sin_traslape_csv)
    df_out, csv_out = cruzar_carbono_con_perdida(df_carbono, df_perdida, id_proyecto, carpeta_salida=carpeta_salida,
                                                  causa_default=causa_default,
                                                  factor_combustion_por_anio=factor_combustion_por_anio,
                                                  causa_por_anio=causa_por_anio)
    if generar_balance:
        generar_balance_stock_vs_perdida(df_carbono, df_out, id_proyecto, carpeta_salida=carpeta_salida)
    return df_out, csv_out


# ==============================================================================
# --- MAPA 3D: CO2e LIBERADO por anillo (el "espejo" del mapa de carbono.py) ---
# ==============================================================================
# core/carbono.py ya tiene un mapa 3D que muestra, por anillo, cuánto CO2e
# ALMACENA hoy. Este es su espejo: mismo terreno, mismo estilo de tarjeta
# (core/reportes_html.py, sin duplicar diseño), pero mostrando cuánto CO2e
# ya se LIBERÓ a la atmósfera por la deforestación confirmada (Hansen) de
# cada anillo -- y cuántas hectáreas de ESE anillo se perdieron, no cuánto
# mide el anillo (por eso aquí 'area_ha' de la tarjeta es hectáreas
# PERDIDAS, no el tamaño total del anillo -- ver docstring de
# _construir_html_perdida_con_tarjetas).
def generar_mapa_3d_perdida_carbono(geojson_path, id_proyecto, zonas_m=None, percentil_cauce=None,
                                     carpeta_salida=None, carpeta_srtm=None):
    """Regenera el mismo terreno 3D de geomatica.py (reusa el .tif de SRTM
    ya descargado, no vuelve a bajar ni a consultar Earth Engine) pero con
    tarjetas de CO2e LIBERADO por anillo en vez de CO2e almacenado.

    Requiere que ya existan, en carpeta_salida (corridos antes, en este
    orden -- si falta alguno se explica en el log y se devuelve None):
      1. core/carbono.py --mapa-3d (o al menos combinar_con_geomatica())
         -> resumen_terreno_y_carbono_{id}.csv (carbono almacenado hoy)
      2. cruzar_carbono_con_perdida() de este mismo módulo (--balance no
         hace falta, esta función corre generar_balance_stock_vs_perdida()
         internamente para no desincronizarse del % que muestra)
         -> co2e_asociado_perdida_{id}.csv (CO2e liberado por año/anillo)

    NO propaga incertidumbre en el CO2e liberado (a diferencia del mapa de
    stock): cruzar_carbono_con_perdida() multiplica una densidad puntual
    (t CO2e/ha del anillo) por hectáreas perdidas, sin arrastrar la
    desviación estándar del dataset satelital -- se documenta en la nota al
    pie de la página en vez de inventar un número de incertidumbre que no
    se calculó."""
    from core import geomatica, reportes_html

    zonas_m = zonas_m if zonas_m is not None else ZONAS_ANALISIS_M
    carpeta_salida = carpeta_salida or os.path.expanduser(f"~/resultados_{id_proyecto.lower()}")
    carpeta_srtm = carpeta_srtm or CARPETA_SRTM

    path_terreno_carbono = os.path.join(carpeta_salida, f"resumen_terreno_y_carbono_{id_proyecto.lower()}.csv")
    path_perdida = os.path.join(carpeta_salida, f"co2e_asociado_perdida_{id_proyecto.lower()}.csv")
    if not os.path.exists(path_terreno_carbono) or not os.path.exists(path_perdida):
        log("No se pudo generar el mapa 3D de pérdida: falta "
            f"{'resumen_terreno_y_carbono_*.csv (corre core.carbono --mapa-3d primero)' if not os.path.exists(path_terreno_carbono) else ''} "
            f"{'co2e_asociado_perdida_*.csv (corre cruzar_carbono_con_perdida()/este módulo primero)' if not os.path.exists(path_perdida) else ''}",
            nivel="WARN")
        return None

    df_terreno_carbono = pd.read_csv(path_terreno_carbono)
    df_terreno_carbono = df_terreno_carbono[~df_terreno_carbono["zona"].astype(str).str.startswith("TOTAL")]
    df_co2e_perdida = pd.read_csv(path_perdida)

    df_balance, _ = generar_balance_stock_vs_perdida(df_terreno_carbono, df_co2e_perdida, id_proyecto,
                                                       carpeta_salida=carpeta_salida)
    anillos_omitidos = df_balance.attrs.get("anillos_omitidos", [])
    if anillos_omitidos:
        log(f"El mapa 3D de pérdida saldrá SIN tarjeta para: {', '.join(anillos_omitidos)} -- no hay dato de "
            "pérdida para ese/esos anillo(s) en co2e_asociado_perdida_*.csv (ver aviso arriba). El TOTAL de la "
            "página también quedará sin ese anillo, y se marca como incompleto en el propio HTML.", nivel="WARN")
    _, zona_de_buffer, nombre_anillo = _anillo_exclusivo_de_carbono(df_terreno_carbono)
    nombres_anillo_reales = set(nombre_anillo.values())
    desglose_causa = _desglose_causa_por_anillo(df_co2e_perdida, nombres_anillo_reales)

    incluir_gedi = df_balance["stock_gedi_tco2e"].notna().any()
    etiqueta_periodo = df_balance["periodo_perdida_evaluado"].iloc[0]

    geom_utm_nucleo, dst_array, meta_utm, utm_crs = geomatica.cargar_dem_utm(geojson_path, zonas_m, carpeta_srtm)
    hidrologia = geomatica.calcular_hidrologia_d8(
        dst_array, meta_utm, geom_utm_nucleo, zonas_m, max(zonas_m), utm_crs,
        percentil_cauce or PERCENTIL_CAUCE_HIDROLOGIA, carpeta_srtm, id_proyecto,
    )

    hover_por_zona = {}
    for buf_m, etiqueta in zona_de_buffer.items():
        anillo = nombre_anillo[buf_m]
        filas = df_balance[df_balance["anillo"] == anillo]
        if filas.empty:
            # Sin dato de pérdida para este anillo (ver anillos_omitidos) -- se avisa explícitamente en vez de
            # dejar el hover genérico "Límite de <zona>" de geomatica.py, que se leería como "no se marcó" en
            # vez de "no hay dato", justo la confusión que causó este mismo aviso en el resto de la página.
            hover_por_zona[etiqueta] = f"{anillo} -- SIN DATO de pérdida (falta en co2e_asociado_perdida_*.csv)"
            continue
        fila = filas.iloc[0]
        hover_txt = f"{anillo} -- CO2e liberado (AGB, {etiqueta_periodo}): {fila['co2e_liberado_agb_t']:,.0f} t"
        if pd.notna(fila.get("pct_stock_agb_liberado")):
            hover_txt += f" ({fila['pct_stock_agb_liberado']:.1f}% del stock)"
        if incluir_gedi and pd.notna(fila.get("co2e_liberado_gedi_t")):
            hover_txt += f"<br>{anillo} -- CO2e liberado (GEDI, {etiqueta_periodo}): {fila['co2e_liberado_gedi_t']:,.0f} t"
        hover_por_zona[etiqueta] = hover_txt
    capas_extra = geomatica.construir_anillos_visuales_3d(hidrologia, hover_por_zona=hover_por_zona)

    titulo_base = f"{id_proyecto} -- CO2e liberado por deforestación confirmada (Hansen {etiqueta_periodo})"
    fig = geomatica.generar_mapa_3d(
        hidrologia, id_proyecto, html_path=None,
        subtitulo="rota, acerca, y pasa el mouse sobre los anillos de color para ver su CO2e liberado",
        utm_crs=utm_crs, titulo_base="Vista 3D interactiva", capas_extra=capas_extra, devolver_fig=True,
    )

    html_completo = _construir_html_perdida_con_tarjetas(fig, df_balance, desglose_causa, zona_de_buffer,
                                                           nombre_anillo, id_proyecto, titulo_base,
                                                           etiqueta_periodo, incluir_gedi, anillos_omitidos)
    html_path = os.path.join(carpeta_salida, f"{id_proyecto.lower()}_3d_perdida_carbono.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_completo)
    log(f"Mapa 3D de CO2e liberado por pérdida (con tarjetas): {html_path}", nivel="OK")
    return html_path


# ==============================================================================
# --- TARJETAS HTML: CO2e liberado, fuera de la escena 3D ---------------------
# ==============================================================================
# Mismo diseño de tarjeta/página que core/carbono.py (vive en
# core/reportes_html.py, un solo lugar) -- aquí solo queda la lógica de QUÉ
# dato va en cada tarjeta, específica de "CO2e liberado".
def _construir_html_perdida_con_tarjetas(fig, df_balance, desglose_causa, zona_de_buffer, nombre_anillo,
                                          id_proyecto, titulo_base, etiqueta_periodo, incluir_gedi,
                                          anillos_omitidos=None):
    """OJO -- a diferencia de la tarjeta de stock (core/carbono.py), aquí el
    número chico arriba del valor principal ('area_ha' en tarjeta_html) NO
    es el tamaño del anillo -- es cuántas hectáreas DE ESE anillo se
    perdieron en el periodo evaluado (justo lo que se pidió: 'cuántas
    hectáreas se perdieron por anillo'). Se calcula sumando las hectáreas
    de desglose_causa (que ya vienen sin traslape, anillo exclusivo) --
    mismo origen de datos que el propio CO2e liberado, no un cálculo nuevo.

    anillos_omitidos: anillos de zona_de_buffer/nombre_anillo que
    generar_balance_stock_vs_perdida() no pudo calcular (sin fila de
    pérdida para ese anillo) -- NO se rellenan con 0 (sería fingir 'sin
    pérdida' donde en realidad es 'sin dato'), simplemente no tienen
    tarjeta. Para que eso no se lea como un total completo por accidente,
    el título de la tarjeta TOTAL y la nota al pie avisan explícitamente
    cuáles faltan.

    Igual que en core/carbono.py: este mapa arma su HTML con fig.to_html()
    directo (para envolverlo en las tarjetas), así que el clic -> Google
    Earth (geomatica.CLICK_ABRE_GOOGLE_EARTH_JS) no llega gratis desde
    geomatica.generar_mapa_3d() y hay que pasarlo aquí a mano."""
    from core import geomatica, reportes_html

    tarjetas = []
    total_ha_perdida = 0.0
    causa_total = {}

    for buf_m in sorted(zona_de_buffer):
        zona = zona_de_buffer[buf_m]
        anillo = nombre_anillo[buf_m]
        filas = df_balance[df_balance["anillo"] == anillo]
        if filas.empty:
            continue
        fila = filas.iloc[0]
        co2e_agb = fila["co2e_liberado_agb_t"]
        pct = fila.get("pct_stock_agb_liberado")

        entradas_causa = desglose_causa.get(anillo, [])
        perdida_ha_anillo = sum(e["perdida_ha"] for e in entradas_causa)
        total_ha_perdida += perdida_ha_anillo

        lineas_secundarias = []
        if incluir_gedi and pd.notna(fila.get("co2e_liberado_gedi_t")):
            lineas_secundarias.append(reportes_html.linea_secundaria_html(
                f"{fila['co2e_liberado_gedi_t']:,.0f} t CO2e", "GEDI L4A"))
        for entrada in entradas_causa:
            if entrada["co2e_agb_t"] <= 0 and entrada["perdida_ha"] <= 0:
                continue
            lineas_secundarias.append(reportes_html.linea_secundaria_html(
                f"{entrada['co2e_agb_t']:,.0f} t CO2e · {entrada['perdida_ha']:,.1f} ha", entrada["categoria_causa"]))
            cat = entrada["categoria_causa"]
            acc = causa_total.setdefault(cat, {"perdida_ha": 0.0, "co2e_agb_t": 0.0})
            acc["perdida_ha"] += entrada["perdida_ha"]
            acc["co2e_agb_t"] += entrada["co2e_agb_t"]

        nota = (f"{pct:.1f}% del carbono almacenado hoy (ESA CCI) ya liberado" if pd.notna(pct)
                else "ESA CCI Biomass")
        tarjetas.append(reportes_html.tarjeta_html(
            zona, reportes_html.COLORES_ZONA_HEX.get(zona, "#666"), f"{co2e_agb:,.0f}", "t CO2e liberado",
            nota_principal=nota, lineas_secundarias=lineas_secundarias, area_ha=perdida_ha_anillo,
        ))

    fila_total = df_balance[df_balance["anillo"] == "TOTAL (suma sin traslape)"].iloc[0]
    lineas_total = []
    if incluir_gedi and pd.notna(fila_total.get("co2e_liberado_gedi_t")):
        lineas_total.append(reportes_html.linea_secundaria_html(
            f"{fila_total['co2e_liberado_gedi_t']:,.0f} t CO2e", "GEDI L4A"))
    for cat, datos in sorted(causa_total.items(), key=lambda kv: kv[1]["co2e_agb_t"], reverse=True):
        lineas_total.append(reportes_html.linea_secundaria_html(
            f"{datos['co2e_agb_t']:,.0f} t CO2e · {datos['perdida_ha']:,.1f} ha", cat))

    pct_total = fila_total.get("pct_stock_agb_liberado")
    anillos_omitidos = anillos_omitidos or []
    nombre_total = (f"Total ({etiqueta_periodo}, anillos sumados)" if not anillos_omitidos
                     else f"Total ({etiqueta_periodo}, INCOMPLETO -- falta {', '.join(anillos_omitidos)})")
    tarjeta_total = reportes_html.tarjeta_html(
        "TOTAL", reportes_html.COLOR_TOTAL_HEX, f"{fila_total['co2e_liberado_agb_t']:,.0f}", "t CO2e liberado",
        nota_principal=(f"{pct_total:.1f}% del carbono almacenado hoy (ESA CCI) ya liberado"
                         if pd.notna(pct_total) else "ESA CCI Biomass"),
        lineas_secundarias=lineas_total, nombre_mostrado=nombre_total,
        area_ha=total_ha_perdida, es_total=True,
    )

    div_mapa = fig.to_html(full_html=False, include_plotlyjs=True, config={"displaylogo": False},
                            div_id="mapa3d", post_script=geomatica.CLICK_ABRE_GOOGLE_EARTH_JS)
    subtitulo = (f"CO2e liberado a la atmósfera por deforestación confirmada -- anillo exclusivo, sí sumable "
                 f"entre tarjetas (Hansen {etiqueta_periodo}, densidad de carbono ESA CCI Above-Ground Biomass"
                 f"{' + GEDI L4A' if incluir_gedi else ''})")
    aviso_faltantes = ""
    if anillos_omitidos:
        aviso_faltantes = (
            f"⚠ FALTA DATO para {', '.join(anillos_omitidos)}: co2e_asociado_perdida_*.csv no trae ninguna fila "
            f"'{etiqueta_periodo}' para ese/esos anillo(s), así que no tienen tarjeta y el TOTAL de esta página NO "
            "los incluye (no es que ese anillo no perdió nada -- es que no hay dato para calcularlo). Revisa que "
            "core/deforestacion.py y cruzar_carbono_con_perdida() hayan corrido con las mismas zonas_m que "
            "core/carbono.py (--zonas 0,500,1000) y que el CSV de pérdida sin traslape realmente incluya ese "
            "anillo, y vuelve a correr --mapa-3d-perdida. "
        )
    nota_pie = (
        aviso_faltantes +
        "Las hectáreas de cada tarjeta son las PERDIDAS dentro de ese anillo en todo el periodo evaluado, no el "
        "tamaño del anillo. El porcentaje compara la pérdida ya confirmada contra el carbono que el anillo "
        "almacena HOY (foto del año del dataset satelital) -- puede superar 100% cuando la pérdida acumulada "
        "incluye años posteriores a esa foto; no es un error de cálculo, es información real (ver docstring de "
        "generar_balance_stock_vs_perdida). El desglose por causa separa tala/desmonte (remoción casi completa de "
        "biomasa asumida) de incendio confirmado (~45% de la biomasa se combuste, factor IPCC) -- nunca se suman "
        "con el mismo peso. GEDI y ESA CCI son dos estimaciones independientes, no se promedian entre sí. A "
        "diferencia del mapa de carbono almacenado, aquí no se muestra incertidumbre (±): el cruce carbono x "
        "pérdida multiplica una densidad puntual por hectáreas perdidas, sin propagar la desviación estándar del "
        "dataset satelital."
    )
    return reportes_html.pagina_html_con_tarjetas(
        f"{id_proyecto} -- CO2e liberado por zona", titulo_base, subtitulo,
        "".join(tarjetas) + tarjeta_total, div_mapa, nota_pie,
    )


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

    df_balance, _ = generar_balance_stock_vs_perdida(df_carbono, df_out, "demo", carpeta_salida="/tmp/mn_demo_carbono_perdida")
    print("\n--- Balance stock vs. CO2e liberado (demo) ---")
    print(df_balance.to_string(index=False))
    return df_out, df_balance


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
    ap.add_argument("--balance", action="store_true",
                     help="Además del cruce año x anillo, genera balance_stock_vs_perdida_*.csv: por anillo, "
                          "cuánto CO2e almacena HOY (stock, carbono.py) vs. cuánto ya se liberó según la pérdida "
                          "confirmada del periodo (% del stock), ESA CCI y GEDI por separado.")
    ap.add_argument("--mapa-3d-perdida", action="store_true",
                     help="Genera el mapa 3D 'espejo' del de core/carbono.py: mismo terreno y mismo estilo de "
                          "tarjeta, pero con el CO2e YA LIBERADO por anillo (y hectáreas perdidas por anillo) en "
                          "vez del CO2e almacenado. Requiere --geojson y --id-proyecto, y que ya existan en "
                          "--carpeta-salida tanto resumen_terreno_y_carbono_*.csv (core.carbono --mapa-3d) como "
                          "co2e_asociado_perdida_*.csv (--carbono-csv/--perdida-sin-traslape-csv de este módulo, "
                          "corridos antes, sin necesidad de --balance).")
    ap.add_argument("--geojson", type=str, default=None, help="Ruta al GeoJSON del polígono núcleo (solo para --mapa-3d-perdida)")
    ap.add_argument("--zonas", type=str, default=None, help="Buffers en metros separados por coma, ej. '0,500,1000' (solo para --mapa-3d-perdida)")
    ap.add_argument("--percentil-cauce", type=float, default=None, help="Percentil D8 para declarar cauce (solo para --mapa-3d-perdida)")
    ap.add_argument("--carpeta-srtm", type=str, default=None, help="Carpeta con el SRTM ya descargado (solo para --mapa-3d-perdida)")
    args = ap.parse_args()

    if args.demo:
        demo()
        return

    if args.mapa_3d_perdida:
        if not (args.geojson and args.id_proyecto):
            ap.error("--mapa-3d-perdida requiere --geojson e --id-proyecto")
        zonas_m = [int(z) for z in args.zonas.split(",")] if args.zonas else None
        generar_mapa_3d_perdida_carbono(
            args.geojson, args.id_proyecto, zonas_m=zonas_m, percentil_cauce=args.percentil_cauce,
            carpeta_salida=args.carpeta_salida, carpeta_srtm=args.carpeta_srtm,
        )
        return

    if not (args.carbono_csv and args.perdida_sin_traslape_csv and args.id_proyecto):
        ap.error("--carbono-csv, --perdida-sin-traslape-csv e --id-proyecto son obligatorios fuera de --demo/--mapa-3d-perdida")

    anios_incendio = [int(a) for a in args.anios_incendio.split(",")] if args.anios_incendio else []
    factor_combustion_por_anio = {a: args.factor_combustion for a in anios_incendio}
    causa_por_anio = {a: f"incendio confirmado (factor de combustion={args.factor_combustion}, "
                          "ver core/validacion_incendios.py)" for a in anios_incendio}

    procesar_carbono_perdida_real(
        carbono_csv=args.carbono_csv, perdida_sin_traslape_csv=args.perdida_sin_traslape_csv,
        id_proyecto=args.id_proyecto, carpeta_salida=args.carpeta_salida, causa_default=args.causa,
        factor_combustion_por_anio=factor_combustion_por_anio or None,
        causa_por_anio=causa_por_anio or None, generar_balance=args.balance,
    )


if __name__ == "__main__":
    main()

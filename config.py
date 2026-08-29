#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Configuracion central del Motor Nacional v9.2 modular. Fuente unica de
constantes -- todos los demas modulos importan de aqui, nunca redefinen."""
import os
from datetime import datetime

DB_PATH = os.environ.get("IRDCLOUD_DB_PATH", "./motor_nacional_poderoso.db")
CHECKPOINT_PATH = os.environ.get("IRDCLOUD_CHECKPOINT_PATH", "./checkpoint_nacional_poderoso.json")
EXPORT_DIR = os.environ.get("IRDCLOUD_EXPORT_DIR", "./export_icm")
VERSION = "EINSTEIN_NACIONAL_v9.2_MODULAR"

HANSEN_DATASET = "UMD/hansen/global_forest_change_2025_v1_13"

CORTE_EUDR = 2020
CORTE_LOSSYEAR = 20
JRC_MIN_VERDE = 15.0
JRC_MIN_VERDE_NICFI = 10.0
AREA_MIN_HANSEN_HA = 0.5
MIN_PIXELES_CONFIANZA_ALTA = 30
MIN_PIXELES_CONFIANZA_MEDIA = 10
MIN_CLUSTER_HANSEN = 3
CAIDA_NDVI_MAX = 25.0
UMBRAL_AGUA = 80
UMBRAL_PENDIENTE = 30
UMBRAL_SILVO = (20, 40)
UMBRAL_SILVO_NDVI = 0.55
UMBRAL_CAFE_TREECOVER = (30, 70)
UMBRAL_CAFE_NDVI = 0.70
UMBRAL_CAFE_NDVI_STD = 0.08
ALTITUD_CAFE = (800, 1500)
ALTITUD_AGUACATE = (1500, 2400)
PESOS_PRIORIZACION = {'riesgo': 0.30, 'carbono': 0.25, 'aptitud': 0.20, 'marginacion': 0.15, 'conectividad': 0.10}

# --- core/geomatica.py: pipeline geomático/hidrológico por zonas (D8) ---
CARPETA_SRTM = os.environ.get("IRDCLOUD_SRTM_DIR", os.path.expanduser("~/srtm_temp/"))
ZONAS_ANALISIS_M = [0, 500, 1000]  # buffers en metros: nucleo, +500m, +1000m
PERCENTIL_CAUCE_HIDROLOGIA = 98  # percentil de acumulacion de flujo D8 para declarar cauce.
# Decision empirica (no teorica): se probo 92 vs 98 en 2 sitios geomorfologicamente
# distintos -- Ramsar 1601 Texolo (barranca angosta) y Cofre de Perote (volcan
# conico) -- validando contra la Red Hidrografica INEGI 1:50,000. En ambos, 98
# redujo drasticamente el "ruido" de cauces falsos en zonas planas/laderas anchas
# (ej. Texolo nucleo: 48.5% -> 73.1% de cauces D8 dentro de 30m de la red oficial).
# Si se agregan mas sitios y el patron cambia, reconsiderar este valor -- son
# solo 2 pruebas, no una ley general (ver core/validacion_hidrologica.py).

# --- core/validacion_hidrologica.py: comparacion D8 vs Red Hidrografica INEGI ---
TOLERANCIA_VALIDACION_HIDRO_M = 30.0  # metros -- igual a la resolucion del SRTM

# --- core/anp_lookup.py: busqueda de ANP por nombre en el shapefile nacional de CONANP ---
ANP_SHAPEFILE_NACIONAL = os.environ.get(
    "IRDCLOUD_ANP_SHAPEFILE",
    os.path.expanduser("~/MotorNacional/ANP_shapefile/232-ANP_ITRF08_19162026.shp"),
)

# --- core/carbono.py: estimacion de biomasa/carbono/CO2e por zona (Earth Engine) ---
# OJO: la banda 'agb' de ESA CCI viene en Mg de BIOMASA seca/ha, NO en
# carbono directo (a diferencia del dataset anterior NASA/ORNL). Por eso
# aqui la conversion es en DOS pasos: biomasa -> carbono (x FRACCION_CARBONO_
# BIOMASA) -> CO2e (x FACTOR_C_A_CO2). Ver el docstring de core/carbono.py.
CARBONO_DATASET_GEE = "ESA/CCI/Above_Ground_Biomass/V6_0"
CARBONO_ANIO = 2022  # anio mas reciente disponible en el dataset (agosto 2026)
FRACCION_CARBONO_BIOMASA = 0.47  # fraccion de carbono en biomasa seca (IPCC, valor tipico)
FACTOR_C_A_CO2 = 44 / 12  # razon molecular CO2/C

# GEDI L4A: mediciones LiDAR DIRECTAS por huella (~25m cada disparo), no un
# raster continuo interpolado como ESA CCI -- sirve para comparar "lo local"
# (mediciones reales, pero dispersas/con huecos) contra "lo satelital"
# (ESA CCI, cobertura completa pero modelada). Puede no haber ninguna huella
# GEDI dentro de un predio chico -- eso se reporta explicito, no se inventa.
GEDI_DATASET_L4A = "LARSE/GEDI/GEDI04_A_002_MONTHLY"
GEDI_ESCALA_M = 25  # resolucion nominal de la huella GEDI

# --- core/cuenca_completa.py: delimitacion de cuenca real entre dos sitios ---
MARGEN_CORREDOR_GRADOS = 0.05  # ~5.5km de margen extra alrededor del bbox que envuelve
# origen + punto de salida, para que la cuenca real delimitada por catchment() no se
# corte artificialmente en el borde del DEM del corredor.
PERCENTIL_SNAP_CAUCE = 95  # percentil de acumulacion de flujo D8 usado por snap_to_mask()
# para ajustar el punto de salida dado (coordenadas aproximadas) a la celda de cauce
# real mas cercana -- sin esto, catchment() puede delimitar una cuenca chueca o vacia.

os.makedirs(EXPORT_DIR, exist_ok=True)


def log(msg, nivel="INFO"):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {nivel}: {msg}")
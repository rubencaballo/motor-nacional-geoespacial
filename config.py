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

os.makedirs(EXPORT_DIR, exist_ok=True)


def log(msg, nivel="INFO"):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {nivel}: {msg}")

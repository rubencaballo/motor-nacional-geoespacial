#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Base de datos (SQLite WAL) y checkpoint por estado."""
import json
import os
import sqlite3
from datetime import datetime, timezone

from config import DB_PATH, CHECKPOINT_PATH


def init_db_nacional():
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute("PRAGMA cache_size=-100000;")
    conn.execute("""
        CREATE TABLE IF NOT EXISTS nacional_poderoso (
            id_predio TEXT PRIMARY KEY, estado TEXT, municipio TEXT, superficie_ha REAL,
            B_treecover REAL, B_ndvi REAL, B_ndvi_std REAL, B_textura REAL, B_altitud REAL, B_bfast INTEGER,
            JRC_pct REAL, ESRI_crops REAL, T_perdida INTEGER,
            D_ha_validada REAL, D_ha_descartada REAL, D_ha_historica_pre2020 REAL, D_ha_2025 REAL, D_jrc REAL, D_pend REAL,
            ndvi_2020 REAL, ndvi_2024 REAL, nbdi_2020 REAL, nbdi_2024 REAL,
            n_pixeles_hansen INTEGER, confianza_hansen TEXT,
            nicfi_disp INTEGER, nicfi_ndvi_ant REAL, nicfi_ndvi_rec REAL, nicfi_caida REAL,
            curp TEXT, rfc TEXT, net_mass REAL, legal_docs TEXT, production_date TEXT,
            E_elegible INTEGER, color TEXT, dictamen_corto TEXT, dictamen_largo TEXT, texto_evidencia TEXT,
            sistema TEXT, requiere_revision INTEGER, status_dds TEXT, score REAL,
            hash_geo TEXT, checksum_integridad TEXT, fecha TEXT, version TEXT
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_estado ON nacional_poderoso(estado)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_color ON nacional_poderoso(color)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_sistema ON nacional_poderoso(sistema)")
    conn.commit()
    conn.close()


def checkpoint_load():
    if os.path.exists(CHECKPOINT_PATH):
        with open(CHECKPOINT_PATH) as f:
            return json.load(f)
    return {"estados_completados": [], "total_predios": 0, "fecha_inicio": datetime.now(timezone.utc).isoformat()}


def checkpoint_save(chk):
    with open(CHECKPOINT_PATH, 'w') as f:
        json.dump(chk, f, indent=2)



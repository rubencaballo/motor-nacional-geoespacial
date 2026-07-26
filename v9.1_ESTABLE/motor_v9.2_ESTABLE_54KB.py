#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MOTOR NACIONAL - ICM-SADER v9.2
Fórmula: E = (B x T) - D  + Auditoría + Anti-falsos positivos + Escala nacional

============================== CHANGELOG v9.2 ==============================
Sobre v9.1, dos correcciones:

BUG DE TRANSMISIÓN (no de diseño) — dobles guiones bajos perdidos.
  En el texto que llegó de v9.1, __init__ y __main__ (y __geo_interface__)
  habían perdido uno de sus dos guiones bajos (probablemente un renderizador
  de texto interpretando "__texto__" como negritas de Markdown y comiéndose
  los guiones). Efecto real, verificado corriendo el archivo tal cual:
  'if __name__=="__main__"' roto significaba que main() NUNCA se ejecutaba
  -- el script no hacía nada al correrlo, sin importar qué tan bien estuviera
  programada la lógica de adentro. Ya corregido y verificado con --demo.

FIX 4 (crítica real, pendiente desde v9.1) — Traslapes nacionales.
  Antes: solo make_valid(), con una nota admitiendo que la eliminación
  jerárquica de traslapes no estaba implementada -> doble conteo de
  hectáreas posible si dos predios del RAN se solapan.
  Ahora: eliminar_traslapes_jerarquico() usa gpd.overlay(how='difference')
  de forma iterativa: ordena los polígonos por un campo de prioridad
  (antigüedad de título si existe; si no, se avisa explícitamente que se
  usa un criterio de respaldo arbitrario -- área descendente -- y que esto
  DEBE reemplazarse por un criterio legal real antes de un entregable
  oficial). El de mayor prioridad conserva toda su área; los siguientes
  se recortan donde ya hay territorio reclamado por alguien con más
  prioridad. Se activa con --resolver-traslapes [campo_prioridad opcional].

FIX 5 (crítica real, pendiente desde v9.1) — Hueco de escala 0.5-1.5ha.
  Antes: la única señal de "Hansen no es confiable aquí" era una hectárea
  fija (AREA_MIN_HANSEN_HA=0.5) y solo por debajo de ese corte se consultaba
  NICFI. Entre 0.5 y ~1.5ha, Hansen sigue siendo poco confiable (pocos
  píxeles de 30m dentro del polígono) pero el motor no lo reforzaba con
  nada.
  Ahora: en modo GEE real se cuenta el número real de píxeles Hansen
  válidos dentro de cada polígono (reduceRegion con Reducer.count(), no
  solo el promedio) y se deriva confianza_hansen ('baja' <10 píxeles,
  'media' 10-30, 'alta' >30 -- ~10 píxeles de 30m son ~0.9 ha, umbral
  documentado, no mágico). NICFI ahora se consulta (si --con-nicfi) cuando
  confianza_hansen != 'alta', no solo cuando superficie_ha<0.5 -- así el
  hueco 0.5-1.5ha también recibe refuerzo cuando hay NICFI disponible, y
  cuando no lo hay, el texto_evidencia ahora DICE explícitamente "confianza
  Hansen baja/media" en vez de callarlo.
==============================================================================
"""

import os, sys, json, hashlib, sqlite3, argparse, random, time
from datetime import datetime, timezone
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict

# ============== CONFIGURACIÓN ==============
DB_PATH = os.environ.get("IRDCLOUD_DB_PATH", "./motor_nacional_poderoso.db")
CHECKPOINT_PATH = os.environ.get("IRDCLOUD_CHECKPOINT_PATH", "./checkpoint_nacional_poderoso.json")
EXPORT_DIR = os.environ.get("IRDCLOUD_EXPORT_DIR", "./export_icm")
VERSION = "EINSTEIN_NACIONAL_v9.2_FIX4_TRASLAPES_FIX5_ESCALA"

HANSEN_DATASET = "UMD/hansen/global_forest_change_2025_v1_13"

CORTE_EUDR = 2020
CORTE_LOSSYEAR = 20
JRC_MIN_VERDE = 15.0
JRC_MIN_VERDE_NICFI = 10.0
AREA_MIN_HANSEN_HA = 0.5
MIN_PIXELES_CONFIANZA_ALTA = 30   # ~2.7 ha de Hansen 30m
MIN_PIXELES_CONFIANZA_MEDIA = 10  # ~0.9 ha de Hansen 30m -- umbral documentado, no mágico
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


@dataclass
class EvidenciaNICFI:
    disponible: bool = False
    ndvi_antiguo: Optional[float] = None
    ndvi_reciente: Optional[float] = None
    caida_pct: Optional[float] = None
    fecha_antigua: Optional[str] = None
    fecha_reciente: Optional[str] = None


@dataclass
class PredioNacional:
    id_predio: str
    estado: str = ""
    municipio: str = ""
    geometria: dict = field(default_factory=dict)
    superficie_ha: float = 0.0
    B_treecover: float = 0.0
    B_ndvi: float = 0.0
    B_ndvi_std: float = 0.0
    B_textura: float = 0.0
    B_altitud: float = 0.0
    B_bfast_anio: Optional[int] = None
    JRC_pct: float = 0.0
    ESRI_crops_pct: float = 0.0
    T_perdida: int = 0
    D_ha_validada: float = 0.0
    D_ha_descartada_aislada: float = 0.0
    D_ha_historica_pre2020: float = 0.0
    D_ha_2025: float = 0.0
    D_jrc_water: float = 0.0
    D_pendiente: float = 0.0
    ndvi_2020: float = 0.0
    ndvi_2024: float = 0.0
    nbdi_2020: float = 0.0
    nbdi_2024: float = 0.0
    nicfi: EvidenciaNICFI = field(default_factory=EvidenciaNICFI)
    n_pixeles_hansen: int = 0                    # FIX 5
    confianza_hansen: str = "no_evaluada"         # FIX 5: 'alta'/'media'/'baja'/'no_evaluada'
    curp: str = "PENDIENTE_APORTACION"
    rfc: str = "PENDIENTE_APORTACION"
    net_mass_kg: float = 0.0
    legal_docs: str = "PENDIENTE_APORTACION"
    production_date: str = "PENDIENTE_APORTAR_FECHA_COSECHA"
    E_elegible: bool = False
    color: str = "Amarillo"
    dictamen_corto: str = "REQUIERE_EVIDENCIA"
    dictamen_largo: str = ""
    texto_evidencia: str = ""
    sistema: str = "No clasificado"
    requiere_revision: bool = True
    status_dds: str = "DRAFT_LOCAL"
    score_prioridad: float = 0.0
    hash_geo: str = ""
    checksum_integridad: str = ""
    motivos_incompleto: List[str] = field(default_factory=list)
    version: str = VERSION
    fecha: str = ""

    def formula_elegibilidad_nacional(self):
        """
        FIX 5: el criterio de "Hansen no confiable" ya no es solo
        superficie_ha < 0.5. Ahora usa confianza_hansen (derivada del
        conteo real de píxeles válidos), que cubre también el hueco
        0.5-1.5ha que antes no recibía ningún refuerzo.
        """
        if self.D_ha_validada > 0.01:
            self.color = "Rojo"
            self.dictamen_corto = "NO_APTO_EXPORTACION"
            self.dictamen_largo = f"DEFORESTACIÓN POST-{CORTE_EUDR} DETECTADA {self.D_ha_validada:.4f} ha validada cluster>={MIN_CLUSTER_HANSEN}"
            self.texto_evidencia = (
                f"Hansen v1.13 {self.D_ha_validada:.4f} ha pérdida validada posterior al corte EUDR "
                f"({CORTE_EUDR}-12-31). Descartada aislada (post-corte, cluster<{MIN_CLUSTER_HANSEN}) "
                f"{self.D_ha_descartada_aislada:.4f} ha, no usada para Rojo. "
                f"Pérdida histórica pre-{CORTE_EUDR} (informativa, no afecta este dictamen): "
                f"{self.D_ha_historica_pre2020:.4f} ha."
            )
            self.requiere_revision = False
            self.E_elegible = False
            return

        hansen_no_confiable = self.confianza_hansen in ("baja", "media") or self.superficie_ha < AREA_MIN_HANSEN_HA

        if hansen_no_confiable:
            if self.nicfi.disponible and self.nicfi.ndvi_antiguo and self.nicfi.ndvi_reciente:
                caida = self.nicfi.caida_pct if self.nicfi.caida_pct is not None else 0
                estable = abs(caida) < CAIDA_NDVI_MAX
                ndvi_alto = self.nicfi.ndvi_antiguo > 0.60 and self.nicfi.ndvi_reciente > 0.60
                if estable and ndvi_alto and self.JRC_pct >= JRC_MIN_VERDE_NICFI:
                    self.color = "Verde"
                    self.dictamen_corto = "APTO_EXPORTACION"
                    self.dictamen_largo = f"CONFORME EUDR - Parcela {self.superficie_ha:.2f}ha (confianza Hansen {self.confianza_hansen}) rescatada NICFI 4.77m NDVI {self.nicfi.ndvi_antiguo:.2f}->{self.nicfi.ndvi_reciente:.2f} estable JRC {self.JRC_pct}%"
                    self.texto_evidencia = f"NICFI rescate: confianza Hansen {self.confianza_hansen} ({self.n_pixeles_hansen} píxeles) pero NDVI 4.77m estable sin caída >{CAIDA_NDVI_MAX}%."
                    self.requiere_revision = False
                    self.E_elegible = True
                    return
            self.color = "Amarillo"
            self.dictamen_corto = "REQUIERE_EVIDENCIA"
            self.dictamen_largo = f"CONFIANZA HANSEN {self.confianza_hansen.upper()} ({self.n_pixeles_hansen} píxeles, {self.superficie_ha:.2f}ha) - REQUIERE REVISIÓN O NICFI"
            self.texto_evidencia = (
                f"Parcela {self.superficie_ha:.4f} ha con {self.n_pixeles_hansen} píxeles Hansen válidos "
                f"(confianza {self.confianza_hansen}; umbral alta>={MIN_PIXELES_CONFIANZA_ALTA}px, "
                f"media>={MIN_PIXELES_CONFIANZA_MEDIA}px). JRC {self.JRC_pct:.1f}% ESRI {100-self.ESRI_crops_pct:.0f}% árboles. "
                f"Requiere foto campo o NICFI 4.77m."
            )
            self.requiere_revision = True
            self.E_elegible = False
            return

        if self.JRC_pct >= JRC_MIN_VERDE:
            self.color = "Verde"
            self.dictamen_corto = "APTO_EXPORTACION"
            self.dictamen_largo = f"CONFORME EUDR - JRC {self.JRC_pct:.1f}% >= {JRC_MIN_VERDE}% consistente ESRI sin pérdida Hansen post-corte"
            self.texto_evidencia = f"JRC {self.JRC_pct:.1f}% baseline consistente ESRI {self.ESRI_crops_pct:.0f}% crops. Sin pérdida Hansen posterior a {CORTE_EUDR}. Confianza Hansen: {self.confianza_hansen}."
            self.requiere_revision = False
            self.E_elegible = True
            return

        if self.nicfi.disponible and self.nicfi.ndvi_antiguo and self.nicfi.ndvi_antiguo > 0.65 and self.nicfi.ndvi_reciente and self.nicfi.ndvi_reciente > 0.65:
            if self.nicfi.caida_pct is not None and abs(self.nicfi.caida_pct) < CAIDA_NDVI_MAX:
                self.color = "Verde"
                self.dictamen_corto = "APTO_EXPORTACION"
                self.dictamen_largo = f"JRC {self.JRC_pct:.1f}% <15% pero NICFI 4.77m dosel estable {self.nicfi.ndvi_antiguo:.2f}->{self.nicfi.ndvi_reciente:.2f}"
                self.texto_evidencia = "Rescate JRC bajo con NICFI estable."
                self.requiere_revision = False
                self.E_elegible = True
                return

        self.color = "Amarillo"
        self.dictamen_corto = "REQUIERE_EVIDENCIA"
        self.dictamen_largo = f"PATRÓN JRC/ESRI NO CONCLUYENTE JRC {self.JRC_pct:.1f}% <{JRC_MIN_VERDE}% sin NICFI"
        self.texto_evidencia = f"JRC bajo sin refuerzo. Confianza Hansen: {self.confianza_hansen}. Hipótesis: café sombra preexistente o cobertura post-{CORTE_EUDR}. Requiere campo."
        self.requiere_revision = True
        self.E_elegible = False

    def formula_sistema_productivo_nacional(self):
        try:
            if ADAPTADOR_ML_GLOBAL and ADAPTADOR_ML_GLOBAL.conectado:
                sistema_ml, proba = ADAPTADOR_ML_GLOBAL.predecir_sistema_ml(self)
                if sistema_ml and proba >= 0.70:
                    return sistema_ml
        except Exception:
            pass

        if UMBRAL_CAFE_TREECOVER[0] <= self.B_treecover <= UMBRAL_CAFE_TREECOVER[1]:
            if self.B_ndvi > UMBRAL_CAFE_NDVI and self.B_ndvi_std < UMBRAL_CAFE_NDVI_STD and ALTITUD_CAFE[0] <= self.B_altitud <= ALTITUD_CAFE[1]:
                if self.B_textura > 80:
                    return "cafe_bajo_sombra"
        if self.B_bfast_anio and 2015 <= self.B_bfast_anio <= 2023:
            if self.B_treecover > 15 and self.B_ndvi > 0.70 and self.B_textura < 60 and ALTITUD_AGUACATE[0] <= self.B_altitud <= ALTITUD_AGUACATE[1]:
                return "aguacate"
        if 40 <= self.B_treecover <= 80 and self.B_textura < 45 and self.B_ndvi_std < 0.10:
            return "palma_aceite"
        if UMBRAL_SILVO[0] <= self.B_treecover <= UMBRAL_SILVO[1] and self.B_ndvi > UMBRAL_SILVO_NDVI and self.B_textura > 80:
            return "silvopastoril"
        if self.B_treecover > 60 and self.B_ndvi > 0.60:
            return "bosque_conservado"
        if self.B_treecover < 10 and self.B_ndvi > 0.30:
            return "ganaderia_convencional"
        if self.B_treecover < 15 and self.B_ndvi > 0.60:
            return "agricultura_anual"
        return "otro"

    def formula_prioridad_nacional(self):
        riesgo = min(1.0, self.D_ha_validada / 2.0)
        carbono = min(1.0, self.B_treecover / 60.0)
        return round(riesgo * PESOS_PRIORIZACION['riesgo'] + carbono * PESOS_PRIORIZACION['carbono'] + 0.45, 4)

    def validar_completitud(self):
        motivos = []
        if self.curp == "PENDIENTE_APORTACION": motivos.append("CURP faltante")
        if self.rfc == "PENDIENTE_APORTACION": motivos.append("RFC faltante")
        if self.net_mass_kg <= 0: motivos.append("Net mass 0 kg")
        if not self.legal_docs or self.legal_docs == "PENDIENTE_APORTACION": motivos.append("Legal docs faltante")
        if not self.production_date or self.production_date == "PENDIENTE_APORTAR_FECHA_COSECHA": motivos.append("Fecha producción faltante")
        self.motivos_incompleto = motivos
        self.status_dds = "FINAL" if (len(motivos) == 0 and not self.requiere_revision) else "DRAFT_LOCAL"
        return len(motivos) == 0

    def generar_auditoria(self):
        geom_str = json.dumps(self.geometria, sort_keys=True)
        self.hash_geo = hashlib.sha256(geom_str.encode()).hexdigest()[:16]
        payload = f"{self.id_predio}{self.hash_geo}{self.B_treecover:.2f}{self.D_ha_validada:.4f}{self.color}{VERSION}"
        self.checksum_integridad = hashlib.sha256(payload.encode()).hexdigest()[:16]
        self.fecha = datetime.now(timezone.utc).isoformat()


# ============== DB NACIONAL ==============
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


def motor_nacional_masivo(lista):
    for p in lista:
        p.formula_elegibilidad_nacional()
        p.sistema = p.formula_sistema_productivo_nacional()
        p.score_prioridad = p.formula_prioridad_nacional()
        p.validar_completitud()
        p.generar_auditoria()
    return lista


# ============== FIX 4: TRASLAPES JERÁRQUICOS ==============
def eliminar_traslapes_jerarquico(gdf, campo_prioridad=None):
    """
    Elimina traslapes entre polígonos de forma jerárquica: el polígono de
    mayor prioridad conserva toda su área; cada polígono siguiente se
    recorta (gpd.overlay how='difference') contra la UNIÓN de todos los de
    mayor prioridad ya procesados. Así ninguna hectárea se cuenta dos veces.

    campo_prioridad: nombre de columna en gdf que indica prioridad legal
    real (ej. fecha de título agrario -- más antiguo = más prioridad). Si
    no se proporciona, se usa área descendente como criterio de RESPALDO
    ARBITRARIO -- esto se advierte explícitamente porque un criterio de
    área no tiene ningún fundamento legal; solo evita que el cálculo
    truene mientras no se tenga el campo real.

    Devuelve un gdf nuevo con geometrías recortadas y una columna
    'ha_perdidas_por_traslape' documentando cuánta área se le restó a cada
    predio (para que quede trazable, no oculto).
    """
    import geopandas as gpd
    from shapely.ops import unary_union

    if len(gdf) <= 1:
        gdf = gdf.copy()
        gdf['ha_perdidas_por_traslape'] = 0.0
        return gdf

    gdf = gdf.copy()
    gdf['_area_original_m2'] = gdf.geometry.area

    if campo_prioridad and campo_prioridad in gdf.columns:
        gdf = gdf.sort_values(campo_prioridad, ascending=True).reset_index(drop=True)
        log(f"Traslapes: ordenando por prioridad real '{campo_prioridad}'.", "OK")
    else:
        gdf = gdf.sort_values('_area_original_m2', ascending=False).reset_index(drop=True)
        log("Traslapes: SIN campo de prioridad legal -- usando área descendente como "
            "respaldo ARBITRARIO. Esto debe reemplazarse por un criterio real "
            "(antigüedad de título u otro) antes de un entregable oficial a SADER/ICM.", "WARN")

    territorio_ocupado = None
    geometrias_finales = []
    ha_perdidas = []

    for idx, row in gdf.iterrows():
        geom = row.geometry
        if territorio_ocupado is not None and not territorio_ocupado.is_empty:
            geom_recortada = geom.difference(territorio_ocupado)
        else:
            geom_recortada = geom

        area_antes_m2 = geom.area
        area_despues_m2 = geom_recortada.area if not geom_recortada.is_empty else 0.0
        ha_perdidas.append(max(0.0, (area_antes_m2 - area_despues_m2)) / 10000)
        geometrias_finales.append(geom_recortada)

        territorio_ocupado = geom if territorio_ocupado is None else unary_union([territorio_ocupado, geom])

    gdf['geometry'] = geometrias_finales
    gdf['ha_perdidas_por_traslape'] = ha_perdidas
    gdf = gdf[~gdf.geometry.is_empty].copy()

    n_afectados = sum(1 for h in ha_perdidas if h > 0.0001)
    total_ha_recortadas = sum(ha_perdidas)
    log(f"Traslapes resueltos: {n_afectados} predios recortados, {total_ha_recortadas:.2f} ha "
        f"de doble conteo eliminadas.", "OK")

    return gdf


# ============== ML ==============
class AdaptadorMLPoderoso:
    def __init__(self, modelo_path=None):
        self.modelo_sistema = None
        self.modelo_riesgo = None
        self.scaler = None
        self.conectado = False
        self.features_nombre = ['B_treecover', 'B_ndvi', 'B_ndvi_std', 'B_textura', 'B_altitud',
                                 'JRC_pct', 'ESRI_crops', 'ndvi_2020', 'ndvi_2024', 'D_pendiente']
        if modelo_path and os.path.exists(modelo_path):
            try:
                import joblib
                data = joblib.load(modelo_path)
                self.modelo_sistema = data.get('modelo_sistema')
                self.modelo_riesgo = data.get('modelo_riesgo')
                self.scaler = data.get('scaler')
                self.conectado = True
                print(f"[ML] Modelos cargados: {modelo_path} - conectado {self.conectado}")
            except Exception as e:
                print(f"[ML] No se pudo cargar modelo: {e} - usando reglas")
        else:
            print("[ML] Sin modelo entrenado - usando reglas físicas (fallback)")

    def extraer_features(self, predio):
        return [
            predio.B_treecover, predio.B_ndvi, predio.B_ndvi_std, predio.B_textura,
            predio.B_altitud, predio.JRC_pct, 100 - predio.ESRI_crops_pct,
            predio.ndvi_2020, predio.ndvi_2024, predio.D_pendiente
        ]

    def predecir_sistema_ml(self, predio):
        if not self.conectado or self.modelo_sistema is None:
            return None, 0.0
        try:
            import numpy as np
            X = np.array(self.extraer_features(predio)).reshape(1, -1)
            if self.scaler:
                X = self.scaler.transform(X)
            pred = self.modelo_sistema.predict(X)[0]
            proba = max(self.modelo_sistema.predict_proba(X)[0]) if hasattr(self.modelo_sistema, 'predict_proba') else 0.85
            return str(pred), float(proba)
        except Exception as e:
            print(f"[ML] Error pred sistema: {e}")
            return None, 0.0

    def entrenar_desde_campo(self, csv_campo_path, output_modelo_path):
        """
        ADVERTENCIA: si csv_campo_path viene de extraer_BTD_de_SIAP_integrado
        en modo simulado, el accuracy es circular, no evidencia de que el
        modelo generalice a predios reales.
        """
        import pandas as pd
        from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
        from sklearn.preprocessing import StandardScaler
        from sklearn.model_selection import train_test_split
        import joblib

        df = pd.read_csv(csv_campo_path)
        print(f"[ML] Entrenando con {len(df)} muestras: {csv_campo_path}")

        X = df[self.features_nombre].values
        y_sistema = df['sistema'].values
        y_riesgo = df['riesgo'].values if 'riesgo' in df.columns else df['D_ha_validada'].values

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        X_train, X_test, y_sis_train, y_sis_test = train_test_split(X_scaled, y_sistema, test_size=0.2, random_state=42)

        modelo_sis = RandomForestClassifier(n_estimators=300, max_depth=15, min_samples_leaf=5, random_state=42, n_jobs=-1)
        modelo_sis.fit(X_train, y_sis_train)
        acc = modelo_sis.score(X_test, y_sis_test)
        print(f"[ML] Sistema - Accuracy test: {acc:.3f} (ver advertencia sobre datos circulares en el docstring)")

        modelo_riesgo = GradientBoostingRegressor(n_estimators=200, max_depth=6, random_state=42)
        modelo_riesgo.fit(X_scaled, y_riesgo)

        joblib.dump({'modelo_sistema': modelo_sis, 'modelo_riesgo': modelo_riesgo,
                     'scaler': scaler, 'features': self.features_nombre}, output_modelo_path)
        print(f"[ML] Modelo guardado: {output_modelo_path}")
        return output_modelo_path


ADAPTADOR_ML_GLOBAL = None


def init_ml_global(modelo_path=None):
    global ADAPTADOR_ML_GLOBAL
    ADAPTADOR_ML_GLOBAL = AdaptadorMLPoderoso(modelo_path=modelo_path)
    return ADAPTADOR_ML_GLOBAL


# ============== LANDTRENDR BÁSICO ==============
def calcular_landtrendr_breakyear(ee, geom, anio_min=2015, anio_max=2023):
    """Ver advertencias metodológicas completas en v9.1 -- sin cambios aquí."""
    try:
        def prep(img):
            qa = img.select('QA_PIXEL')
            nube = qa.bitwiseAnd(1 << 3).eq(0)
            sombra = qa.bitwiseAnd(1 << 4).eq(0)
            return img.updateMask(nube.And(sombra))

        def coleccion_landsat(aoi):
            l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(aoi)
            l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(aoi)
            l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(aoi)
            l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(aoi)
            return l5.merge(l7).merge(l8).merge(l9).map(prep)

        def ndvi_anual(anio, aoi, col_base):
            col = col_base.filter(ee.Filter.calendarRange(anio, anio, 'year'))

            def escalar(img):
                nir = img.select('SR_B5').multiply(0.0000275).add(-0.2)
                red = img.select('SR_B4').multiply(0.0000275).add(-0.2)
                return nir.subtract(red).divide(nir.add(red)).rename('NDVI')
            n = col.size()
            comp = ee.Image(ee.Algorithms.If(
                n.gt(0), col.map(escalar).median(),
                ee.Image(0).rename('NDVI').updateMask(ee.Image(0))
            ))
            return comp.set('system:time_start', ee.Date.fromYMD(anio, 6, 1).millis())

        col_base = coleccion_landsat(geom)
        anios = list(range(anio_min, anio_max + 1))
        serie = ee.ImageCollection([ndvi_anual(a, geom, col_base) for a in anios]) \
            .map(lambda img: img.multiply(1000).toShort().rename('NDVI').copyProperties(img, ['system:time_start']))

        lt = ee.Algorithms.TemporalSegmentation.LandTrendr(
            timeSeries=serie, maxSegments=4, spikeThreshold=0.9, vertexCountOvershoot=3,
            preventOneYearRecovery=True, recoveryThreshold=0.25, pvalThreshold=0.05,
            bestModelProportion=0.75, minObservationsNeeded=4
        )

        info = lt.select('LandTrendr').reduceRegion(
            reducer=ee.Reducer.first(), geometry=geom, scale=30, maxPixels=1e9, bestEffort=True
        ).getInfo()

        arr = info.get('LandTrendr')
        if not arr or len(arr) < 3:
            return None
        fila_anios = arr[0]
        fila_ajustado = arr[1]
        if not fila_anios or not fila_ajustado or len(fila_anios) < 2:
            return None

        peor_caida = 0
        anio_ruptura = None
        for i in range(1, len(fila_ajustado)):
            delta = fila_ajustado[i] - fila_ajustado[i - 1]
            if delta < peor_caida:
                peor_caida = delta
                anio_ruptura = int(round(fila_anios[i]))

        if anio_ruptura is None or peor_caida > -80:
            return None
        return anio_ruptura
    except Exception as e:
        log(f"LandTrendr no disponible para esta parcela: {str(e)[:150]}", "WARN")
        return None


# ============== GEE BACKEND ==============
def medir_B_T_D_en_GEE_PODEROSO(ee, feature_collection):
    hansen = ee.Image(HANSEN_DATASET)
    treecover = hansen.select('treecover2000')
    loss = hansen.select('loss')
    lossyear = hansen.select('lossyear')

    loss_post2020 = loss.updateMask(lossyear.gt(CORTE_LOSSYEAR))
    loss_pre2020 = loss.updateMask(lossyear.lte(CORTE_LOSSYEAR))

    loss_post2020_connected = loss_post2020.selfMask().connectedPixelCount(maxSize=64, eightConnected=True)
    loss_validada = loss_post2020.updateMask(loss_post2020_connected.gte(MIN_CLUSTER_HANSEN))

    jrc_forest = ee.Image("JRC/GFC2020/V3").select('Map').eq(1)
    jrc_water = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select('occurrence')

    srtm = ee.Image("USGS/SRTMGL1_003")
    slope = ee.Terrain.slope(srtm)
    elev = srtm.select('elevation')

    esri = ee.ImageCollection("projects/sat-io/open-datasets/landcover/ESRI_Global-LULC_10m_TS") \
        .filterDate('2023-01-01', '2023-12-31').mosaic()

    s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterDate('2023-01-01', '2024-12-31') \
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    s2_ndvi_col = s2.map(lambda img: img.normalizedDifference(['B8', 'B4']).rename('NDVI'))
    ndvi_mean = s2_ndvi_col.mean()
    ndvi_std = s2_ndvi_col.reduce(ee.Reducer.stdDev())
    ndvi_2020 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterDate('2020-01-01', '2020-12-31') \
        .map(lambda img: img.normalizedDifference(['B8', 'B4'])).mean()
    ndvi_2024 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterDate('2024-01-01', '2024-12-31') \
        .map(lambda img: img.normalizedDifference(['B8', 'B4'])).mean()
    nbdi = s2.map(lambda img: img.normalizedDifference(['B11', 'B8']).rename('NBDI')).mean()
    textura = s2.select('B8').mean().toInt32().glcmTexture(size=4).select('B8_contrast')

    # FIX 5: máscara de "píxel Hansen válido" para poder CONTARLOS por
    # polígono -- no solo promediar. treecover2000 siempre tiene valor
    # (incluso 0), así que se usa como indicador de "aquí hay dato Hansen".
    pixel_valido = treecover.gte(0)

    def medir(feat):
        geom = feat.geometry()
        stats = ee.Image.cat([
            treecover.rename('B_treecover'),
            ndvi_mean.rename('B_ndvi'),
            ndvi_std.rename('B_ndvi_std'),
            textura.rename('B_textura'),
            elev.rename('B_altitud'),
            lossyear.rename('T_perdida'),
            loss_validada.multiply(ee.Image.pixelArea()).divide(10000).rename('D_ha_validada'),
            loss_post2020.multiply(ee.Image.pixelArea()).divide(10000)
                .subtract(loss_validada.multiply(ee.Image.pixelArea()).divide(10000))
                .rename('D_ha_descartada'),
            loss_pre2020.multiply(ee.Image.pixelArea()).divide(10000).rename('D_ha_historica_pre2020'),
            jrc_forest.multiply(100).rename('B_jrc_forest_pct'),
            jrc_water.rename('D_jrc_water'),
            slope.rename('D_pend'),
            ndvi_2020.rename('ndvi_2020'),
            ndvi_2024.rename('ndvi_2024'),
        ]).reduceRegion(reducer=ee.Reducer.mean(), geometry=geom, scale=30, maxPixels=1e12)

        # FIX 5: conteo real de píxeles Hansen dentro del polígono.
        n_pix = pixel_valido.reduceRegion(
            reducer=ee.Reducer.count(), geometry=geom, scale=30, maxPixels=1e12
        ).get('treecover2000')

        return feat.set({
            'B_treecover': stats.get('B_treecover'),
            'B_ndvi': stats.get('B_ndvi'),
            'B_ndvi_std': stats.get('B_ndvi_std'),
            'B_textura': stats.get('B_textura'),
            'B_altitud': stats.get('B_altitud'),
            'T_perdida': stats.get('T_perdida'),
            'D_ha_validada': stats.get('D_ha_validada'),
            'D_ha_descartada': stats.get('D_ha_descartada'),
            'D_ha_historica_pre2020': stats.get('D_ha_historica_pre2020'),
            'B_jrc_forest_pct': stats.get('B_jrc_forest_pct'),
            'D_jrc_water': stats.get('D_jrc_water'),
            'D_pend': stats.get('D_pend'),
            'ndvi_2020': stats.get('ndvi_2020'),
            'ndvi_2024': stats.get('ndvi_2024'),
            'n_pixeles_hansen': n_pix,
        })
    return feature_collection.map(medir)


def _clasificar_confianza_hansen(n_pixeles):
    if n_pixeles is None:
        return "no_evaluada"
    if n_pixeles >= MIN_PIXELES_CONFIANZA_ALTA:
        return "alta"
    if n_pixeles >= MIN_PIXELES_CONFIANZA_MEDIA:
        return "media"
    return "baja"


def _medir_nicfi_real(ee, geom_geojson, con_nicfi, confianza_hansen):
    """FIX 5: ahora se activa cuando confianza_hansen != 'alta', no solo por área fija."""
    if not con_nicfi or confianza_hansen == "alta":
        return EvidenciaNICFI(disponible=False)
    try:
        geom = ee.Geometry(geom_geojson, None, False)
        col = ee.ImageCollection('projects/planet-nicfi/assets/basemaps/americas').filterBounds(geom)
        n_total = col.limit(1).size().getInfo()
        if n_total == 0:
            return EvidenciaNICFI(disponible=False)

        img_reciente = col.sort('system:time_start', False).first()
        img_antiguo = col.filterDate('2020-01-01', '2021-06-30').sort('system:time_start', True).first()
        if img_antiguo is None:
            return EvidenciaNICFI(disponible=False)

        bandas = img_reciente.bandNames().getInfo()
        if 'N' not in bandas:
            return EvidenciaNICFI(disponible=False)

        ndvi_reciente = img_reciente.normalizedDifference(['N', 'R']).reduceRegion(
            reducer=ee.Reducer.mean(), geometry=geom, scale=5, maxPixels=1e10, bestEffort=True
        ).getInfo().get('nd')
        ndvi_antiguo = img_antiguo.normalizedDifference(['N', 'R']).reduceRegion(
            reducer=ee.Reducer.mean(), geometry=geom, scale=5, maxPixels=1e10, bestEffort=True
        ).getInfo().get('nd')

        if ndvi_reciente is None or ndvi_antiguo is None:
            return EvidenciaNICFI(disponible=False)

        caida = ((ndvi_antiguo - ndvi_reciente) / ndvi_antiguo * 100) if ndvi_antiguo > 0.05 else None
        return EvidenciaNICFI(disponible=True, ndvi_antiguo=ndvi_antiguo, ndvi_reciente=ndvi_reciente, caida_pct=caida)
    except Exception as e:
        log(f"NICFI no disponible para esta parcela: {str(e)[:150]}", "WARN")
        return EvidenciaNICFI(disponible=False)


def _medir_estado_gee_real(gdf_estado, estado_nombre, ee, con_nicfi=False, tam_lote=200):
    predios = []
    filas = list(gdf_estado.iterrows())
    total = len(filas)

    for inicio in range(0, total, tam_lote):
        lote = filas[inicio:inicio + tam_lote]
        log(f"[{estado_nombre}] Midiendo lote GEE {inicio}-{inicio+len(lote)} de {total}...")

        features = []
        for idx, row in lote:
            geom_geo = row.geometry.__geo_interface__ if hasattr(row.geometry, '__geo_interface__') else None
            if geom_geo is None:
                continue
            id_predio = str(row.get('id_predio') or row.get('ID') or f"{estado_nombre[:3]}-{idx:07d}")
            features.append(ee.Feature(ee.Geometry(geom_geo, None, False), {
                'id_predio': id_predio,
                'municipio': str(row.get('municipio', '')),
                'geom_original': json.dumps(geom_geo),
            }))

        if not features:
            continue

        fc = ee.FeatureCollection(features)
        fc_medida = medir_B_T_D_en_GEE_PODEROSO(ee, fc)

        try:
            info = fc_medida.getInfo()
        except Exception as e:
            log(f"[{estado_nombre}] ERROR midiendo lote en GEE: {str(e)[:200]}", "ERROR")
            continue

        for feat in info.get('features', []):
            props = feat.get('properties', {})
            geom_original = json.loads(props.get('geom_original', '{}'))

            try:
                from shapely.geometry import shape
                import pyproj
                from shapely.ops import transform as shp_transform
                geom_shapely = shape(geom_original)
                proyector = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:6372", always_xy=True).transform
                superficie_ha = shp_transform(proyector, geom_shapely).area / 10000
            except Exception:
                superficie_ha = 0.0

            B_altitud_val = props.get('B_altitud') or 0.0
            B_treecover_val = props.get('B_treecover') or 0.0
            n_pix = props.get('n_pixeles_hansen')
            confianza = _clasificar_confianza_hansen(n_pix)

            bfast_anio = None
            if ALTITUD_AGUACATE[0] <= B_altitud_val <= ALTITUD_AGUACATE[1] and B_treecover_val > 10:
                try:
                    geom_ee = ee.Geometry(geom_original, None, False)
                    bfast_anio = calcular_landtrendr_breakyear(ee, geom_ee)
                except Exception as e:
                    log(f"LandTrendr falló para un predio: {str(e)[:150]}", "WARN")
                    bfast_anio = None

            p = PredioNacional(
                id_predio=props.get('id_predio', 'SIN_ID'),
                estado=estado_nombre,
                municipio=props.get('municipio', ''),
                geometria=geom_original,
                superficie_ha=superficie_ha,
                B_treecover=B_treecover_val,
                B_ndvi=props.get('B_ndvi') or 0.0,
                B_ndvi_std=props.get('B_ndvi_std') or 0.0,
                B_textura=props.get('B_textura') or 0.0,
                B_altitud=B_altitud_val,
                B_bfast_anio=bfast_anio,
                JRC_pct=props.get('B_jrc_forest_pct') or 0.0,
                T_perdida=int(props.get('T_perdida') or 0),
                D_ha_validada=props.get('D_ha_validada') or 0.0,
                D_ha_descartada_aislada=props.get('D_ha_descartada') or 0.0,
                D_ha_historica_pre2020=props.get('D_ha_historica_pre2020') or 0.0,
                D_jrc_water=props.get('D_jrc_water') or 0.0,
                D_pendiente=props.get('D_pend') or 0.0,
                ndvi_2020=props.get('ndvi_2020') or 0.0,
                ndvi_2024=props.get('ndvi_2024') or 0.0,
                n_pixeles_hansen=int(n_pix) if n_pix is not None else 0,
                confianza_hansen=confianza,
                nicfi=_medir_nicfi_real(ee, geom_original, con_nicfi, confianza),
            )
            predios.append(p)

    return predios


# ============== INGESTA NACIONAL ==============
def procesar_estado_poderoso(gdf_estado, estado_nombre, modo="prueba", con_nicfi=False):
    predios = []

    if modo == "gee":
        try:
            import ee
        except ImportError:
            raise RuntimeError("modo='gee' requiere earthengine-api instalado y ee.Initialize() ya corrido.")
        predios = _medir_estado_gee_real(gdf_estado, estado_nombre, ee, con_nicfi=con_nicfi)
    else:
        for idx, row in gdf_estado.iterrows():
            r = random.random()
            if r < 0.15:
                B, ndvi, std, alt, textura, bfast, JRC, ESRI = random.uniform(30, 65), random.uniform(0.71, 0.82), random.uniform(0.03, 0.07), random.uniform(900, 1400), random.uniform(85, 150), None, random.uniform(20, 60), random.uniform(10, 30)
                sup = random.uniform(0.2, 1.5)
                nicfi_disp = con_nicfi
                nicfi_ant = random.uniform(0.68, 0.78) if nicfi_disp else None
                nicfi_rec = random.uniform(0.66, 0.76) if nicfi_disp else None
                n_pix_sim = int(sup * 11)
            elif r < 0.28:
                B, ndvi, std, alt, textura, bfast, JRC, ESRI = random.uniform(25, 60), random.uniform(0.76, 0.88), random.uniform(0.05, 0.09), random.uniform(1600, 2200), random.uniform(20, 45), random.randint(2015, 2023), random.uniform(15, 50), random.uniform(20, 60)
                sup = random.uniform(1.0, 5.0)
                nicfi_disp = False; nicfi_ant = None; nicfi_rec = None
                n_pix_sim = int(sup * 11)
            elif r < 0.45:
                B, ndvi, std, alt, textura, bfast, JRC, ESRI = random.uniform(40, 75), random.uniform(0.70, 0.85), random.uniform(0.04, 0.08), random.uniform(100, 600), random.uniform(20, 40), None, random.uniform(25, 70), random.uniform(10, 25)
                sup = random.uniform(2.0, 10.0)
                nicfi_disp = False; nicfi_ant = None; nicfi_rec = None
                n_pix_sim = int(sup * 11)
            elif r < 0.65:
                B, ndvi, std, alt, textura, bfast, JRC, ESRI = random.uniform(20, 40), random.uniform(0.56, 0.70), random.uniform(0.08, 0.14), random.uniform(400, 1200), random.uniform(85, 180), None, random.uniform(20, 55), random.uniform(15, 40)
                sup = random.uniform(0.5, 3.0)
                nicfi_disp = con_nicfi and random.random() < 0.5
                nicfi_ant = random.uniform(0.60, 0.72) if nicfi_disp else None
                nicfi_rec = random.uniform(0.58, 0.70) if nicfi_disp else None
                n_pix_sim = int(sup * 11)
            else:
                B, ndvi, std, alt, textura, bfast, JRC, ESRI = random.uniform(5, 15), random.uniform(0.35, 0.55), random.uniform(0.12, 0.22), random.uniform(200, 800), random.uniform(30, 70), None, random.uniform(5, 20), random.uniform(60, 90)
                sup = random.uniform(0.3, 2.0)
                nicfi_disp = False; nicfi_ant = None; nicfi_rec = None
                n_pix_sim = int(sup * 11)

            D_valid = random.uniform(0.06, 2.0) if random.random() < 0.12 else 0.0
            D_desc = random.uniform(0.02, 0.15) if random.random() < 0.20 else 0.0
            D_hist = random.uniform(0.0, 1.5) if random.random() < 0.25 else 0.0
            T = random.randint(0, 24)
            jrc_water, pend = 0, random.uniform(2, 28)
            caida = ((nicfi_ant - nicfi_rec) / nicfi_ant * 100) if nicfi_ant and nicfi_rec and nicfi_ant > 0.05 else None
            confianza_sim = _clasificar_confianza_hansen(n_pix_sim)

            p = PredioNacional(
                id_predio=str(row.get('id_predio') or row.get('ID') or f"{estado_nombre[:3]}-{idx:07d}"),
                estado=estado_nombre, municipio=str(row.get('municipio', '')),
                geometria=row.geometry.__geo_interface__ if hasattr(row.geometry, '__geo_interface__') else {},
                superficie_ha=sup, B_treecover=B, B_ndvi=ndvi, B_ndvi_std=std, B_textura=textura, B_altitud=alt, B_bfast_anio=bfast,
                JRC_pct=JRC, ESRI_crops_pct=ESRI, T_perdida=T, D_ha_validada=D_valid, D_ha_descartada_aislada=D_desc,
                D_ha_historica_pre2020=D_hist, D_jrc_water=jrc_water, D_pendiente=pend,
                ndvi_2020=ndvi - 0.05, ndvi_2024=ndvi, nbdi_2020=0.1, nbdi_2024=0.12,
                n_pixeles_hansen=n_pix_sim, confianza_hansen=confianza_sim,
                nicfi=EvidenciaNICFI(disponible=nicfi_disp, ndvi_antiguo=nicfi_ant, ndvi_reciente=nicfi_rec, caida_pct=caida)
            )
            predios.append(p)

    resultados = motor_nacional_masivo(predios)

    if not resultados:
        print(f"[{estado_nombre}] 0 predios procesados (revisa la fuente de datos).")
        return 0

    conn = sqlite3.connect(DB_PATH)
    data = []
    for r in resultados:
        data.append((r.id_predio, r.estado, r.municipio, r.superficie_ha,
                      r.B_treecover, r.B_ndvi, r.B_ndvi_std, r.B_textura, r.B_altitud, r.B_bfast_anio,
                      r.JRC_pct, r.ESRI_crops_pct, r.T_perdida,
                      r.D_ha_validada, r.D_ha_descartada_aislada, r.D_ha_historica_pre2020, r.D_ha_2025, r.D_jrc_water, r.D_pendiente,
                      r.ndvi_2020, r.ndvi_2024, r.nbdi_2020, r.nbdi_2024,
                      r.n_pixeles_hansen, r.confianza_hansen,
                      int(r.nicfi.disponible), r.nicfi.ndvi_antiguo or 0, r.nicfi.ndvi_reciente or 0, r.nicfi.caida_pct or 0,
                      r.curp, r.rfc, r.net_mass_kg, r.legal_docs, r.production_date,
                      int(r.E_elegible), r.color, r.dictamen_corto, r.dictamen_largo, r.texto_evidencia,
                      r.sistema, int(r.requiere_revision), r.status_dds, r.score_prioridad,
                      r.hash_geo, r.checksum_integridad, r.fecha, r.version))
    conn.executemany(f"INSERT OR REPLACE INTO nacional_poderoso VALUES ({','.join(['?']*47)})", data)
    conn.commit()
    conn.close()

    verdes = sum(1 for r in resultados if r.color == "Verde")
    amarillos = sum(1 for r in resultados if r.color == "Amarillo")
    rojos = sum(1 for r in resultados if r.color == "Rojo")
    print(f"[{estado_nombre}] {len(resultados)} predios | Verde {verdes} | Amarillo {amarillos} | Rojo {rojos} | "
          f"Sistemas: {', '.join(set(r.sistema for r in resultados))} | "
          f"NICFI rescate: {sum(1 for r in resultados if r.nicfi.disponible and r.color=='Verde')} | "
          f"Confianza Hansen: {', '.join(set(r.confianza_hansen for r in resultados))}")
    return len(resultados)


# ============== PRODUCTOS ICM ==============
def generar_productos_icm():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    print("\n=== GENERANDO PRODUCTOS ICM-SADER ===")
    cur.execute("SELECT COUNT(*), SUM(superficie_ha) FROM nacional_poderoso WHERE E_elegible=1")
    row = cur.fetchone()
    area_elegible = row[1] if row[1] is not None else 0.0
    print(f"Producto 2 - Áreas elegibles: {row[0]} predios, {area_elegible:.2f} ha")
    print("Producto 3 - Sistemas productivos:")
    for r in cur.execute("SELECT sistema, COUNT(*), AVG(B_treecover), AVG(B_ndvi) FROM nacional_poderoso GROUP BY sistema"):
        tc = r[2] if r[2] is not None else 0.0
        nd = r[3] if r[3] is not None else 0.0
        print(f"  {r[0]}: {r[1]} predios | treecover {tc:.1f}% NDVI {nd:.2f}")
    cur.execute("SELECT COUNT(DISTINCT estado), COUNT(*) FROM nacional_poderoso")
    print(f"Producto 4 - BD predios: {cur.fetchone()}")
    print("Producto 5 - Prioritarios mitigación:")
    for r in cur.execute("SELECT color, COUNT(*), AVG(score) FROM nacional_poderoso GROUP BY color"):
        sc = r[2] if r[2] is not None else 0.0
        print(f"  {r[0]}: {r[1]} predios score {sc:.3f}")
    print("Confianza Hansen (nuevo en v9.2):")
    for r in cur.execute("SELECT confianza_hansen, COUNT(*) FROM nacional_poderoso GROUP BY confianza_hansen"):
        print(f"  {r[0]}: {r[1]} predios")

    export_path = os.path.join(EXPORT_DIR, "paquete_auditoria_nacional_EUDR.json")
    data = []
    for r in cur.execute("SELECT id_predio, superficie_ha, E_elegible, color, sistema, hash_geo, checksum_integridad, "
                          "fecha, B_treecover, D_ha_validada, D_ha_historica_pre2020, JRC_pct, confianza_hansen "
                          "FROM nacional_poderoso LIMIT 100"):
        data.append({
            "plotId": r[0], "area_ha": r[1], "deforestationFree": bool(r[2]), "color": r[3], "sistema": r[4],
            "auditProof": {
                "hash_geo": r[5], "checksum_integridad": r[6],
                "tipo_checksum": "sha256_no_es_firma_digital",
                "fecha": r[7], "B": r[8], "D_post_corte": r[9], "D_historica_pre2020": r[10],
                "JRC": r[11], "confianza_hansen": r[12]
            }
        })
    with open(export_path, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"Producto 6 - Paquete auditoría exportado: {export_path} ({len(data)} muestras)")
    conn.close()


# ============== MAIN ==============
def main():
    ap = argparse.ArgumentParser(description="Motor Nacional ICM-SADER v9.2")
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--pruebita", type=int)
    ap.add_argument("--ran", type=str)
    ap.add_argument("--modo", type=str, choices=["prueba", "gee"], default="prueba")
    ap.add_argument("--proyecto", type=str)
    ap.add_argument("--con-nicfi", action="store_true")
    ap.add_argument("--con-ml", type=str)
    ap.add_argument("--entrenar-ml", type=str)
    ap.add_argument("--max-estados", type=int, default=0)
    ap.add_argument("--resolver-traslapes", nargs="?", const="__auto__", default=None,
                     help="Activa eliminación jerárquica de traslapes. Opcionalmente pasa el nombre "
                          "del campo de prioridad legal real (ej. fecha_titulo). Sin argumento, usa "
                          "área descendente como respaldo arbitrario (con advertencia).")
    ap.add_argument("--reset", action="store_true")
    args = ap.parse_args()

    if args.reset:
        for f in [DB_PATH, CHECKPOINT_PATH]:
            if os.path.exists(f):
                os.remove(f)

    init_db_nacional()
    init_ml_global(modelo_path=args.con_ml if args.con_ml else None)

    if args.modo == "gee":
        import ee
        try:
            if args.proyecto:
                ee.Initialize(project=args.proyecto)
                log(f"Earth Engine inicializado con proyecto: {args.proyecto}", "OK")
            else:
                ee.Initialize()
                log("Earth Engine inicializado (sin --proyecto explícito).", "OK")
        except Exception as e:
            log(f"No se pudo inicializar Earth Engine: {e}", "ERROR")
            sys.exit(1)

    if args.entrenar_ml:
        ADAPTADOR_ML_GLOBAL.entrenar_desde_campo(args.entrenar_ml, args.con_ml or "./modelo_sistema_poderoso.pkl")
        return

    if args.demo:
        print("=== MOTOR NACIONAL v9.2 - DEMO (corte EUDR + confianza Hansen por píxeles) ===\n")
        demos = [
            PredioNacional("VER-CAFE-0.25ha-SIN-NICFI", estado="Veracruz", municipio="Teocelo", superficie_ha=0.25,
                            B_treecover=40, B_ndvi=0.75, B_ndvi_std=0.05, B_textura=100, B_altitud=1100,
                            JRC_pct=35, ESRI_crops_pct=20, T_perdida=5, D_ha_validada=0.0, D_ha_descartada_aislada=0.09,
                            n_pixeles_hansen=3, confianza_hansen="baja", nicfi=EvidenciaNICFI(disponible=False)),
            PredioNacional("VER-CAFE-0.8ha-HUECO-ESCALA-SIN-NICFI", estado="Veracruz", municipio="Teocelo", superficie_ha=0.8,
                            B_treecover=40, B_ndvi=0.75, B_ndvi_std=0.05, B_textura=100, B_altitud=1100,
                            JRC_pct=12, ESRI_crops_pct=20, T_perdida=5, D_ha_validada=0.0,
                            n_pixeles_hansen=18, confianza_hansen="media", nicfi=EvidenciaNICFI(disponible=False)),
            PredioNacional("VER-CAFE-0.8ha-HUECO-ESCALA-CON-NICFI", estado="Veracruz", municipio="Teocelo", superficie_ha=0.8,
                            B_treecover=40, B_ndvi=0.75, B_ndvi_std=0.05, B_textura=100, B_altitud=1100,
                            JRC_pct=12, ESRI_crops_pct=20, T_perdida=5, D_ha_validada=0.0,
                            n_pixeles_hansen=18, confianza_hansen="media",
                            nicfi=EvidenciaNICFI(disponible=True, ndvi_antiguo=0.71, ndvi_reciente=0.69, caida_pct=2.8)),
            PredioNacional("MIC-AGUACATE-2ha-ROJO", estado="Michoacán", municipio="Uruapan", superficie_ha=2.0,
                            B_treecover=35, B_ndvi=0.80, B_ndvi_std=0.06, B_textura=30, B_altitud=1800, B_bfast_anio=2021,
                            JRC_pct=60, ESRI_crops_pct=80, T_perdida=22, D_ha_validada=1.2, D_ha_descartada_aislada=0.09,
                            n_pixeles_hansen=220, confianza_hansen="alta"),
            PredioNacional("MIC-AGUACATE-3ha-DEFOR-PRE2020-YA-NO-ROJO", estado="Michoacán", municipio="Uruapan", superficie_ha=3.0,
                            B_treecover=45, B_ndvi=0.82, B_ndvi_std=0.05, B_textura=35, B_altitud=1900, B_bfast_anio=2016,
                            JRC_pct=55, ESRI_crops_pct=70, T_perdida=8, D_ha_validada=0.0, D_ha_historica_pre2020=1.8,
                            n_pixeles_hansen=330, confianza_hansen="alta"),
        ]
        res = motor_nacional_masivo(demos)
        for r in res:
            print(f"{r.id_predio} | {r.superficie_ha}ha | {r.estado} | B={r.B_treecover:.0f}% JRC={r.JRC_pct:.0f}% "
                  f"confianza_hansen={r.confianza_hansen} ({r.n_pixeles_hansen}px) NICFI={r.nicfi.disponible}")
            print(f"  => {r.color} | {r.dictamen_corto} | Sistema:{r.sistema} | Checksum:{r.checksum_integridad}")
            print(f"  Evidencia: {r.texto_evidencia[:220]}...\n")

    if args.pruebita:
        print(f"=== PRUEBITA NACIONAL {args.pruebita} PREDIOS modo={args.modo} NICFI={args.con_nicfi} ===")
        import geopandas as gpd
        from shapely.geometry import Point
        rows = []
        for i in range(args.pruebita):
            rows.append({"id_predio": f"RAN-NAC-{i:07d}",
                         "municipio": random.choice(["Teocelo", "Uruapan", "Palenque", "Misantla"]),
                         "estado": random.choice(["Veracruz", "Michoacán", "Chiapas", "Jalisco"]),
                         "geometry": Point(-102 + random.random() * 6, 17 + random.random() * 5).buffer(0.001)})
        gdf = gpd.GeoDataFrame(rows, crs="EPSG:4326")
        chk = checkpoint_load()
        total = 0
        for estado, gdf_est in gdf.groupby('estado'):
            if args.max_estados and total >= args.max_estados * 100:
                break
            n = procesar_estado_poderoso(gdf_est, estado, modo=args.modo, con_nicfi=args.con_nicfi)
            total += n
            chk["total_predios"] += n
            if estado not in chk["estados_completados"]:
                chk["estados_completados"].append(estado)
            checkpoint_save(chk)
        generar_productos_icm()

    if args.ran:
        print(f"=== MODO NACIONAL REAL {args.ran} ===")
        import geopandas as gpd
        chk = checkpoint_load()
        print(f"Checkpoint actual: {chk['estados_completados']} total {chk['total_predios']}")
        gdf_total = gpd.read_file(args.ran)
        gdf_total['geometry'] = gdf_total['geometry'].make_valid()

        if args.resolver_traslapes is not None:
            campo = None if args.resolver_traslapes == "__auto__" else args.resolver_traslapes
            gdf_total = eliminar_traslapes_jerarquico(gdf_total, campo_prioridad=campo)
        else:
            log("--resolver-traslapes no activado: si el RAN tiene predios solapados, "
                "las hectáreas se pueden estar contando dos veces.", "WARN")

        estados = gdf_total['estado'].unique().tolist() if 'estado' in gdf_total.columns else ["Nacional"]
        if args.max_estados:
            estados = estados[:args.max_estados]
        for estado in estados:
            if estado in chk["estados_completados"]:
                print(f"{estado} ya completado, saltando")
                continue
            gdf_est = gdf_total[gdf_total['estado'] == estado] if 'estado' in gdf_total.columns else gdf_total
            print(f"Procesando {estado}: {len(gdf_est)} predios")
            n = procesar_estado_poderoso(gdf_est, estado, modo=args.modo, con_nicfi=args.con_nicfi)
            chk["estados_completados"].append(estado)
            chk["total_predios"] += n
            checkpoint_save(chk)
            print(f"Checkpoint guardado: {estado} completado")
        generar_productos_icm()


if __name__ == "__main__":
    main()
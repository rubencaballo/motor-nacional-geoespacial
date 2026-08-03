#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PredioNacional: unica definicion. Depende de config (constantes) y de
modelos.evidencia (EvidenciaNICFI) -- ambos importados explicitamente, a
diferencia de la version rota anterior que no importaba nada."""
import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import List, Optional

from config import (
    CORTE_EUDR, MIN_CLUSTER_HANSEN, AREA_MIN_HANSEN_HA, CAIDA_NDVI_MAX,
    JRC_MIN_VERDE, JRC_MIN_VERDE_NICFI, UMBRAL_CAFE_TREECOVER, UMBRAL_CAFE_NDVI,
    UMBRAL_CAFE_NDVI_STD, ALTITUD_CAFE, ALTITUD_AGUACATE, UMBRAL_SILVO,
    UMBRAL_SILVO_NDVI, PESOS_PRIORIZACION, MIN_PIXELES_CONFIANZA_ALTA,
    MIN_PIXELES_CONFIANZA_MEDIA, VERSION,
)
from modelos.evidencia import EvidenciaNICFI

# ADAPTADOR_ML_GLOBAL se inyecta desde core.motor al arrancar main() para
# evitar import circular (ml.modelo tambien podria necesitar PredioNacional
# en el futuro). Si nadie lo inyecta, sigue en None y el motor usa reglas.
ADAPTADOR_ML_GLOBAL = None


def set_adaptador_ml(adaptador):
    global ADAPTADOR_ML_GLOBAL
    ADAPTADOR_ML_GLOBAL = adaptador


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
    hash_resultado: str = ""
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
        payload_contenido = f"{self.id_predio}{self.hash_geo}{self.B_treecover:.2f}{self.D_ha_validada:.4f}{self.color}"
        self.hash_resultado = hashlib.sha256(payload_contenido.encode()).hexdigest()[:16]
        self.fecha = datetime.now(timezone.utc).isoformat()



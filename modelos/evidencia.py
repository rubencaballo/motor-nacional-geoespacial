from dataclasses import dataclass
from typing import Optional, Dict

@dataclass
class EvidenciaNICFI:
    ndvi_antiguo: float = 0.0
    ndvi_reciente: float = 0.0
    caida_pct: float = 0.0
    anio_perdida: int = 0
    hansen_loss_year: int = 0
    area_ha: float = 0.0
    ndvi_serie: Optional[Dict[int, float]] = None
    nbdi_serie: Optional[Dict[int, float]] = None
    jrc_cobertura_pct: float = 0.0
    jrc_bosque_2020: Optional[bool] = None
    esri_crops_pct: float = 0.0
    esri_trees_pct: float = 0.0
    esri_landcover_2020: Optional[int] = None
    es_urbano_esri: bool = False
    hansen_perdida_ha: float = 0.0
    nbdi_2020: Optional[float] = None
    planet_ndvi_2020: Optional[float] = None
    planet_ndvi_reciente: Optional[float] = None
    planet_caida_pct: float = 0.0
    planet_textura: str = "heterogenea"
    planet_disponible: bool = False

    def __post_init__(self):
        if self.hansen_loss_year and not self.anio_perdida:
            self.anio_perdida = self.hansen_loss_year
        if self.anio_perdida and not self.hansen_loss_year:
            self.hansen_loss_year = self.anio_perdida

    def _ndvi_estable(self):
        if self.ndvi_serie and len(self.ndvi_serie)>=3:
            vals=list(self.ndvi_serie.values())
            return (max(vals)-min(vals))<0.15 and min(vals)>0.4
        return self.caida_pct < 0.2

    def _nbdi_negativo(self):
        if self.nbdi_serie and len(self.nbdi_serie)>=2:
            return all(v<0 for v in self.nbdi_serie.values())
        if self.nbdi_2020 is not None:
            return self.nbdi_2020 < 0
        return False

    def _contradiccion_jrc_esri(self):
        return (self.jrc_cobertura_pct < 15 and self.esri_trees_pct > 80) or (self.jrc_cobertura_pct >= 15 and self.esri_trees_pct == 100 and self.esri_crops_pct == 0)

    def es_exportable_eudr(self):
        B = (self.jrc_cobertura_pct >= 15) if self.jrc_cobertura_pct > 0 else (self.ndvi_antiguo > 0.6)
        D_hansen = self.hansen_perdida_ha > 0.01 or self.anio_perdida > 2020
        D_nicfi = self.caida_pct > 0.2
        return B and not (D_hansen or D_nicfi)

    def dictamen_completo(self):
        B_jrc = self.jrc_cobertura_pct >= 15
        B_esri = self.esri_trees_pct > 50 or self.esri_crops_pct == 0
        D_hansen = self.hansen_perdida_ha > 0.01 or self.anio_perdida > 2020
        ndvi_estable = self._ndvi_estable()
        nbdi_neg = self._nbdi_negativo()
        contradiccion = self._contradiccion_jrc_esri()

        if self.area_ha > 0 and self.area_ha < 0.5:
            if contradiccion and not D_hansen and ndvi_estable and nbdi_neg:
                if self.planet_disponible and self.planet_ndvi_2020 and self.planet_ndvi_2020>0.6:
                    return "VERDE - APTO_EXPORTACION - Planet 4.7m confirma agroforestal"
                return "AMARILLO - REQUIERE PLANET 4.7m"

        if B_jrc and B_esri and not D_hansen and ndvi_estable and nbdi_neg:
            return "VERDE - APTO_EXPORTACION - CONFORME EUDR"
        if contradiccion and not D_hansen and ndvi_estable:
            return "VERDE - APTO_EXPORTACION - FAQ v5 Sec 2.3"
        if D_hansen:
            return "ROJO - NO APTO - Deforestacion Hansen >2020"
        return "AMARILLO - REQUIERE_EVIDENCIA_DE_CAMPO"

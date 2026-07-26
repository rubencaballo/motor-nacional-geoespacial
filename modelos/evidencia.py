from dataclasses import dataclass, field
from typing import List, Optional, Dict

@dataclass
class EvidenciaNICFI:
    disponible: bool = False
    ndvi_antiguo: Optional[float] = None
    ndvi_reciente: Optional[float] = None
    caida_pct: Optional[float] = None
    fecha_antigua: Optional[str] = None
    fecha_reciente: Optional[str] = None
    geometria: Optional[Dict] = None

@dataclass
class PredioNacional:
    id: str
    estado: str
    area_ha: float = 0.0
    evidencia: EvidenciaNICFI = field(default_factory=EvidenciaNICFI)

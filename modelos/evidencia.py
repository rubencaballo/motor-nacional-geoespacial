class EvidenciaNICFI:
    disponible: bool = False
    ndvi_antiguo: Optional[float] = None
    ndvi_reciente: Optional[float] = None
    caida_pct: Optional[float] = None
    fecha_antigua: Optional[str] = None
    fecha_reciente: Optional[str] = None


@dataclass

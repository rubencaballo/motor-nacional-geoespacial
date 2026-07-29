#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EvidenciaNICFI: unica definicion en todo el proyecto (antes habia dos
versiones incompatibles con campos distintos -- ese era el bug real)."""
from dataclasses import dataclass
from typing import Optional


@dataclass
class EvidenciaNICFI:
    disponible: bool = False
    ndvi_antiguo: Optional[float] = None
    ndvi_reciente: Optional[float] = None
    caida_pct: Optional[float] = None
    fecha_antigua: Optional[str] = None
    fecha_reciente: Optional[str] = None

"""
salamandra_biodiversidad.py
============================
Modulo nucleo de biodiversidad + carbono para la plataforma Salamandra.
Aplicable a cualquier poligono: ANP, predio privado, finca cafetalera, etc.

Este archivo es el resultado de auditar ~60 scripts historicos del proyecto
Ramsar 1601 Texolo (versiones V19 a V28). Incluye UNICAMENTE las tecnicas que
sobrevivieron esa revision de integridad cientifica. Cada funcion documenta
que bug real corrige, para que quede como referencia al extender el modulo.

INCLUIDO (validado):
  - Area SIEMPRE calculada desde la geometria real del poligono, nunca
    hardcodeada. (Bug encontrado: un poligono real de ~503 ha circulaba
    etiquetado como "29 ha" en nombres de archivo y se uso ese numero falso
    en todos los calculos posteriores durante varias versiones.)
  - Validacion de clase taxonomica contra familia, con tabla de referencia.
    (Bug mas grave encontrado: columnas "CLASE"/"CLASE_LIMPIA" con aves
    reales etiquetadas como Mammalia/Amphibia/Reptilia. Ningun intento de
    "correccion" en 4 versiones sucesivas lo resolvio de raiz porque
    nunca se cruzo contra una tabla de autoridad familia->clase.)
  - Rarefaccion real por remuestreo Monte Carlo, para comparar riqueza de
    especies entre zonas con distinto esfuerzo de muestreo.
  - Correccion de area real 3D via DEM (pendiente del terreno).
  - Modelo especie-area de Wilson-MacArthur (S = c * A^z), con z de
    literatura por default; z empirico solo si hay evidencia suficiente.
  - Comparacion honesta de modelos: validacion cruzada K-fold + AIC,
    no solo ajuste en la muestra de entrenamiento.
  - Estimacion de carbono con doble fuente independiente (biomasa satelital
    tipo ESA CCI + lidar tipo GEDI L4A) reportando SIEMPRE su incertidumbre.

DELIBERADAMENTE NO INCLUIDO (por hallazgos de la auditoria):
  - Formulas de "presion humana" ajustadas hasta dar el resultado esperado.
  - Datos de comparacion nacional/regional simulados (np.random como si
    fueran datos reales).
  - Exclusion de registros "inconvenientes" para mejorar una metrica.
  - Ajuste de z con menos de 5 sitios independientes (con 3 puntos el
    exponente sale inestable e incluso negativo, ecologicamente absurdo).
  - Cualquier "correccion" de clase que no se verifique familia por familia.

Dependencias: numpy, pandas (requeridas). statsmodels (para AIC/OLS,
opcional pero recomendada). rasterio (opcional, solo para area 3D real).
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
import pandas as pd

try:
    import statsmodels.api as sm
    _HAS_STATSMODELS = True
except ImportError:  # pragma: no cover
    _HAS_STATSMODELS = False

try:
    import rasterio
    _HAS_RASTERIO = True
except ImportError:  # pragma: no cover
    _HAS_RASTERIO = False


# ---------------------------------------------------------------------------
# 1. AREA DESDE GEOMETRIA REAL (nunca hardcodear hectareas)
# ---------------------------------------------------------------------------

def compute_polygon_area_ha(coordinates: Sequence[Sequence[float]]) -> float:
    """
    Calcula el area real (en hectareas) de un poligono a partir de sus
    coordenadas geograficas (lon, lat), usando una proyeccion local
    equirrectangular centrada en la latitud media del poligono y la
    formula del shoelace.

    Parametros
    ----------
    coordinates : lista de (lon, lat) del anillo exterior del poligono
        (formato GeoJSON: coordinates[0] de un Polygon, o cada anillo de
        un MultiPolygon).

    Por que existe esta funcion
    ----------------------------
    En la auditoria del proyecto Texolo, un poligono geojson/kml cuyo
    NOMBRE decia "NUCLEO_29HA" resulto medir realmente 504.24 ha al
    calcularlo desde sus propias coordenadas. Ese "29" nunca fue el area
    real de ningun poligono: fue un numero mal tecleado en algun momento
    que se copio como SUP_HA=29.0 en cada resumen posterior, durante mas
    de 5 versiones del proyecto, sin que nadie lo recalculara del archivo
    geografico real. La regla de este modulo es: el area SIEMPRE se
    deriva de la geometria, nunca se lee de un nombre de archivo ni de
    un valor fijo en el codigo.
    """
    coords = list(coordinates)
    if len(coords) < 3:
        raise ValueError("Un poligono necesita al menos 3 vertices.")

    lats = [c[1] for c in coords]
    lat0 = sum(lats) / len(lats)
    R = 6371000.0  # radio medio de la Tierra, metros
    m_per_deg_lat = R * (math.pi / 180)
    m_per_deg_lon = R * math.cos(math.radians(lat0)) * (math.pi / 180)

    pts = [(lon * m_per_deg_lon, lat * m_per_deg_lat) for lon, lat in coords]
    area = 0.0
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        area += x1 * y2 - x2 * y1
    area = abs(area) / 2.0
    return area / 10000.0  # m2 -> ha


def area_ha_from_geojson_feature(feature: dict) -> float:
    """
    Extrae el area en hectareas de una Feature GeoJSON de tipo Polygon o
    MultiPolygon, sumando todos los sub-poligonos si aplica.
    """
    geom = feature["geometry"]
    gtype = geom["type"]
    if gtype == "Polygon":
        rings = [geom["coordinates"][0]]
    elif gtype == "MultiPolygon":
        rings = [poly[0] for poly in geom["coordinates"]]
    else:
        raise ValueError(f"Tipo de geometria no soportado: {gtype}")
    return sum(compute_polygon_area_ha(ring) for ring in rings)


# ---------------------------------------------------------------------------
# 2. VALIDACION DE CLASE TAXONOMICA CONTRA FAMILIA
# ---------------------------------------------------------------------------

# Tabla de referencia minima. Amplia esta tabla segun las familias que vayan
# apareciendo en cada sitio nuevo (ANP o finca) -- nunca confies en la
# columna CLASE de un CSV sin cruzarla contra algo asi.
FAMILY_CLASS_MAP: dict[str, str] = {
    # Aves (ejemplos vistos en el proyecto Texolo)
    "Accipitridae": "Aves", "Falconidae": "Aves", "Psittacidae": "Aves",
    "Turdidae": "Aves", "Tyrannidae": "Aves", "Icteridae": "Aves",
    "Parulidae": "Aves", "Cathartidae": "Aves", "Thraupidae": "Aves",
    "Passerellidae": "Aves", "Picidae": "Aves", "Vireonidae": "Aves",
    "Columbidae": "Aves", "Fringillidae": "Aves", "Troglodytidae": "Aves",
    "Corvidae": "Aves", "Trochilidae": "Aves", "Apodidae": "Aves",
    # Reptilia
    "Phrynosomatidae": "Reptilia", "Colubridae": "Reptilia",
    "Dactyloidae": "Reptilia", "Viperidae": "Reptilia",
    # Amphibia
    "Plethodontidae": "Amphibia", "Craugastoridae": "Amphibia",
    "Bufonidae": "Amphibia", "Hylidae": "Amphibia", "Ranidae": "Amphibia",
    # Mammalia
    "Felidae": "Mammalia", "Mustelidae": "Mammalia",
    "Phyllostomidae": "Mammalia", "Cervidae": "Mammalia",
    "Procyonidae": "Mammalia", "Didelphidae": "Mammalia",
    "Sciuridae": "Mammalia", "Canidae": "Mammalia",
    # Grupos que NO son fauna objetivo pero aparecen mezclados en
    # descargas GBIF/iNaturalist sin filtrar -- deben salir de cualquier
    # conteo de Aves/Mammalia/Amphibia/Reptilia.
    "Asteraceae": "Plantae", "Melastomataceae": "Plantae",
    "Bromeliaceae": "Plantae", "Orchidaceae": "Plantae",
    "Zingiberaceae": "Plantae", "Piperaceae": "Plantae",
    "Convolvulaceae": "Plantae", "Gesneriaceae": "Plantae",
    "Nymphalidae": "Insecta", "Apidae": "Insecta", "Pieridae": "Insecta",
    "Noctuidae": "Insecta", "Apatelodidae": "Insecta",
    "Chrysomelidae": "Insecta",
    "Salticidae": "Arachnida", "Araneidae": "Arachnida",
    "Hymenogastraceae": "Fungi",

    # -----------------------------------------------------------------
    # Ampliacion 2026-09-02: 188 familias reales encontradas al validar
    # el CSV real de biodiversidad de V19.4 (proyecto Texolo, 4954
    # registros). Clasificadas a mano una por una. Los hallazgos mas
    # importantes: solo aparecieron 2 familias nuevas de anfibios reales
    # (Ambystomatidae, Centrolenidae -- salamandras de montaña y ranas
    # de cristal, ambas de interes de conservacion real), 2 de reptiles
    # (Gekkonidae, Scincidae) y 1 de mamiferos (Geomyidae). Todo lo
    # demas en esta lista es planta, hongo, insecto, arácnido, molusco,
    # crustaceo o miriapodo -- ninguno de estos cuenta como Aves,
    # Mammalia, Amphibia o Reptilia aunque haya llegado etiquetado asi
    # en el CSV original.
    # -----------------------------------------------------------------

    # Aves (22 familias nuevas confirmadas reales)
    "Alcedinidae": "Aves", "Anatidae": "Aves", "Ardeidae": "Aves",
    "Caprimulgidae": "Aves", "Cardinalidae": "Aves", "Cinclidae": "Aves",
    "Cotingidae": "Aves", "Cracidae": "Aves", "Cuculidae": "Aves",
    "Furnariidae": "Aves", "Hirundinidae": "Aves", "Mimidae": "Aves",
    "Momotidae": "Aves", "Passeridae": "Aves", "Polioptilidae": "Aves",
    "Ptilogonatidae": "Aves", "Rallidae": "Aves", "Ramphastidae": "Aves",
    "Regulidae": "Aves", "Scolopacidae": "Aves", "Strigidae": "Aves",
    "Trogonidae": "Aves",

    # Mammalia (1 familia nueva confirmada real)
    "Geomyidae": "Mammalia",  # tuzas / pocket gophers

    # Amphibia (2 familias nuevas confirmadas reales -- hallazgo notable)
    "Ambystomatidae": "Amphibia",  # salamandras de montaña
    "Centrolenidae": "Amphibia",   # ranas de cristal

    # Reptilia (2 familias nuevas confirmadas reales)
    "Gekkonidae": "Reptilia",   # gecos
    "Scincidae": "Reptilia",    # escincos / lagartijas escamosas

    # Plantae
    "Acanthaceae": "Plantae", "Altingiaceae": "Plantae",
    "Anemiaceae": "Plantae", "Annonaceae": "Plantae",
    "Apocynaceae": "Plantae", "Araceae": "Plantae",
    "Araliaceae": "Plantae", "Arecaceae": "Plantae",
    "Aristolochiaceae": "Plantae", "Asparagaceae": "Plantae",
    "Asphodelaceae": "Plantae", "Balsaminaceae": "Plantae",
    "Bignoniaceae": "Plantae", "Cactaceae": "Plantae",
    "Calceolariaceae": "Plantae", "Campanulaceae": "Plantae",
    "Cannaceae": "Plantae", "Commelinaceae": "Plantae",
    "Cucurbitaceae": "Plantae", "Ericaceae": "Plantae",
    "Euphorbiaceae": "Plantae", "Fabaceae": "Plantae",
    "Heliconiaceae": "Plantae", "Iridaceae": "Plantae",
    "Lamiaceae": "Plantae", "Lentibulariaceae": "Plantae",
    "Lythraceae": "Plantae", "Malvaceae": "Plantae",
    "Marantaceae": "Plantae", "Moraceae": "Plantae",
    "Myrtaceae": "Plantae", "Onagraceae": "Plantae",
    "Orobanchaceae": "Plantae", "Papaveraceae": "Plantae",
    "Passifloraceae": "Plantae", "Phytolaccaceae": "Plantae",
    "Pinaceae": "Plantae", "Platanaceae": "Plantae",
    "Poaceae": "Plantae", "Polygonaceae": "Plantae",
    "Primulaceae": "Plantae", "Pteridaceae": "Plantae",
    "Rosaceae": "Plantae", "Rubiaceae": "Plantae",
    "Rutaceae": "Plantae", "Siparunaceae": "Plantae",
    "Solanaceae": "Plantae", "Urticaceae": "Plantae",
    "Violaceae": "Plantae",
    # Briofitas/hepaticas -- se agrupan como Plantae para fines de este
    # modulo (no son fauna objetivo), aunque filogeneticamente no son
    # plantas vasculares.
    "Dumortieraceae": "Plantae", "Eurhynchidae": "Plantae",

    # Fungi
    "Auriculariaceae": "Fungi", "Bolbitiaceae": "Fungi",
    "Boletaceae": "Fungi", "Cerrenaceae": "Fungi",
    "Clavicipitaceae": "Fungi", "Dacrymycetaceae": "Fungi",
    "Hygrophoraceae": "Fungi", "Hypocreaceae": "Fungi",
    "Laetiporaceae": "Fungi", "Mycenaceae": "Fungi",
    "Nectriaceae": "Fungi", "Omphalotaceae": "Fungi",
    "Polyporaceae": "Fungi", "Psathyrellaceae": "Fungi",
    "Schizophyllaceae": "Fungi", "Strophariaceae": "Fungi",

    # Insecta (ordenes varios: Coleoptera, Lepidoptera, Hymenoptera,
    # Diptera, Hemiptera, Odonata, Orthoptera, Blattodea, Mantodea,
    # Phasmatodea, Megaloptera)
    "Acrididae": "Insecta", "Alydidae": "Insecta",
    "Anostostomatidae": "Insecta", "Asilidae": "Insecta",
    "Bibionidae": "Insecta", "Bombyliidae": "Insecta",
    "Buprestidae": "Insecta", "Calliphoridae": "Insecta",
    "Calopterygidae": "Insecta", "Cantharidae": "Insecta",
    "Castniidae": "Insecta", "Cerambycidae": "Insecta",
    "Cercopidae": "Insecta", "Cicadellidae": "Insecta",
    "Cicadidae": "Insecta", "Coccinellidae": "Insecta",
    "Coenagrionidae": "Insecta", "Coreidae": "Insecta",
    "Corydalidae": "Insecta", "Crabronidae": "Insecta",
    "Crambidae": "Insecta", "Culicidae": "Insecta",
    "Curculionidae": "Insecta", "Cynipidae": "Insecta",
    "Diapheromeridae": "Insecta", "Dryophthoridae": "Insecta",
    "Ectobiidae": "Insecta", "Elateridae": "Insecta",
    "Endomychidae": "Insecta", "Erebidae": "Insecta",
    "Erotylidae": "Insecta", "Eumenidae": "Insecta",
    "Formicidae": "Insecta", "Geometridae": "Insecta",
    "Gerridae": "Insecta", "Halictidae": "Insecta",
    "Hesperiidae": "Insecta", "Hippoboscidae": "Insecta",
    "Ichneumonidae": "Insecta", "Lampyridae": "Insecta",
    "Lestidae": "Insecta", "Leucospidae": "Insecta",
    "Libellulidae": "Insecta", "Limacodidae": "Insecta",
    "Limoniidae": "Insecta", "Lycaenidae": "Insecta",
    "Lycidae": "Insecta", "Mantidae": "Insecta",
    "Megachilidae": "Insecta", "Megalopodidae": "Insecta",
    "Megalopygidae": "Insecta", "Meloidae": "Insecta",
    "Membracidae": "Insecta", "Miridae": "Insecta",
    "Mutillidae": "Insecta", "Notodontidae": "Insecta",
    "Papilionidae": "Insecta", "Passalidae": "Insecta",
    "Pentatomidae": "Insecta", "Platystictidae": "Insecta",
    "Platystomatidae": "Insecta", "Pompilidae": "Insecta",
    "Reduviidae": "Insecta", "Rhipiceridae": "Insecta",
    "Riodinidae": "Insecta", "Romaleidae": "Insecta",
    "Saturniidae": "Insecta", "Scarabaeidae": "Insecta",
    "Scoliidae": "Insecta", "Sphecidae": "Insecta",
    "Sphingidae": "Insecta", "Staphylinidae": "Insecta",
    "Stratiomyidae": "Insecta", "Syrphidae": "Insecta",
    "Tachinidae": "Insecta", "Tenebrionidae": "Insecta",
    "Tetrigidae": "Insecta", "Tettigoniidae": "Insecta",
    "Thyrididae": "Insecta", "Tipulidae": "Insecta",
    "Vespidae": "Insecta",

    # Arachnida
    "Hersiliidae": "Arachnida", "Lycosidae": "Arachnida",
    "Pisauridae": "Arachnida", "Tetragnathidae": "Arachnida",
    "Theridiidae": "Arachnida", "Trechaleidae": "Arachnida",

    # Mollusca (caracoles terrestres)
    "Bulimulidae": "Mollusca", "Helicidae": "Mollusca",

    # Crustacea (isopodos terrestres -- cochinillas de humedad)
    "Armadillidiidae": "Crustacea", "Porcellionidae": "Crustacea",

    # Myriapoda (ciempies, milpies)
    "Aphelidesmidae": "Myriapoda", "Scolopendridae": "Myriapoda",

    # Platyhelminthes (planarias terrestres)
    "Geoplanidae": "Platyhelminthes",
}


@dataclass
class ValidacionClaseReporte:
    total_registros: int
    registros_corregidos: int
    detalle_por_clase_declarada: pd.DataFrame
    familias_desconocidas: list[str]

    @property
    def porcentaje_corregido(self) -> float:
        if self.total_registros == 0:
            return 0.0
        return 100.0 * self.registros_corregidos / self.total_registros


def validate_taxonomic_class(
    df: pd.DataFrame,
    class_col: str = "CLASE",
    family_col: str = "FAMILIA",
    family_class_map: Optional[dict[str, str]] = None,
) -> tuple[pd.DataFrame, ValidacionClaseReporte]:
    """
    Cruza la columna de clase declarada contra la familia real, usando una
    tabla de autoridad (FAMILY_CLASS_MAP por default). Devuelve el
    dataframe con una columna nueva 'CLASE_VALIDADA' (la que manda) y un
    reporte de cuantos registros no coincidian.

    Por que existe esta funcion
    ----------------------------
    Es la funcion mas importante del modulo. En la auditoria del proyecto
    Texolo se encontro que, durante 4 versiones sucesivas (V24 a V27B),
    la columna CLASE/CLASE_LIMPIA tenia aves reales (Cathartes aura,
    Amazona albifrons, Streptoprocne zonaris...) etiquetadas como
    Mammalia y Amphibia. Cada intento de "correccion" (columnas con
    nombres como CLASE_CORREGIDA, CLASE_LIMPIA, "AUDITADO" en el nombre
    del archivo) heredo el mismo error sin corregirlo, porque nunca se
    verifico contra una tabla real de familia->clase. Cualquier
    breakdown de biodiversidad por clase taxonomica en Salamandra debe
    pasar por esta funcion antes de mostrarse o reportarse.
    """
    if family_class_map is None:
        family_class_map = FAMILY_CLASS_MAP

    out = df.copy()
    out["CLASE_VALIDADA"] = out[family_col].map(family_class_map)

    desconocidas = sorted(
        out.loc[out["CLASE_VALIDADA"].isna(), family_col].dropna().unique().tolist()
    )
    if desconocidas:
        warnings.warn(
            f"{len(desconocidas)} familias no estan en la tabla de referencia "
            f"y no se pudieron validar (se conserva su CLASE original sin "
            f"garantia): {desconocidas[:15]}{'...' if len(desconocidas) > 15 else ''}"
        )
    out["CLASE_VALIDADA"] = out["CLASE_VALIDADA"].fillna(out[class_col])

    mismatch = out["CLASE_VALIDADA"] != out[class_col]
    detalle = (
        out.assign(_mismatch=mismatch)
        .groupby(class_col)["_mismatch"]
        .agg(total="count", corregidos="sum")
        .assign(pct_corregido=lambda d: 100 * d["corregidos"] / d["total"])
        .reset_index()
    )

    reporte = ValidacionClaseReporte(
        total_registros=len(out),
        registros_corregidos=int(mismatch.sum()),
        detalle_por_clase_declarada=detalle,
        familias_desconocidas=desconocidas,
    )
    if reporte.porcentaje_corregido > 10:
        warnings.warn(
            f"ALERTA: {reporte.porcentaje_corregido:.1f}% de los registros "
            f"tenian una clase declarada que no coincide con su familia. "
            f"Esto fue exactamente el patron del bug historico Texolo -- "
            f"revisa la fuente de datos antes de reportar cifras por clase."
        )
    return out, reporte


# ---------------------------------------------------------------------------
# 3. RAREFACCION REAL (Monte Carlo) -- controla por esfuerzo de muestreo
# ---------------------------------------------------------------------------

def rarefy_richness(
    df: pd.DataFrame,
    level_col: str,
    species_col: str,
    n_target: int,
    n_iter: int = 1000,
    random_state: Optional[int] = None,
) -> pd.DataFrame:
    """
    Rarefaccion real por remuestreo sin reemplazo: para cada nivel/zona,
    remuestrea n_target registros al azar (n_iter veces) y cuenta especies
    unicas cada vez. Devuelve la media y desviacion estandar de la riqueza
    rarefecida -- comparable entre zonas con distinto numero de registros.

    Por que existe esta funcion
    ----------------------------
    Comparar riqueza cruda entre zonas con esfuerzo de muestreo muy
    distinto es enganoso: en el proyecto Texolo el nucleo tenia 3,312
    registros contra 537-864 en los buffers, asi que su riqueza cruda
    (568 spp) estaba inflada solo por tener mas observaciones, no
    necesariamente por ser mas biodiverso. Rarefectando a un mismo N la
    diferencia real resulto mucho mas chica (198 vs 158-165 spp) pero
    seguia siendo positiva -- ese es el numero defendible, no el crudo.
    Zonas con menos registros que n_target se excluyen (no se puede
    rarefectar hacia arriba) y se reportan aparte.
    """
    rng = np.random.default_rng(random_state)
    filas = []
    excluidos = []
    for nivel, grupo in df.groupby(level_col):
        n_disponibles = len(grupo)
        if n_disponibles < n_target:
            excluidos.append((nivel, n_disponibles))
            continue
        especies = grupo[species_col].to_numpy()
        conteos = np.empty(n_iter)
        for i in range(n_iter):
            muestra = rng.choice(especies, size=n_target, replace=False)
            conteos[i] = len(np.unique(muestra))
        filas.append(
            {
                level_col: nivel,
                "N_disponibles": n_disponibles,
                "N_target": n_target,
                "S_raw": grupo[species_col].nunique(),
                "S_rarefecida_media": conteos.mean(),
                "S_rarefecida_sd": conteos.std(ddof=1),
                "n_iter": n_iter,
            }
        )
    if excluidos:
        warnings.warn(
            f"Zonas excluidas de la rarefaccion por tener menos de "
            f"{n_target} registros (no se puede rarefectar hacia arriba): "
            f"{excluidos}"
        )
    return pd.DataFrame(filas)


# ---------------------------------------------------------------------------
# 4. AREA REAL 3D VIA DEM (correccion por pendiente del terreno)
# ---------------------------------------------------------------------------

def compute_area3d_ha(polygon_ring: Sequence[Sequence[float]], dem_path: str,
                       cell_size_m: float = 30.0) -> dict:
    """
    Estima el area real de superficie (siguiendo la pendiente del
    terreno) integrando un DEM sobre el poligono, y la compara con el
    area en planta (proyectada).

    Requiere rasterio. Si no esta instalado, lanza un error explicito en
    vez de fallar en silencio o inventar un factor de correccion.

    Devuelve un dict con: area_plana_ha, area_3d_ha, factor_relieve,
    altitud_min_m, altitud_max_m, dH_m (rango altitudinal).
    """
    if not _HAS_RASTERIO:
        raise ImportError(
            "compute_area3d_ha requiere 'rasterio' (pip install rasterio). "
            "No se aplica ningun factor de correccion aproximado por "
            "default -- si se necesita una estimacion sin DEM, calculala "
            "explicitamente y marcala como NO VERIFICADA en el reporte."
        )

    area_plana_ha = compute_polygon_area_ha(polygon_ring)

    with rasterio.open(dem_path) as src:
        from rasterio.mask import mask
        from shapely.geometry import shape

        poligono = shape({"type": "Polygon", "coordinates": [polygon_ring]})
        elevacion, _ = mask(src, [poligono], crop=True, nodata=np.nan)
        elevacion = elevacion[0].astype(float)
        elevacion[elevacion == src.nodata] = np.nan

        dy, dx = np.gradient(elevacion, cell_size_m)
        pendiente_rad = np.arctan(np.sqrt(dx ** 2 + dy ** 2))
        factor_pixel = 1.0 / np.cos(pendiente_rad)
        factor_pixel = np.where(np.isnan(elevacion), np.nan, factor_pixel)

        factor_relieve = np.nanmean(factor_pixel)
        area_3d_ha = area_plana_ha * factor_relieve

        alt_min = float(np.nanmin(elevacion))
        alt_max = float(np.nanmax(elevacion))

    return {
        "area_plana_ha": area_plana_ha,
        "area_3d_ha": area_3d_ha,
        "factor_relieve": float(factor_relieve),
        "altitud_min_m": alt_min,
        "altitud_max_m": alt_max,
        "dH_m": alt_max - alt_min,
    }


# ---------------------------------------------------------------------------
# 5. MODELO ESPECIE-AREA DE WILSON-MACARTHUR (S = c * A^z)
# ---------------------------------------------------------------------------

@dataclass
class WilsonFit:
    c: float
    z: float
    z_fuente: str  # "literatura" o "empirico"
    n_sitios: int


def fit_species_area(
    sites: pd.DataFrame,
    area_col: str,
    richness_col: str,
    z_prior: float = 0.26,
    allow_empirical_z: bool = False,
    min_sitios_para_z_empirico: int = 5,
) -> WilsonFit:
    """
    Ajusta S = c * A^z. Por default usa z de literatura (0.26 es el valor
    citado para bosque mesofilo de montana fragmentado en Veracruz,
    Villaseñor et al. 2007 -- VERIFICAR esta cita antes de citarla tu
    mismo en un dictamen, no se confirmo su exactitud en esta auditoria).

    Solo permite ajustar z de forma empirica (regresion log-log) si hay
    al menos `min_sitios_para_z_empirico` sitios independientes.

    Por que existe esta restriccion
    ---------------------------------
    En el proyecto Texolo se intento ajustar z empiricamente con
    unicamente 3 puntos (nucleo, buffer 500m, buffer 1000m). El
    resultado fue z ~= -0.48, es decir, "menos area = mas especies" --
    ecologicamente sin sentido, un artefacto puro de sobreajustar una
    recta con 3 puntos. Con datos rarefectados el problema fue menor
    pero z seguia saliendo negativo (~-0.10). Regla de este modulo:
    si no hay al menos 5 sitios/poligonos independientes, no se ofrece
    ajuste empirico, se usa z de literatura y se documenta esa eleccion.
    """
    n_sitios = len(sites)
    if allow_empirical_z:
        if n_sitios < min_sitios_para_z_empirico:
            warnings.warn(
                f"Se pidio z empirico pero solo hay {n_sitios} sitios "
                f"(minimo {min_sitios_para_z_empirico}). Se usa z de "
                f"literatura ({z_prior}) en su lugar -- un ajuste con "
                f"tan pocos puntos produce exponentes inestables o "
                f"negativos, como paso en la auditoria del proyecto "
                f"Texolo (z ~= -0.48 con solo 3 puntos)."
            )
        else:
            log_A = np.log(sites[area_col].to_numpy())
            log_S = np.log(sites[richness_col].to_numpy())
            z_emp, log_c = np.polyfit(log_A, log_S, 1)
            if z_emp <= 0:
                warnings.warn(
                    f"z empirico salio <= 0 ({z_emp:.3f}), lo cual no es "
                    f"biologicamente plausible para una relacion "
                    f"especie-area. Se usa z de literatura ({z_prior}) "
                    f"en su lugar. Revisa los datos de entrada."
                )
            else:
                return WilsonFit(c=math.exp(log_c), z=z_emp,
                                  z_fuente="empirico", n_sitios=n_sitios)

    # z de literatura: c se despeja para cada sitio y se promedia
    c_por_sitio = sites[richness_col] / (sites[area_col] ** z_prior)
    return WilsonFit(c=float(c_por_sitio.mean()), z=z_prior,
                      z_fuente="literatura", n_sitios=n_sitios)


# ---------------------------------------------------------------------------
# 6. COMPARACION HONESTA DE MODELOS: K-FOLD + AIC
# ---------------------------------------------------------------------------

def compare_models_kfold(
    df: pd.DataFrame,
    y_col: str,
    modelos: dict[str, list[str]],
    k: int = 5,
    random_state: Optional[int] = None,
) -> pd.DataFrame:
    """
    Compara modelos candidatos (cada uno definido por su lista de columnas
    predictoras) con validacion cruzada K-fold (MAE fuera de muestra) y
    AIC (penaliza complejidad). Requiere statsmodels.

    `modelos` es un dict {nombre: [columnas_predictoras]}, por ejemplo:
        {"Lineal": ["AREA_HA"],
         "Wilson": ["LOG_AREA_HA"],
         "Wilson_Humboldt_3D": ["LOG_AREA_HA", "DH_M"]}

    Por que existe esta funcion
    ----------------------------
    Es la unica forma honesta de justificar que un modelo mas complejo
    (como la extension "Wilson-Humboldt 3D" con termino de altitud) vale
    la pena: si no reduce el error fuera de muestra y no mejora el AIC
    respecto a un modelo mas simple, no se adopta, sin importar cuanto
    esfuerzo se invirtio en construirlo. En la auditoria del proyecto
    Texolo, esta prueba se corrio una sola vez (V27) y el modelo simple
    (Wilson clasico) le gano al extendido en ambos criterios -- ese
    resultado se conservo tal cual, no se manipulo para justificar la
    version mas elaborada.
    """
    if not _HAS_STATSMODELS:
        raise ImportError("compare_models_kfold requiere 'statsmodels'.")

    n = len(df)
    rng = np.random.default_rng(random_state)
    idx = rng.permutation(n)
    folds = np.array_split(idx, k)

    resultados = []
    for nombre, cols in modelos.items():
        X_full = sm.add_constant(df[cols])
        y_full = df[y_col]
        modelo_full = sm.OLS(y_full, X_full).fit()

        errores_abs = []
        for i in range(k):
            test_idx = folds[i]
            train_idx = np.concatenate([folds[j] for j in range(k) if j != i])
            X_train = sm.add_constant(df.iloc[train_idx][cols])
            y_train = df.iloc[train_idx][y_col]
            X_test = sm.add_constant(df.iloc[test_idx][cols], has_constant="add")
            y_test = df.iloc[test_idx][y_col]

            modelo = sm.OLS(y_train, X_train).fit()
            pred = modelo.predict(X_test)
            errores_abs.extend(np.abs(pred.to_numpy() - y_test.to_numpy()))

        resultados.append(
            {
                "modelo": nombre,
                "n_predictores": len(cols),
                "kfold_MAE": float(np.mean(errores_abs)),
                "AIC": float(modelo_full.aic),
                "R2": float(modelo_full.rsquared),
            }
        )

    out = pd.DataFrame(resultados).sort_values("kfold_MAE").reset_index(drop=True)
    out["ganador_por_MAE"] = out["modelo"] == out.loc[out["kfold_MAE"].idxmin(), "modelo"]
    out["ganador_por_AIC"] = out["modelo"] == out.loc[out["AIC"].idxmin(), "modelo"]
    return out


# ---------------------------------------------------------------------------
# 7. CARBONO CON DOBLE FUENTE INDEPENDIENTE (siempre con incertidumbre)
# ---------------------------------------------------------------------------

@dataclass
class EstimacionCarbono:
    fuente_a_nombre: str
    fuente_a_valor_t: float
    fuente_a_incertidumbre_t: float
    fuente_b_nombre: str
    fuente_b_valor_t: float
    fuente_b_incertidumbre_t: float
    combinado_valor_t: float
    combinado_incertidumbre_t: float
    discrepancia_relativa_pct: float
    alerta_discrepancia: bool


def estimate_carbon_dual_source(
    fuente_a_nombre: str, fuente_a_valor_t: float, fuente_a_incertidumbre_t: float,
    fuente_b_nombre: str, fuente_b_valor_t: float, fuente_b_incertidumbre_t: float,
    umbral_alerta_pct: float = 25.0,
) -> EstimacionCarbono:
    """
    Combina dos estimaciones independientes de carbono (p.ej. un mapa
    satelital de biomasa tipo ESA CCI y huellas lidar tipo GEDI L4A) por
    promedio ponderado por varianza inversa, y reporta la incertidumbre
    combinada -- nunca solo el numero central.

    Por que existe esta funcion
    ----------------------------
    Tener dos fuentes independientes que coinciden en orden de magnitud
    es la validacion cruzada mas fuerte que se encontro en todo el
    proyecto Texolo. Pero el margen de error de un producto satelital
    global (ESA CCI) puede ser ~80% del valor central, mientras que
    lidar con miles de huellas puede bajar a ~10% -- reportar solo el
    numero central sin ese margen es enganoso para un dictamen tecnico.
    """
    var_a = fuente_a_incertidumbre_t ** 2
    var_b = fuente_b_incertidumbre_t ** 2
    peso_a = 1.0 / var_a if var_a > 0 else 0.0
    peso_b = 1.0 / var_b if var_b > 0 else 0.0

    if peso_a + peso_b == 0:
        combinado_valor = (fuente_a_valor_t + fuente_b_valor_t) / 2
        combinado_incert = max(fuente_a_incertidumbre_t, fuente_b_incertidumbre_t)
    else:
        combinado_valor = (peso_a * fuente_a_valor_t + peso_b * fuente_b_valor_t) / (peso_a + peso_b)
        combinado_incert = math.sqrt(1.0 / (peso_a + peso_b))

    discrepancia_pct = 100 * abs(fuente_a_valor_t - fuente_b_valor_t) / max(
        combinado_valor, 1e-9
    )

    return EstimacionCarbono(
        fuente_a_nombre=fuente_a_nombre,
        fuente_a_valor_t=fuente_a_valor_t,
        fuente_a_incertidumbre_t=fuente_a_incertidumbre_t,
        fuente_b_nombre=fuente_b_nombre,
        fuente_b_valor_t=fuente_b_valor_t,
        fuente_b_incertidumbre_t=fuente_b_incertidumbre_t,
        combinado_valor_t=combinado_valor,
        combinado_incertidumbre_t=combinado_incert,
        discrepancia_relativa_pct=discrepancia_pct,
        alerta_discrepancia=discrepancia_pct > umbral_alerta_pct,
    )


# ---------------------------------------------------------------------------
# AUTO-TEST -- corre `python salamandra_biodiversidad.py` para verificar
# que el modulo funciona antes de integrarlo a Salamandra.
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=== Auto-test salamandra_biodiversidad.py ===\n")

    # 1. Area desde geometria (cuadrado de ~1km x 1km ~= 100 ha)
    lat0, lon0 = 19.40, -97.00
    d = 0.0045  # ~500m en grados
    cuadrado = [
        (lon0 - d, lat0 - d), (lon0 + d, lat0 - d),
        (lon0 + d, lat0 + d), (lon0 - d, lat0 + d), (lon0 - d, lat0 - d),
    ]
    area = compute_polygon_area_ha(cuadrado)
    print(f"1. Area de poligono de prueba (~100 ha esperadas): {area:.2f} ha")

    # 2. Validacion de clase
    df_test = pd.DataFrame({
        "CLASE": ["Aves", "Mammalia", "Amphibia", "Reptilia"],
        "FAMILIA": ["Parulidae", "Tyrannidae", "Plethodontidae", "Colubridae"],
        "NOMBRE_CIENTIFICO": ["Setophaga virens", "Cathartes aura (mal etiquetado)",
                               "Bolitoglossa platydactyla", "Lampropeltis polyzona"],
    })
    df_validado, reporte = validate_taxonomic_class(df_test)
    print(f"\n2. Validacion de clase: {reporte.registros_corregidos}/"
          f"{reporte.total_registros} registros corregidos "
          f"({reporte.porcentaje_corregido:.0f}%)")
    print(df_validado[["FAMILIA", "CLASE", "CLASE_VALIDADA"]])

    # 3. Rarefaccion
    rng = np.random.default_rng(42)
    especies_pool = [f"spp_{i}" for i in range(50)]
    filas = []
    for nivel, n in [("NUCLEO", 800), ("BUFFER", 200)]:
        for _ in range(n):
            filas.append({"NIVEL": nivel, "ESPECIE": rng.choice(especies_pool)})
    df_rare = pd.DataFrame(filas)
    resultado_rare = rarefy_richness(df_rare, "NIVEL", "ESPECIE", n_target=150, n_iter=200, random_state=1)
    print(f"\n3. Rarefaccion (n_target=150):")
    print(resultado_rare[["NIVEL", "S_raw", "S_rarefecida_media", "S_rarefecida_sd"]])

    # 5. Wilson-MacArthur (con y sin permitir z empirico, pocos sitios)
    sitios = pd.DataFrame({"AREA_HA": [29, 145, 312], "RIQUEZA": [568, 171, 204]})
    fit_lit = fit_species_area(sitios, "AREA_HA", "RIQUEZA")
    print(f"\n5. Wilson (z literatura): c={fit_lit.c:.2f}, z={fit_lit.z} ({fit_lit.z_fuente})")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        fit_emp = fit_species_area(sitios, "AREA_HA", "RIQUEZA", allow_empirical_z=True)
        print(f"   Wilson (z empirico pedido, solo 3 sitios): c={fit_emp.c:.2f}, "
              f"z={fit_emp.z} ({fit_emp.z_fuente}) -- advertencias: {len(w)}")

    # 6. Comparacion de modelos K-fold + AIC
    rng2 = np.random.default_rng(7)
    n_obs = 40
    area_sim = rng2.uniform(5, 500, n_obs)
    dh_sim = rng2.uniform(0, 300, n_obs)
    ruido = rng2.normal(0, 0.15, n_obs)
    log_area = np.log(area_sim)
    log_riqueza = 3.0 + 0.26 * log_area + ruido
    df_modelos = pd.DataFrame({
        "AREA_HA": area_sim,
        "LOG_AREA_HA": log_area,
        "DH_M": dh_sim,
        "LOG_RIQUEZA": log_riqueza,
    })
    if _HAS_STATSMODELS:
        comparacion = compare_models_kfold(
            df_modelos, y_col="LOG_RIQUEZA",
            modelos={
                "Lineal": ["AREA_HA"],
                "Wilson": ["LOG_AREA_HA"],
                "Wilson_Humboldt_3D": ["LOG_AREA_HA", "DH_M"],
            },
            k=5, random_state=1,
        )
        print("\n6. Comparacion de modelos (K-fold MAE + AIC):")
        print(comparacion)
    else:
        print("\n6. (statsmodels no instalado, se omite el test de compare_models_kfold)")

    # 7. Carbono doble fuente
    carbono = estimate_carbon_dual_source(
        "ESA_CCI", 108763, 69433,
        "GEDI_L4A", 123078, 8929,
    )
    print(f"\n7. Carbono combinado: {carbono.combinado_valor_t:.0f} t "
          f"+/- {carbono.combinado_incertidumbre_t:.0f} t "
          f"(discrepancia {carbono.discrepancia_relativa_pct:.1f}%, "
          f"alerta={carbono.alerta_discrepancia})")

    print("\n=== Auto-test completado sin errores ===")

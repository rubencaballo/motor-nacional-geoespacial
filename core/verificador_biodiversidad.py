#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verificador_biodiversidad.py
=============================

Módulo generalizado (cualquier polígono: ANP, finca, predio) con TRES
técnicas, cada una construida para corregir un problema concreto y
verificable en datos de biodiversidad/carbono:

  2. validar_clase_por_familia()
     Cruza cada registro contra una tabla de autoridad FAMILIA -> CLASE,
     en vez de confiar ciegamente en la columna "CLASE" de un CSV
     (esa columna suele venir mal etiquetada en descargas de GBIF/iNat
     cuando el proveedor original clasificó mal el registro).

  3. rarefaccion_monte_carlo()
     Curva de rarefacción por remuestreo sin reemplazo (Monte Carlo),
     para poder comparar riqueza entre zonas/años con esfuerzo de
     muestreo desigual, en vez de comparar conteos crudos de registros.
     Incluye un estimador Chao1 opcional como ancla de riqueza mínima
     esperada (extra, no sustituye la curva).

  4. area_3d_dem()
     Área real de superficie (corregida por relieve) calculada célula
     por célula a partir de un DEM real, usando rasterio. Si rasterio
     no está instalado, o si no se entrega un DEM, la función AVISA
     con un error explícito -- nunca inventa un factor de corrección
     (p. ej. "x1.19") sin datos de elevación detrás.

Diseño explícito:
------------------
- Ninguna superficie se recibe como número hardcodeado: todo polígono
  se lee de un archivo vectorial real (GeoJSON/Shapefile/GPKG) y su
  área se calcula desde la geometría, en una proyección métrica
  adecuada (UTM automático según el centroide). Esto es lo que evita
  el tipo de bug "el archivo se llama _29 pero el polígono mide 504 ha":
  el nombre del archivo nunca es la fuente del área.

Dependencias:
-------------
  obligatorias:  pandas, numpy, geopandas, shapely, pyproj
  opcional:      rasterio  (solo necesaria para area_3d_dem)

Uso como librería:
-------------------
    from verificador_biodiversidad import (
        validar_clase_por_familia,
        rarefaccion_monte_carlo,
        area_3d_dem,
    )

Uso como CLI:
-------------
    python verificador_biodiversidad.py validar-clase \
        --csv registros.csv --col-familia FAMILIA --col-clase CLASE \
        --tabla-familia tabla_familia_clase.csv --out registros_validados.csv

    python verificador_biodiversidad.py rarefaccion \
        --csv registros.csv --col-especie NOMBRE_CIENTIFICO \
        --col-zona NIVEL --n-iter 999 --out curva_rarefaccion.csv

    python verificador_biodiversidad.py area-3d \
        --poligono poligono.geojson --dem dem.tif --out area_3d.json
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

try:
    import geopandas as gpd
    from shapely.geometry.base import BaseGeometry
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "Este módulo requiere geopandas y shapely. Instala con:\n"
        "  pip install geopandas shapely pyproj --break-system-packages"
    ) from exc


# ============================================================================
# 2. VALIDACIÓN DE CLASE CONTRA FAMILIA
# ============================================================================

# Tabla de arranque (NO es la tabla de 240 familias mencionada en tus
# dictámenes -- esa nunca me fue entregada). Cubre familias comunes de
# fauna de bosque templado/mesófilo de montaña en México como PUNTO DE
# PARTIDA para que no arranques de una tabla vacía. Para uso en producción,
# reemplázala por tu tabla de autoridad real vía --tabla-familia (CSV con
# columnas: familia,clase), o mejor aún, derívala del GBIF Backbone
# Taxonomy o del catálogo de CONABIO.
TABLA_FAMILIA_CLASE_DEFECTO: dict[str, str] = {
    # Aves (una muestra representativa, no exhaustiva)
    "Trochilidae": "Aves", "Turdidae": "Aves", "Parulidae": "Aves",
    "Tyrannidae": "Aves", "Corvidae": "Aves", "Columbidae": "Aves",
    "Accipitridae": "Aves", "Picidae": "Aves", "Icteridae": "Aves",
    "Fringillidae": "Aves", "Trogonidae": "Aves", "Strigidae": "Aves",
    "Cardinalidae": "Aves", "Emberizidae": "Aves", "Vireonidae": "Aves",
    # Mammalia
    "Felidae": "Mammalia", "Canidae": "Mammalia", "Cervidae": "Mammalia",
    "Mustelidae": "Mammalia", "Procyonidae": "Mammalia",
    "Sciuridae": "Mammalia", "Vespertilionidae": "Mammalia",
    "Phyllostomidae": "Mammalia", "Didelphidae": "Mammalia",
    "Leporidae": "Mammalia",
    # Reptilia
    "Colubridae": "Reptilia", "Viperidae": "Reptilia",
    "Phrynosomatidae": "Reptilia", "Anguidae": "Reptilia",
    "Teiidae": "Reptilia", "Kinosternidae": "Reptilia",
    # Amphibia
    "Plethodontidae": "Amphibia", "Hylidae": "Amphibia",
    "Ranidae": "Amphibia", "Bufonidae": "Amphibia",
    "Craugastoridae": "Amphibia", "Ambystomatidae": "Amphibia",
    # Insecta
    "Nymphalidae": "Insecta", "Papilionidae": "Insecta",
    "Coccinellidae": "Insecta", "Formicidae": "Insecta",
    "Scarabaeidae": "Insecta", "Cerambycidae": "Insecta",
    "Apidae": "Insecta", "Libellulidae": "Insecta",
    # Arachnida
    "Salticidae": "Arachnida", "Araneidae": "Arachnida",
    "Theraphosidae": "Arachnida",
}


@dataclass
class ResultadoValidacionClase:
    df_validado: pd.DataFrame
    n_total: int
    n_corregidos: int
    pct_corregidos: float
    familias_no_encontradas: list[str]
    resumen_correcciones: pd.DataFrame  # familia, clase_original, clase_validada, n

    def __str__(self) -> str:
        return (
            f"Validación CLASE vs FAMILIA: {self.n_corregidos}/{self.n_total} "
            f"registros ({self.pct_corregidos:.1f}%) tenían CLASE incorrecta "
            f"y fueron corregidos a CLASE_VALIDADA. "
            f"{len(self.familias_no_encontradas)} familias no están en la "
            f"tabla de autoridad y quedaron sin validar (revísalas a mano)."
        )


def cargar_tabla_familia_clase(path: Optional[str] = None) -> dict[str, str]:
    """
    Carga una tabla de autoridad FAMILIA -> CLASE desde un CSV
    (columnas: familia, clase). Si no se da ruta, usa la tabla de
    arranque incluida arriba (NO apta para producción por sí sola).
    """
    if path is None:
        warnings.warn(
            "Usando la tabla FAMILIA->CLASE de arranque incluida en este "
            "script (~50 familias de muestra). Para un dictamen real, "
            "entrega tu propia tabla de autoridad con --tabla-familia.",
            stacklevel=2,
        )
        return dict(TABLA_FAMILIA_CLASE_DEFECTO)

    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    if "familia" not in cols or "clase" not in cols:
        raise ValueError(
            f"{path} debe tener columnas 'familia' y 'clase'. "
            f"Columnas encontradas: {list(df.columns)}"
        )
    tabla = dict(zip(df[cols["familia"]].astype(str), df[cols["clase"]].astype(str)))
    return tabla


def validar_clase_por_familia(
    df: pd.DataFrame,
    tabla_familia_clase: dict[str, str],
    col_familia: str = "FAMILIA",
    col_clase: str = "CLASE",
) -> ResultadoValidacionClase:
    """
    Para cada registro, busca su FAMILIA en la tabla de autoridad y
    determina la CLASE correcta (CLASE_VALIDADA). Si la FAMILIA no está
    en la tabla, se conserva la CLASE original y se marca para revisión
    manual (no se descarta el registro, no se asume nada).

    Esto es lo que atrapa el caso "ave etiquetada como mamífero":
    si CLASE dice "Mammalia" pero la FAMILIA corresponde a un ave
    (p. ej. Trochilidae), CLASE_VALIDADA corrige a "Aves" y el registro
    queda marcado en el reporte de correcciones.
    """
    if col_familia not in df.columns:
        raise KeyError(f"No existe la columna '{col_familia}' en el DataFrame.")
    if col_clase not in df.columns:
        raise KeyError(f"No existe la columna '{col_clase}' en el DataFrame.")

    df = df.copy()
    familia_a_clase = df[col_familia].astype(str).map(tabla_familia_clase)

    familias_no_encontradas = sorted(
        df.loc[familia_a_clase.isna(), col_familia].astype(str).unique().tolist()
    )

    # Si la familia está en la tabla, esa es la verdad; si no, se conserva
    # la clase original tal cual (no hay de dónde más sacarla).
    df["CLASE_VALIDADA"] = familia_a_clase.where(
        familia_a_clase.notna(), df[col_clase].astype(str)
    )

    clase_original = df[col_clase].astype(str)
    difiere = df["CLASE_VALIDADA"] != clase_original
    # Solo cuenta como "corregido" cuando SÍ había una familia conocida
    # (si la familia era desconocida, no hubo corrección real, solo
    # se dejó el valor original).
    corregido_real = difiere & familia_a_clase.notna()

    n_total = len(df)
    n_corregidos = int(corregido_real.sum())
    pct = (n_corregidos / n_total * 100.0) if n_total else 0.0

    resumen = (
        df.loc[corregido_real, [col_familia, col_clase, "CLASE_VALIDADA"]]
        .rename(columns={col_clase: "clase_original"})
        .value_counts()
        .reset_index(name="n")
        .rename(columns={col_familia: "familia"})
        .sort_values("n", ascending=False)
        .reset_index(drop=True)
    )

    return ResultadoValidacionClase(
        df_validado=df,
        n_total=n_total,
        n_corregidos=n_corregidos,
        pct_corregidos=pct,
        familias_no_encontradas=familias_no_encontradas,
        resumen_correcciones=resumen,
    )


# ============================================================================
# 3. RAREFACCIÓN MONTE CARLO
# ============================================================================

def rarefaccion_monte_carlo(
    especies: Sequence[str],
    n_iter: int = 999,
    tamanos_muestra: Optional[Sequence[int]] = None,
    seed: Optional[int] = None,
) -> pd.DataFrame:
    """
    Curva de rarefacción individual-based por remuestreo Monte Carlo
    (sin reemplazo). Para cada tamaño de muestra n, remuestrea `n_iter`
    veces sin reemplazo del conjunto completo de registros y calcula la
    riqueza de especies esperada (media) con intervalo del 95%.

    Esto es lo que permite comparar riqueza entre dos zonas con esfuerzo
    de muestreo distinto: en vez de comparar "298 spp en núcleo" contra
    "200 spp en buffer1000" directamente (donde el núcleo tuvo ~7x más
    registros), se compara la riqueza esperada AL MISMO tamaño de
    muestra en ambas curvas.

    Parámetros
    ----------
    especies : secuencia de identificadores de especie, uno por registro
               (p. ej. df["NOMBRE_CIENTIFICO"].tolist())
    n_iter   : número de remuestreos Monte Carlo por tamaño de muestra
    tamanos_muestra : tamaños de muestra a evaluar; por defecto usa un
               espaciado log-ish desde 1 hasta N total (máx. 30 puntos)
    seed     : semilla para reproducibilidad

    Retorna
    -------
    DataFrame con columnas: n_muestra, riqueza_media, riqueza_ci_low,
    riqueza_ci_high, riqueza_std
    """
    especies_arr = np.asarray(especies)
    n_total = len(especies_arr)
    if n_total == 0:
        raise ValueError("La lista de especies está vacía.")

    rng = np.random.default_rng(seed)

    if tamanos_muestra is None:
        max_puntos = min(30, n_total)
        tamanos_muestra = sorted(
            set(np.unique(np.linspace(1, n_total, max_puntos, dtype=int)))
        )

    filas = []
    for n in tamanos_muestra:
        n = int(n)
        if n < 1 or n > n_total:
            continue
        riquezas = np.empty(n_iter, dtype=int)
        for i in range(n_iter):
            muestra = rng.choice(especies_arr, size=n, replace=False)
            riquezas[i] = len(np.unique(muestra))
        filas.append(
            {
                "n_muestra": n,
                "riqueza_media": float(np.mean(riquezas)),
                "riqueza_ci_low": float(np.percentile(riquezas, 2.5)),
                "riqueza_ci_high": float(np.percentile(riquezas, 97.5)),
                "riqueza_std": float(np.std(riquezas)),
            }
        )

    return pd.DataFrame(filas)


def chao1_riqueza(especies: Sequence[str]) -> dict:
    """
    Estimador Chao1 de riqueza mínima esperada (extra, complementa la
    curva de rarefacción -- no es uno de los tres puntos pedidos, pero
    es prácticamente gratis calcularlo con los mismos datos y da una
    segunda ancla de referencia).

        Chao1 = S_obs + (F1^2) / (2 * F2)     si F2 > 0
        Chao1 = S_obs + F1*(F1-1) / 2         si F2 == 0 (corrección de sesgo)

    donde F1 = # especies representadas por exactamente 1 registro
    (singletons) y F2 = # especies con exactamente 2 registros (doubletons).
    """
    conteos = pd.Series(especies).value_counts()
    s_obs = len(conteos)
    f1 = int((conteos == 1).sum())
    f2 = int((conteos == 2).sum())

    if f2 > 0:
        chao1 = s_obs + (f1 ** 2) / (2 * f2)
    else:
        chao1 = s_obs + f1 * (f1 - 1) / 2

    return {
        "s_observada": s_obs,
        "singletons_f1": f1,
        "doubletons_f2": f2,
        "chao1_estimado": round(float(chao1), 1),
        "especies_no_detectadas_estimadas": round(float(chao1) - s_obs, 1),
    }


# ============================================================================
# 4. ÁREA 3D VÍA DEM
# ============================================================================

@dataclass
class ResultadoArea3D:
    area_plana_ha: float
    area_3d_ha: float
    factor_correccion: float
    pendiente_media_grados: float
    pendiente_max_grados: float
    n_celdas_validas: int
    fuente_dem: str
    cobertura_dem_pct: float = 100.0

    def __str__(self) -> str:
        cobertura_txt = "" if self.cobertura_dem_pct >= 99.95 else f" | Cobertura DEM: {self.cobertura_dem_pct:.1f}%"
        return (
            f"Área plana: {self.area_plana_ha:.2f} ha | "
            f"Área 3D real: {self.area_3d_ha:.2f} ha | "
            f"Factor: x{self.factor_correccion:.3f} | "
            f"Pendiente media: {self.pendiente_media_grados:.1f}° "
            f"(máx {self.pendiente_max_grados:.1f}°) | "
            f"DEM: {self.fuente_dem}{cobertura_txt}"
        )


def _utm_epsg_para_geom(geom_wgs84_centroid_lon: float, geom_wgs84_centroid_lat: float) -> int:
    zona = int((geom_wgs84_centroid_lon + 180) / 6) + 1
    hemisferio = 32600 if geom_wgs84_centroid_lat >= 0 else 32700
    return hemisferio + zona


def area_3d_dem(
    poligono_path: str,
    dem_path: str,
) -> ResultadoArea3D:
    """
    Calcula el área real de superficie (corregida por relieve) de un
    polígono, célula por célula, a partir de un DEM real usando rasterio.

    Método: para cada píxel del DEM dentro del polígono, se estima la
    pendiente local por gradiente (dz/dx, dz/dy) y el área 3D del
    píxel = área_plana_del_píxel / cos(pendiente). Se suman todos los
    píxeles válidos dentro del polígono.

    Si rasterio no está instalado, o si no se entrega DEM, esta función
    LANZA UN ERROR explícito. Nunca devuelve un factor de corrección
    inventado (p. ej. "x1.19" sin datos de elevación detrás).

    Parámetros
    ----------
    poligono_path : ruta a un archivo vectorial (GeoJSON, Shapefile, GPKG...)
                     con el/los polígono(s) de interés
    dem_path       : ruta a un raster DEM (GeoTIFF u otro formato soportado
                     por rasterio) que cubra el polígono

    Retorna
    -------
    ResultadoArea3D
    """
    try:
        import rasterio
        from rasterio.windows import from_bounds as _window_from_bounds, Window as _Window
        from rasterio.features import rasterize as _rasterize
    except ImportError as exc:
        raise ImportError(
            "area_3d_dem() requiere 'rasterio' y no está instalado. "
            "Instala con:\n"
            "  pip install rasterio --break-system-packages\n"
            "Este módulo NO inventa un factor de corrección de superficie "
            "sin datos de elevación reales -- si no puedes instalar "
            "rasterio, usa solo el área plana del polígono y decláralo "
            "explícitamente como área en planta, sin corrección 3D."
        ) from exc

    gdf = gpd.read_file(poligono_path)
    if gdf.empty:
        raise ValueError(f"{poligono_path} no contiene geometrías.")

    # Área plana: reproyecta a UTM automático según el centroide (WGS84)
    gdf_wgs84 = gdf.to_crs(4326)
    centroide = gdf_wgs84.geometry.union_all().centroid if hasattr(
        gdf_wgs84.geometry, "union_all"
    ) else gdf_wgs84.unary_union.centroid
    epsg_utm = _utm_epsg_para_geom(centroide.x, centroide.y)
    gdf_utm = gdf.to_crs(epsg_utm)
    area_plana_m2 = float(gdf_utm.geometry.union_all().area) if hasattr(
        gdf_utm.geometry, "union_all"
    ) else float(gdf_utm.unary_union.area)
    area_plana_ha = area_plana_m2 / 10_000.0

    with rasterio.open(dem_path) as src:
        if src.crs is not None and src.crs.is_geographic:
            raise ValueError(
                f"El DEM ({dem_path}) está en un CRS geográfico (grados), "
                "no proyectado. Reproyéctalo a un CRS métrico (p. ej. la "
                "UTM de la zona) antes de calcular pendiente/área real."
            )

        # Reproyecta el polígono al CRS del DEM
        gdf_dem_crs = gdf.to_crs(src.crs)
        geom_dem_crs = (
            gdf_dem_crs.geometry.union_all() if hasattr(gdf_dem_crs.geometry, "union_all")
            else gdf_dem_crs.unary_union
        )

        px_x = abs(src.transform.a)
        px_y = abs(src.transform.e)
        area_pixel_plana_m2 = px_x * px_y

        # BUG CORREGIDO (detectado al probar con un polígono no rectangular
        # real -- las pruebas originales de esta función solo usaban una
        # caja/rectángulo, donde este problema no se manifiesta):
        # recortar exactamente al polígono (rasterio.mask, crop=True) deja en
        # NaN toda celda fuera del polígono pero dentro de su bounding box --
        # aunque el DEM SÍ tenga elevación real ahí. np.gradient() no es
        # NaN-aware: CUALQUIER celda pegada al borde real del polígono (no
        # solo las esquinas) queda con un vecino NaN y su pendiente sale mal
        # (o se pierde por completo). Para un polígono no rectangular eso es
        # todo el perímetro, no un caso raro -- y sistemáticamente hacía ver
        # el terreno más plano de lo real (factor de relieve subestimado; en
        # la prueba que lo detectó, x1.006 en vez de x1.035 con 15° de
        # pendiente real conocida).
        # Un primer intento de arreglo (rellenar el hueco con el promedio
        # global de elevación del polígono) generó el problema contrario:
        # con pendiente regional fuerte, ese relleno-constante crea un salto
        # artificial justo en el borde real, e infla la pendiente ahí (se
        # detectó con el mismo caso de prueba: pendiente_max subía a 83.7°,
        # un valor que no existe en el DEM real).
        # Corrección final: se lee una VENTANA del raster con un margen de 2
        # píxeles alrededor del polígono (datos de elevación REALES del DEM,
        # no inventados) para que el gradiente en el borde se calcule contra
        # vecinos verdaderos -- y solo se rasteriza el polígono EXACTO
        # (sin margen) para decidir qué celdas cuentan en la suma de área.
        # El margen nunca se cuenta como área propia, solo evita contaminar
        # la pendiente de las celdas reales del borde.
        margen_px = 2
        minx, miny, maxx, maxy = geom_dem_crs.bounds
        minx -= margen_px * px_x
        maxx += margen_px * px_x
        miny -= margen_px * px_y
        maxy += margen_px * px_y
        ventana = _window_from_bounds(minx, miny, maxx, maxy, transform=src.transform)
        ventana = ventana.round_offsets().round_lengths()
        ventana = ventana.intersection(_Window(0, 0, src.width, src.height))
        if ventana.width <= 0 or ventana.height <= 0:
            raise ValueError(
                f"El polígono no se superpone con el DEM ({dem_path}). "
                "Verifica que ambos cubran la misma zona."
            )

        elev = src.read(1, window=ventana).astype(float)
        transform = src.window_transform(ventana)
        nodata = src.nodata
        if nodata is not None:
            elev = np.where(elev == nodata, np.nan, elev)

    # Máscara EXACTA del polígono (sin margen) sobre la malla de la ventana leída.
    poligono_mask = _rasterize(
        [(geom_dem_crs, 1)], out_shape=elev.shape, transform=transform,
        fill=0, default_value=1, dtype=np.uint8,
    ).astype(bool)
    n_celdas_poligono = int(np.sum(poligono_mask))
    validos = poligono_mask & ~np.isnan(elev)
    n_celdas_validas = int(np.sum(validos))
    if n_celdas_validas == 0:
        raise ValueError(
            "No hay celdas válidas del DEM dentro del polígono "
            "(revisa que el DEM cubra la zona y que no sea todo nodata)."
        )

    # COBERTURA REAL DEL DEM DENTRO DEL POLÍGONO -- detectado con datos reales
    # (corrida de Texolo, 02/09/2026): si el DEM tiene huecos (vacíos SRTM ya
    # excluidos como nodata en vez de contaminar la elevación -- ver
    # core/geomatica.cargar_dem_utm()) dentro del polígono, area_3d_ha se
    # calcula solo sobre las celdas que SÍ tienen dato, mientras que
    # area_plana_ha (líneas arriba) es el área EXACTA del polígono vectorial
    # completo -- comparar una contra otra cuando falta cobertura real
    # entonces puede dar factor_correccion < 1.0, que es matemáticamente
    # imposible (1/cos(pendiente) siempre es >=1) y sería una alarma real de
    # bug SI la cobertura fuera completa. Con cobertura incompleta, ese <1.0
    # no es un bug de este código -- es que hay un hueco real de dato que no
    # se debe tapar con un número. Igual que en core/geomatica.py (ver su
    # REGLA DE ORO), NO se aplica ningún piso artificial que esconda esto --
    # se lanza un error explícito en vez de devolver un factor engañoso.
    cobertura_pct = (n_celdas_validas / n_celdas_poligono * 100.0) if n_celdas_poligono > 0 else 0.0
    UMBRAL_COBERTURA_MINIMA_PCT = 95.0
    if cobertura_pct < UMBRAL_COBERTURA_MINIMA_PCT:
        raise ValueError(
            f"El DEM ({dem_path}) solo cubre {cobertura_pct:.1f}% del polígono con datos de "
            f"elevación válidos ({n_celdas_validas}/{n_celdas_poligono} celdas) -- por debajo del "
            f"mínimo aceptable ({UMBRAL_COBERTURA_MINIMA_PCT:.0f}%). Esto normalmente significa que "
            "el DEM tiene huecos reales (vacíos SRTM) dentro de esta zona, o que no cubre el polígono "
            "por completo. No se calcula un factor de relieve con cobertura tan incompleta -- revisa "
            "el DEM (ej. cuántas celdas 'SRTM crudo: ... celdas fuera de rango' imprimió "
            "cargar_dem_utm()/preparar_dem_utm.py para este sitio) antes de usarlo en el dictamen."
        )

    # Última red de seguridad: si TODAVÍA queda algún NaN dentro de la
    # ventana con margen (nodata real del DEM, ej. un vacío SRTM, o el
    # polígono está pegado al borde del raster) se rellena solo para el
    # cálculo del gradiente -- con el promedio LOCAL de esta ventana chica
    # (2 píxeles de margen), no un promedio global del sitio, así que el
    # riesgo de un salto artificial como el descrito arriba es mínimo.
    elev_para_gradiente = elev
    if np.any(np.isnan(elev)):
        relleno_local = np.nanmean(elev) if np.any(~np.isnan(elev)) else 0.0
        elev_para_gradiente = np.where(np.isnan(elev), relleno_local, elev)

    # Gradiente en Y y X (metros por metro)
    dzdy, dzdx = np.gradient(elev_para_gradiente, px_y, px_x)
    pendiente_rad = np.arctan(np.sqrt(dzdx ** 2 + dzdy ** 2))
    pendiente_grados = np.where(validos, np.degrees(pendiente_rad), np.nan)

    area_3d_por_celda_m2 = area_pixel_plana_m2 / np.cos(pendiente_rad)
    area_3d_m2 = float(np.nansum(np.where(validos, area_3d_por_celda_m2, np.nan)))
    area_3d_ha = area_3d_m2 / 10_000.0

    factor = area_3d_ha / area_plana_ha if area_plana_ha > 0 else float("nan")

    return ResultadoArea3D(
        area_plana_ha=round(area_plana_ha, 2),
        area_3d_ha=round(area_3d_ha, 2),
        factor_correccion=round(factor, 4),
        pendiente_media_grados=round(float(np.nanmean(np.where(validos, pendiente_grados, np.nan))), 2),
        pendiente_max_grados=round(float(np.nanmax(np.where(validos, pendiente_grados, np.nan))), 2),
        n_celdas_validas=n_celdas_validas,
        fuente_dem=str(Path(dem_path).name),
        cobertura_dem_pct=round(cobertura_pct, 2),
    )


# ============================================================================
# CLI
# ============================================================================

def _cmd_validar_clase(args: argparse.Namespace) -> None:
    df = pd.read_csv(args.csv)
    tabla = cargar_tabla_familia_clase(args.tabla_familia)
    resultado = validar_clase_por_familia(
        df, tabla, col_familia=args.col_familia, col_clase=args.col_clase
    )
    print(resultado)
    if not resultado.resumen_correcciones.empty:
        print("\nCorrecciones aplicadas (top 20):")
        print(resultado.resumen_correcciones.head(20).to_string(index=False))
    if resultado.familias_no_encontradas:
        print(
            f"\n[AVISO] {len(resultado.familias_no_encontradas)} familias sin "
            f"entrada en la tabla de autoridad, sin validar:\n  "
            + ", ".join(resultado.familias_no_encontradas[:30])
            + (" ..." if len(resultado.familias_no_encontradas) > 30 else "")
        )
    if args.out:
        resultado.df_validado.to_csv(args.out, index=False)
        print(f"\nGuardado: {args.out}")


def _cmd_rarefaccion(args: argparse.Namespace) -> None:
    df = pd.read_csv(args.csv)
    if args.col_zona and args.col_zona in df.columns:
        grupos = df.groupby(args.col_zona)
    else:
        grupos = [(None, df)]

    resultados = []
    for zona, sub in grupos:
        curva = rarefaccion_monte_carlo(
            sub[args.col_especie].astype(str).tolist(),
            n_iter=args.n_iter,
            seed=args.seed,
        )
        curva["zona"] = zona if zona is not None else "TODO"
        resultados.append(curva)

        chao1 = chao1_riqueza(sub[args.col_especie].astype(str).tolist())
        print(f"\n--- Zona: {zona or 'TODO'} (n={len(sub)} registros) ---")
        print(f"  Chao1 estimado: {chao1['chao1_estimado']} spp "
              f"(observadas: {chao1['s_observada']}, "
              f"no detectadas est.: {chao1['especies_no_detectadas_estimadas']})")
        print(f"  Riqueza al máximo n muestreado: "
              f"{curva.iloc[-1]['riqueza_media']:.1f} "
              f"[{curva.iloc[-1]['riqueza_ci_low']:.0f}-{curva.iloc[-1]['riqueza_ci_high']:.0f}]")

    tabla_final = pd.concat(resultados, ignore_index=True)
    if args.out:
        tabla_final.to_csv(args.out, index=False)
        print(f"\nGuardado: {args.out}")


def _cmd_area_3d(args: argparse.Namespace) -> None:
    resultado = area_3d_dem(args.poligono, args.dem)
    print(resultado)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(resultado.__dict__, f, indent=2, ensure_ascii=False)
        print(f"\nGuardado: {args.out}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Verificador de biodiversidad/carbono generalizado "
                     "(polígono ANP, finca o predio) — puntos 2, 3 y 4."
    )
    sub = parser.add_subparsers(dest="comando", required=True)

    p2 = sub.add_parser("validar-clase", help="Validación CLASE vs FAMILIA")
    p2.add_argument("--csv", required=True)
    p2.add_argument("--col-familia", default="FAMILIA")
    p2.add_argument("--col-clase", default="CLASE")
    p2.add_argument("--tabla-familia", default=None,
                     help="CSV con columnas familia,clase (tabla de autoridad)")
    p2.add_argument("--out", default=None)
    p2.set_defaults(func=_cmd_validar_clase)

    p3 = sub.add_parser("rarefaccion", help="Rarefacción Monte Carlo")
    p3.add_argument("--csv", required=True)
    p3.add_argument("--col-especie", default="NOMBRE_CIENTIFICO")
    p3.add_argument("--col-zona", default=None,
                     help="Columna para calcular una curva por zona (opcional)")
    p3.add_argument("--n-iter", type=int, default=999)
    p3.add_argument("--seed", type=int, default=None)
    p3.add_argument("--out", default=None)
    p3.set_defaults(func=_cmd_rarefaccion)

    p4 = sub.add_parser("area-3d", help="Área 3D real vía DEM")
    p4.add_argument("--poligono", required=True)
    p4.add_argument("--dem", required=True)
    p4.add_argument("--out", default=None)
    p4.set_defaults(func=_cmd_area_3d)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

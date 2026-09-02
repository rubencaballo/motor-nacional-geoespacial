#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Diagnóstico puntual, DESECHABLE (no modifica core/cuenca_completa.py):
por qué el catchment() delimitado desde el punto de salida dado (Cascada de
Texolo, 19.4213,-97.0295) sale tan chico (2.94 ha con --percentil-snap 85,
29.85 ha con el default 95) a pesar de que el ajuste (snap) cae muy cerca
del punto dado (56m).

Reusa el DEM del corredor ya descargado en caché (no vuelve a bajar nada) y
corre el mismo pipeline D8 de pysheds hasta accumulation() -- pero en vez de
delimitar la cuenca completa, imprime:
  1) La acumulación de flujo EXACTA en el pixel del punto dado.
  2) La acumulación MÁXIMA dentro de una ventana de ~300m alrededor (busca
     si hay un cauce más fuerte muy cerca que quizá el snap no encontró).
  3) En qué percentil, respecto a TODO el corredor, caen esos dos valores --
     para confirmar si el problema es que un cauce real pero "chico a
     escala de corredor" (típico de una barranca angosta como Texolo) queda
     por debajo de cualquier percentil global razonable.

Ajusta las 4 constantes de abajo si tu geojson/coordenadas/id difieren."""
import os

import numpy as np
import pyproj
import geopandas as gpd
import rasterio
from pysheds.grid import Grid

from core.cuenca_completa import calcular_bbox_corredor, cargar_dem_corredor

# --- AJUSTA AQUÍ SI ES NECESARIO ---
GEOJSON_ORIGEN = "Cofre_de_Perote.geojson"  # ruta relativa, igual que en tus comandos anteriores
PUNTO_SALIDA_LATLON = (19.4213, -97.0295)  # (lat, lon) de la Cascada de Texolo
ID_CORREDOR = "Perote_a_Texolo"
CARPETA_SRTM = os.path.expanduser("~/srtm_temp")
RADIO_VENTANA_PIXELES = 10  # ~300m a resolución SRTM de 30m
# ------------------------------------

gdf_origen = gpd.read_file(GEOJSON_ORIGEN)
geom_wgs84_origen = (gdf_origen.geometry.union_all() if hasattr(gdf_origen.geometry, "union_all")
                      else gdf_origen.geometry.unary_union)
lat_salida, lon_salida = PUNTO_SALIDA_LATLON
punto_salida_wgs84 = (lon_salida, lat_salida)
bbox = calcular_bbox_corredor(geom_wgs84_origen, punto_salida_wgs84)

print(f"Reusando DEM del corredor en caché (bbox: {tuple(round(b, 4) for b in bbox)})...")
dst_array, meta_utm, utm_crs = cargar_dem_corredor(bbox, ID_CORREDOR, CARPETA_SRTM)

proyector = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True).transform
x_utm, y_utm = proyector(lon_salida, lat_salida)

# --- Mismo pipeline D8 que delimitar_cuenca(), pero solo hasta accumulation() ---
valid = ~np.isnan(dst_array) & (dst_array > -1000)
z_filled = np.where(valid, dst_array, np.nanmean(dst_array[valid]) if np.any(valid) else 0.0)
meta_pysheds = {
    "driver": "GTiff", "height": z_filled.shape[0], "width": z_filled.shape[1],
    "count": 1, "dtype": "float32", "crs": utm_crs, "transform": meta_utm["transform"], "nodata": -9999,
}
temp_path = os.path.join(CARPETA_SRTM, f"temp_diagnostico_{ID_CORREDOR}.tif")
with rasterio.open(temp_path, "w", **meta_pysheds) as dst:
    dst.write(np.where(valid, z_filled, -9999).astype(np.float32), 1)

print("Corriendo D8 (fill_pits -> fill_depressions -> resolve_flats -> flowdir -> accumulation)...")
grid = Grid.from_raster(temp_path)
dem = grid.read_raster(temp_path)
pit_filled = grid.fill_pits(dem)
flooded = grid.fill_depressions(pit_filled)
inflated = grid.resolve_flats(flooded)
dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
fdir = grid.flowdir(inflated, dirmap=dirmap)
acc = grid.accumulation(fdir, dirmap=dirmap)
acc_arr = np.asarray(acc)

transform = meta_utm["transform"]
rows, cols = acc_arr.shape
col_f, row_f = ~transform * (x_utm, y_utm)
col_i, row_i = int(round(col_f)), int(round(row_f))

print(f"\nPunto dado (UTM): ({x_utm:.1f}, {y_utm:.1f}) -> pixel (fila={row_i}, col={col_i}) de {rows}x{cols}")
acc_punto = float(acc_arr[row_i, col_i])
print(f"Acumulación de flujo EXACTA en ese pixel: {acc_punto:.1f} celdas contribuyentes")

proj_a_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)
validos = acc_arr[valid]

# --- Barrido a radios crecientes (no solo un radio fijo) -- para ver a qué
# distancia del punto dado empieza a aparecer un cauce de verdad importante,
# en vez de conformarnos con "el mejor dentro de 300m", que puede seguir
# siendo un tributario chico si el punto dado está cerca de la cabecera. ---
print("\n--- Barrido de acumulación máxima a radios crecientes ---")
radios_m = [150, 300, 600, 1000, 1500, 2500, 4000]
mejor_global = (acc_punto, row_i, col_i, 0.0)
for radio_m in radios_m:
    radio_px = max(1, int(round(radio_m / 30.0)))
    r0, r1 = max(0, row_i - radio_px), min(rows, row_i + radio_px + 1)
    c0, c1 = max(0, col_i - radio_px), min(cols, col_i + radio_px + 1)
    ventana = acc_arr[r0:r1, c0:c1]
    idx_local = np.unravel_index(np.argmax(ventana), ventana.shape)
    row_max, col_max = r0 + idx_local[0], c0 + idx_local[1]
    acc_max_local = float(acc_arr[row_max, col_max])
    x_max, y_max = transform * (col_max + 0.5, row_max + 0.5)
    dist_max = float(np.hypot(x_max - x_utm, y_max - y_utm))
    lon_max, lat_max = proj_a_wgs84.transform(x_max, y_max)
    pctl_max = float((validos < acc_max_local).mean() * 100)
    print(f"  Radio ~{radio_m}m: máximo = {acc_max_local:.0f} celdas (percentil {pctl_max:.1f}), "
          f"a {dist_max:.0f}m del punto dado -- lat={lat_max:.6f}, lon={lon_max:.6f}")
    if acc_max_local > mejor_global[0]:
        mejor_global = (acc_max_local, row_max, col_max, dist_max)

acc_mejor, row_mejor, col_mejor, dist_mejor = mejor_global
x_mejor, y_mejor = transform * (col_mejor + 0.5, row_mejor + 0.5)
lon_mejor, lat_mejor = proj_a_wgs84.transform(x_mejor, y_mejor)
print(f"\nEl cauce MÁS FUERTE encontrado en todo el barrido: {acc_mejor:.0f} celdas, "
      f"a {dist_mejor:.0f}m del punto dado -- lat={lat_mejor:.6f}, lon={lon_mejor:.6f}")
print("  -> revísalo en Google Earth/Maps: si ESE sí coincide con el cauce visible de la barranca,")
print("     probablemente sea mejor punto de salida que las coordenadas originales.")

acc_max_local, dist_max = acc_mejor, dist_mejor  # para que el resumen final de abajo use el mejor del barrido

pctl_punto = float((validos < acc_punto).mean() * 100)
pctl_max_local = float((validos < acc_max_local).mean() * 100)
print(f"\nPercentil del punto dado dentro de TODO el corredor: {pctl_punto:.1f}")
print(f"Percentil del máximo local (~{RADIO_VENTANA_PIXELES * 30}m) dentro de TODO el corredor: {pctl_max_local:.1f}")
print("\n(Si estos percentiles están muy por debajo de 85-95, confirma que Texolo es un cauce real pero")
print(" 'chico' comparado con todo el corredor -- típico de una barranca angosta mal resuelta por SRTM de 30m")
print(" a esta escala. Con ese dato decidimos el siguiente paso: bajar mucho el percentil, o cambiar la")
print(" estrategia de snap para que busque el máximo LOCAL en vez de comparar contra el corredor entero.)")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Todo lo que habla con Earth Engine: medicion Hansen/JRC/Sentinel-2,
LandTrendr basico, NICFI real, y clasificacion de confianza por pixeles."""
import json

from config import (
    HANSEN_DATASET, MIN_CLUSTER_HANSEN, CORTE_LOSSYEAR,
    MIN_PIXELES_CONFIANZA_ALTA, MIN_PIXELES_CONFIANZA_MEDIA,
    AREA_MIN_HANSEN_HA, ALTITUD_AGUACATE, log,
)
from modelos.evidencia import EvidenciaNICFI
from modelos.predio import PredioNacional


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



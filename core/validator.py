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


def eliminar_traslapes_jerarquico(gdf, campo_prioridad=None):
    """
    Elimina traslapes entre polígonos de forma jerárquica: el polígono de
    mayor prioridad conserva toda su área; cada polígono siguiente se
    recorta (gpd.overlay how='difference') contra la UNIÓN de todos los de
    mayor prioridad ya procesados. Así ninguna hectárea se cuenta dos veces.

    campo_prioridad: nombre de columna en gdf que indica prioridad legal
    real (ej. fecha de título agrario -- más antiguo = más prioridad). Si
    no se proporciona, se usa área descendente como criterio de RESPALDO
    ARBITRARIO -- esto se advierte explícitamente porque un criterio de
    área no tiene ningún fundamento legal; solo evita que el cálculo
    truene mientras no se tenga el campo real.

    Devuelve un gdf nuevo con geometrías recortadas y una columna
    'ha_perdidas_por_traslape' documentando cuánta área se le restó a cada
    predio (para que quede trazable, no oculto).
    """
    import geopandas as gpd
    from shapely.ops import unary_union

    if len(gdf) <= 1:
        gdf = gdf.copy()
        gdf['ha_perdidas_por_traslape'] = 0.0
        return gdf

    gdf = gdf.copy()
    gdf['_area_original_m2'] = gdf.geometry.area

    if campo_prioridad and campo_prioridad in gdf.columns:
        gdf = gdf.sort_values(campo_prioridad, ascending=True).reset_index(drop=True)
        log(f"Traslapes: ordenando por prioridad real '{campo_prioridad}'.", "OK")
    else:
        gdf = gdf.sort_values('_area_original_m2', ascending=False).reset_index(drop=True)
        log("Traslapes: SIN campo de prioridad legal -- usando área descendente como "
            "respaldo ARBITRARIO. Esto debe reemplazarse por un criterio real "
            "(antigüedad de título u otro) antes de un entregable oficial a SADER/ICM.", "WARN")

    territorio_ocupado = None
    geometrias_finales = []
    ha_perdidas = []

    for idx, row in gdf.iterrows():
        geom = row.geometry
        if territorio_ocupado is not None and not territorio_ocupado.is_empty:
            geom_recortada = geom.difference(territorio_ocupado)
        else:
            geom_recortada = geom

        area_antes_m2 = geom.area
        area_despues_m2 = geom_recortada.area if not geom_recortada.is_empty else 0.0
        ha_perdidas.append(max(0.0, (area_antes_m2 - area_despues_m2)) / 10000)
        geometrias_finales.append(geom_recortada)

        territorio_ocupado = geom if territorio_ocupado is None else unary_union([territorio_ocupado, geom])

    gdf['geometry'] = geometrias_finales
    gdf['ha_perdidas_por_traslape'] = ha_perdidas
    gdf = gdf[~gdf.geometry.is_empty].copy()

    n_afectados = sum(1 for h in ha_perdidas if h > 0.0001)
    total_ha_recortadas = sum(ha_perdidas)
    log(f"Traslapes resueltos: {n_afectados} predios recortados, {total_ha_recortadas:.2f} ha "
        f"de doble conteo eliminadas.", "OK")

    return gdf


# ============== ML ==============

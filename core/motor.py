def motor_nacional_masivo(lista):
    for p in lista:
        p.formula_elegibilidad_nacional()
        p.sistema = p.formula_sistema_productivo_nacional()
        p.score_prioridad = p.formula_prioridad_nacional()
        p.validar_completitud()
        p.generar_auditoria()
    return lista


# ============== FIX 4: TRASLAPES JERÁRQUICOS ==============
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
def procesar_estado_poderoso(gdf_estado, estado_nombre, modo="prueba", con_nicfi=False):
    predios = []

    if modo == "gee":
        try:
            import ee
        except ImportError:
            raise RuntimeError("modo='gee' requiere earthengine-api instalado y ee.Initialize() ya corrido.")
        predios = _medir_estado_gee_real(gdf_estado, estado_nombre, ee, con_nicfi=con_nicfi)
    else:
        for idx, row in gdf_estado.iterrows():
            r = random.random()
            if r < 0.15:
                B, ndvi, std, alt, textura, bfast, JRC, ESRI = random.uniform(30, 65), random.uniform(0.71, 0.82), random.uniform(0.03, 0.07), random.uniform(900, 1400), random.uniform(85, 150), None, random.uniform(20, 60), random.uniform(10, 30)
                sup = random.uniform(0.2, 1.5)
                nicfi_disp = con_nicfi
                nicfi_ant = random.uniform(0.68, 0.78) if nicfi_disp else None
                nicfi_rec = random.uniform(0.66, 0.76) if nicfi_disp else None
                n_pix_sim = int(sup * 11)
            elif r < 0.28:
                B, ndvi, std, alt, textura, bfast, JRC, ESRI = random.uniform(25, 60), random.uniform(0.76, 0.88), random.uniform(0.05, 0.09), random.uniform(1600, 2200), random.uniform(20, 45), random.randint(2015, 2023), random.uniform(15, 50), random.uniform(20, 60)
                sup = random.uniform(1.0, 5.0)
                nicfi_disp = False; nicfi_ant = None; nicfi_rec = None
                n_pix_sim = int(sup * 11)
            elif r < 0.45:
                B, ndvi, std, alt, textura, bfast, JRC, ESRI = random.uniform(40, 75), random.uniform(0.70, 0.85), random.uniform(0.04, 0.08), random.uniform(100, 600), random.uniform(20, 40), None, random.uniform(25, 70), random.uniform(10, 25)
                sup = random.uniform(2.0, 10.0)
                nicfi_disp = False; nicfi_ant = None; nicfi_rec = None
                n_pix_sim = int(sup * 11)
            elif r < 0.65:
                B, ndvi, std, alt, textura, bfast, JRC, ESRI = random.uniform(20, 40), random.uniform(0.56, 0.70), random.uniform(0.08, 0.14), random.uniform(400, 1200), random.uniform(85, 180), None, random.uniform(20, 55), random.uniform(15, 40)
                sup = random.uniform(0.5, 3.0)
                nicfi_disp = con_nicfi and random.random() < 0.5
                nicfi_ant = random.uniform(0.60, 0.72) if nicfi_disp else None
                nicfi_rec = random.uniform(0.58, 0.70) if nicfi_disp else None
                n_pix_sim = int(sup * 11)
            else:
                B, ndvi, std, alt, textura, bfast, JRC, ESRI = random.uniform(5, 15), random.uniform(0.35, 0.55), random.uniform(0.12, 0.22), random.uniform(200, 800), random.uniform(30, 70), None, random.uniform(5, 20), random.uniform(60, 90)
                sup = random.uniform(0.3, 2.0)
                nicfi_disp = False; nicfi_ant = None; nicfi_rec = None
                n_pix_sim = int(sup * 11)

            D_valid = random.uniform(0.06, 2.0) if random.random() < 0.12 else 0.0
            D_desc = random.uniform(0.02, 0.15) if random.random() < 0.20 else 0.0
            D_hist = random.uniform(0.0, 1.5) if random.random() < 0.25 else 0.0
            T = random.randint(0, 24)
            jrc_water, pend = 0, random.uniform(2, 28)
            caida = ((nicfi_ant - nicfi_rec) / nicfi_ant * 100) if nicfi_ant and nicfi_rec and nicfi_ant > 0.05 else None
            confianza_sim = _clasificar_confianza_hansen(n_pix_sim)

            p = PredioNacional(
                id_predio=str(row.get('id_predio') or row.get('ID') or f"{estado_nombre[:3]}-{idx:07d}"),
                estado=estado_nombre, municipio=str(row.get('municipio', '')),
                geometria=row.geometry.__geo_interface__ if hasattr(row.geometry, '__geo_interface__') else {},
                superficie_ha=sup, B_treecover=B, B_ndvi=ndvi, B_ndvi_std=std, B_textura=textura, B_altitud=alt, B_bfast_anio=bfast,
                JRC_pct=JRC, ESRI_crops_pct=ESRI, T_perdida=T, D_ha_validada=D_valid, D_ha_descartada_aislada=D_desc,
                D_ha_historica_pre2020=D_hist, D_jrc_water=jrc_water, D_pendiente=pend,
                ndvi_2020=ndvi - 0.05, ndvi_2024=ndvi, nbdi_2020=0.1, nbdi_2024=0.12,
                n_pixeles_hansen=n_pix_sim, confianza_hansen=confianza_sim,
                nicfi=EvidenciaNICFI(disponible=nicfi_disp, ndvi_antiguo=nicfi_ant, ndvi_reciente=nicfi_rec, caida_pct=caida)
            )
            predios.append(p)

    resultados = motor_nacional_masivo(predios)

    if not resultados:
        print(f"[{estado_nombre}] 0 predios procesados (revisa la fuente de datos).")
        return 0

    conn = sqlite3.connect(DB_PATH)
    data = []
    for r in resultados:
        data.append((r.id_predio, r.estado, r.municipio, r.superficie_ha,
                      r.B_treecover, r.B_ndvi, r.B_ndvi_std, r.B_textura, r.B_altitud, r.B_bfast_anio,
                      r.JRC_pct, r.ESRI_crops_pct, r.T_perdida,
                      r.D_ha_validada, r.D_ha_descartada_aislada, r.D_ha_historica_pre2020, r.D_ha_2025, r.D_jrc_water, r.D_pendiente,
                      r.ndvi_2020, r.ndvi_2024, r.nbdi_2020, r.nbdi_2024,
                      r.n_pixeles_hansen, r.confianza_hansen,
                      int(r.nicfi.disponible), r.nicfi.ndvi_antiguo or 0, r.nicfi.ndvi_reciente or 0, r.nicfi.caida_pct or 0,
                      r.curp, r.rfc, r.net_mass_kg, r.legal_docs, r.production_date,
                      int(r.E_elegible), r.color, r.dictamen_corto, r.dictamen_largo, r.texto_evidencia,
                      r.sistema, int(r.requiere_revision), r.status_dds, r.score_prioridad,
                      r.hash_geo, r.checksum_integridad, r.fecha, r.version))
    conn.executemany(f"INSERT OR REPLACE INTO nacional_poderoso VALUES ({','.join(['?']*47)})", data)
    conn.commit()
    conn.close()

    verdes = sum(1 for r in resultados if r.color == "Verde")
    amarillos = sum(1 for r in resultados if r.color == "Amarillo")
    rojos = sum(1 for r in resultados if r.color == "Rojo")
    print(f"[{estado_nombre}] {len(resultados)} predios | Verde {verdes} | Amarillo {amarillos} | Rojo {rojos} | "
          f"Sistemas: {', '.join(set(r.sistema for r in resultados))} | "
          f"NICFI rescate: {sum(1 for r in resultados if r.nicfi.disponible and r.color=='Verde')} | "
          f"Confianza Hansen: {', '.join(set(r.confianza_hansen for r in resultados))}")
    return len(resultados)


# ============== PRODUCTOS ICM ==============
def generar_productos_icm():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    print("\n=== GENERANDO PRODUCTOS ICM-SADER ===")
    cur.execute("SELECT COUNT(*), SUM(superficie_ha) FROM nacional_poderoso WHERE E_elegible=1")
    row = cur.fetchone()
    area_elegible = row[1] if row[1] is not None else 0.0
    print(f"Producto 2 - Áreas elegibles: {row[0]} predios, {area_elegible:.2f} ha")
    print("Producto 3 - Sistemas productivos:")
    for r in cur.execute("SELECT sistema, COUNT(*), AVG(B_treecover), AVG(B_ndvi) FROM nacional_poderoso GROUP BY sistema"):
        tc = r[2] if r[2] is not None else 0.0
        nd = r[3] if r[3] is not None else 0.0
        print(f"  {r[0]}: {r[1]} predios | treecover {tc:.1f}% NDVI {nd:.2f}")
    cur.execute("SELECT COUNT(DISTINCT estado), COUNT(*) FROM nacional_poderoso")
    print(f"Producto 4 - BD predios: {cur.fetchone()}")
    print("Producto 5 - Prioritarios mitigación:")
    for r in cur.execute("SELECT color, COUNT(*), AVG(score) FROM nacional_poderoso GROUP BY color"):
        sc = r[2] if r[2] is not None else 0.0
        print(f"  {r[0]}: {r[1]} predios score {sc:.3f}")
    print("Confianza Hansen (nuevo en v9.2):")
    for r in cur.execute("SELECT confianza_hansen, COUNT(*) FROM nacional_poderoso GROUP BY confianza_hansen"):
        print(f"  {r[0]}: {r[1]} predios")

    export_path = os.path.join(EXPORT_DIR, "paquete_auditoria_nacional_EUDR.json")
    data = []
    for r in cur.execute("SELECT id_predio, superficie_ha, E_elegible, color, sistema, hash_geo, checksum_integridad, "
                          "fecha, B_treecover, D_ha_validada, D_ha_historica_pre2020, JRC_pct, confianza_hansen "
                          "FROM nacional_poderoso LIMIT 100"):
        data.append({
            "plotId": r[0], "area_ha": r[1], "deforestationFree": bool(r[2]), "color": r[3], "sistema": r[4],
            "auditProof": {
                "hash_geo": r[5], "checksum_integridad": r[6],
                "tipo_checksum": "sha256_no_es_firma_digital",
                "fecha": r[7], "B": r[8], "D_post_corte": r[9], "D_historica_pre2020": r[10],
                "JRC": r[11], "confianza_hansen": r[12]
            }
        })
    with open(export_path, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"Producto 6 - Paquete auditoría exportado: {export_path} ({len(data)} muestras)")
    conn.close()


# ============== MAIN ==============
def main():
    ap = argparse.ArgumentParser(description="Motor Nacional ICM-SADER v9.2")
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--pruebita", type=int)
    ap.add_argument("--ran", type=str)
    ap.add_argument("--modo", type=str, choices=["prueba", "gee"], default="prueba")
    ap.add_argument("--proyecto", type=str)
    ap.add_argument("--con-nicfi", action="store_true")
    ap.add_argument("--con-ml", type=str)
    ap.add_argument("--entrenar-ml", type=str)
    ap.add_argument("--max-estados", type=int, default=0)
    ap.add_argument("--resolver-traslapes", nargs="?", const="__auto__", default=None,
                     help="Activa eliminación jerárquica de traslapes. Opcionalmente pasa el nombre "
                          "del campo de prioridad legal real (ej. fecha_titulo). Sin argumento, usa "
                          "área descendente como respaldo arbitrario (con advertencia).")
    ap.add_argument("--reset", action="store_true")
    args = ap.parse_args()

    if args.reset:
        for f in [DB_PATH, CHECKPOINT_PATH]:
            if os.path.exists(f):
                os.remove(f)

    init_db_nacional()
    init_ml_global(modelo_path=args.con_ml if args.con_ml else None)

    if args.modo == "gee":
        import ee
        try:
            if args.proyecto:
                ee.Initialize(project=args.proyecto)
                log(f"Earth Engine inicializado con proyecto: {args.proyecto}", "OK")
            else:
                ee.Initialize()
                log("Earth Engine inicializado (sin --proyecto explícito).", "OK")
        except Exception as e:
            log(f"No se pudo inicializar Earth Engine: {e}", "ERROR")
            sys.exit(1)

    if args.entrenar_ml:
        ADAPTADOR_ML_GLOBAL.entrenar_desde_campo(args.entrenar_ml, args.con_ml or "./modelo_sistema_poderoso.pkl")
        return

    if args.demo:
        print("=== MOTOR NACIONAL v9.2 - DEMO (corte EUDR + confianza Hansen por píxeles) ===\n")
        demos = [
            PredioNacional("VER-CAFE-0.25ha-SIN-NICFI", estado="Veracruz", municipio="Teocelo", superficie_ha=0.25,
                            B_treecover=40, B_ndvi=0.75, B_ndvi_std=0.05, B_textura=100, B_altitud=1100,
                            JRC_pct=35, ESRI_crops_pct=20, T_perdida=5, D_ha_validada=0.0, D_ha_descartada_aislada=0.09,
                            n_pixeles_hansen=3, confianza_hansen="baja", nicfi=EvidenciaNICFI(disponible=False)),
            PredioNacional("VER-CAFE-0.8ha-HUECO-ESCALA-SIN-NICFI", estado="Veracruz", municipio="Teocelo", superficie_ha=0.8,
                            B_treecover=40, B_ndvi=0.75, B_ndvi_std=0.05, B_textura=100, B_altitud=1100,
                            JRC_pct=12, ESRI_crops_pct=20, T_perdida=5, D_ha_validada=0.0,
                            n_pixeles_hansen=18, confianza_hansen="media", nicfi=EvidenciaNICFI(disponible=False)),
            PredioNacional("VER-CAFE-0.8ha-HUECO-ESCALA-CON-NICFI", estado="Veracruz", municipio="Teocelo", superficie_ha=0.8,
                            B_treecover=40, B_ndvi=0.75, B_ndvi_std=0.05, B_textura=100, B_altitud=1100,
                            JRC_pct=12, ESRI_crops_pct=20, T_perdida=5, D_ha_validada=0.0,
                            n_pixeles_hansen=18, confianza_hansen="media",
                            nicfi=EvidenciaNICFI(disponible=True, ndvi_antiguo=0.71, ndvi_reciente=0.69, caida_pct=2.8)),
            PredioNacional("MIC-AGUACATE-2ha-ROJO", estado="Michoacán", municipio="Uruapan", superficie_ha=2.0,
                            B_treecover=35, B_ndvi=0.80, B_ndvi_std=0.06, B_textura=30, B_altitud=1800, B_bfast_anio=2021,
                            JRC_pct=60, ESRI_crops_pct=80, T_perdida=22, D_ha_validada=1.2, D_ha_descartada_aislada=0.09,
                            n_pixeles_hansen=220, confianza_hansen="alta"),
            PredioNacional("MIC-AGUACATE-3ha-DEFOR-PRE2020-YA-NO-ROJO", estado="Michoacán", municipio="Uruapan", superficie_ha=3.0,
                            B_treecover=45, B_ndvi=0.82, B_ndvi_std=0.05, B_textura=35, B_altitud=1900, B_bfast_anio=2016,
                            JRC_pct=55, ESRI_crops_pct=70, T_perdida=8, D_ha_validada=0.0, D_ha_historica_pre2020=1.8,
                            n_pixeles_hansen=330, confianza_hansen="alta"),
        ]
        res = motor_nacional_masivo(demos)
        for r in res:
            print(f"{r.id_predio} | {r.superficie_ha}ha | {r.estado} | B={r.B_treecover:.0f}% JRC={r.JRC_pct:.0f}% "
                  f"confianza_hansen={r.confianza_hansen} ({r.n_pixeles_hansen}px) NICFI={r.nicfi.disponible}")
            print(f"  => {r.color} | {r.dictamen_corto} | Sistema:{r.sistema} | Checksum:{r.checksum_integridad}")
            print(f"  Evidencia: {r.texto_evidencia[:220]}...\n")

    if args.pruebita:
        print(f"=== PRUEBITA NACIONAL {args.pruebita} PREDIOS modo={args.modo} NICFI={args.con_nicfi} ===")
        import geopandas as gpd
        from shapely.geometry import Point
        rows = []
        for i in range(args.pruebita):
            rows.append({"id_predio": f"RAN-NAC-{i:07d}",
                         "municipio": random.choice(["Teocelo", "Uruapan", "Palenque", "Misantla"]),
                         "estado": random.choice(["Veracruz", "Michoacán", "Chiapas", "Jalisco"]),
                         "geometry": Point(-102 + random.random() * 6, 17 + random.random() * 5).buffer(0.001)})
        gdf = gpd.GeoDataFrame(rows, crs="EPSG:4326")
        chk = checkpoint_load()
        total = 0
        for estado, gdf_est in gdf.groupby('estado'):
            if args.max_estados and total >= args.max_estados * 100:
                break
            n = procesar_estado_poderoso(gdf_est, estado, modo=args.modo, con_nicfi=args.con_nicfi)
            total += n
            chk["total_predios"] += n
            if estado not in chk["estados_completados"]:
                chk["estados_completados"].append(estado)
            checkpoint_save(chk)
        generar_productos_icm()

    if args.ran:
        print(f"=== MODO NACIONAL REAL {args.ran} ===")
        import geopandas as gpd
        chk = checkpoint_load()
        print(f"Checkpoint actual: {chk['estados_completados']} total {chk['total_predios']}")
        gdf_total = gpd.read_file(args.ran)
        gdf_total['geometry'] = gdf_total['geometry'].make_valid()

        if args.resolver_traslapes is not None:
            campo = None if args.resolver_traslapes == "__auto__" else args.resolver_traslapes
            gdf_total = eliminar_traslapes_jerarquico(gdf_total, campo_prioridad=campo)
        else:
            log("--resolver-traslapes no activado: si el RAN tiene predios solapados, "
                "las hectáreas se pueden estar contando dos veces.", "WARN")

        estados = gdf_total['estado'].unique().tolist() if 'estado' in gdf_total.columns else ["Nacional"]
        if args.max_estados:
            estados = estados[:args.max_estados]
        for estado in estados:
            if estado in chk["estados_completados"]:
                print(f"{estado} ya completado, saltando")
                continue
            gdf_est = gdf_total[gdf_total['estado'] == estado] if 'estado' in gdf_total.columns else gdf_total
            print(f"Procesando {estado}: {len(gdf_est)} predios")
            n = procesar_estado_poderoso(gdf_est, estado, modo=args.modo, con_nicfi=args.con_nicfi)
            chk["estados_completados"].append(estado)
            chk["total_predios"] += n
            checkpoint_save(chk)
            print(f"Checkpoint guardado: {estado} completado")
        generar_productos_icm()


if __name__ == "__main__":
    main()

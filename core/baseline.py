class AdaptadorMLPoderoso:
    def __init__(self, modelo_path=None):
        self.modelo_sistema = None
        self.modelo_riesgo = None
        self.scaler = None
        self.conectado = False
        self.features_nombre = ['B_treecover', 'B_ndvi', 'B_ndvi_std', 'B_textura', 'B_altitud',
                                 'JRC_pct', 'ESRI_crops', 'ndvi_2020', 'ndvi_2024', 'D_pendiente']
        if modelo_path and os.path.exists(modelo_path):
            try:
                import joblib
                data = joblib.load(modelo_path)
                self.modelo_sistema = data.get('modelo_sistema')
                self.modelo_riesgo = data.get('modelo_riesgo')
                self.scaler = data.get('scaler')
                self.conectado = True
                print(f"[ML] Modelos cargados: {modelo_path} - conectado {self.conectado}")
            except Exception as e:
                print(f"[ML] No se pudo cargar modelo: {e} - usando reglas")
        else:
            print("[ML] Sin modelo entrenado - usando reglas físicas (fallback)")

    def extraer_features(self, predio):
        return [
            predio.B_treecover, predio.B_ndvi, predio.B_ndvi_std, predio.B_textura,
            predio.B_altitud, predio.JRC_pct, 100 - predio.ESRI_crops_pct,
            predio.ndvi_2020, predio.ndvi_2024, predio.D_pendiente
        ]

    def predecir_sistema_ml(self, predio):
        if not self.conectado or self.modelo_sistema is None:
            return None, 0.0
        try:
            import numpy as np
            X = np.array(self.extraer_features(predio)).reshape(1, -1)
            if self.scaler:
                X = self.scaler.transform(X)
            pred = self.modelo_sistema.predict(X)[0]
            proba = max(self.modelo_sistema.predict_proba(X)[0]) if hasattr(self.modelo_sistema, 'predict_proba') else 0.85
            return str(pred), float(proba)
        except Exception as e:
            print(f"[ML] Error pred sistema: {e}")
            return None, 0.0

    def entrenar_desde_campo(self, csv_campo_path, output_modelo_path):
        """
        ADVERTENCIA: si csv_campo_path viene de extraer_BTD_de_SIAP_integrado
        en modo simulado, el accuracy es circular, no evidencia de que el
        modelo generalice a predios reales.
        """
        import pandas as pd
        from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
        from sklearn.preprocessing import StandardScaler
        from sklearn.model_selection import train_test_split
        import joblib

        df = pd.read_csv(csv_campo_path)
        print(f"[ML] Entrenando con {len(df)} muestras: {csv_campo_path}")

        X = df[self.features_nombre].values
        y_sistema = df['sistema'].values
        y_riesgo = df['riesgo'].values if 'riesgo' in df.columns else df['D_ha_validada'].values

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        X_train, X_test, y_sis_train, y_sis_test = train_test_split(X_scaled, y_sistema, test_size=0.2, random_state=42)

        modelo_sis = RandomForestClassifier(n_estimators=300, max_depth=15, min_samples_leaf=5, random_state=42, n_jobs=-1)
        modelo_sis.fit(X_train, y_sis_train)
        acc = modelo_sis.score(X_test, y_sis_test)
        print(f"[ML] Sistema - Accuracy test: {acc:.3f} (ver advertencia sobre datos circulares en el docstring)")

        modelo_riesgo = GradientBoostingRegressor(n_estimators=200, max_depth=6, random_state=42)
        modelo_riesgo.fit(X_scaled, y_riesgo)

        joblib.dump({'modelo_sistema': modelo_sis, 'modelo_riesgo': modelo_riesgo,
                     'scaler': scaler, 'features': self.features_nombre}, output_modelo_path)
        print(f"[ML] Modelo guardado: {output_modelo_path}")
        return output_modelo_path


ADAPTADOR_ML_GLOBAL = None


def init_ml_global(modelo_path=None):
    global ADAPTADOR_ML_GLOBAL
    ADAPTADOR_ML_GLOBAL = AdaptadorMLPoderoso(modelo_path=modelo_path)
    return ADAPTADOR_ML_GLOBAL


# ============== LANDTRENDR BÁSICO ==============
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

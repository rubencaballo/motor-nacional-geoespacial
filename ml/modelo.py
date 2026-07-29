#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AdaptadorMLPoderoso: unica definicion (antes habia dos, con nombres
identicos y contenido incompatible, en core/baseline.py y ml/modelo.py)."""
import os
from typing import Optional, Tuple

from modelos.predio import PredioNacional


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



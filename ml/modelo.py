from typing import List
import numpy as np

class AdaptadorMLPoderoso:
    def __init__(self):
        self.scaler = None
        self.model = None

    def extraer_features(self, predio):
        # Lee de EvidenciaNICFI real
        ev = getattr(predio, 'evidencia', None)
        if ev:
            return [
                float(getattr(ev, 'ndvi_antiguo', 0)),
                float(getattr(ev, 'ndvi_reciente', 0)),
                float(getattr(ev, 'caida_pct', 0)),
                float(getattr(ev, 'hansen_loss_year', 0)),
                float(getattr(predio, 'area_ha', 0)),
            ] + [0.0]*5
        return [0.0]*10

    def predecir_sistema(self, predio):
        try:
            X = np.array(self.extraer_features(predio)).reshape(1, -1)
            if self.scaler:
                X = self.scaler.transform(X)
            # Aquí va tu modelo entrenado
            return "milpa_rozatumba" if X[0][2] > 0.2 else "acahual"
        except Exception as e:
            print(f"Error MLP: {e}")
            return "desconocido"

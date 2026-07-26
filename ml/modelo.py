from typing import List
import numpy as np

class AdaptadorMLPoderoso:
    def __init__(self):
        self.scaler = None
        self.model = None

    def extraer_features(self, predio):
        # TODO: tu logica de features Hansen + NDVI
        return [0.0]*10

    def predecir_sistema(self, predio):
        try:
            X = np.array(self.extraer_features(predio)).reshape(1, -1)
            if self.scaler:
                X = self.scaler.transform(X)
            return "sistema_predicho"
        except Exception as e:
            print(f"Error MLP: {e}")
            return "desconocido"

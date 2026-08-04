# motor-nacional-geoespacial
Motor Nacional Geoespacial EUDR - Fórmula E=(B x T)-D - SADER ICM
# Motor Nacional Geoespacial — EUDR

Motor de evaluación de riesgo de deforestación por parcela, alineado al Reglamento UE 2023/1115 (EUDR), desarrollado como piloto propio en Teocelo y Xico, Veracruz. Proyecto vinculado a SADER/AGRICULTURA e Iniciativa Climática de México (ICM).

## ¿Qué hace?

Para cada parcela, el motor evalúa si existe evidencia de pérdida de cobertura forestal posterior al corte temporal del EUDR (31 de diciembre de 2020) y emite un dictamen:

- 🟢 **Verde** — Apto para exportación
- 🟡 **Amarillo** — Requiere evidencia adicional (foto de campo o imagen Planet NICFI)
- 🔴 **Rojo** — No apto para exportación

El criterio central se resume en la fórmula `E = (B × T) − D`:

- **B** — cobertura base de bosque/sombra (Hansen Global Forest Change, JRC Forest Cover, ESRI Land Cover)
- **T** — umbral/tipo de sistema productivo evaluado
- **D** — pérdida de cobertura detectada posterior al corte EUDR

La confianza del dictamen depende del número de píxeles Hansen válidos dentro de la parcela (`MIN_PIXELES_CONFIANZA_ALTA = 30`, `MIN_PIXELES_CONFIANZA_MEDIA = 10`, definidos en `config.py`). En casos de confianza baja o media, el motor contrasta contra JRC, ESRI y, si está disponible, series NDVI de Planet NICFI (4.77m) antes de resolver.

Cada resultado incluye evidencia textual trazable y dos identificadores de integridad:

- `checksum_integridad` — incluye la versión del motor; sirve para auditoría en producción (saber con qué versión se generó un registro).
- `hash_resultado` — no incluye versión; sirve para verificar que dos versiones del motor producen el mismo resultado lógico (pruebas de regresión).

## Sistemas productivos cubiertos

Actualmente implementados: **café bajo sombra** y **aguacate**. La arquitectura está pensada para añadir sistemas adicionales (palma, ganadería) como nuevos módulos, sin modificar el núcleo de decisión.

## Arquitectura
config.py Constantes: corte EUDR, umbrales de confianza, VERSION
core/
motor.py Orquestador: motor_nacional_masivo(), CLI (--demo)
geo.py Utilidades geométricas
gee_backend.py Integración con Google Earth Engine (Hansen, JRC, Sentinel-2, LandTrendr)
db.py Persistencia y checkpoints
modelos/
predio.py Dataclass PredioNacional (resultado por parcela + auditoría)
evidencia.py EvidenciaNICFI (contraste con imágenes Planet NICFI)
ml/
modelo.py AdaptadorMLPoderoso (interfaz para modelo entrenado; hoy corre en fallback de reglas físicas)
## Cómo correrlo

Desde la raíz del repositorio:

```bash
python3 -m core.motor --demo
```

Esto corre 5 casos de prueba con valores sintéticos (no consulta Earth Engine en vivo) y muestra, por cada uno, el dictamen, la evidencia y ambos identificadores de integridad.

## Estado actual — honesto, sin inflar

- ✅ Lógica de decisión modular, sin clases duplicadas, cubierta por pruebas de regresión (`hash_resultado`) entre la versión monolítica original y la modular.
- ✅ `core/gee_backend.py` contiene llamadas reales a la API de Earth Engine (Hansen GFC, JRC, Landsat, LandTrendr) — no es un mock.
- ⚠️ El `--demo` usa datos sintéticos, no una consulta en vivo a Earth Engine; el pipeline de extracción real aún no se ha verificado end-to-end.
- ⚠️ El módulo de ML (`AdaptadorMLPoderoso`) no tiene un modelo entrenado con datos reales; hoy opera con reglas físicas de respaldo.
- ⚠️ La arquitectura está diseñada para escalar a volúmenes grandes de predios, pero no se ha probado en producción a esa escala.

## Contexto

Desarrollado por Rubén Viccon Morales como parte de IRD CLOUD Engine, en el marco de un piloto de trazabilidad EUDR con productores de café de Teocelo, Veracruz (ver también GEODICTUM, sello de trazabilidad entregado a productores en 2026).

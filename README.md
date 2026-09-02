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

    config.py           Constantes: corte EUDR, umbrales de confianza, VERSION
    core/
      motor.py           Orquestador: motor_nacional_masivo(), CLI (--demo)
      geo.py             Utilidades geométricas
      gee_backend.py     Integración con Google Earth Engine (Hansen, JRC, Sentinel-2, LandTrendr)
      db.py              Persistencia y checkpoints
    modelos/
      predio.py          Dataclass PredioNacional (resultado por parcela + auditoría)
      evidencia.py       EvidenciaNICFI (contraste con imágenes Planet NICFI)
    ml/
      modelo.py          AdaptadorMLPoderoso (interfaz para modelo entrenado; hoy corre en fallback de reglas físicas)

## Cómo correrlo

Desde la raíz del repositorio:

    python3 -m core.motor --demo

Esto corre 5 casos de prueba con valores sintéticos (no consulta Earth Engine en vivo) y muestra, por cada uno, el dictamen, la evidencia y ambos identificadores de integridad.

## Estado actual — honesto, sin inflar

- ✅ Lógica de decisión modular, sin clases duplicadas, cubierta por pruebas de regresión (`hash_resultado`) entre la versión monolítica original y la modular.
- ✅ `core/gee_backend.py` contiene llamadas reales a la API de Earth Engine (Hansen GFC, JRC, Landsat, LandTrendr) — no es un mock.
- ⚠️ El `--demo` usa datos sintéticos, no una consulta en vivo a Earth Engine; el pipeline de extracción real aún no se ha verificado end-to-end.
- ⚠️ El módulo de ML (`AdaptadorMLPoderoso`) no tiene un modelo entrenado con datos reales; hoy opera con reglas físicas de respaldo.
- ⚠️ La arquitectura está diseñada para escalar a volúmenes grandes de predios, pero no se ha probado en producción a esa escala.

## Plataforma de monitoreo de ANPs (Salamandra / IRD Cloud Engine)

Este mismo `core/` también contiene un segundo pipeline, independiente del motor
EUDR de arriba: monitoreo 3D de Áreas Naturales Protegidas (terreno, hidrología D8,
carbono/CO2e almacenado y liberado, deforestación Hansen, validación de incendios),
piloteado en Cofre de Perote, Veracruz.

- `core/analizar_sitio.py` — orquestador (terreno + carbono + validación hidrológica)
- `core/geomatica.py` — DEM/hidrología D8 y el mapa 3D base
- `core/carbono.py` — biomasa/CO2e almacenado por zona (ESA CCI + GEDI L4A)
- `core/deforestacion.py` — pérdida Hansen por año + mapa 3D
- `core/carbono_perdida.py` — cruce carbono × pérdida, balance, mapa 3D de CO2e liberado
- `core/validacion_incendios.py` — validación de causa (incendio vs. tala) por dNBR
- `core/agua_superficial.py` — agua superficial visible por año (JRC Global Surface Water), experimental

**Los comandos exactos para correr este pipeline completo, en el orden correcto,
están en [`RUNBOOK.md`](RUNBOOK.md)** — no en este README. Para correrlo de un
tirón sobre un sitio nuevo, edita las variables al principio de
[`run_pipeline.sh`](run_pipeline.sh) y corre `bash run_pipeline.sh`.

## Dictamen integrado CO2 + Biodiversidad (V19.4 Híbrido Corregido)

Un tercer pipeline, complementario al de arriba: combina el CO2e por zona (mismo
dataset ESA CCI + GEDI del punto anterior) con un inventario de biodiversidad por
zona (GBIF) en un solo dictamen técnico-científico para CONANP/SEDEMA. Pensado
para cualquier ANP o predio, sea o no sitio Ramsar.

- `core/salamandra_biodiversidad.py` — `validate_taxonomic_class()`: valida la
  columna `CLASE` de cada registro de biodiversidad contra su `FAMILIA` (tabla de
  autoridad de 240 familias) y corrige a `CLASE_VALIDADA`. Existe porque los
  scripts históricos de descarga GBIF (V16/V16.5) escribían en `CLASE` la clase
  taxonómica *pedida* a la API, no la que GBIF *realmente* devolvía por registro
  — un bug que corrompió el `CLASE` del CSV histórico de Texolo durante 4
  versiones seguidas (hasta 34% de discrepancia visto en corridas reales).
  `V19_4_HIBRIDO_CORREGIDO.py` ya trae `CLASE` correcta desde el origen (toma
  `o.get("class")` de la respuesta de GBIF, nunca la clase pedida), pero de
  todas formas pasa por esta validación como segunda red de seguridad.
- `core/V19_4_HIBRIDO_CORREGIDO.py` — orquestador de este pipeline: lee (o
  descarga automáticamente de GBIF) biodiversidad por zona, cruza con el CSV de
  carbono, y genera el dictamen (PDF + TXT), 4 mapas HD y los CSV de resumen. Si
  no se pasa `--csv-biodiversidad`, descarga de GBIF en vivo (todas las clases
  taxonómicas, anillo exclusivo núcleo/500m/1000m sin traslape entre zonas, con
  reintento automático ante caídas 503 del servidor de GBIF).

**Importante sobre el Sitio Ramsar 1601:** el texto de "sitio Ramsar confirmado"
en el dictamen no sale por default — requiere pasar explícitamente
`--confirmar-ramsar-1601-texolo`, y esa bandera solo debe usarse para Cascadas de
Texolo. Ver el porqué (un bug real que sí llegó a generarse) y los comandos
exactos en [`RUNBOOK.md`](RUNBOOK.md).

## Contexto

Desarrollado por Rubén Viccon Morales como parte de IRD CLOUD Engine, en el marco de un piloto de trazabilidad EUDR con productores de café de Teocelo, Veracruz (ver también GEODICTUM, sello de trazabilidad entregado a productores en 2026).

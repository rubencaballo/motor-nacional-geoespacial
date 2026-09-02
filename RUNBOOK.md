# RUNBOOK — cómo correr el pipeline completo de un sitio

Este archivo existe por una razón concreta: los comandos exactos (banderas, orden,
qué CSV alimenta a qué mapa) no deben depender de que alguien —yo incluido— se
acuerde de memoria. Cada comando de este documento está copiado directamente de la
definición real de `argparse` en el código (no reconstruido de memoria), verificado
el 2026-08-31 (sección de biodiversidad verificada el 2026-09-02). Si algún día un
comando de aquí deja de coincidir con el código, el código manda — pero avisa para
corregir este archivo también.

Para un sitio nuevo (otro ANP, un predio, el sitio Ramsar que sigue después de
Cofre de Perote), solo cambian las variables de la sección "Parámetros del sitio";
el orden y las banderas son los mismos.

## Orden y dependencias (por qué importa la secuencia)

1. **Terreno + Carbono** (`core.analizar_sitio`) — siempre primero. Genera el DEM/malla
   3D y, si se da `--proyecto-gee`, el CO2e almacenado por zona.
2. **Deforestación, CSVs oficiales** (`core.deforestacion`, SIN `--mapa-3d`) — genera
   `deforestacion_resumen_sin_traslape_*.csv` y `deforestacion_historial_anual_*.csv`.
   Este paso **no genera el mapa 3D**, pero es obligatorio antes de los pasos 4-6
   porque ellos leen estos dos CSV.
3. **Deforestación, mapa 3D** (`core.deforestacion --mapa-3d`) — el mapa visual con
   las manchas por año. Si el paso 2 ya se corrió, este mapa cita automáticamente el
   total oficial en el subtítulo (la nota de 723 vs 771 ha, ver más abajo).
4. **Cruce carbono × pérdida + balance** (`core.carbono_perdida`, sin `--mapa-3d-perdida`)
   — junta el CO2e almacenado (paso 1) con la pérdida oficial (paso 2) y genera
   `co2e_asociado_perdida_*.csv` (+ `balance_stock_vs_perdida_*.csv` con `--balance`).
5. **Mapa 3D de CO2e liberado** (`core.carbono_perdida --mapa-3d-perdida`) — el mapa
   "espejo" del de carbono almacenado, pero con lo ya liberado. Necesita que el
   paso 1 y el paso 4 ya hayan corrido (lee sus CSV de `--carpeta-salida`).
6. **Validación histórica de incendios** (`core.validacion_incendios --historial --mapa-3d`)
   — opcional pero recomendado; si se le da `--historial-csv-existente` apuntando al
   CSV del paso 2, fusiona causa probable con la pérdida ya calculada.

Los pasos 2 y 3 son **el mismo script, corrido dos veces con banderas distintas** —
es el detalle más fácil de olvidar de todo el pipeline: `--mapa-3d` por sí solo
**no** genera `deforestacion_resumen_sin_traslape_*.csv` (verificado leyendo
`generar_mapa_3d_deforestacion()`, que no llama a `procesar_sitio_real()`). Sin el
paso 2 corrido antes, el paso 4 no tiene con qué cruzar y los pasos 5-6 no pueden
citar el total oficial.

El **dictamen integrado CO2 + Biodiversidad** (`V19_4_HIBRIDO_CORREGIDO.py`) es un
pipeline aparte, no uno de estos 6 pasos — ver su propia sección más abajo.

## Requisito previo: Earth Engine

Todos los pasos que consultan datos reales (todo excepto `--demo`) necesitan
`earthengine-api` instalado y una cuenta de Earth Engine ya autenticada una vez
(`earthengine authenticate`, o el flujo equivalente ya usado). El proyecto de
Google Cloud se pasa con `--proyecto-gee` en cada comando; si no se pasa, se usa
`ee.Initialize()` sin proyecto explícito.

## Paso 1 — Terreno + Carbono (+ validación hidrológica opcional)

```bash
python3 -m core.analizar_sitio \
  --geojson Cofre_de_Perote.geojson \
  --id-proyecto Cofre_de_Perote \
  --proyecto-gee ee-rvicconmorales \
  --zonas 0,500,1000 \
  --carpeta-salida ~/resultados_cofre_de_perote
```

Alternativa: `--anp "Pico de Orizaba"` en vez de `--geojson` si el sitio está en el
shapefile nacional de CONANP (usa `core.anp_lookup` para resolverlo). Para incluir
también la validación de cauces D8 contra INEGI, agrega
`--shapefile-inegi ruta/a/red_hidrografica.gpkg`.

`--zonas` es opcional; si se omite, usa el default de `config.py`
(`ZONAS_ANALISIS_M = [0, 500, 1000]`).

**Validación hidrológica (D8 vs INEGI) — qué dice el % y cómo leerlo.** El CSV
(`validacion_hidrologica_*.csv`) y el mapa 3D (con `core.validacion_hidrologica
--mapa-3d`) ahora incluyen, además de las zonas nucleo/buffer_500m/buffer_1000m
(ya exclusivas por diseño -- cada píxel de cauce D8 pertenece a una sola zona,
sin traslape), una fila/línea **TOTAL calculada sobre la unión real de los
puntos**, nunca como promedio simple de los % por zona (ese promedio simple sí
sería un error: con los datos reales de Cofre de Perote da 50.9%, cuando el
total correcto, ponderado, es 48.9%). El mapa 3D también agrega en el propio
subtítulo el aviso honesto de por qué el % puede verse bajo: INEGI (escala
1:50,000) no digitaliza cauces efímeros/pequeños que el D8 (SRTM 30m, más fino)
sí detecta -- un % de 40-55% no es necesariamente un error del modelo.

**Comando exacto para correr SOLO la validación hidrológica** (en vez de pasar
por `analizar_sitio.py`) -- así es como Ruben lo corre en la práctica, CSV +
mapa 3D en un solo comando, reusando el SRTM ya descargado:

```bash
python3 -m core.validacion_hidrologica \
  --geojson Cofre_de_Perote.geojson \
  --shapefile-inegi conjunto_de_datos/cnit50k.gpkg \
  --id-proyecto Cofre_de_Perote \
  --mapa-3d
```

Sin `--proyecto-gee` porque este módulo no toca Earth Engine (solo SRTM local +
el shapefile de INEGI). Sin `--zonas`/`--percentil-cauce`/`--carpeta-salida`
usa los defaults de `config.py` (`ZONAS_ANALISIS_M=[0,500,1000]`,
`PERCENTIL_CAUCE_HIDROLOGIA=98` -- de ahí el `_p98` en el nombre del CSV/HTML
de salida). `--layer` default ya es `corriente_ag_l` (la capa lineal de
cauces del CNIT50k), así que tampoco hace falta pasarlo si el `.gpkg` usa ese
nombre de capa.

## Paso 2 — Deforestación: CSVs oficiales (SIN `--mapa-3d`)

```bash
python3 -m core.deforestacion \
  --geojson Cofre_de_Perote.geojson \
  --id-proyecto Cofre_de_Perote \
  --proyecto-gee ee-rvicconmorales \
  --zonas 0,500,1000 \
  --anio-inicio 2010 --anio-fin 2025 \
  --carpeta-salida ~/resultados_cofre_de_perote
```

`--anio-inicio` default: `2010` (`config.DEFORESTACION_ANIO_INICIO_DEFAULT`).
`--anio-fin` default: año actual − 1. Genera `deforestacion_resumen_*.csv`,
`deforestacion_historial_anual_*.csv` y `deforestacion_resumen_sin_traslape_*.csv`.

## Paso 3 — Deforestación: mapa 3D

```bash
python3 -m core.deforestacion \
  --geojson Cofre_de_Perote.geojson \
  --id-proyecto Cofre_de_Perote \
  --proyecto-gee ee-rvicconmorales \
  --zonas 0,500,1000 \
  --anio-inicio 2010 --anio-fin 2025 \
  --carpeta-salida ~/resultados_cofre_de_perote \
  --mapa-3d
```

Mismos parámetros que el paso 2 más `--mapa-3d`. Genera
`cofre_de_perote_3d_deforestacion.html` y `deforestacion_desglose_anual_visual_*.csv`.

## Paso 4 — Cruce carbono × pérdida + balance

```bash
python3 -m core.carbono_perdida \
  --carbono-csv ~/resultados_cofre_de_perote/resumen_terreno_y_carbono_cofre_de_perote.csv \
  --perdida-sin-traslape-csv ~/resultados_cofre_de_perote/deforestacion_resumen_sin_traslape_cofre_de_perote.csv \
  --id-proyecto Cofre_de_Perote \
  --carpeta-salida ~/resultados_cofre_de_perote \
  --anios-incendio 2025 \
  --balance
```

`--anios-incendio` solo si hay años con incendio confirmado (ver paso 6); esos años
usan `--factor-combustion` (default `0.45`) en vez de remoción completa. Si no hay
ningún incendio confirmado, omite esa bandera por completo. `--balance` es opcional
pero recomendado (genera `balance_stock_vs_perdida_*.csv` como archivo aparte; el
mapa del paso 5 calcula su propio balance internamente de todas formas).

## Paso 5 — Mapa 3D de CO2e liberado (espejo del de carbono almacenado)

```bash
python3 -m core.carbono_perdida \
  --mapa-3d-perdida \
  --geojson Cofre_de_Perote.geojson \
  --id-proyecto Cofre_de_Perote \
  --zonas 0,500,1000 \
  --carpeta-salida ~/resultados_cofre_de_perote
```

Requiere que los pasos 1 y 4 ya hayan corrido (lee `resumen_terreno_y_carbono_*.csv`
y `co2e_asociado_perdida_*.csv` de `--carpeta-salida`). Genera
`cofre_de_perote_3d_perdida_carbono.html`.

## Paso 6 — Validación histórica de incendios

```bash
python3 -m core.validacion_incendios \
  --historial \
  --geojson Cofre_de_Perote.geojson \
  --id-proyecto Cofre_de_Perote \
  --anio-inicio 2010 --anio-fin 2025 \
  --zonas 0,500,1000 \
  --eventos-confirmados 2025:2025-04-17 \
  --historial-csv-existente ~/resultados_cofre_de_perote/deforestacion_historial_anual_cofre_de_perote.csv \
  --proyecto-gee ee-rvicconmorales \
  --carpeta-salida ~/resultados_cofre_de_perote \
  --mapa-3d
```

`--eventos-confirmados` solo si hay fecha real conocida de algún incendio (formato
`año:YYYY-MM-DD`, varios separados por coma); los años sin fecha conocida se
evalúan con el screening genérico (dNBR + `VALIDACION_INCENDIO_UMBRAL_DNBR_QUEMADO
= 0.10`). Si no hay ningún evento confirmado, omite esa bandera.

Para validar un solo evento puntual en vez del historial completo, existe también
el modo sin `--historial`: `python3 -m core.validacion_incendios --geojson ... --id-proyecto ... --fecha-evento YYYY-MM-DD --anio-hansen 2025 ...` (ver `python3 -m core.validacion_incendios --help` para el resto de banderas; ese modo aún no tiene la nota de 723-vs-771 aplicada, solo el modo `--historial`).

## Atajo: `run_pipeline.sh`

`run_pipeline.sh` (raíz del repo) corre los 6 pasos de arriba en orden con un solo
comando. Edita las variables al principio del archivo (o expórtalas como variables
de entorno antes de llamarlo) y corre:

```bash
bash run_pipeline.sh
```

Está pensado para copiarlo y ajustar solo las variables de arriba cuando se analice
un sitio nuevo — el resto del script no debería necesitar tocarse.

## La nota de 723 vs 771 hectáreas — por qué existen dos cifras

Dos métodos, ambos válidos, para lo mismo (pérdida de bosque Hansen 2010-2025 en
Cofre de Perote), y por eso dan números distintos:

- **771.1 ha ("oficial")** — `deforestacion_resumen_sin_traslape_*.csv` (paso 2),
  consulta vectorial directa a Earth Engine sobre el polígono exacto. Es la cifra
  que alimenta el mapa de CO2e liberado (paso 5).
- **723.0 ha ("visual")** — `deforestacion_desglose_anual_visual_*.csv` (paso 3) y
  el mapa histórico de incendios (paso 6): ambos remuestrean el raster de Hansen a
  la malla local del terreno (rasterio, `Resampling.nearest`) para poder dibujarlo
  encima del DEM 3D; ese remuestreo pierde/gana algunos píxeles de borde.

Ya está documentado en el propio docstring de `generar_desglose_anual_visual()`
como comportamiento esperado, no un error de doble conteo. Los mapas de
deforestación (paso 3) y de incendio histórico (paso 6) ya citan automáticamente
el total oficial en su subtítulo cuando el CSV del paso 2 existe, así que el
usuario final ve ambas cifras explicadas en el mismo mapa.

## Para un sitio nuevo (ej. el sitio Ramsar)

Lo único que cambia entre sitios son las variables de la sección "Parámetros del
sitio" en `run_pipeline.sh` (o las banderas `--geojson`/`--id-proyecto`/
`--carpeta-salida` en cada comando de arriba): ningún archivo de `core/` tiene
valores de Cofre de Perote incrustados en la lógica (verificado con
`grep -ri "perote\|cofre" core/`, solo aparece en comentarios/docstrings). La zona
UTM y la descarga del SRTM se detectan automáticamente del centroide del polígono.

Caveat real, no de código sino de método: Hansen (deforestación) y ESA CCI/GEDI
(biomasa aérea) están calibrados para bosque. Un sitio Ramsar (humedal, laguna,
vegetación no leñosa) puede dar señal casi nula o poco significativa en estos
mismos datasets, porque ahí el carbono relevante suele estar en el suelo/sedimento
("blue carbon"), no en biomasa leñosa aérea. El pipeline correrá sin errores sobre
ese polígono, pero los números de carbono/deforestación pueden no ser
representativos del sistema real hasta no revisar si aplica otro dataset (o al
menos documentar la limitación en el reporte de ese sitio).

## Dictamen integrado CO2 + Biodiversidad (V19.4 Híbrido Corregido)

Pipeline aparte de los 6 pasos de arriba (no pasa por `run_pipeline.sh`). Combina
el CO2e por zona (mismo CSV de carbono del paso 1) con un inventario de
biodiversidad por zona (GBIF) en un solo dictamen PDF+TXT para CONANP/SEDEMA. Se
corre directo, no como `python3 -m core.X`, porque importa a su módulo hermano
`salamandra_biodiversidad.py` por la ruta del propio archivo (`sys.path.insert`),
no como paquete:

```bash
python3 core/V19_4_HIBRIDO_CORREGIDO.py \
  --geojson Cofre_de_Perote.geojson \
  --csv-carbono ~/resultados_cofre_de_perote/resumen_terreno_y_carbono_cofre_de_perote.csv
```

Sin `--csv-biodiversidad`, descarga de GBIF en vivo: todas las clases
taxonómicas, anillo exclusivo por zona (núcleo/500m/1000m sin traslape entre
ellas -- ver `zonas_anillo_exclusivo()`), hasta 5,000 regs/zona por default
(ajustable con `--gbif-max-regs-por-zona`), y guarda el CSV descargado en la
carpeta de salida. Si ya tienes un CSV de biodiversidad de una corrida anterior
(o un CSV histórico), pásalo con `--csv-biodiversidad ruta.csv` para saltarte la
descarga.

El geojson pasado en `--geojson` **debe tener las 3 zonas** (núcleo + buffer 500m
+ buffer 1000m, como features separadas con `NIVEL`/`tipo` conteniendo
"Nucleo"/"DESGLOSE_29", "500" y "1000" respectivamente) para que la descarga de
GBIF consulte los 3 anillos -- un geojson de un solo polígono (por ejemplo el
límite oficial WDPA/CONANP tal cual, sin buffers) solo le da geometría al
núcleo, y las zonas de buffer se quedan en 0 registros sin avisar con error
(se ve en el log como `[GBIF] BUFFER_500m: geometría vacía, se omite.` sí avisa,
pero conviene revisar el `V19_4_RESUMEN_BIO.csv` de salida para confirmar que
los 3 niveles tienen registros, no solo el núcleo). Si no tienes un geojson de 3
zonas, se construye buffereando el polígono oficial en la proyección UTM
correspondiente (no en Web Mercator/EPSG:3857 -- infla el área ~13% a esta
latitud) y verificando el área resultante contra el CSV de carbono.

Reintenta automáticamente (hasta 6 veces, espera creciente 10s/20s/.../60s) ante
errores 503 "Backend fetch failed" de GBIF -- son caídas cortas del servidor de
GBIF, no errores del polígono ni de la consulta. Si tras 6 reintentos una zona
sigue en 0 registros, GBIF está caído de verdad: espera unos minutos y vuelve a
correr el mismo comando.

**Sitio Ramsar 1601 (Cascadas de Texolo) -- bandera obligatoria y a propósito.**
`--sitio-ramsar` NO tiene default `"1601"` (default: `None`) precisamente porque
ese default causó, el 2026-09-02, que una corrida de OTRO sitio (Cofre de
Perote, que no es Ramsar) generara un dictamen que afirmaba falsamente ser el
"Sitio Ramsar 1601... Cascadas de Texolo y su entorno" con criterios (ii);(iv)
confirmados -- se detectó antes de enviarse a CONANP/SEDEMA y se corrigió ese
mismo día. Para Texolo, pasa explícitamente:

```bash
python3 core/V19_4_HIBRIDO_CORREGIDO.py \
  --confirmar-ramsar-1601-texolo \
  --geojson Buffer_500_1000m_Ramsar_1601_Texolo.geojson \
  --csv-carbono ~/resultados_ramsar_1601_texolo/resumen_terreno_y_carbono_ramsar_1601_texolo.csv
```

Para cualquier otro sitio (Ramsar o no), omite ambas banderas -- el dictamen
simplemente no menciona Ramsar (texto: "no se reporta como sitio Ramsar en esta
corrida"). `--sitio-ramsar <ID>` sin `--confirmar-ramsar-1601-texolo` imprime ese
ID tal cual pero marcado como "NO verificado por este script" -- úsalo solo si
de verdad tienes otro sitio Ramsar distinto de Texolo y sabes que sus criterios
no están verificados aquí. Pasar `--sitio-ramsar 1601` junto con
`--confirmar-ramsar-1601-texolo` para un ID distinto de 1601 es un error fatal
a propósito (`sys.exit`), no un valor que el script intente adivinar.

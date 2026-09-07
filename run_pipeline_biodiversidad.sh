#!/usr/bin/env bash
# run_pipeline_biodiversidad.sh -- corre el modulo PUBLICO de biodiversidad
# (core.biodiversidad) y, opcionalmente, el DICTAMEN tecnico CONANP/SEDEMA
# (core/V19_4_HIBRIDO_CORREGIDO.py, version con el fix Ramsar del 2026-09-04),
# ambos con descarga GBIF EN VIVO -- por eso este script NO sirve dentro del
# sandbox de esta sesion (sin salida a internet a api.gbif.org): correlo en tu
# maquina Ubuntu, igual que ayer corriste V19_4_HIBRIDO_CORREGIDO.py a mano.
#
# Es un script APARTE de run_pipeline.sh (a proposito, ver RUNBOOK.md): el
# modulo publico de biodiversidad no depende del CSV de carbono, pero el
# dictamen CONANP si lo necesita -- corre run_pipeline.sh primero si vas a
# generar el dictamen para un sitio nuevo.
#
# Uso:
#   1. Edita la seccion "Parametros del sitio" (o exportalas como variables
#      de entorno antes de llamar al script -- override sin tocar el archivo).
#   2. Corre desde la raiz del repo:  bash run_pipeline_biodiversidad.sh

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# ============================================================
# --- Parametros del sitio -- EDITA ESTO para cada ANP/predio ---
# ============================================================
ID_PROYECTO="${ID_PROYECTO:-Cofre_de_Perote}"
GEOJSON="${GEOJSON:-Cofre_de_Perote.geojson}"
ZONAS="${ZONAS:-0,500,1000}"
CARPETA_SALIDA="${CARPETA_SALIDA:-$HOME/resultados_$(echo "$ID_PROYECTO" | tr '[:upper:]' '[:lower:]')}"

# --- [1/2] Modulo publico Salamandra -- core.biodiversidad ---
CORRER_MODULO_PUBLICO="${CORRER_MODULO_PUBLICO:-1}"
GBIF_MAX_REGS_PUBLICO="${GBIF_MAX_REGS_PUBLICO:-2000}"
SIN_IUCN="${SIN_IUCN:-0}"                 # 1 para saltar consultas IUCN (mas rapido, solo NOM-059)
SIN_MAPA_3D="${SIN_MAPA_3D:-0}"           # 1 para solo CSV + resumen, sin mapa 3D
TABLA_NOM059_LOCAL="${TABLA_NOM059_LOCAL:-}"     # opcional, override/adicion por sitio
TABLA_ENDEMISMO_LOCAL="${TABLA_ENDEMISMO_LOCAL:-}"  # opcional, override/adicion por sitio

# --- [2/2] Dictamen tecnico CONANP/SEDEMA -- V19_4_HIBRIDO_CORREGIDO.py ---
# Apagado por default: es el documento formal, no algo que quieras regenerar
# en automatico sin revisar los parametros de Ramsar cada vez.
CORRER_DICTAMEN_CONANP="${CORRER_DICTAMEN_CONANP:-0}"
CSV_CARBONO="${CSV_CARBONO:-$CARPETA_SALIDA/resumen_terreno_y_carbono_$(echo "$ID_PROYECTO" | tr '[:upper:]' '[:lower:]').csv}"
GBIF_MAX_REGS_DICTAMEN="${GBIF_MAX_REGS_DICTAMEN:-5000}"
# Ramsar -- OBLIGATORIO decidir a mano, no hay default seguro (ver bug del
# 2026-09-02 en RUNBOOK.md). Deja los dos en 0/vacio para un sitio que NO es
# Ramsar (ej. Cofre de Perote): el dictamen simplemente no menciona Ramsar.
CONFIRMAR_RAMSAR_1601_TEXOLO="${CONFIRMAR_RAMSAR_1601_TEXOLO:-0}"  # 1 SOLO para Texolo (sitio 1601)
SITIO_RAMSAR_ID="${SITIO_RAMSAR_ID:-}"     # ID de OTRO sitio Ramsar no verificado por este script (opcional)
DEM_RUTA="${DEM_RUTA:-}"                   # opcional, GeoTIFF real para area 3D via verificador_biodiversidad
N_ITER_RAREFACCION="${N_ITER_RAREFACCION:-200}"  # sube a 999+ para el dictamen final
NDWI="${NDWI:-}"                           # opcional
SPP_ALTO_RIESGO="${SPP_ALTO_RIESGO:-}"     # opcional

id_lower="$(echo "$ID_PROYECTO" | tr '[:upper:]' '[:lower:]')"
mkdir -p "$CARPETA_SALIDA"

echo "=================================================================="
echo " Sitio: $ID_PROYECTO"
echo " Geojson: $GEOJSON"
echo " Carpeta de salida: $CARPETA_SALIDA"
echo "=================================================================="

if [ "$CORRER_MODULO_PUBLICO" = "1" ]; then
    echo
    echo "=== [1/2] Modulo publico Salamandra -- core.biodiversidad (GBIF en vivo) ==="
    ARGS_BIO=(--geojson "$GEOJSON" --id-proyecto "$ID_PROYECTO" --zonas "$ZONAS"
              --carpeta-salida "$CARPETA_SALIDA" --gbif-max-regs-por-zona "$GBIF_MAX_REGS_PUBLICO")
    [ "$SIN_IUCN" = "1" ] && ARGS_BIO+=(--sin-iucn)
    [ "$SIN_MAPA_3D" = "1" ] && ARGS_BIO+=(--sin-mapa-3d)
    [ -n "$TABLA_NOM059_LOCAL" ] && ARGS_BIO+=(--tabla-nom059-local "$TABLA_NOM059_LOCAL")
    [ -n "$TABLA_ENDEMISMO_LOCAL" ] && ARGS_BIO+=(--tabla-endemismo-local "$TABLA_ENDEMISMO_LOCAL")
    python3 -m core.biodiversidad "${ARGS_BIO[@]}"
else
    echo
    echo "=== [1/2] Modulo publico Salamandra: SALTADO (CORRER_MODULO_PUBLICO=0) ==="
fi

if [ "$CORRER_DICTAMEN_CONANP" = "1" ]; then
    if [ ! -f "$CSV_CARBONO" ]; then
        echo "ERROR: no encuentro $CSV_CARBONO -- corre run_pipeline.sh (o core.analizar_sitio) para este sitio primero." >&2
        exit 1
    fi
    if [ "$CONFIRMAR_RAMSAR_1601_TEXOLO" = "1" ] && [ -n "$SITIO_RAMSAR_ID" ] && [ "$SITIO_RAMSAR_ID" != "1601" ]; then
        echo "ERROR: CONFIRMAR_RAMSAR_1601_TEXOLO=1 junto con SITIO_RAMSAR_ID distinto de 1601 no tiene sentido -- son excluyentes." >&2
        exit 1
    fi
    echo
    echo "=== [2/2] Dictamen tecnico CONANP/SEDEMA -- V19_4_HIBRIDO_CORREGIDO.py (GBIF en vivo) ==="
    ARGS_DICT=(--geojson "$GEOJSON" --csv-carbono "$CSV_CARBONO" --out-dir "$CARPETA_SALIDA"
               --gbif-max-regs-por-zona "$GBIF_MAX_REGS_DICTAMEN" --n-iter-rarefaccion "$N_ITER_RAREFACCION")
    if [ "$CONFIRMAR_RAMSAR_1601_TEXOLO" = "1" ]; then
        ARGS_DICT+=(--confirmar-ramsar-1601-texolo)
    elif [ -n "$SITIO_RAMSAR_ID" ]; then
        ARGS_DICT+=(--sitio-ramsar "$SITIO_RAMSAR_ID")
    fi
    [ -n "$DEM_RUTA" ] && ARGS_DICT+=(--dem "$DEM_RUTA")
    [ -n "$NDWI" ] && ARGS_DICT+=(--ndwi "$NDWI")
    [ -n "$SPP_ALTO_RIESGO" ] && ARGS_DICT+=(--spp-alto-riesgo "$SPP_ALTO_RIESGO")
    python3 core/V19_4_HIBRIDO_CORREGIDO.py "${ARGS_DICT[@]}"
else
    echo
    echo "=== [2/2] Dictamen tecnico CONANP/SEDEMA: SALTADO (CORRER_DICTAMEN_CONANP=0) ==="
fi

echo
echo "=================================================================="
echo " LISTO. Resultados en: $CARPETA_SALIDA"
echo "=================================================================="

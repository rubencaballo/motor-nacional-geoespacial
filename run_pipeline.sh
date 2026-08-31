#!/usr/bin/env bash
# run_pipeline.sh -- corre el pipeline COMPLETO de un sitio (ANP o poligono propio)
# con un solo comando, en el orden correcto. Ver RUNBOOK.md para el detalle y el
# porque de cada paso -- cada bandera de aqui esta sacada directamente del argparse
# real de cada modulo (core/analizar_sitio.py, core/deforestacion.py,
# core/carbono_perdida.py, core/validacion_incendios.py), no reconstruida de memoria.
#
# Uso:
#   1. Edita las variables de la seccion "Parametros del sitio" (o expportalas como
#      variables de entorno antes de llamar al script -- override sin tocar el archivo).
#   2. Corre desde la raiz del repo:  bash run_pipeline.sh
#
# Para un sitio nuevo (ej. el sitio Ramsar despues de Cofre de Perote) solo hace
# falta cambiar esta seccion -- el resto del script no deberia necesitar tocarse.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# ============================================================
# --- Parametros del sitio -- EDITA ESTO para cada ANP/predio ---
# ============================================================
ID_PROYECTO="${ID_PROYECTO:-Cofre_de_Perote}"
GEOJSON="${GEOJSON:-Cofre_de_Perote.geojson}"
PROYECTO_GEE="${PROYECTO_GEE:-ee-rvicconmorales}"
ZONAS="${ZONAS:-0,500,1000}"
ANIO_INICIO="${ANIO_INICIO:-2001}"  # 2001 es el minimo real de Hansen (lossyear codigo 1) -- adoptado
                                     # como ventana OFICIAL de Cofre de Perote tras la prueba de
                                     # correlacion perdida-agua (ver RUNBOOK.md, Paso 7c). Antes era
                                     # 2010 -- si trabajas otro sitio y quieres esa ventana mas corta,
                                     # exporta ANIO_INICIO=2010 antes de llamar a este script.
ANIO_FIN="${ANIO_FIN:-2025}"
CARPETA_SALIDA="${CARPETA_SALIDA:-$HOME/resultados_$(echo "$ID_PROYECTO" | tr '[:upper:]' '[:lower:]')}"

# --- Opcionales -- deja vacio ("") para saltar ese paso/bandera ---
SHAPEFILE_INEGI="${SHAPEFILE_INEGI:-}"               # valida cauces D8 vs INEGI dentro del paso 1
ANIOS_INCENDIO="${ANIOS_INCENDIO:-}"                 # ej. "2025" -- anios con incendio YA confirmado (paso 4)
EVENTOS_CONFIRMADOS="${EVENTOS_CONFIRMADOS:-}"       # ej. "2025:2025-04-17" -- fecha real conocida (paso 6)
CORRER_VALIDACION_INCENDIOS="${CORRER_VALIDACION_INCENDIOS:-1}"  # 0 para saltar el paso 6 completo
CORRER_AGUA_SUPERFICIAL="${CORRER_AGUA_SUPERFICIAL:-0}"  # 1 para correr el paso 7 (experimental, ver RUNBOOK.md)
CORRER_CORREDOR_DESCENDENTE="${CORRER_CORREDOR_DESCENDENTE:-0}"  # 1 para correr el paso 7d (experimental,
                                                                  # ver RUNBOOK.md) -- requiere FECHA_EVENTO_INCENDIO
FECHA_EVENTO_INCENDIO="${FECHA_EVENTO_INCENDIO:-}"    # ej. "2025-04-17" -- fecha real del incendio a evaluar
ANIO_HANSEN_INCENDIO="${ANIO_HANSEN_INCENDIO:-}"      # ej. "2025" -- año Hansen de ese mismo evento

id_lower="$(echo "$ID_PROYECTO" | tr '[:upper:]' '[:lower:]')"
mkdir -p "$CARPETA_SALIDA"

echo "=================================================================="
echo " Sitio: $ID_PROYECTO"
echo " Geojson: $GEOJSON"
echo " Carpeta de salida: $CARPETA_SALIDA"
echo "=================================================================="

echo
echo "=== [1/6] Terreno + Carbono (+ validacion hidrologica si hay INEGI) -- analizar_sitio.py ==="
ARGS_SITIO=(--geojson "$GEOJSON" --id-proyecto "$ID_PROYECTO" --proyecto-gee "$PROYECTO_GEE"
            --zonas "$ZONAS" --carpeta-salida "$CARPETA_SALIDA")
if [ -n "$SHAPEFILE_INEGI" ]; then
    ARGS_SITIO+=(--shapefile-inegi "$SHAPEFILE_INEGI")
fi
python3 -m core.analizar_sitio "${ARGS_SITIO[@]}"

echo
echo "=== [2/6] Deforestacion -- CSVs oficiales (SIN --mapa-3d) ==="
python3 -m core.deforestacion --geojson "$GEOJSON" --id-proyecto "$ID_PROYECTO" \
    --proyecto-gee "$PROYECTO_GEE" --zonas "$ZONAS" \
    --anio-inicio "$ANIO_INICIO" --anio-fin "$ANIO_FIN" --carpeta-salida "$CARPETA_SALIDA"

echo
echo "=== [3/6] Deforestacion -- mapa 3D (--mapa-3d) ==="
python3 -m core.deforestacion --geojson "$GEOJSON" --id-proyecto "$ID_PROYECTO" \
    --proyecto-gee "$PROYECTO_GEE" --zonas "$ZONAS" \
    --anio-inicio "$ANIO_INICIO" --anio-fin "$ANIO_FIN" --carpeta-salida "$CARPETA_SALIDA" --mapa-3d

echo
echo "=== [4/6] Cruce carbono x perdida + balance -- carbono_perdida.py ==="
ARGS_CRUCE=(--carbono-csv "$CARPETA_SALIDA/resumen_terreno_y_carbono_${id_lower}.csv"
            --perdida-sin-traslape-csv "$CARPETA_SALIDA/deforestacion_resumen_sin_traslape_${id_lower}.csv"
            --id-proyecto "$ID_PROYECTO" --carpeta-salida "$CARPETA_SALIDA" --balance)
if [ -n "$ANIOS_INCENDIO" ]; then
    ARGS_CRUCE+=(--anios-incendio "$ANIOS_INCENDIO")
fi
python3 -m core.carbono_perdida "${ARGS_CRUCE[@]}"

echo
echo "=== [5/6] Mapa 3D de CO2e liberado -- carbono_perdida.py --mapa-3d-perdida ==="
python3 -m core.carbono_perdida --mapa-3d-perdida --geojson "$GEOJSON" --id-proyecto "$ID_PROYECTO" \
    --zonas "$ZONAS" --carpeta-salida "$CARPETA_SALIDA"

if [ "$CORRER_VALIDACION_INCENDIOS" = "1" ]; then
    echo
    echo "=== [6/6] Validacion historica de incendios -- validacion_incendios.py --historial --mapa-3d ==="
    ARGS_INCENDIO=(--historial --geojson "$GEOJSON" --id-proyecto "$ID_PROYECTO"
                   --anio-inicio "$ANIO_INICIO" --anio-fin "$ANIO_FIN" --zonas "$ZONAS"
                   --historial-csv-existente "$CARPETA_SALIDA/deforestacion_historial_anual_${id_lower}.csv"
                   --proyecto-gee "$PROYECTO_GEE" --carpeta-salida "$CARPETA_SALIDA" --mapa-3d)
    if [ -n "$EVENTOS_CONFIRMADOS" ]; then
        ARGS_INCENDIO+=(--eventos-confirmados "$EVENTOS_CONFIRMADOS")
    fi
    python3 -m core.validacion_incendios "${ARGS_INCENDIO[@]}"
else
    echo
    echo "=== [6/6] Validacion historica de incendios: SALTADA (CORRER_VALIDACION_INCENDIOS=0) ==="
fi

if [ "$CORRER_AGUA_SUPERFICIAL" = "1" ]; then
    echo
    echo "=== [7 - opcional/experimental] Agua superficial visible por año -- agua_superficial.py ==="
    python3 -m core.agua_superficial --geojson "$GEOJSON" --id-proyecto "$ID_PROYECTO" \
        --proyecto-gee "$PROYECTO_GEE" --zonas "$ZONAS" --carpeta-salida "$CARPETA_SALIDA"
fi

if [ "$CORRER_CORREDOR_DESCENDENTE" = "1" ]; then
    if [ -z "$FECHA_EVENTO_INCENDIO" ] || [ -z "$ANIO_HANSEN_INCENDIO" ]; then
        echo "ERROR: CORRER_CORREDOR_DESCENDENTE=1 requiere FECHA_EVENTO_INCENDIO y ANIO_HANSEN_INCENDIO." >&2
        exit 1
    fi
    echo
    echo "=== [7d - opcional/experimental] Corredor hidrológico descendente (D8) -- corredor_descendente.py ==="
    python3 -m core.corredor_descendente --geojson "$GEOJSON" --id-proyecto "$ID_PROYECTO" \
        --proyecto-gee "$PROYECTO_GEE" --zonas "$ZONAS" --carpeta-salida "$CARPETA_SALIDA" \
        --fecha-evento "$FECHA_EVENTO_INCENDIO" --anio-hansen "$ANIO_HANSEN_INCENDIO" --mapa-3d
fi

echo
echo "=================================================================="
echo " LISTO. Resultados en: $CARPETA_SALIDA"
echo "=================================================================="

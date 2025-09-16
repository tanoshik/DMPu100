#!/usr/bin/env bash
set -euo pipefail

# --- settings ---
QUERY="data/query_profile_seed101.csv"

# タイムスタンプ付き出力ディレクトリ
ts=$(date +%Y%m%d%H%M%S)
OUTDIR="output/smoke_${ts}"
mkdir -p "$OUTDIR"

# サイズと表示名（name列の末尾に入るコア名）
declare -A LABELS=(
  [1000000]="S1M"
  [2000000]="S2M"
  [4000000]="S4M"
  [25000000]="S25M"
)

SIZES=(1000000 2000000 4000000 25000000)

echo "[INFO] Output dir: $OUTDIR"
for S in "${SIZES[@]}"; do
  DB="data/virtual_db_u100_S${S}_seed101.rds"
  TAG="${LABELS[$S]}_seed101"
  OUT="${OUTDIR}/scores_${TAG}.csv"

  if [[ ! -f "$DB" ]]; then
    echo "[WARN] missing DB: $DB (skip)"
    continue
  fi

  echo "[RUN] ${TAG}"
  Rscript scripts/matcher_onepass.R \
    --db="$DB" \
    --query="$QUERY" \
    --out="$OUT" \
    --report=all \
    --n_cap=200 \
    --score_min=0
done

echo "[DONE] See: ${OUTDIR}/time.log  (ヘッダ: LOAD_DB_SEC,LOAD_Q_SEC,COMP_SEC,TOTAL_SEC,PEAK_MiB,name)"

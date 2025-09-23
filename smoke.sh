#!/usr/bin/env bash
# smoke.sh - Phase-2 minimal smoke (ratio=0/1) with snake_case CLI.
# ASCII only. Git Bash / bash compatible.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")"/.. && pwd)"
cd "$ROOT_DIR"

R_BIN="${R_BIN:-Rscript}"

DB1K="data/virtual_db_u100_S1000_seed101.rds"
QCSV="data/query_profile_seed101.csv"
IDFCSV="data/gf_idf_mask.csv"

OUTDIR="output/smoke_phase2_1k"
BASE_DIR="$OUTDIR/baseline"
R0_DIR="$OUTDIR/ratio0"
R1_DIR="$OUTDIR/ratio1"
MASKED_R0="data/_masked_r0.rds"
MASKED_R1="data/_masked_r1.rds"

mkdir -p "$BASE_DIR" "$R0_DIR" "$R1_DIR"

echo "[SMK] Using:"
echo "  DB : $DB1K"
echo "  Q  : $QCSV"
echo "  IDF: $IDFCSV"
echo "  OUT: $OUTDIR"

# 0) Baseline (original DB) - score_min omitted (no threshold)
$R_BIN scripts/matcher_onepass.R \
  --db "$DB1K" \
  --query "$QCSV" \
  --out "$BASE_DIR/scores.csv" \
  --report all \
  --n_cap 1000 \
  --any_code 9999 \
  --h2a_on 0 \
  --debug 0

# 1) ratio=0 -> masked DB should be identical to baseline
$R_BIN scripts/cli/dmp_apply_idf_mask.R \
  --in_db "$DB1K" \
  --out_db "$MASKED_R0" \
  --idf_csv "$IDFCSV" \
  --ratio 0 \
  --seed 0 \
  --any_code 9999 \
  --debug 0

$R_BIN scripts/matcher_onepass.R \
  --db "$MASKED_R0" \
  --query "$QCSV" \
  --out "$R0_DIR/scores.csv" \
  --report all \
  --n_cap 1000 \
  --any_code 9999 \
  --h2a_on 0 \
  --debug 0

# Assert equality: baseline vs ratio0 scores
if diff -q "$BASE_DIR/scores.csv" "$R0_DIR/scores.csv" > /dev/null ; then
  echo "[OK] ratio=0 scores identical to baseline"
else
  echo "[ERR] ratio=0 scores differ from baseline" >&2
  exit 1
fi

# 2) ratio=1 -> masked DB; histogram should shift to higher scores on average
$R_BIN scripts/cli/dmp_apply_idf_mask.R \
  --in_db "$DB1K" \
  --out_db "$MASKED_R1" \
  --idf_csv "$IDFCSV" \
  --ratio 1 \
  --seed 0 \
  --any_code 9999 \
  --debug 0

$R_BIN scripts/matcher_onepass.R \
  --db "$MASKED_R1" \
  --query "$QCSV" \
  --out "$R1_DIR/scores.csv" \
  --report all \
  --n_cap 1000 \
  --any_code 9999 \
  --h2a_on 0 \
  --debug 0

# Simple check: mean(hist score) ratio1 >= baseline
$R_BIN - <<'RS'
b <- read.csv("output/smoke_phase2_1k/baseline/hist.csv",  stringsAsFactors=FALSE)
r <- read.csv("output/smoke_phase2_1k/ratio1/hist.csv",     stringsAsFactors=FALSE)
nm_b <- if ("Score" %in% names(b)) "Score" else names(b)[1]
nm_r <- if ("Score" %in% names(r)) "Score" else names(r)[1]
mb <- sum(as.integer(b[[nm_b]]) * as.integer(b[["Count"]])) / sum(as.integer(b[["Count"]]))
mr <- sum(as.integer(r[[nm_r]]) * as.integer(r[["Count"]])) / sum(as.integer(r[["Count"]]))
cat(sprintf("[CHK] mean score baseline=%.3f, ratio1=%.3f\n", mb, mr))
if (is.na(mb) || is.na(mr)) stop("NA in histogram means")
if (mr + 1e-9 < mb) stop("ratio1 mean score did not increase")
RS

echo "[DONE] smoke_phase2_1k passed"

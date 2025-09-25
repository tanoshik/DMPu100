#!/usr/bin/env bash
# smoke.sh - Phase-2 minimal smoke (ratio=0/1) with snake_case CLI.
# ASCII only. Git Bash / bash compatible.

set -euo pipefail
set -x

# --- Hardened R launcher: run Rscript.exe with clean env/vanilla ---
_RS_EXE=""
if [ -z "${_RS_EXE}" ] && [ -x "/c/Program Files/R/R-4.5.1/bin/Rscript.exe" ]; then
  _RS_EXE="/c/Program Files/R/R-4.5.1/bin/Rscript.exe"
elif [ -z "${_RS_EXE}" ] && [ -x "/c/Program Files/R/R-4.5.0/bin/Rscript.exe" ]; then
  _RS_EXE="/c/Program Files/R/R-4.5.0/bin/Rscript.exe"
else
  _RS_EXE="$(command -v Rscript || true)"
fi
if [ -z "${_RS_EXE}" ]; then
  echo "[ERR] Rscript not found" >&2; exit 1
fi

runR() {
  # 環境を空にしない。vanilla のみ指定（ユーザ/サイトプロファイルは無効化される）
  DMP_DBG="${DMP_DBG:-}" \
  MSYS2_ARG_CONV_EXCL="*" \
  CYGWIN="nodosfilewarning" \
  "$_RS_EXE" --vanilla "$@"
}

echo "[SMK] R_EXE=${_RS_EXE}"

# --- repo root をこのファイルの場所に固定 ---
ROOT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd -P)"
cd "$ROOT_DIR"
echo "[SMK] PWD=$PWD"

# サニティ
test -f scripts/matcher_onepass.R          || { echo "[ERR] scripts/matcher_onepass.R not found"; exit 1; }
test -f scripts/cli/dmp_apply_idf_mask.R   || { echo "[ERR] scripts/cli/dmp_apply_idf_mask.R not found"; exit 1; }

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

# 0) Baseline
DMP_DBG=1 runR scripts/matcher_onepass.R \
  --db "$DB1K" \
  --query "$QCSV" \
  --out "$BASE_DIR/scores.csv" \
  --report all \
  --n_cap 1000 \
  --any_code 9999 \
  --h2a_on 0 \
  --debug 0 \
|| { echo "[ERR] matcher_onepass baseline crashed"; test -f output/_trace_matcher.log && sed -n '1,200p' output/_trace_matcher.log; exit 1; }

# 1) ratio=0
runR scripts/cli/dmp_apply_idf_mask.R \
  --in_db "$DB1K" \
  --out_db "$MASKED_R0" \
  --idf_csv "$IDFCSV" \
  --ratio 0 \
  --seed 0 \
  --any_code 9999 \
  --debug 0

runR scripts/matcher_onepass.R \
  --db "$MASKED_R0" \
  --query "$QCSV" \
  --out "$R0_DIR/scores.csv" \
  --report all \
  --n_cap 1000 \
  --any_code 9999 \
  --h2a_on 0 \
  --debug 0

# 1’) histogram も出す（後段チェック用）
runR scripts/matcher_onepass.R \
  --db "$DB1K" \
  --query "$QCSV" \
  --out "$BASE_DIR/hist.csv" \
  --report hist \
  --n_cap 0 \
  --any_code 9999 \
  --h2a_on 0 \
  --debug 0
runR scripts/matcher_onepass.R \
  --db "$MASKED_R0" \
  --query "$QCSV" \
  --out "$R0_DIR/hist.csv" \
  --report hist \
  --n_cap 0 \
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

# 2) ratio=1
runR scripts/cli/dmp_apply_idf_mask.R \
  --in_db "$DB1K" \
  --out_db "$MASKED_R1" \
  --idf_csv "$IDFCSV" \
  --ratio 1 \
  --seed 0 \
  --any_code 9999 \
  --debug 0

runR scripts/matcher_onepass.R \
  --db "$MASKED_R1" \
  --query "$QCSV" \
  --out "$R1_DIR/scores.csv" \
  --report all \
  --n_cap 1000 \
  --any_code 9999 \
  --h2a_on 0 \
  --debug 0

# Simple check: mean(hist score) ratio1 >= baseline
runR - <<'RS'
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

#!/usr/bin/env bash
# smoke_r0.5.sh - Phase-2 ratio=0.5 only (snake_case). ASCII only.

set -euo pipefail

# --- R launcher (smoke.sh と同等の実績ルート) ---
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
  # smoke.sh と同じ最小環境・vanilla 起動
  MSYS2_ARG_CONV_EXCL="*" \
  CYGWIN="nodosfilewarning" \
  PATH="$(dirname "$_RS_EXE"):/usr/bin:/bin" \
  "$_RS_EXE" --vanilla "$@"
}

# --- 位置合わせ ---
ROOT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
cd "$ROOT_DIR"
echo "[R05] PWD=$PWD"
test -f scripts/cli/dmp_apply_idf_mask.R
test -f scripts/matcher_onepass.R

# --- 固定入出力 ---
DB="data/virtual_db_u100_S1000_seed101.rds"
Q="data/query_profile_seed101.csv"
IDF="data/gf_idf_mask.csv"
OUT="output/smoke_r0p5_1k"
MASKED="data/_masked_r0p5.rds"
RATIO="0.5"
SEED="0"
ANY="9999"

mkdir -p "$OUT"

echo "[R05] Using:"
echo "  DB   : $DB"
echo "  Q    : $Q"
echo "  IDF  : $IDF"
echo "  RATIO: $RATIO"
echo "  SEED : $SEED"
echo "  OUT  : $OUT"
echo "  R_EXE: $_RS_EXE"

# --- 事前存在チェック ---
[ -s "$DB" ]  || { echo "[ERR] missing $DB"  >&2; exit 1; }
[ -s "$Q" ]   || { echo "[ERR] missing $Q"   >&2; exit 1; }
[ -s "$IDF" ] || { echo "[ERR] missing $IDF" >&2; exit 1; }

# --- STEP1: マスク適用 (debug=1 で要約を残す) ---
echo "[R05] STEP1: mask -> $MASKED"
set +e
runR scripts/cli/dmp_apply_idf_mask.R \
  --in_db "$DB" \
  --out_db "$MASKED" \
  --idf_csv "$IDF" \
  --ratio "$RATIO" \
  --seed "$SEED" \
  --any_code "$ANY" \
  --debug 1 | tee "$OUT/_log_mask.txt"
rc=${PIPESTATUS[0]}
set -e
if [ $rc -ne 0 ]; then
  echo "[ERR] mask step failed (rc=$rc). See $OUT/_log_mask.txt" >&2
  exit $rc
fi
if [ ! -s "$MASKED" ]; then
  echo "[ERR] masked DB not written: $MASKED" >&2
  sed -n '1,200p' "$OUT/_log_mask.txt" || true
  exit 1
fi

# --- STEP2: matcher (report=all) ---
echo "[R05] STEP2: matcher_onepass (report=all)"
runR scripts/matcher_onepass.R \
  --db "$MASKED" \
  --query "$Q" \
  --out "$OUT/scores.csv" \
  --report all \
  --n_cap 1000 \
  --any_code "$ANY" \
  --h2a_on 0 \
  --debug 0

# --- 最終確認 ---
ls -l "$OUT"/scores.csv "$OUT"/hist.csv "$OUT"/opts_scores.json 2>/dev/null || true
echo "[DONE] smoke_r0p5_1k finished"

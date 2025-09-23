#!/usr/bin/env bash
set -euo pipefail

# test/regression_phase1.sh
# 最小回帰：fix_query -> onepass(strict) -> emit_detail(raw) の3成果物を
# ゴールデン(test/golden/phase1)とsha256で比較（ファイル名はベース名で比較）。

R_BIN="${R_BIN:-Rscript}"
DB1K="data/virtual_db_u100_S1000_seed101.rds"
QCSV="data/query_profile_seed101.csv"
GOLD_DIR="test/golden/phase1"

if [[ ! -f "$DB1K" || ! -f "$QCSV" ]]; then
  echo "[ERR] DB or query CSV not found"; exit 1
fi
mkdir -p "test/tmp/regress_phase1"
RUN_DIR="test/tmp/regress_phase1/$(date +%Y%m%d%H%M%S)/run"
mkdir -p "$RUN_DIR"

# 1) lenient ingest -> RDS
$R_BIN scripts/fix_query.R \
  --db "$DB1K" \
  --query "$QCSV" \
  --out_rds "$RUN_DIR/query_fixed.rds" \
  --no_audit 1 >/dev/null
echo "[OK] wrote $RUN_DIR/query_fixed.rds"

# 2) onepass(strict)
$R_BIN scripts/matcher_onepass.R \
  --db "$DB1K" \
  --query "$RUN_DIR/query_fixed.rds" \
  --out "$RUN_DIR/scores.csv" \
  --report all --n_cap 1000 --score_min 0 --debug 0 >/dev/null

# 3) emit_detail(raw) for all rows
$R_BIN scripts/emit/emit_detail.R \
  --opt_path "$RUN_DIR/opts_scores.json" \
  --scores_path "$RUN_DIR/scores.csv" \
  --out_dir "$RUN_DIR" \
  --any_code 9999 >/dev/null

# 4) manifest（ベース名で固定）
pushd "$RUN_DIR" >/dev/null
sha256sum scores.csv hist.csv raw_detail.csv > manifest.sha256
popd >/dev/null

# 5) ゴールデン側 manifest が無ければ一度だけ作る（初回整備用）
if [[ ! -f "$GOLD_DIR/manifest.sha256" ]]; then
  echo "[WARN] $GOLD_DIR/manifest.sha256 not found; creating it from current run"
  mkdir -p "$GOLD_DIR"
  cp "$RUN_DIR"/{scores.csv,hist.csv,raw_detail.csv} "$GOLD_DIR"/
  ( cd "$GOLD_DIR" && sha256sum scores.csv hist.csv raw_detail.csv > manifest.sha256 )
fi

# 6) 比較（ベース名で統一）
echo "[INFO] comparing manifests (3 files)"
diff -u "$GOLD_DIR/manifest.sha256" "$RUN_DIR/manifest.sha256" && {
  echo "[PASS] regression matched."
} || {
  echo "[FAIL] regression mismatch."; exit 2
}

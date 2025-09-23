#!/usr/bin/env bash
set -euo pipefail

# smoke.sh (Phase-1 strict + edges)
R_BIN="${R_BIN:-Rscript}"
OUTROOT="output/smoke_$(date +%Y%m%d%H%M%S)"
DB1K="data/virtual_db_u100_S1000_seed101.rds"
QCSV="data/query_profile_seed101.csv"
FREQ="data/freq_table.rds"

[[ -f "$DB1K" && -f "$QCSV" && -f "$FREQ" ]] || { echo "[ERR] DB/QCSV/FREQ missing"; exit 1; }

run_case () {
  local NAME="$1"; local QSRC="$2"; local EDGE_MODE="${3:-0}"
  local OUTD="$OUTROOT/$NAME"
  mkdir -p "$OUTD"

  echo "== RUN: $NAME =="

  # 1) lenient ingest (align L, ANY補完) -> RDS
  $R_BIN scripts/fix_query.R \
    --db "$DB1K" \
    --query "$QSRC" \
    --out_rds "$OUTD/query_fixed.rds" \
    --no_audit 1 >/dev/null
  echo "[OK] fix_query wrote $OUTD/query_fixed.rds"

  # 2) onepass(strict)
  $R_BIN scripts/matcher_onepass.R \
    --db "$DB1K" \
    --query "$OUTD/query_fixed.rds" \
    --out "$OUTD/scores.csv" \
    --report all --n_cap 1000 --score_min 0 --debug 0 >/dev/null

  # 3) emit raw_detail (全行)
  $R_BIN scripts/emit/emit_detail.R \
    --opt_path "$OUTD/opts_scores.json" \
    --scores_path "$OUTD/scores.csv" \
    --out_dir "$OUTD" \
    --any_code 9999 >/dev/null
  echo "[OK] emit_detail raw for 20 samples at $OUTD/_raw20"

  # 3.1) 20件を散らして抽出（上位/中位/下位）
  mkdir -p "$OUTD/_raw20"
  awk -F, 'NR==1||NR==2||NR==51||NR==101||NR==151||NR==201||NR==251||NR==301||NR==351||NR==401||NR==451||NR==501||NR==551||NR==601||NR==651||NR==701||NR==751||NR==801||NR==851||NR==1001{print}' "$OUTD/scores.csv" > "$OUTD/scores_20.csv"
  cut -d, -f1 "$OUTD/scores_20.csv" | sed '1d' > "$OUTD/ids_20.csv"
  while read -r sid; do
    awk -F, -v s="$sid" 'NR==1 || $NF==s' "$OUTD/raw_detail.csv" > "$OUTD/_raw20/${sid}.csv"
  done < "$OUTD/ids_20.csv"

  # 4) 計測（外付け）
  $R_BIN scripts/measure_wrap.R --out_dir "$OUTD" >/dev/null

  # 5) 簡易チェック
  [[ -s "$OUTD/scores.csv" && -s "$OUTD/hist.csv" && -s "$OUTD/raw_detail.csv" ]] || { echo "[ERR] missing outputs"; exit 2; }
  echo "[OK] $NAME"
}

# -----------------------
# main
# -----------------------
mkdir -p "$OUTROOT"

# case 1: 1k_normal_off
run_case "1k_normal_off" "$QCSV"

# case 2: edge_sparse_query（1行だけ＋他ANY）: 元CSVから1行だけ残す
SPQ="output/tmp/edge_sparse_query.csv"
mkdir -p output/tmp
(head -n 1 "$QCSV"; sed -n '2p' "$QCSV") > "$SPQ"
run_case "edge_sparse_query" "$SPQ"

# case 3: edge_typo_locus_NEG（strict onepass へ直接CSVを渡すので FIX しない）
echo "== RUN: edge_typo_locus_NEG =="
NEGCSV="output/tmp/edge_typo.csv"
cp "$QCSV" "$NEGCSV"
# 2列目の Locus 名を typo
awk -F, 'BEGIN{OFS=","} NR==1{print;next} NR==2{$1=$1"_X"; print; next} {print}' "$QCSV" > "$NEGCSV"

set +e
$R_BIN scripts/matcher_onepass.R \
  --db "$DB1K" \
  --query "$NEGCSV" \
  --out "$OUTROOT/edge_typo_locus_NEG/scores.csv" \
  --report top --n_cap 10 --score_min 0 --debug 0 >/dev/null
rc=$?
set -e
if [[ $rc -ne 0 ]]; then
  echo "[OK] NEGATIVE as expected (loci coverage)"
else
  echo "[ERR] NEG expected but succeeded"; exit 3
fi

echo "[DONE] smoke(min strict+edges) => $OUTROOT"

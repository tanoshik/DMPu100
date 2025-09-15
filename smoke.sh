#!/usr/bin/env bash
set -euo pipefail

# ==== config ====
DB="data/virtual_db_u100_S1000000_seed101.rds"
QCSV="data/query_profile_seed101.csv"
OUTDIR="output/smoke_$(date +%Y%m%d%H%M%S)"
mkdir -p "$OUTDIR"

# 小道具: CSV行数（ヘッダ除く）を数える
csv_rows () { awk 'END{print (NR>0?NR-1:0)}' "$1"; }

echo "== smoke start =="
echo "DB  : $DB"
echo "QCSV: $QCSV"
echo "OUT : $OUTDIR"
echo

# 1) report=all（TopN, n_cap=200）  ← 基本動作
echo "[1/4] report=all (TopN, n_cap=200)"
Rscript scripts/matcher_onepass.R \
  --db "$DB" --query "$QCSV" \
  --out "$OUTDIR/scores_topn_all.csv" \
  --report all \
  --n_cap 200

S1="$OUTDIR/scores_topn_all.csv"
H1="${S1/$OUTDIR\/scores/$OUTDIR\/hist}"
echo "  scores: $S1  rows=$(csv_rows "$S1")"
echo "  hist  : $H1  rows=$(csv_rows "$H1")"
echo

# 2) report=top（Threshold, score_min=30, n_cap=200） ← T指定のスコア出力のみ
echo "[2/4] report=top (Threshold, score_min=30, n_cap=200)"
Rscript scripts/matcher_onepass.R \
  --db "$DB" --query "$QCSV" \
  --out "$OUTDIR/scores_t30_top.csv" \
  --report top \
  --score_min 30 \
  --n_cap 200

S2="$OUTDIR/scores_t30_top.csv"
echo "  scores: $S2  rows=$(csv_rows "$S2")"
echo

# 3) report=hist（Hist-only） ← heapスキップ動作の確認
echo "[3/4] report=hist (Hist-only fast path)"
Rscript scripts/matcher_onepass.R \
  --db "$DB" --query "$QCSV" \
  --out "$OUTDIR/scores_histonly.csv" \
  --report hist \
  --n_cap 200

H3="$OUTDIR/hist_histonly.csv"
# 実名は scores_* → hist_* に置換されます（上と同じ命名規則）
H3="${H3/$OUTDIR\/hist_histonly.csv/$(ls "$OUTDIR"/hist_*.csv | xargs -n1 basename | grep -E '^hist(_t[0-9]+)?_?histonly\.csv$' || true)}"
# フォールバック（単一ファイル想定）
if [ -z "$H3" ]; then H3="$(ls "$OUTDIR"/hist_*.csv | head -n1)"; else H3="$OUTDIR/$H3"; fi
echo "  hist  : $H3  rows=$(csv_rows "$H3")"
echo

# 4) report=all + n_cap=0（Histは出す／scoresは0件） ← 研究用途
echo "[4/4] report=all (n_cap=0 → scores=0件, histは出力)"
Rscript scripts/matcher_onepass.R \
  --db "$DB" --query "$QCSV" \
  --out "$OUTDIR/scores_all_n0.csv" \
  --report all \
  --n_cap 0

S4="$OUTDIR/scores_all_n0.csv"
H4="${S4/$OUTDIR\/scores/$OUTDIR\/hist}"
echo "  scores: $S4  rows=$(csv_rows "$S4")  (期待値=0)"
echo "  hist  : $H4  rows=$(csv_rows "$H4")"
echo

echo "== smoke done =="

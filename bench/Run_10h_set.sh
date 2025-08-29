#!/usr/bin/env bash
set -u

# ==== user knobs ====
DB_2M="data/virtual_db_u100_S2000000_seed123.rds"
DB_1M="data/virtual_db_u100_S1000000_seed123.rds"
DB_100K="data/virtual_db_u100_S100000_seed123.rds"
QMODE="csv"               # csv | db_first
BS_MAIN=32000             # main block_size
TAG_PREFIX="_AUTO$(date +%Y%m%d_%H%M)"  # unique tag root
RSC="Rscript bench/bench_run_match_fast_chunks_u100.R"

log() { echo "[$(date +%F\ %T)] $*"; }

run() {
  log "RUN: $*"
  eval "$@" || log "WARN: command returned non-zero"
}

# ---- 2) workers 微調整：80k × (6,8,10,12) ----
for wk in 6 8 10 12; do
  run "$RSC --db=$DB_2M --query=$QMODE --use_cpp=TRUE --block_size=$BS_MAIN \
      --plan=80000x${wk} \
      --tag=${TAG_PREFIX}_CPP_2M_bs${BS_MAIN}_80k_w${wk}"
done

# ---- 3) block_size スイープ：100k×6, bs={16k,32k,48k} ----
for bs in 16000 32000 48000; do
  run "$RSC --db=$DB_2M --query=$QMODE --use_cpp=TRUE --block_size=$bs \
      --plan=100000x6 \
      --tag=${TAG_PREFIX}_CPP_2M_bs${bs}_100k6"
done

# ---- 4) Soak：2M / 50k×10 を連投（回数は必要に応じ調整）----
for i in $(seq 1 30); do
  run "$RSC --db=$DB_2M --query=$QMODE --use_cpp=TRUE --block_size=$BS_MAIN \
      --plan=50000x10 \
      --tag=${TAG_PREFIX}_CPP_2M_bs${BS_MAIN}_soak_${i}"
done

# ---- 5) ごく一部だけ整合性チェック（Rcpp, 100k と 2M の2本）----
run "$RSC --db=$DB_100K --query=$QMODE --use_cpp=TRUE --block_size=16000 \
    --plan=100000x6 \
    --consistency=true --consistency_rows=1000 \
    --tag=${TAG_PREFIX}_CONSIST_CPP_100k"
run "$RSC --db=$DB_2M --query=$QMODE --use_cpp=TRUE --block_size=$BS_MAIN \
    --plan=50000x10 \
    --consistency=true --consistency_rows=1000 \
    --tag=${TAG_PREFIX}_CONSIST_CPP_2M_w10_bs${BS_MAIN}"

# ---- 6) 参考：R版ベースライン（単発）----
run "$RSC --db=$DB_1M --query=$QMODE --use_cpp=FALSE --block_size=$BS_MAIN \
    --plan=100000x6 \
    --tag=${TAG_PREFIX}_R_1M_bs${BS_MAIN}_100k6"

log "ALL DONE."

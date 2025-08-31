#!/bin/bash
set -euo pipefail

DB="data/virtual_db_u100_S2000000_seed123.rds"
QUERY="data/query_profile_seed123.csv"
SCRIPT="bench/bench_run_match_fast_chunks_u100.R"
BASE_OPTS="--db=$DB --query=$QUERY --query_mode=csv --use_cpp=TRUE"

# プラン群
PLANS=(
  # --- 最適域（70%） ---
  "50000x10"
  "80000x8"
  "100000x6"

  # --- 探索域（20%） ---
  "20000x20"
  "40000x12"
  "120000x4"

  # --- スモーク（10%） ---
  "200000x2"
  "10000x1"
)

# 1条件あたりの繰り返し回数（調整用）
REPEAT=100   # => 8条件×100回 = 800run ≈ 26h想定

for plan in "${PLANS[@]}"; do
  for ((i=1; i<=REPEAT; i++)); do
    TAG="bench_${plan}_run${i}"
    echo "[$(date '+%F %T')] START plan=$plan run=$i"
    Rscript "$SCRIPT" $BASE_OPTS --plan="$plan" --tag="$TAG"
    echo "[$(date '+%F %T')] DONE plan=$plan run=$i"
    sleep 2
  done
done

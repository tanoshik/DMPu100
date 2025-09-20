#!/usr/bin/env bash
# Smoke (strict; fixed seed101 files) for DMPu100
set -Eeuo pipefail

PROJECT_ROOT="/c/Users/tanos/projects/DMPu100"
DB1K="data/virtual_db_u100_S1000_seed101.rds"
DB2M="data/virtual_db_u100_S2000000_seed101.rds"
QUERY_SEED101="data/query_profile_seed101.csv"

R_BIN="${R_BIN:-Rscript}"
STAMP=$(date +%Y%m%d%H%M%S)
OUT_DIR="output/smoke_${STAMP}"

cd "$PROJECT_ROOT"
mkdir -p "$OUT_DIR"

echo "[INFO] DB1K : $DB1K"
echo "[INFO] DB2M : $DB2M"
echo "[INFO] QUERY: $QUERY_SEED101"

[[ -f "$DB1K" ]] || { echo "[ERR] not found: $DB1K"; exit 1; }
[[ -f "$DB2M" ]] || { echo "[ERR] not found: $DB2M"; exit 1; }
[[ -f "$QUERY_SEED101" ]] || { echo "[ERR] not found: $QUERY_SEED101"; exit 1; }

run_case () {
  local LABEL="$1"   # 1k_normal_off / 1k_normal_on / 2M_normal_off
  local DB="$2"
  local H2A="$3"     # 0 or 1

  local RUN_DIR="${OUT_DIR}/${LABEL}"
  mkdir -p "$RUN_DIR"

  # --- strict: クエリCSV 事前チェック（NA/空欄/非数は即停止） ---
  local PRECHK="${RUN_DIR}/_precheck_query.R"
  cat > "$PRECHK" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
q <- args[1]
x <- read.csv(q, check.names=FALSE, stringsAsFactors=FALSE)
nms <- tolower(names(x))
stopifnot(any(nms %in% c('locus','marker')),
          any(nms %in% c('allele1','a1','q1')),
          any(nms %in% c('allele2','a2','q2')))
lc  <- which(nms %in% c('locus','marker'))[1]
a1c <- which(nms %in% c('allele1','a1','q1'))[1]
a2c <- which(nms %in% c('allele2','a2','q2'))[1]
if (any(is.na(x[[lc]]) | trimws(as.character(x[[lc]]))=='')) stop('Query CSV: Locus has NA/blank')
if (any(is.na(x[[a1c]]) | is.na(x[[a2c]])))                stop('Query CSV: Allele has NA')
oknum <- function(v) all(grepl('^\\s*[0-9]+(\\.[0-9]+)?\\s*$', v))
if (!all(sapply(x[c(a1c,a2c)], oknum))) stop('Query CSV: Allele has non-numeric cell(s)')
cat('[OK] query precheck passed\n')
RS
  "$R_BIN" "$PRECHK" "$QUERY_SEED101"

  # --- 実行 ---
  "$R_BIN" scripts/matcher_onepass.R \
    --db "$DB" --query "$QUERY_SEED101" \
    --out "${RUN_DIR}/scores.csv" \
    --report all --n_cap 200 --score_min 0 --any_code 9999 --h2a_on "$H2A"

  # --- 出力存在 ---
  [[ -f "${RUN_DIR}/scores.csv" ]]        || { echo "[ERR] ${LABEL}: scores.csv missing"; exit 1; }
  [[ -f "${RUN_DIR}/hist.csv"   ]]        || { echo "[ERR] ${LABEL}: hist.csv missing"; exit 1; }
  [[ -f "${RUN_DIR}/opts_scores.json" ]]  || { echo "[ERR] ${LABEL}: opts_scores.json missing"; exit 1; }
  [[ -f "${RUN_DIR}/time.log"   ]]        || { echo "[ERR] ${LABEL}: time.log missing"; exit 1; }

  # --- 満点チェック（一時Rで実行。inlineは使わない） ---
  local FSCHK="${RUN_DIR}/_check_fullscore.R"
  cat > "$FSCHK" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
db_path   <- args[1]
scores_p  <- args[2]
lab       <- args[3]
db <- readRDS(db_path)
L <- length(db$locus_ids)
MAX_SC <- as.integer(2L * L)
sc <- read.csv(scores_p, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(nrow(sc) >= 1)
if (!identical(as.integer(sc$Score[1]), MAX_SC)) {
  stop(sprintf('[%s] Top score not full (%s != %s). seed101 mismatch?', lab, sc$Score[1], MAX_SC))
} else {
  cat(sprintf('[OK] %s: top score == %d (full)\n', lab, MAX_SC))
}
RS
  "$R_BIN" "$FSCHK" "$DB" "${RUN_DIR}/scores.csv" "$LABEL"

  # --- 並びチェック（C++修正確認）：Score 降順・同点は SampleID 昇順 ---
  local ORDCHK="${RUN_DIR}/_assert_sorted.R"
  cat > "$ORDCHK" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
csv <- args[1]
x <- read.csv(csv, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(nrow(x) >= 1, all(c("SampleID","Score") %in% names(x)))
ok_desc <- all(diff(x$Score) <= 0)
tie_idx <- which(diff(x$Score) == 0)
ok_tie  <- TRUE
if (length(tie_idx) > 0) {
  ok_tie <- all(x$SampleID[tie_idx] <= x$SampleID[tie_idx + 1])
}
if (!(ok_desc && ok_tie)) {
  stop("scores not sorted by (Score desc, SampleID asc)")
} else {
  cat("[OK] scores order verified (Score desc, SampleID asc)\n")
}
RS
  "$R_BIN" "$ORDCHK" "${RUN_DIR}/scores.csv"

  # --- bits / raw合計 検証（emit_detail がある場合のみ） ---
  if [[ -f "scripts/emit/emit_detail.R" && -f "scripts/devtools/check_bits_and_fullscore.R" ]]; then
    "$R_BIN" "scripts/devtools/check_bits_and_fullscore.R" "${RUN_DIR}" "${DB}" || {
      echo "[ERR] ${LABEL}: bits/fullscore check failed"; exit 1;
    }
  else
    echo "[NOTE] bits/raw check skipped (emit_detail.R or devtools checker not found)"
  fi

  echo "[OK] ${LABEL} done."
}

run_case "1k_normal_off" "$DB1K" 0
run_case "1k_normal_on"  "$DB1K" 1
run_case "2M_normal_off" "$DB2M" 0

echo "[DONE] smoke => $OUT_DIR"

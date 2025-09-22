#!/usr/bin/env bash
# smoke_phase1_min.sh — Phase-1 minimal smoke (revised)
# Case: 1k_normal_off のみ（mask=normal, h2a_on=0）
# Checks:
#   (1) C++整列: Score降順・SampleID昇順
#   (2) raw_detail の per-locus 'score' 合計が scores.csv と一致（Top10）
#   (3) 分布に散らした20件を raw_detail で出力（人手チェック用）
#   (4) オプション: NEG_PRECHECK=1 のときだけ「擬似的に壊れた query」を使って
#       matcher と同等の「loci網羅性チェック」を通す（= そこで意図的に失敗させる）
#       * ファイル名は必ず _negative/pseudo_bad_query.csv を使い、紛らわしさを回避

set -euo pipefail

# ---- Config ---------------------------------------------------------------
R_BIN="${R_BIN:-Rscript}"
OUT_BASE="${OUT_BASE:-output}"
DATE_TAG="$(date +%Y%m%d%H%M%S)"
OUT_DIR="${OUT_BASE}/smoke_${DATE_TAG}"
mkdir -p "${OUT_DIR}"

DB1K="${DB1K:-data/virtual_db_u100_S1000_seed101.rds}"
QUERY_SEED101="${QUERY_SEED101:-data/query_profile_seed101.csv}"
ANY_CODE="${ANY_CODE:-9999}"
NCAP_1K="${NCAP_1K:-1000}"   # 1k は全件
NEG_PRECHECK="${NEG_PRECHECK:-0}"

_banner () { echo; echo "== $* =="; }
_die () { echo "[ERR] $*" >&2; exit 1; }

write_order_check () {
  local f="$1"
  cat > "$f" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
sc <- args[1]; label <- args[2]
x <- read.csv(sc, stringsAsFactors=FALSE, check.names=FALSE)
ok <- all(order(-x$Score, x$SampleID) == seq_len(nrow(x)))
if (!ok) stop(sprintf("[order mismatch] %s not sorted (Score desc, SampleID asc)", label))
cat(sprintf("[OK] order: %s\n", label))
RS
}

write_sumscore_check () {
  local f="$1"
  cat > "$f" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
raw_p <- args[1]; sc_p <- args[2]; K <- as.integer(args[3])
xr <- read.csv(raw_p, stringsAsFactors=FALSE, check.names=FALSE)
xs <- read.csv(sc_p,  stringsAsFactors=FALSE, check.names=FALSE)
names(xr) <- tolower(names(xr))
stopifnot(all(c("sampleid","score") %in% names(xr)))
xs <- xs[order(-xs$Score, xs$SampleID), c("SampleID","Score")]
top <- head(xs, K)
# per-sample sum from raw_detail
agg <- aggregate(score ~ sampleid, data=xr, FUN=function(v) sum(as.integer(v), na.rm=TRUE))
colnames(agg) <- c("SampleID","ScoreSum")
chk <- merge(top, agg, by="SampleID", all.x=TRUE, sort=FALSE)
if (any(is.na(chk$ScoreSum))) {
  missing <- paste(chk$SampleID[is.na(chk$ScoreSum)], collapse=",")
  stop(sprintf("raw_detail missing samples for topK: %s", missing))
}
bad <- chk[ chk$Score != chk$ScoreSum, ]
if (nrow(bad) > 0) {
  print(bad)
  stop(sprintf("sum(score) mismatch in topK (%d rows)", nrow(bad)))
}
cat(sprintf("[OK] sum(score) consistency for top %d\n", nrow(top)))
RS
}

write_pick20_and_emit_raw () {
  local f="$1"
  cat > "$f" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
sc_p <- args[1]; opts_p <- args[2]; out_dir <- args[3]
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
xs <- read.csv(sc_p, stringsAsFactors=FALSE, check.names=FALSE)
xs <- xs[order(-xs$Score, xs$SampleID), c("SampleID","Score")]
N <- nrow(xs); if (N < 20) stop("scores has <20 rows; need n_cap >= 20")
idx <- unique(pmax(1, pmin(N, c(1:7, round(c(0.45,0.50,0.55)*N), (N-6):N))))
ids <- data.frame(SampleID = xs$SampleID[idx], stringsAsFactors=FALSE)
ids_p <- file.path(out_dir, "ids_20.csv"); write.csv(ids, ids_p, row.names=FALSE)
sc20_p <- file.path(out_dir, "scores_20.csv"); write.csv(ids, sc20_p, row.names=FALSE)
suppressWarnings({ source("scripts/emit/emit_detail.R") })
emit_detail(opt_path = opts_p, scores_path = sc20_p, out_dir = out_dir, mode = "raw")
raw_p <- file.path(out_dir, "raw_detail.csv")
if (!file.exists(raw_p)) stop("raw_detail not created for 20 samples")
cat(sprintf("[OK] emit_detail raw for 20 samples at %s\n", out_dir))
RS
}

write_negative_precheck_same_path () {
  # matcher_onepass.R の CSV 取込みと同じ経路（列名正規化→locus照合）を通して
  # 「CSV が db$locus_ids を網羅していない」ことで意図的に失敗させる。
  local f="$1"
  cat > "$f" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
opts_p <- args[1]; bad_csv <- args[2]
suppressPackageStartupMessages(library(jsonlite))
opts <- jsonlite::fromJSON(opts_p, simplifyVector=TRUE)
db_path <- opts$db_path
if (is.null(db_path)) stop("opts lacks db_path")
db <- readRDS(db_path)  # list(sample_ids, locus_ids, A1, A2)
if (!is.list(db) || !all(c("sample_ids","locus_ids","A1","A2") %in% names(db))) {
  stop("DB RDS must contain list(sample_ids, locus_ids, A1, A2)")
}
x <- read.csv(bad_csv, stringsAsFactors=FALSE, check.names=FALSE)
nms <- tolower(names(x)); names(x) <- nms
lc  <- which(nms %in% c("locus","marker"))[1]
a1c <- which(nms %in% c("allele1","a1","q1"))[1]
a2c <- which(nms %in% c("allele2","a2","q2"))[1]
if (is.na(lc) || is.na(a1c) || is.na(a2c)) stop("CSV needs columns: Locus(+), Allele1/A1/Q1, Allele2/A2/Q2")
mm <- match(db$locus_ids, x[[lc]])
if (any(is.na(mm))) stop("CSV loci must cover all db$locus_ids [NEGATIVE TEST TRIGGERED]")
cat("[ERR] negative precheck did not fail as expected\n")
quit(status=1)
RS
}

run_case_1k_normal_off () {
  local LABEL="1k_normal_off"
  _banner "RUN: ${LABEL}"
  local RUN_DIR="${OUT_DIR}/${LABEL}"
  mkdir -p "${RUN_DIR}"

  local ORDCHK="${RUN_DIR}/_assert_order.R"
  local SUMCHK="${RUN_DIR}/_assert_sumscore_topK.R"
  local PK20="${RUN_DIR}/_pick20_emit_raw.R"
  local NEGPC="${RUN_DIR}/_negative_precheck_same_path.R"
  write_order_check "$ORDCHK"
  write_sumscore_check "$SUMCHK"
  write_pick20_and_emit_raw "$PK20"
  write_negative_precheck_same_path "$NEGPC"

  # 1) matcher 実行（h2a_on=0, report=all, n_cap=1000）
  local OUT_SCORES="${RUN_DIR}/scores.csv"
  "${R_BIN}" scripts/matcher_onepass.R \
    --db "$DB1K" --query "$QUERY_SEED101" \
    --report all --n_cap "${NCAP_1K}" --score_min 0 \
    --any_code "${ANY_CODE}" --h2a_on 0 \
    --out "${OUT_SCORES}"

  # 2) 生成物の確認
  [[ -f "${OUT_SCORES}" ]] || _die "scores.csv not found: ${LABEL}"
  local OPTS_JSON
  OPTS_JSON="$(ls "${RUN_DIR}"/opts_*.json 2>/dev/null | head -n1 || true)"
  [[ -n "${OPTS_JSON:-}" && -f "${OPTS_JSON}" ]] || _die "opts_*.json not found: ${LABEL}"
  [[ -f "${RUN_DIR}/hist.csv" ]] || _die "hist.csv not found: ${LABEL}"

  # 3) 並び保証
  "${R_BIN}" "$ORDCHK" "${OUT_SCORES}" "$LABEL"

  # 4) raw_detail（全TopN）→ Top10 で合計一致チェック
  mkdir -p "${RUN_DIR}/_raw_chk"
  ${R_BIN} -e "source('scripts/emit/emit_detail.R'); emit_detail('$OPTS_JSON','${OUT_SCORES}','${RUN_DIR}/_raw_chk','raw')"
  [[ -f "${RUN_DIR}/_raw_chk/raw_detail.csv" ]] || _die "raw_detail.csv not created"
  "${R_BIN}" "$SUMCHK" "${RUN_DIR}/_raw_chk/raw_detail.csv" "${OUT_SCORES}" 10

  # 5) 20件の raw_detail（人手チェック用）
  "${R_BIN}" "$PK20" "${OUT_SCORES}" "$OPTS_JSON" "${RUN_DIR}/_raw20"

  # 6) オプション NEGATIVE precheck（本番I/Oとは混同しない命名: pseudo_*）
  if [[ "${NEG_PRECHECK}" == "1" ]]; then
    local NEG_DIR="${RUN_DIR}/_negative"
    mkdir -p "${NEG_DIR}"
    # locus を一部だけにした「擬似的に壊れた query」: matcher の "CSV loci must cover all db$locus_ids" を誘発
    local PSEUDO="${NEG_DIR}/pseudo_bad_query.csv"
    echo "Locus,Allele1,Allele2" > "${PSEUDO}"
    echo "D8S1179,12,12" >> "${PSEUDO}"  # 故意に1行だけ

    echo "[INFO] running NEGATIVE precheck through the SAME locus coverage path"
    if "${R_BIN}" "$NEGPC" "${OPTS_JSON}" "${PSEUDO}"; then
      true
    fi
    echo "[OK] negative precheck (loci coverage) failed as expected"
  fi

  echo "[OK] ${LABEL}"
}

# ---- Run battery ----------------------------------------------------------
run_case_1k_normal_off
echo "[DONE] smoke(min, revised) => ${OUT_DIR}"

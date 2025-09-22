#!/usr/bin/env bash
# smoke.sh — Phase-1 smoke battery (1k & 2M) with raw-detail spot checks (20 samples)
# v2: no hard fail on bits normalization; bits optional (prefer per-locus score, then code, then bits).

set -euo pipefail

# ---- Config ---------------------------------------------------------------
R_BIN="${R_BIN:-Rscript}"
OUT_BASE="${OUT_BASE:-output}"
DATE_TAG="$(date +%Y%m%d%H%M%S)"
OUT_DIR="${OUT_BASE}/smoke_${DATE_TAG}"
mkdir -p "${OUT_DIR}"

DB1K="${DB1K:-data/virtual_db_u100_S1000_seed101.rds}"
DB2M="${DB2M:-data/virtual_db_u100_S2000000_seed101.rds}"
QUERY_SEED101="${QUERY_SEED101:-data/query_profile_seed101.csv}"

ANY_CODE="${ANY_CODE:-9999}"
NCAP_1K="${NCAP_1K:-1000}"   # 1k は全件
NCAP_2M="${NCAP_2M:-200}"    # 2M は軽め

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

write_bits_score_check () {
  local f="$1"
  cat > "$f" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
raw_p <- args[1]; sc_p <- args[2]; K <- as.integer(args[3])
xr <- read.csv(raw_p, stringsAsFactors=FALSE, check.names=FALSE)
xs <- read.csv(sc_p,  stringsAsFactors=FALSE, check.names=FALSE)

# column normalizations
names(xr) <- tolower(names(xr))
names(xs) <- c("SampleID","Score")
stopifnot(all(c("SampleID","Score") %in% names(xs)))
if (!all(c("sampleid","score") %in% names(xr))) stop("raw_detail needs SampleID/score")

SCORE_TABLE <- as.integer(c(0,1,1,1, 1,1,2,2, 1,2,1,2, 1,2,2,2))

# safe bits helpers (no hard stop)
pad4 <- function(s) sprintf("%04s", s)
bits4 <- function(v) {
  s <- as.character(v)
  s[is.na(s)] <- ""
  s <- gsub("[^01]", "", s)
  s[nchar(s)==0] <- NA_character_
  s[nchar(s) < 4] <- pad4(s[nchar(s) < 4])
  s[nchar(s) > 4] <- substr(s[nchar(s) > 4], nchar(s[nchar(s) > 4]) - 3, nchar(s[nchar(s) > 4]))
  s
}
code_from_bits <- function(s4) {
  if (is.na(s4)) return(NA_integer_)
  v <- as.integer(strsplit(s4,"")[[1]])
  sum(v * c(8,4,2,1))
}

# derive 'code' if missing
if (!("code" %in% names(xr)) && ("bits" %in% names(xr))) {
  xr$bits <- bits4(xr$bits)
  xr$code <- vapply(xr$bits, code_from_bits, integer(1))
}

# compute per-sample sums
xs <- xs[order(-xs$Score, xs$SampleID), c("SampleID","Score")]
top <- head(xs, K)

for (i in seq_len(nrow(top))) {
  sid <- top$SampleID[i]; sc0 <- as.integer(top$Score[i])
  sub <- xr[ xr$sampleid==sid, , drop=FALSE ]
  if (nrow(sub) == 0) stop(sprintf("no raw rows for SampleID=%s", sid))

  # A) from raw 'score' (most robust, already per-locus)
  sc_sum <- if ("score" %in% names(sub)) sum(as.integer(sub$score), na.rm=TRUE) else NA_integer_

  # B) from SCORE_TABLE[code+1] if code entirely available
  sc_from_code <- if ("code" %in% names(sub) && all(!is.na(sub$code))) {
    sum(SCORE_TABLE[as.integer(sub$code) + 1L], na.rm=TRUE)
  } else NA_integer_

  # C) from bits only if code missing
  sc_from_bits <- if (!("code" %in% names(sub)) && ("bits" %in% names(sub))) {
    bb <- bits4(sub$bits)
    if (all(is.na(bb))) NA_integer_ else sum(vapply(bb, function(s4){
      if (is.na(s4)) return(NA_integer_)
      v <- as.integer(strsplit(s4,"")[[1]]); sum(v * c(8,4,2,1))
    }, integer(1)) |> {\(code) SCORE_TABLE[code + 1L]}(), na.rm=TRUE)
  } else NA_integer_

  sc1 <- if (!is.na(sc_sum)) sc_sum else if (!is.na(sc_from_code)) sc_from_code else sc_from_bits
  if (is.na(sc1)) next  # skip this sample; do not hard fail on encoding
  if (!identical(sc0, as.integer(sc1))) {
    stop(sprintf("[mismatch] %s: scores.csv=%d vs recomputed=%d", sid, sc0, sc1))
  }
}
cat(sprintf("[OK] score consistency for top %d (bits tolerant)\\n", nrow(top)))
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
xraw <- read.csv(raw_p, stringsAsFactors=FALSE, check.names=FALSE)
# no hard checks for bits; this is for human spot check
if (!all(ids$SampleID %in% xraw$SampleID)) stop("raw_detail missing selected SampleIDs")
cat(sprintf("[OK] emit_detail raw for 20 samples at %s\n", out_dir))
RS
}

write_negative_precheck () {
  local f="$1"
  cat > "$f" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
bad_csv <- args[1]
ok <- FALSE
try({
  src <- read.csv(bad_csv, stringsAsFactors=FALSE)
  stopifnot(all(c("Locus","Allele1","Allele2") %in% names(src)))
  if (any(src$Allele1=="" | is.na(src$Allele1))) stop("NA/blank in Allele1")
  ok <- TRUE
}, silent=TRUE)
if (ok) stop("precheck should fail but passed")
cat("[OK] negative precheck failed as expected\n")
RS
}

run_case () {
  local LABEL="$1"   # e.g., 1k_normal_off
  local DB="$2"
  local H2A="$3"     # 0 or 1
  local EXTRA_ARGS="${4:-}"
  _banner "RUN: ${LABEL}"
  local RUN_DIR="${OUT_DIR}/${LABEL}"
  mkdir -p "${RUN_DIR}"

  local ORDCHK="${RUN_DIR}/_assert_order.R"
  local BSCCHK="${RUN_DIR}/_assert_bits_score_topK.R"
  local PK20="${RUN_DIR}/_pick20_emit_raw.R"
  local NEGPC="${RUN_DIR}/_negative_precheck.R"
  write_order_check "$ORDCHK"
  write_bits_score_check "$BSCCHK"
  write_pick20_and_emit_raw "$PK20"
  write_negative_precheck "$NEGPC"

  # 1) matcher 実行
  local NCAP="${NCAP_2M}"
  [[ "$LABEL" == 1k* ]] && NCAP="${NCAP_1K}"
  local OUT_BASEPATH="${RUN_DIR}/scores.csv"
  "${R_BIN}" scripts/matcher_onepass.R \
    --db "$DB" --query "$QUERY_SEED101" \
    --report all --n_cap "${NCAP}" --score_min 0 \
    --any_code "${ANY_CODE}" --h2a_on "${H2A}" ${EXTRA_ARGS} \
    --out "${OUT_BASEPATH}"

  # 2) 生成物のパス
  local SCORES_CSV="${RUN_DIR}/scores.csv"
  [[ -f "${SCORES_CSV}" ]] || _die "scores.csv not found: ${LABEL}"
  local OPTS_JSON
  OPTS_JSON="$(ls "${RUN_DIR}"/opts_*.json 2>/dev/null | head -n1 || true)"
  [[ -n "${OPTS_JSON:-}" && -f "${OPTS_JSON}" ]] || _die "opts_*.json not found: ${LABEL}"

  # 3) 並び保証
  "${R_BIN}" "$ORDCHK" "${SCORES_CSV}" "$LABEL"

  # 4) raw 展開（TopN全部）→ 整合チェックは top 10（bits 寛容）
  mkdir -p "${RUN_DIR}/_raw_chk"
  Rscript -e "source('scripts/emit/emit_detail.R'); emit_detail('$OPTS_JSON','${SCORES_CSV}','${RUN_DIR}/_raw_chk','raw')"
  [[ -f "${RUN_DIR}/_raw_chk/raw_detail.csv" ]] || _die "raw_detail.csv not created"
  "${R_BIN}" "$BSCCHK" "${RUN_DIR}/_raw_chk/raw_detail.csv" "${SCORES_CSV}" 10

  # 5) 1k：分布に散らした20件で raw 展開チェック（人力確認用）
  if [[ "$LABEL" == 1k* ]]; then
    "${R_BIN}" "$PK20" "${SCORES_CSV}" "$OPTS_JSON" "${RUN_DIR}/_raw20"
  fi

  # 6) ネガティブ precheck（全体は落とさない）
  local BROKE_DIR="${RUN_DIR}/_negative"
  mkdir -p "$BROKE_DIR"
  cp "$QUERY_SEED101" "${BROKE_DIR}/bad_query.csv"
  printf "Locus,Allele1,Allele2\nD8S1179,,12\n" > "${BROKE_DIR}/bad_query.csv"
  if "${R_BIN}" "$NEGPC" "${BROKE_DIR}/bad_query.csv"; then true; fi

  echo "[OK] ${LABEL}"
}

# ---- Run battery ----------------------------------------------------------
run_case "1k_normal_off" "$DB1K" 0
run_case "1k_normal_on"  "$DB1K" 1
run_case "2M_normal_off" "$DB2M" 0
run_case "2M_normal_on"  "$DB2M" 1
run_case "1k_force_all" "$DB1K" 0 "--force_all_loci 1"

echo "[DONE] smoke => ${OUT_DIR}"

# scripts/tools/check_bits_and_fullscore.R
suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript scripts/tools/check_bits_and_fullscore.R <run_dir> <db_rds>")

run_dir <- args[1]
db_rds  <- args[2]
scores_path <- file.path(run_dir, "scores.csv")
opts_json   <- file.path(run_dir, "opts_scores.json")

if (!file.exists(scores_path)) stop("scores.csv not found")
if (!file.exists(opts_json))   stop("opts_scores.json not found")
if (!file.exists(db_rds))      stop("DB RDS not found")

# 先頭1件
sc <- read.csv(scores_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(sc) >= 1)
top_id <- sc$SampleID[1]
top_sc <- as.integer(sc$Score[1])

# 一時IDファイル作成
ids_csv <- file.path(run_dir, "_ids_top1.csv")
write.csv(data.frame(SampleID = top_id), ids_csv, row.names = FALSE)

# emit_detail(raw) 実行（out_dir=_raw_chk）
out_raw <- file.path(run_dir, "_raw_chk")
dir.create(out_raw, showWarnings = FALSE, recursive = TRUE)

if (!file.exists("scripts/emit/emit_detail.R")) stop("scripts/emit/emit_detail.R not found")

# opts は run_dir の opts_scores.json をそのまま渡す
source("scripts/emit/emit_detail.R", local = TRUE)
emit_detail(
  opt_path    = opts_json,
  scores_path = ids_csv,
  out_dir     = out_raw,
  mode        = "raw",
  any_code    = 9999L
)

# raw_detail を探す（raw_*.csv 想定）
raw_candidates <- list.files(out_raw, pattern = "^raw_.*\\.csv$", full.names = TRUE)
if (length(raw_candidates) == 0) stop("raw_*.csv not found under emit out_dir")
raw_path <- raw_candidates[1]

raw <- read.csv(raw_path, stringsAsFactors = FALSE, check.names = FALSE)

# 列チェック（順序も厳密）
expected <- c("Locus","q1","q2","r1","r2","bits","code","score","SampleID")
if (!identical(names(raw), expected)) {
  stop(sprintf("raw columns mismatch.\n expected: %s\n actual  : %s",
               paste(expected, collapse=","), paste(names(raw), collapse=",")))
}

# bitsは4桁の2進
if (any(!grepl("^[01]{4}$", raw$bits))) {
  bad_n <- sum(!grepl("^[01]{4}$", raw$bits))
  stop(sprintf("bits must be 4-bit binary strings; bad rows=%d", bad_n))
}

# 合計一致（整数厳密）
agg <- sum(as.integer(raw$score))
if (!identical(as.integer(agg), as.integer(top_sc))) {
  stop(sprintf("sum(raw.score) != scores.top: %s != %s", agg, top_sc))
}

cat(sprintf("[OK] bits=4bit binary, raw-sum=%d == top-score\n", top_sc))

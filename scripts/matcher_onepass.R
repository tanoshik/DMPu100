# scripts/matcher_onepass.R
# One-pass matcher CLI wrapper for dmp_match_cpp (schema-strict)
# ASCII-only; JP comments allowed for dev.

suppressPackageStartupMessages({
  library(optparse)
  library(Rcpp)
  library(jsonlite)
})

# Rcpp core
sourceCpp("src/dmp_match.cpp")  # exports: dmp_match_cpp, dmp_hist_cpp

# SCORE_TABLE (DMP spec; index: 0..15)
SCORE_TABLE <- as.integer(c(0,1,1,1, 1,1,2,2, 1,2,1,2, 1,2,2,2))

# CLI options
option_list <- list(
  make_option("--db",      type="character", help="path to virtual_db_u100_*.rds"),
  make_option("--query",   type="character", help="path to query (RDS or CSV)"),
  make_option("--out",     type="character", default="output/scores.csv"),
  make_option("--mode",    type="character", default="topn"),      # topn|threshold
  make_option("--report",  type="character", default="all"),       # all|top|hist  (★追加)
  make_option("--top_n",   type="integer",   default=200L),
  make_option("--threshold", type="integer", default=0L),
  make_option("--display_limit", type="integer", default=200L),
  make_option("--any_code", type="integer",  default=9999L),
  make_option("--idf_csv", type="character", default=NA),
  make_option("--force_all_loci", type="integer", default=0L),
  make_option("--debug",   type="integer",   default=1L),
  make_option("--path",    type="character", default=NA)  # optional
)
opt <- parse_args(OptionParser(option_list=option_list))

# guards
if (is.null(opt$db) || is.null(opt$query)) stop("--db and --query are required")
if (!is.na(opt$path) && !dir.exists(opt$path)) stop(sprintf("--path not found: %s", opt$path))
if (!file.exists(opt$db))    stop(sprintf("db not found: %s", opt$db))
if (!file.exists(opt$query)) stop(sprintf("query not found: %s", opt$query))

t0_total <- proc.time()[["elapsed"]]

# ---- load DB (schema-strict) ----
t0 <- proc.time()[["elapsed"]]
db <- readRDS(opt$db)  # expected: list(sample_ids, locus_ids, A1, A2)
if (!all(c("sample_ids","locus_ids","A1","A2") %in% names(db))) {
  stop("DB rds must contain sample_ids, locus_ids, A1, A2")
}
to_mat <- function(X) {
  if (is.matrix(X)) {
    storage.mode(X) <- "integer"; return(X)
  }
  if (is.list(X)) {
    stopifnot(length(X) == length(db$locus_ids))
    S <- length(X[[1]]); L <- length(X)
    M <- matrix(0L, nrow=S, ncol=L)
    for (j in seq_len(L)) M[,j] <- as.integer(X[[j]])
    colnames(M) <- db$locus_ids
    return(M)
  }
  stop("A1/A2 must be matrix(int) or list-of-int")
}
A1m <- to_mat(db$A1)
A2m <- to_mat(db$A2)
stopifnot(nrow(A1m) == length(db$sample_ids), ncol(A1m) == length(db$locus_ids))
stopifnot(nrow(A2m) == length(db$sample_ids), ncol(A2m) == length(db$locus_ids))
t1 <- proc.time()[["elapsed"]]
load_db_sec <- t1 - t0

# ---- load query (zip現物の実装を厳密踏襲) ----
t0 <- proc.time()[["elapsed"]]
if (grepl("\\.rds$", opt$query, ignore.case = TRUE)) {
  q <- readRDS(opt$query)     # list(q1,q2) aligned to locus_ids
  if (!all(c("q1","q2") %in% names(q))) stop("Query RDS must provide q1, q2")
  if (!(length(q$q1) == length(db$locus_ids) && length(q$q2) == length(db$locus_ids))) {
    stop("q1/q2 length must equal length(db$locus_ids)")
  }
  # RDSは×100済み前提
} else {
  x <- read.csv(opt$query, stringsAsFactors = FALSE, check.names = FALSE)
  nms <- tolower(names(x)); names(x) <- nms
  has_locus <- "locus" %in% names(x)
  ANY <- as.integer(opt$any_code)
  
  norm_scale100 <- function(v) {
    u <- suppressWarnings(as.numeric(v))
    out <- integer(length(u))
    na <- is.na(u)
    out[na | (u <= 0)] <- ANY             # NA/非数/<=0 は ANY
    small <- (!na) & (u > 0) & (u < 100)  # 正の <100 は *100
    out[small] <- as.integer(round(u[small] * 100))
    big <- (!na) & (u >= 100)
    out[big] <- as.integer(round(u[big]))
    out
  }
  sort_pair <- function(a1, a2, ANY) {
    # ANY を右に寄せ、かつ昇順
    swap <- (a1 == ANY & a2 != ANY) | ((a1 != ANY & a2 != ANY) & (a1 > a2))
    a1s <- ifelse(swap, a2, a1)
    a2s <- ifelse(swap, a1, a2)
    list(a1=a1s, a2=a2s)
  }
  
  if (has_locus) {
    key_db <- toupper(trimws(db$locus_ids))
    key_q  <- toupper(trimws(as.character(x$locus)))
    m <- match(key_db, key_q)
    q1 <- rep.int(ANY, length(db$locus_ids))
    q2 <- rep.int(ANY, length(db$locus_ids))
    a1col <- if ("allele1" %in% names(x)) "allele1" else "allele_1"
    a2col <- if ("allele2" %in% names(x)) "allele2" else "allele_2"
    hit <- !is.na(m)
    if (!all(hit)) {
      miss <- db$locus_ids[!hit]
      warning(sprintf("[WARN] %d loci missing in query CSV; filled with ANY_CODE (%s): %s",
                      sum(!hit), as.character(ANY), paste(miss, collapse = ",")))
    }
    if (!all(c(a1col, a2col) %in% names(x))) {
      stop("Query CSV must have Allele1/Allele2 columns")
    }
    a1v <- norm_scale100(x[[a1col]][m])
    a2v <- norm_scale100(x[[a2col]][m])
    ord <- sort_pair(a1v, a2v, ANY)
    q1 <- ord$a1; q2 <- ord$a2
    q <- list(q1 = q1, q2 = q2)
  } else if (all(c("allele1","allele2") %in% names(x)) &&
             nrow(x) == length(db$locus_ids)) {
    a1v <- norm_scale100(x$allele1)
    a2v <- norm_scale100(x$allele2)
    ord <- sort_pair(a1v, a2v, ANY)
    q <- list(q1 = ord$a1, q2 = ord$a2)
  } else {
    stop("Query CSV must include Locus, Allele1, Allele2 (or be pre-aligned with matching length).")
  }
}
t1 <- proc.time()[["elapsed"]]
load_q_sec <- t1 - t0

if (opt$debug) {
  n_any <- sum(q$q1 == as.integer(opt$any_code) | q$q2 == as.integer(opt$any_code))
  n_scaled <- sum(q$q1 < 100L | q$q2 < 100L)  # 正常化後は 0 が期待値
  cat(sprintf("[DBG] query norm: ANY_used=%d, (<100 after norm)=%d\n", n_any, n_scaled))
}

# ---- idf mask bits ----
if (!is.na(opt$idf_csv) && file.exists(opt$idf_csv)) {
  idf <- read.csv(opt$idf_csv, stringsAsFactors = FALSE)
  j <- match(db$locus_ids, idf$locus)
  enabled <- ifelse(is.na(j), FALSE, idf$enabled[j] != 0)
} else if (as.integer(opt$force_all_loci) != 0L) {
  enabled <- rep(TRUE, length(db$locus_ids))
} else {
  enabled <- rep(TRUE, length(db$locus_ids))  # default: all enabled
}
idf_mask_bits <- 0L
for (j in seq_along(db$locus_ids)) {
  if (enabled[j]) idf_mask_bits <- bitwOr(idf_mask_bits, bitwShiftL(1L, (j-1L) %% 32L))
}

# ---- run core ----
opts_cpp <- list(
  any_code = as.integer(opt$any_code),
  mode = as.character(opt$mode),
  top_n = as.integer(opt$top_n),
  threshold = as.integer(opt$threshold),
  display_limit = as.integer(opt$display_limit),
  force_all_loci = as.logical(opt$force_all_loci),
  h2a_Q = FALSE, h2a_R = FALSE
)

t0 <- proc.time()[["elapsed"]]
res <- dmp_match_cpp(
  A1 = A1m, A2 = A2m,
  q1 = q$q1, q2 = q$q2,
  score_table = SCORE_TABLE,
  idf_mask_bits = as.integer(idf_mask_bits),
  opts = opts_cpp,
  sample_ids = db$sample_ids
)
t1 <- proc.time()[["elapsed"]]
comp_sec <- t1 - t0

# ---- outputs ----
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)

rep_mode <- tolower(opt$report)
if (!rep_mode %in% c("all","top","hist")) stop("--report must be all|top|hist")

# scores.csv
if (rep_mode %in% c("all","top")) {
  write.csv(res, opt$out, row.names = FALSE)
  if (opt$debug) cat(sprintf("[OK] wrote %s (%d rows)\n", opt$out, nrow(res)))
}

# histogram (Score,Count; 0..MAX_SC 全行) to hist_*.csv
# 例: out=scores_t18A.csv → hist_t18A.csv （threshold>0 のとき tNN 前置）
base <- basename(opt$out); dirp <- dirname(opt$out)
hist_name <- sub("^scores", "hist", base)
if (as.integer(opt$threshold) > 0L) {
  hist_name <- sub("^hist", sprintf("hist_t%d", as.integer(opt$threshold)), hist_name)
}
hist_path <- file.path(dirp, hist_name)

if (rep_mode %in% c("all","hist")) {
  t0h <- proc.time()[["elapsed"]]
  hdf <- dmp_hist_cpp(
    A1 = A1m, A2 = A2m,
    q1 = q$q1, q2 = q$q2,
    score_table = SCORE_TABLE,
    idf_mask_bits = as.integer(idf_mask_bits),
    opts = list(any_code = as.integer(opt$any_code), force_all_loci = as.logical(opt$force_all_loci))
  )
  t1h <- proc.time()[["elapsed"]]
  write.csv(hdf, hist_path, row.names = FALSE)
  if (rep_mode == "all") comp_sec <- comp_sec + (t1h - t0h)
  if (opt$debug) cat(sprintf("[OK] wrote %s (%d rows)\n", hist_path, nrow(hdf)))
}

# ---- opts.json & time.log ----
opts_list <- list(
  version = "0.1.3",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"),
  args = as.list(opt),
  run = list(S = nrow(A1m), L = ncol(A1m))
)
json_path <- file.path(dirname(opt$out), "opts.json")
writeLines(jsonlite::toJSON(opts_list, auto_unbox = TRUE, pretty = TRUE), json_path)
if (opt$debug) cat(sprintf("[OK] wrote %s\n", json_path))

header <- "LOAD_DB_SEC,LOAD_Q_SEC,COMP_SEC,TOTAL_SEC,PEAK_MiB"
total_sec <- proc.time()[["elapsed"]] - t0_total
line <- sprintf("%.3f,%.3f,%.3f,%.3f,%.0f", load_db_sec, load_q_sec, comp_sec, total_sec, NA_real_)
logp <- file.path(dirname(opt$out), "time.log")
if (!file.exists(logp)) writeLines(header, logp)
cat(paste0(line, "\n"), file = logp, append = TRUE)
if (opt$debug) cat("[DBG] done\n")

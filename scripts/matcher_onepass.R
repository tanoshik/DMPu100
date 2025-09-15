# scripts/matcher_onepass.R
# One-pass matcher CLI wrapper for dmp_match_cpp/dmp_hist_cpp (schema-strict)
# ASCII-only; JP comments allowed for dev.

suppressPackageStartupMessages({
  library(optparse)
  library(Rcpp)
  library(jsonlite)
})

# Rcpp core
sourceCpp("src/dmp_match.cpp")  # exports: dmp_match_cpp, dmp_hist_cpp

# SCORE_TABLE (DMP spec; index: 0..15)
SCORE_TABLE <- as.integer(c(
  0,1,1,1,
  1,1,2,2,
  1,2,1,2,
  1,2,2,2
))

# CLI options
option_list <- list(
  make_option("--db",      type="character", help="path to virtual_db_u100_*.rds"),
  make_option("--query",   type="character", help="path to query (RDS or CSV)"),
  make_option("--out",     type="character", default="output/scores.csv"),
  make_option("--report",  type="character", default="all"),       # all|top|hist
  make_option("--n_cap",   type="integer",   default=200L),        # scores max rows (TopN or Threshold cap)
  make_option("--score_min", type="integer", default=NA),          # Threshold if specified; else TopN
  make_option("--any_code",  type="integer", default=9999L),
  make_option("--idf_csv",   type="character", default=NA),
  make_option("--force_all_loci", type="integer", default=0L),
  make_option("--debug",   type="integer",   default=1L),
  make_option("--path",    type="character", default=NA)           # optional
)
opt <- parse_args(OptionParser(option_list=option_list))

# guards
if (is.null(opt$db) || is.null(opt$query)) stop("--db and --query are required")
if (!is.na(opt$path) && !dir.exists(opt$path)) stop(sprintf("--path not found: %s", opt$path))
if (!file.exists(opt$db))    stop(sprintf("db not found: %s", opt$db))
if (!file.exists(opt$query)) stop(sprintf("query not found: %s", opt$query))
if (as.integer(opt$n_cap) < 0L) stop("--n_cap must be >= 0")

t0_total <- proc.time()[["elapsed"]]

# ---- load DB (schema-strict) ----
t0 <- proc.time()[["elapsed"]]
db <- readRDS(opt$db)
# schema guards
for (nm in c("sample_ids","locus_ids","A1","A2")) {
  if (!nm %in% names(db)) stop(sprintf("db.rds missing '%s'", nm))
}
# normalize A1/A2 to integer matrices SxL
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

# ---- load query (zip-spec exact behavior) ----
t0 <- proc.time()[["elapsed"]]
if (grepl("\\.rds$", opt$query, ignore.case = TRUE)) {
  q <- readRDS(opt$query)     # list(q1,q2) aligned to locus_ids
  if (!all(c("q1","q2") %in% names(q))) stop("Query RDS must provide q1, q2")
  if (!(length(q$q1) == length(db$locus_ids) && length(q$q2) == length(db$locus_ids))) {
    stop("q1/q2 length must equal length(db$locus_ids)")
  }
  # RDS values are pre-scaled x100
} else {
  x <- read.csv(opt$query, stringsAsFactors = FALSE, check.names = FALSE)
  # normalize column names (tolower) and pick locus / allele1 / allele2
  nms <- tolower(names(x))
  names(x) <- nms
  ANY <- as.integer(opt$any_code)
  norm_scale100 <- function(v){
    as.integer(round(as.numeric(v) * 100))
  }
  sort_pair <- function(a1, a2, ANY) {
    a1 <- as.integer(a1); a2 <- as.integer(a2)
    swap <- (a1 == ANY & a2 != ANY) | (a1 > a2 & a2 != ANY)
    a1n <- ifelse(swap, a2, a1)
    a2n <- ifelse(swap, a1, a2)
    list(a1 = a1n, a2 = a2n)
  }
  if (all(c("locus","allele1","allele2") %in% nms)) {
    m <- match(db$locus_ids, x$locus)
    hit <- !is.na(m)
    if (!all(hit)) {
      miss <- db$locus_ids[!hit]
      stop(sprintf("Query CSV missing loci: %d (ANY=%d) %s", sum(!hit), as.character(ANY), paste(miss, collapse=",")))
    }
    a1v <- norm_scale100(x$allele1[m])
    a2v <- norm_scale100(x$allele2[m])
    ord <- sort_pair(a1v, a2v, ANY)
    q1 <- ord$a1; q2 <- ord$a2
    q <- list(q1 = q1, q2 = q2)
  } else if (all(c("allele1","allele2") %in% nms) &&
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
  n_scaled <- sum(q$q1 < 100L | q$q2 < 100L)  # after normalization, expect 0
  cat(sprintf("[DBG] query norm: ANY_used=%d, (<100 after norm)=%d\n", n_any, n_scaled))
}

# ---- idf mask bits ----
idf_mask_bits <- 0L
if (!is.na(opt$idf_csv) && file.exists(opt$idf_csv)) {
  idf <- read.csv(opt$idf_csv, stringsAsFactors = FALSE)
  loc <- tolower(db$locus_ids); names(loc) <- db$locus_ids
  w <- match(tolower(idf$locus), loc)
  w <- w[!is.na(w)]
  if (length(w) > 0) {
    for (j in w) {
      idf_mask_bits <- bitwOr(idf_mask_bits, bitwShiftL(1L, as.integer(j - 1L)))
    }
  }
}
if (as.integer(opt$force_all_loci) == 1L) {
  idf_mask_bits <- bitwShiftL(1L, length(db$locus_ids)) - 1L
}

# ---- run core ----
rep_mode <- tolower(opt$report)
compute_scores <- (rep_mode %in% c("all","top")) && (as.integer(opt$n_cap) > 0L)
comp_sec <- 0

res <- data.frame(SampleID=character(0), Score=integer(0))
if (compute_scores) {
  t0 <- proc.time()[["elapsed"]]
  opts_cpp <- list(
    any_code = as.integer(opt$any_code),
    score_min = if (is.na(opt$score_min)) NA_integer_ else as.integer(opt$score_min),
    n_cap = as.integer(opt$n_cap),
    force_all_loci = as.logical(opt$force_all_loci),
    h2a_Q = FALSE, h2a_R = FALSE
  )
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
}

# ---- scores output (sorted: Score desc > SampleID asc) ----
if (rep_mode %in% c("all","top")) {
  if (nrow(res) > 0) {
    res$Score <- as.integer(res$Score)
    res$SampleID <- as.character(res$SampleID)
    ord <- order(-res$Score, res$SampleID)
    res <- res[ord, , drop = FALSE]
  }
  dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
  write.csv(res, opt$out, row.names = FALSE)
  if (opt$debug) cat(sprintf("[OK] wrote %s (%d rows)\n", opt$out, nrow(res)))
}

# ---- histogram (Score,Count; 0..MAX_SC full) to hist_*.csv ----
base <- basename(opt$out); dirp <- dirname(opt$out)
hist_name <- sub("^scores", "hist", base)
if (!is.na(opt$score_min) && as.integer(opt$score_min) > 0L) {
  hist_name <- sub("^hist", sprintf("hist_t%d", as.integer(opt$score_min)), hist_name)
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

# ---- opts_<core>.json & time.log ----
# NOTE: rename from "<stem>_opts.json" to "opts_<core>.json"
# core := basename(out) without leading "scores_" and without extension
stem <- tools::file_path_sans_ext(basename(opt$out))
core <- sub("^scores_", "", stem)
json_name <- paste0("opts_", core, ".json")
json_path <- file.path(dirname(opt$out), json_name)

opts_list <- list(
  version = "0.1.4",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"),
  args = as.list(opt),
  run = list(S = nrow(A1m), L = ncol(A1m))
)

writeLines(jsonlite::toJSON(opts_list, auto_unbox = TRUE, pretty = TRUE), json_path)
if (opt$debug) cat(sprintf("[OK] wrote %s\n", json_path))

header <- "LOAD_DB_SEC,LOAD_Q_SEC,COMP_SEC,TOTAL_SEC,PEAK_MiB"
total_sec <- proc.time()[["elapsed"]] - t0_total
line <- sprintf("%.3f,%.3f,%.3f,%.3f,%.0f", load_db_sec, load_q_sec, comp_sec, total_sec, NA_real_)
logp <- file.path(dirname(opt$out), "time.log")
if (!file.exists(logp)) writeLines(header, logp)
cat(paste0(line, "\n"), file = logp, append = TRUE)
if (opt$debug) cat("[DBG] done\n")

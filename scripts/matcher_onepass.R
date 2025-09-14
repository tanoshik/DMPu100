# scripts/matcher_onepass.R
# One-pass matcher CLI wrapper for dmp_match_cpp (schema-strict)
# ASCII-only (public dist); JP comments allowed here for dev.

suppressPackageStartupMessages({
  library(optparse)
  library(Rcpp)
  library(jsonlite)
})

# Rcpp core
sourceCpp("src/dmp_match.cpp")

# SCORE_TABLE (DMP spec; index: 0..15)
SCORE_TABLE <- as.integer(c(0,1,1,1, 1,1,2,2, 1,2,1,2, 1,2,2,2))

# CLI options
option_list <- list(
  make_option("--db",      type="character", help="path to virtual_db_u100_*.rds"),
  make_option("--query",   type="character", help="path to query (RDS or CSV)"),
  make_option("--out",     type="character", default="output/scores.csv"),
  make_option("--mode",    type="character", default="topn"),      # topn|threshold
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

# entrance guard
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

# allow A1/A2 list-of-columns -> convert to SxL integer matrix
to_mat <- function(X) {
  if (is.matrix(X)) {
    storage.mode(X) <- "integer"
    return(X)
  }
  if (is.list(X)) {
    stopifnot(length(X) == length(db$locus_ids))
    S <- length(X[[1]])
    L <- length(X)
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

# ---- load query (align by locus_ids; normalize: 0/NA->ANY, <100 -> *100, sort asc & ANY-right) ----
t0 <- proc.time()[["elapsed"]]
if (grepl("\\.rds$", opt$query, ignore.case = TRUE)) {
    q <- readRDS(opt$query)     # list(q1,q2) aligned to locus_ids
    if (!all(c("q1","q2") %in% names(q))) stop("Query RDS must provide q1, q2")
    if (!(length(q$q1) == length(db$locus_ids) && length(q$q2) == length(db$locus_ids))) {
        stop("q1/q2 length must equal length(db$locus_ids)")
    }
    # RDSは生成段階で×100済み前提なので無加工で通す
  } else {
    x <- read.csv(opt$query, stringsAsFactors = FALSE, check.names = FALSE)
    nms <- tolower(names(x)); names(x) <- nms
    has_locus <- "locus" %in% names(x)
    ANY <- as.integer(opt$any_code)
    norm_scale100 <- function(v) {
      u <- suppressWarnings(as.numeric(v))
      out <- integer(length(u))
      na <- is.na(u)
      # NA/非数/<=0 は ANY
      out[na | (u <= 0)] <- ANY
      # 正の <100 は *100 して丸め
      small <- (!na) & (u > 0) & (u < 100)
      out[small] <- as.integer(round(u[small] * 100))
      # 100 以上は整数化（すでに *100 済み想定）
      big <- (!na) & (u >= 100)
      out[big] <- as.integer(round(u[big]))
      out
    }
    sort_pair <- function(a1, a2, ANY) {
      # ANY を右に寄せ、かつ昇順 (ascending_alleles & any_right_aligned)
      swap <- (a1 == ANY & a2 != ANY) | ((a1 != ANY & a2 != ANY) & (a1 > a2))
      a1s <- ifelse(swap, a2, a1)
      a2s <- ifelse(swap, a1, a2)
      list(a1=a1s, a2=a2s)
    }
          if (has_locus) {
            key_db <- toupper(trimws(db$locus_ids))
            key_q  <- toupper(trimws(as.character(x$locus)))  
            m <- match(key_db, key_q)       # index in query rows for each db locus
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
    n_any <- sum(q$q1 == ANY | q$q2 == ANY)
    n_scaled <- sum(q$q1 < 100L | q$q2 < 100L)  # 正常化後は 0 が期待値
    cat(sprintf("[DBG] query norm: ANY_used=%d, (<100 after norm)=%d\n", n_any, n_scaled))
  }
  
# ---- IDF mask ----
idf_mask_bits <- as.integer(0x1FFFFF) # 21 loci enabled by default
if (!is.na(opt$idf_csv) && nzchar(opt$idf_csv)) {
  idf <- read.csv(opt$idf_csv, stringsAsFactors = FALSE)
  bits <- as.integer(idf$enable[1:21]); bits[is.na(bits)] <- 1L
  val <- 0L
  for (i in seq_len(21L)) if (bits[i] != 0L) val <- bitwOr(val, bitwShiftL(1L, i-1L))
  idf_mask_bits <- as.integer(val)
}

if (opt$debug)
  cat(sprintf("[DBG] SxL=%dx%d idf=0x%X force_all=%d\n",
              nrow(A1m), ncol(A1m), as.integer(idf_mask_bits), as.integer(opt$force_all_loci)))

# ---- run core (no homo_mask) ----
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
write.csv(res, opt$out, row.names = FALSE)
cat(sprintf("[OK] wrote %s (%d rows)\n", opt$out, nrow(res)))

# opts.json (with version)
opts_list <- list(
  version = "0.1.0",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"),
  args = as.list(opt),
  run = list(S = nrow(A1m), L = ncol(A1m))
)
json_path <- file.path(dirname(opt$out), "opts.json")
writeLines(jsonlite::toJSON(opts_list, auto_unbox = TRUE, pretty = TRUE), json_path)
cat(sprintf("[OK] wrote %s\n", json_path))

# time log
get_mem_mib <- function() {
  mib <- NA_real_
  if (requireNamespace("pryr", quietly = TRUE)) mib <- as.numeric(pryr::mem_used()) / (1024^2)
  mib
}
peak_mib <- get_mem_mib()
total_sec <- proc.time()[["elapsed"]] - t0_total
timelog <- sprintf("LOAD_DB_SEC,LOAD_Q_SEC,COMP_SEC,TOTAL_SEC,PEAK_MiB\n%.3f,%.3f,%.3f,%.3f,%.1f",
                   load_db_sec, load_q_sec, comp_sec, total_sec, ifelse(is.na(peak_mib), -1, peak_mib))
time_path <- file.path(dirname(opt$out), "time.log")
writeLines(timelog, time_path)
cat(sprintf("[OK] wrote %s\n", time_path))

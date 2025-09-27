#!/usr/bin/env Rscript
# scripts/matcher_onepass.R
# One-pass matcher CLI wrapper for dmp_match_cpp/dmp_hist_cpp (strict)
# - strict I/O
# - snake_case CLI
# - TRACE LOG to output/_trace_matcher.log

if (identical(Sys.getenv("DMP_DBG"), "1")) options(dmp_dbg=TRUE)

# ---- TRACE (always on) ----
TRACE <- "output/_trace_matcher.log"
try({ dir.create(dirname(TRACE), recursive=TRUE, showWarnings=FALSE) }, silent=TRUE)
cat(sprintf("[TRACE] start %s\n", format(Sys.time(), "%F %T")), file=TRACE, append=TRUE)
.on_end <- function(){ cat(sprintf("[TRACE] end   %s\n", format(Sys.time(), "%F %T")), file=TRACE, append=TRUE) }
.on_err <- function(e){
  cat(sprintf("[TRACE] ERROR: %s\n", conditionMessage(e)), file=TRACE, append=TRUE)
  q(status=1, save="no")
}
options(error = .on_err)
on.exit(.on_end(), add=TRUE)

suppressPackageStartupMessages({
  library(optparse);  cat("[TRACE] loaded: optparse\n", file=TRACE, append=TRUE)
  library(Rcpp);      cat("[TRACE] loaded: Rcpp\n",     file=TRACE, append=TRUE)
  library(jsonlite);  cat("[TRACE] loaded: jsonlite\n", file=TRACE, append=TRUE)
  library(tools);     cat("[TRACE] loaded: tools\n",    file=TRACE, append=TRUE)
})

# ---- SCORE_TABLE (index: 0..15) ----
SCORE_TABLE <- as.integer(c(
  0,1,1,1,
  1,1,2,2,
  1,2,1,2,
  1,2,2,2
))

# ---- CLI options (snake_case) ----
option_list <- list(
  make_option(c("-d","--db"),              type="character", help="path to virtual_db_u100_*.rds"),
  make_option(c("-q","--query"),           type="character", help="path to query (RDS or CSV)"),
  make_option(c("-o","--out"),             type="character", default="output/scores.csv"),
  make_option(c("-r","--report"),          type="character", default="all"),   # all|top|hist
  make_option(c("-N","--n_cap"),           type="integer",   default=200L),
  make_option(c("-t","--score_min"),       type="integer",   default=0L),      # >=0
  make_option(c("-a","--any_code"),        type="integer",   default=9999L),
  make_option(c("-M","--idf_csv"),         type="character", default=NA),      # optional
  make_option(c("-F","--force_all_loci"),  type="integer",   default=0L),      # 1: ignore idf mask
  make_option(c("-H","--h2a_on"),          type="integer",   default=0L),      # 1: DB homo -> RIGHT ANY
  make_option(c("-v","--debug"),           type="integer",   default=1L),
  make_option(c("-C","--path"),            type="character", default=NA)       # optional setwd
)
opt <- parse_args(OptionParser(option_list=option_list))
cat(sprintf("[TRACE] argv parsed; report=%s n_cap=%s\n", opt$report, as.character(opt$n_cap)), file=TRACE, append=TRUE)

# ---- guards ----
if (is.null(opt$db) || is.null(opt$query)) stop("--db and --query are required")
if (!is.na(opt$path)) {
  if (!dir.exists(opt$path)) stop(sprintf("--path not found: %s", opt$path))
  setwd(normalizePath(opt$path, winslash="/"))
  cat(sprintf("[TRACE] setwd=%s\n", getwd()), file=TRACE, append=TRUE)
}
if (!file.exists(opt$db))    stop(sprintf("db not found: %s", opt$db))
if (!file.exists(opt$query)) stop(sprintf("query not found: %s", opt$query))
if (as.integer(opt$n_cap) < 0L) stop("--n_cap must be >= 0")
opt$score_min <- as.integer(opt$score_min)
if (is.na(opt$score_min) || opt$score_min < 0L) stop("--score_min must be an integer >= 0")

# ---- helpers ----
ensure_csv  <- function(p) if (grepl("\\.csv$",  p, ignore.case=TRUE)) p else paste0(p, ".csv")
ensure_json <- function(p) if (grepl("\\.json$", p, ignore.case=TRUE)) p else paste0(p, ".json")

# ---- compile C++ core ----
cat("[TRACE] before sourceCpp\n", file=TRACE, append=TRUE)
tryCatch({
  sourceCpp("src/dmp_match.cpp")  # exports: dmp_match_cpp, dmp_hist_cpp
  cat("[TRACE] after  sourceCpp\n", file=TRACE, append=TRUE)
}, error = .on_err)

t0_total <- proc.time()[["elapsed"]]

# ---- load DB ----
cat("[TRACE] before readRDS(db)\n", file=TRACE, append=TRUE)
db <- readRDS(opt$db)  # list(sample_ids, locus_ids, A1, A2)
cat("[TRACE] after  readRDS(db)\n", file=TRACE, append=TRUE)

if (!is.list(db) || !all(c("sample_ids","locus_ids","A1","A2") %in% names(db))) {
  stop("DB RDS must contain list(sample_ids, locus_ids, A1, A2)")
}
sample_ids <- as.character(db$sample_ids)
sample_ids[is.na(sample_ids)] <- ""  # stabilize comparator
locus_ids  <- as.character(db$locus_ids)
S <- length(sample_ids); L <- length(locus_ids)
cat(sprintf("[TRACE] S=%d L=%d\n", S, L), file=TRACE, append=TRUE)

A1m <- as.matrix(db$A1); storage.mode(A1m) <- "integer"
A2m <- as.matrix(db$A2); storage.mode(A2m) <- "integer"
if (nrow(A1m)!=S || ncol(A1m)!=L || nrow(A2m)!=S || ncol(A2m)!=L) stop("DB A1/A2 shape mismatch")

# ---- h2a (DB homo -> RIGHT ANY) + timing ----
t_h2a <- NA_real_
if (as.integer(opt$h2a_on) == 1L) {
  t0_h2a <- proc.time()[["elapsed"]]
  ANY <- as.integer(opt$any_code)
  mask <- (A1m == A2m) & (A1m != ANY)
  A2m[mask] <- ANY
  t_h2a <- proc.time()[["elapsed"]] - t0_h2a
  cat(sprintf("[TRACE] h2a_on masked=%d (H2A_SEC=%.3f)\n", sum(mask, na.rm=TRUE), t_h2a), file=TRACE, append=TRUE)
} else {
  t_h2a <- 0.0
}


# ---- load query (RDS or CSV, strict) ----
to_pair <- function(a1, a2, ANY) {
  a1 <- as.integer(a1); a2 <- as.integer(a2)
  swap <- ((a1 == ANY & a2 != ANY) | (a1 > a2 & a2 != ANY))
  list(a1 = ifelse(swap, a2, a1), a2 = ifelse(swap, a1, a2))
}
ANY <- as.integer(opt$any_code)

cat("[TRACE] before read query\n", file=TRACE, append=TRUE)
q <- NULL
if (grepl("\\.(rds|RDS)$", opt$query)) {
  qx <- readRDS(opt$query)
  if (!is.list(qx) || !all(c("q1","q2") %in% names(qx))) stop("Query RDS must provide list(q1, q2)")
  if (!(length(qx$q1)==L && length(qx$q2)==L)) stop("query length mismatch")
  qp <- to_pair(qx$q1, qx$q2, ANY); q <- list(q1 = qp$a1, q2 = qp$a2)
} else {
  x <- read.csv(opt$query, stringsAsFactors=FALSE, check.names=FALSE)
  nms <- tolower(names(x)); names(x) <- nms
  lc  <- which(nms %in% c("locus","marker"))[1]
  a1c <- which(nms %in% c("allele1","a1","q1"))[1]
  a2c <- which(nms %in% c("allele2","a2","q2"))[1]
  if (is.na(lc) || is.na(a1c) || is.na(a2c)) stop("CSV needs columns: Locus(+), Allele1/A1/Q1, Allele2/A2/Q2")
  mm <- match(locus_ids, x[[lc]])
  if (any(is.na(mm))) stop("CSV loci must cover all db$locus_ids")
  as_i100 <- function(v) {
    num <- suppressWarnings(as.numeric(v))
    if (any(is.na(num))) stop("query CSV contains non-numeric allele(s) under strict mode")
    as.integer(round(num*100))
  }
  qp <- to_pair(as_i100(x[[a1c]][mm]), as_i100(x[[a2c]][mm]), ANY)
  q <- list(q1 = qp$a1, q2 = qp$a2)
}
cat("[TRACE] after  read query\n", file=TRACE, append=TRUE)

# ---- idf mask bits (uint32: bit j=1 -> enabled) ----
make_all_ones <- function(L) {
  bits <- 0L
  maxj <- min(L-1L, 31L)  # defensive: limit to 32 bits
  for (j in 0:maxj) bits <- bitwOr(bits, bitwShiftL(1L, j))
  bits
}
idf_mask_bits <- make_all_ones(L)
mask_name <- "normal"

if (!is.na(opt$idf_csv) && nzchar(opt$idf_csv) && file.exists(opt$idf_csv)) {
  mask_name <- "gf_idf_mask"
  dfm <- read.csv(opt$idf_csv, stringsAsFactors=FALSE, check.names=FALSE)
  nms <- tolower(names(dfm)); names(dfm) <- nms
  lc <- which(nms %in% c("locus","marker","locus_id"))[1]
  bc <- which(nms %in% c("bit","mask","enabled","enable","use"))[1]
  if (!is.na(lc) && !is.na(bc)) {
    idf_mask_bits <- 0L
    mm <- match(locus_ids, dfm[[lc]])
    maxj <- min(L-1L, 31L)
    for (j in 0:maxj) {
      v <- 1L
      if (!is.na(mm[j+1])) {
        vv <- suppressWarnings(as.integer(dfm[[bc]][ mm[j+1] ]))
        if (!is.na(vv)) v <- ifelse(vv != 0L, 1L, 0L)
      }
      if (v == 1L) idf_mask_bits <- bitwOr(idf_mask_bits, bitwShiftL(1L, j))
    }
  } else {
    warning("idf_csv lacks recognizable columns; falling back to all-ones")
    idf_mask_bits <- make_all_ones(L)
    mask_name <- "normal"
  }
}
if (as.integer(opt$force_all_loci) == 1L) {
  idf_mask_bits <- make_all_ones(L)
  mask_name <- "normal"
}
cat(sprintf("[TRACE] mask_name=%s bits=0x%08X\n", mask_name, as.integer(idf_mask_bits)), file=TRACE, append=TRUE)

# ---- run core ----
rep_mode <- tolower(opt$report)
compute_scores <- rep_mode %in% c("all","top")
compute_hist   <- rep_mode %in% c("all","hist")

# ---- resolve outputs (tag or path) ----
ext <- tolower(tools::file_ext(opt$out))
if (nzchar(ext)) {
  # path mode (backward compatible)
  dir.create(dirname(opt$out), recursive=TRUE, showWarnings=FALSE)
  scores_csv <- ensure_csv(opt$out)
  out_tag    <- tools::file_path_sans_ext(basename(scores_csv))
  out_dir    <- dirname(scores_csv)
} else {
  # tag mode (short)
  out_tag <- opt$out
  out_dir <- file.path("output", out_tag)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  scores_csv <- file.path(out_dir, sprintf("scores_%s.csv", out_tag))
}
hist_csv  <- file.path(out_dir, sprintf("hist_%s.csv", out_tag))
opts_json <- file.path(out_dir, sprintf("opts_scores_%s.json", out_tag))

if (compute_scores) {
  cat("[TRACE] before dmp_match_cpp\n", file=TRACE, append=TRUE)
  res <- dmp_match_cpp(
    A1 = A1m, A2 = A2m, q1 = q$q1, q2 = q$q2,
    score_table   = SCORE_TABLE,
    idf_mask_bits = as.integer(idf_mask_bits),
    opts = list(
      any_code       = as.integer(opt$any_code),
      score_min      = as.integer(opt$score_min),
      n_cap          = as.integer(opt$n_cap),
      force_all_loci = as.logical(opt$force_all_loci)
    ),
    sample_ids = sample_ids
  )
  cat("[TRACE] after  dmp_match_cpp\n", file=TRACE, append=TRUE)
  write.csv(res, scores_csv, row.names=FALSE)
  if (as.integer(opt$debug)==1L) cat(sprintf("[OK] wrote %s (%d rows)\n", scores_csv, nrow(res)))
}

if (compute_hist) {
  cat("[TRACE] before dmp_hist_cpp\n", file=TRACE, append=TRUE)
  hdf <- dmp_hist_cpp(
    A1 = A1m, A2 = A2m, q1 = q$q1, q2 = q$q2,
    score_table   = SCORE_TABLE,
    idf_mask_bits = as.integer(idf_mask_bits),
    opts = list(any_code = as.integer(opt$any_code),
                force_all_loci = as.logical(opt$force_all_loci))
  )
  cat("[TRACE] after  dmp_hist_cpp\n", file=TRACE, append=TRUE)
  write.csv(hdf, hist_csv, row.names=FALSE)
  if (as.integer(opt$debug)==1L) cat(sprintf("[OK] wrote %s (%d rows)\n", hist_csv, nrow(hdf)))
}

opts_list <- list(
  db_path        = normalizePath(opt$db, winslash="/"),
  query_path     = normalizePath(opt$query, winslash="/"),
  mask_name      = mask_name,
  any_code       = as.integer(opt$any_code),
  force_all_loci = as.logical(opt$force_all_loci),
  h2a_on         = as.logical(opt$h2a_on)
)
writeLines(jsonlite::toJSON(opts_list, auto_unbox=TRUE, pretty=TRUE), opts_json)
if (as.integer(opt$debug)==1L) cat(sprintf("[OK] wrote %s\n", opts_json))

# ---- time.log ----
core_name <- sprintf("GF/%s/0.0/%s/%d", mask_name, if (opt$h2a_on==1L) "on" else "off", S)
total_sec <- proc.time()[["elapsed"]] - t0_total
header <- "LOAD_DB_SEC,LOAD_Q_SEC,H2A_SEC,COMP_SEC,TOTAL_SEC,PEAK_MiB,name"
line   <- sprintf("NA,NA,%.3f,NA,%.3f,NA,%s", t_h2a, total_sec, core_name)
logp   <- file.path(dirname(scores_csv), sprintf("time_%s.log", out_tag))
if (!file.exists(logp)) writeLines(header, logp)
cat(paste0(line, "\n"), file=logp, append=TRUE)
cat("[TRACE] done\n", file=TRACE, append=TRUE)

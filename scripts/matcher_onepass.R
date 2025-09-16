# scripts/matcher_onepass.R
# One-pass matcher CLI wrapper for dmp_match_cpp/dmp_hist_cpp (schema-strict)
# ASCII-only; JP comments allowed for dev.

suppressPackageStartupMessages({
  library(optparse)
  library(Rcpp)
  library(jsonlite)
})

# Ensure peakRAM (install if missing; otherwise PEAK_MiB will be NA)
if (!requireNamespace("peakRAM", quietly=TRUE)) {
  tryCatch(install.packages("peakRAM", repos="https://cran.rstudio.com"),
           error=function(e) { message("[WARN] peakRAM install failed; PEAK_MiB will be NA") })
}

# Rcpp core
sourceCpp("src/dmp_match.cpp")  # exports: dmp_match_cpp, dmp_hist_cpp

# SCORE_TABLE (DMP spec; index: 0..15)
SCORE_TABLE <- as.integer(c(
  0,1,1,1,
  1,1,2,2,
  1,2,1,2,
  1,2,2,2
))

# ---- CLI options ----
option_list <- list(
  make_option("--db",      type="character", help="path to virtual_db_u100_*.rds"),
  make_option("--query",   type="character", help="path to query (RDS or CSV)"),
  make_option("--out",     type="character", default="output/scores.csv"),
  make_option("--report",  type="character", default="all"),       # all|top|hist
  make_option("--n_cap",   type="integer",   default=200L),        # TopN/Threshold cap (>=0 で有効; 0=無制限)
  make_option("--score_min", type="integer", default=0L),          # 閾値 (0=無制限/閾値無効)
  make_option("--any_code",  type="integer", default=9999L),
  make_option("--idf_csv",   type="character", default=NA),
  make_option("--force_all_loci", type="integer", default=0L),
  make_option("--debug",   type="integer",   default=1L),
  make_option("--path",    type="character", default=NA)           # optional
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---- normalize/guards ----
if (is.null(opt$db) || is.null(opt$query)) stop("--db and --query are required")
if (!is.na(opt$path) && !dir.exists(opt$path)) stop(sprintf("--path not found: %s", opt$path))
if (!file.exists(opt$db))    stop(sprintf("db not found: %s", opt$db))
if (!file.exists(opt$query)) stop(sprintf("query not found: %s", opt$query))
if (as.integer(opt$n_cap) < 0L) stop("--n_cap must be >= 0")
opt$score_min <- as.integer(opt$score_min)
if (is.na(opt$score_min) || opt$score_min < 0L) stop("--score_min must be an integer >= 0")

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
    storage.mode(X) <- "integer"
    return(X)
  }
  if (is.data.frame(X)) X <- as.matrix(X)
  if (is.vector(X) && length(X) == length(db$sample_ids) * length(db$locus_ids)) {
    X <- matrix(X, nrow=length(db$sample_ids), ncol=length(db$locus_ids), byrow=FALSE)
  }
  storage.mode(X) <- "integer"
  X
}
A1m <- to_mat(db$A1)
A2m <- to_mat(db$A2)
stopifnot(nrow(A1m) == length(db$sample_ids), ncol(A1m) == length(db$locus_ids))
stopifnot(nrow(A2m) == length(db$sample_ids), ncol(A2m) == length(db$locus_ids))
t1 <- proc.time()[["elapsed"]]
load_db_sec <- t1 - t0

# ---- load query (RDS or CSV; zip-spec exact behavior) ----
t0 <- proc.time()[["elapsed"]]
if (grepl("\\.rds$", opt$query, ignore.case = TRUE)) {
  q <- readRDS(opt$query)     # list(q1,q2) aligned to locus_ids
  if (!all(c("q1","q2") %in% names(q))) stop("Query RDS must provide q1, q2")
  if (!(length(q$q1) == length(db$locus_ids) && length(q$q2) == length(db$locus_ids))) {
    stop("q1/q2 length must equal length(db$locus_ids)")
  }
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
  # columns: Locus, Allele1, Allele2 (case-insensitive)
  lc <- which(nms %in% c("locus","marker"))[1]
  a1 <- which(nms %in% c("allele1","a1","q1"))[1]
  a2 <- which(nms %in% c("allele2","a2","q2"))[1]
  if (is.na(lc) || is.na(a1) || is.na(a2)) stop("CSV needs columns: Locus, Allele1, Allele2")
  # align to db$locus_ids order
  mm <- match(db$locus_ids, x[[lc]])
  if (any(is.na(mm))) stop("CSV loci must cover db$locus_ids")
  xp <- sort_pair(norm_scale100(x[[a1]][mm]), norm_scale100(x[[a2]][mm]), ANY)
  q <- list(q1 = xp$a1, q2 = xp$a2)
}
t1 <- proc.time()[["elapsed"]]
load_q_sec <- t1 - t0

# ---- idf mask (from kit CSV) -> 21bit (1=enabled) ----
idf_mask_bits <- 0L
if (!is.na(opt$idf_csv) && file.exists(opt$idf_csv)) {
  k <- read.csv(opt$idf_csv, stringsAsFactors = FALSE)
  # expected columns: locus, enabled
  nms <- tolower(names(k)); names(k) <- nms
  if (!all(c("locus","enabled") %in% names(k))) stop("idf_csv needs locus,enabled")
  k <- k[k$enabled == 1, , drop = FALSE]
  w <- match(db$locus_ids, k$locus)
  w <- which(!is.na(w))
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
compute_scores <- (rep_mode %in% c("all","top")) && (as.integer(opt$n_cap) >= 0L)
comp_sec <- 0
peak_mib <- NA_real_

res <- data.frame(SampleID=character(0), Score=integer(0))
if (compute_scores) {
  if (requireNamespace("peakRAM", quietly=TRUE)) {
    pm <- peakRAM::peakRAM({
      opts_cpp <- list(
        any_code = as.integer(opt$any_code),
        score_min = as.integer(opt$score_min),    # 0=無制限
        n_cap = as.integer(opt$n_cap),            # 0=無制限
        force_all_loci = as.logical(opt$force_all_loci),
        h2a_Q = FALSE, h2a_R = FALSE
      )
      res <<- dmp_match_cpp(
        A1 = A1m, A2 = A2m,
        q1 = q$q1, q2 = q$q2,
        score_table = SCORE_TABLE,
        idf_mask_bits = as.integer(idf_mask_bits),
        opts = opts_cpp,
        sample_ids = db$sample_ids
      )
    })
    comp_sec <- pm$Elapsed_Time_sec
    peak_mib <- pm$Peak_RAM_Used_MiB
  } else {
    t0 <- proc.time()[["elapsed"]]
    opts_cpp <- list(
      any_code = as.integer(opt$any_code),
      score_min = as.integer(opt$score_min),
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
    comp_sec <- proc.time()[["elapsed"]] - t0
    peak_mib <- NA_real_
  }
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
if (as.integer(opt$score_min) > 0L) {
  hist_name <- sub("^hist", sprintf("hist_t%d", as.integer(opt$score_min)), hist_name)
}
hist_path <- file.path(dirp, hist_name)
if (rep_mode %in% c("all","hist")) {
  if (requireNamespace("peakRAM", quietly=TRUE)) {
    pmh <- peakRAM::peakRAM({
      hdf <<- dmp_hist_cpp(
        A1 = A1m, A2 = A2m,
        q1 = q$q1, q2 = q$q2,
        score_table = SCORE_TABLE,
        idf_mask_bits = as.integer(idf_mask_bits),
        opts = list(any_code = as.integer(opt$any_code), force_all_loci = as.logical(opt$force_all_loci))
      )
    })
    # comp_sec accumulation
    if (rep_mode == "hist") comp_sec <- pmh$Elapsed_Time_sec else comp_sec <- comp_sec + pmh$Elapsed_Time_sec
    # peak is max of phases
    if (is.na(peak_mib)) peak_mib <- pmh$Peak_RAM_Used_MiB else peak_mib <- max(peak_mib, pmh$Peak_RAM_Used_MiB, na.rm=TRUE)
  } else {
    t0h <- proc.time()[["elapsed"]]
    hdf <- dmp_hist_cpp(
      A1 = A1m, A2 = A2m,
      q1 = q$q1, q2 = q$q2,
      score_table = SCORE_TABLE,
      idf_mask_bits = as.integer(idf_mask_bits),
      opts = list(any_code = as.integer(opt$any_code), force_all_loci = as.logical(opt$force_all_loci))
    )
    t1h <- proc.time()[["elapsed"]]
    if (rep_mode == "hist") comp_sec <- (t1h - t0h) else comp_sec <- comp_sec + (t1h - t0h)
    # peak_mib remains as is (NA) when peakRAM unavailable
  }
  write.csv(hdf, hist_path, row.names = FALSE)
  if (opt$debug) cat(sprintf("[OK] wrote %s (%d rows)\n", hist_path, nrow(hdf)))
}

# ---- opts_<core>.json & time.log ----
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

# time.log: name 列は右端へ
header <- "LOAD_DB_SEC,LOAD_Q_SEC,COMP_SEC,TOTAL_SEC,PEAK_MiB,name"
total_sec <- proc.time()[["elapsed"]] - t0_total
line <- sprintf("%.3f,%.3f,%.3f,%.3f,%.1f,%s",
                load_db_sec, load_q_sec, comp_sec, total_sec,
                ifelse(is.na(peak_mib), NA_real_, peak_mib), core)
logp <- file.path(dirname(opt$out), "time.log")
if (!file.exists(logp)) writeLines(header, logp)
cat(paste0(line, "\n"), file = logp, append = TRUE)
if (opt$debug) cat("[DBG] done\n")

# scripts/matcher_onepass.R
# One-pass matcher CLI wrapper for dmp_match_cpp/dmp_hist_cpp (schema-strict)
# ASCII-only; JP comments allowed for dev.

suppressPackageStartupMessages({
  library(optparse)
  library(Rcpp)
  library(jsonlite)
  library(tools)
})

# Ensure peakRAM (install if missing; otherwise PEAK_MiB will be NA)
if (!requireNamespace("peakRAM", quietly=TRUE)) {
  tryCatch(install.packages("peakRAM", repos="https://cran.rstudio.com"),
           error=function(e) { message("[WARN] peakRAM install failed; PEAK_MiB will be NA") })
}

# SCORE_TABLE (index: 0..15)
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
  make_option("--n_cap",   type="integer",   default=200L),        # TopN cap
  make_option("--score_min", type="integer", default=0L),          # >=0
  make_option("--any_code",  type="integer", default=9999L),
  make_option("--idf_csv",   type="character", default=NA),        # optional
  make_option("--force_all_loci", type="integer", default=0L),     # 1=ignore idf mask, use all loci
  make_option("--h2a_on", type="integer", default=0L),             # 1=DB側のホモは右any
  make_option("--debug",   type="integer",   default=1L),
  make_option("--path",    type="character", default=NA)           # optional setwd
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---- normalize/guards ----
if (is.null(opt$db) || is.null(opt$query)) stop("--db and --query are required")
if (!is.na(opt$path)) {
  if (!dir.exists(opt$path)) stop(sprintf("--path not found: %s", opt$path))
  setwd(normalizePath(opt$path, winslash = "/"))
}
if (!file.exists(opt$db))    stop(sprintf("db not found: %s", opt$db))
if (!file.exists(opt$query)) stop(sprintf("query not found: %s", opt$query))
if (as.integer(opt$n_cap) < 0L) stop("--n_cap must be >= 0")
opt$score_min <- as.integer(opt$score_min)
if (is.na(opt$score_min) || opt$score_min < 0L) stop("--score_min must be an integer >= 0")

# ---- compile core ----
sourceCpp("src/dmp_match.cpp")  # exports: dmp_match_cpp, dmp_hist_cpp

t0_total <- proc.time()[["elapsed"]]

# ---- load DB (schema-strict) ----
t0 <- proc.time()[["elapsed"]]
db <- readRDS(opt$db)  # list(sample_ids, locus_ids, A1, A2)
if (!is.list(db) || !all(c("sample_ids","locus_ids","A1","A2") %in% names(db))) {
  stop("DB RDS must contain list(sample_ids, locus_ids, A1, A2)")
}
S <- length(db$sample_ids)
L <- length(db$locus_ids)
if (L <= 0L || S <= 0L) stop("DB has empty sample_ids or locus_ids")

to_mat <- function(X) {
  if (is.matrix(X)) {
    storage.mode(X) <- "integer"; return(X)
  }
  if (is.data.frame(X)) X <- as.matrix(X)
  if (is.vector(X) && length(X) == S * L) {
    X <- matrix(X, nrow=S, ncol=L, byrow=FALSE)
  }
  storage.mode(X) <- "integer"
  X
}
A1m <- to_mat(db$A1)
A2m <- to_mat(db$A2)
stopifnot(nrow(A1m) == S, ncol(A1m) == L)
stopifnot(nrow(A2m) == S, ncol(A2m) == L)
t1 <- proc.time()[["elapsed"]]
load_db_sec <- t1 - t0

# ---- h2a (DB homo -> RIGHT ANY) after standardization assumption ----
if (as.integer(opt$h2a_on) == 1L) {
  ANY <- as.integer(opt$any_code)
  mask_homo <- (A1m == A2m) & (A1m != ANY)
  A2m[mask_homo] <- ANY
  if (as.integer(opt$debug) == 1L) {
    cat(sprintf("[DBG] h2a_on: masked %d homozygous cells -> ANY\n", sum(mask_homo, na.rm=TRUE)))
  }
}

# ---- load query (RDS or CSV) ----
t0 <- proc.time()[["elapsed"]]
ANY <- as.integer(opt$any_code)
norm_scale100 <- function(v) {
  # accept numeric or character; "any" -> ANY
  if (is.character(v)) {
    vi <- ifelse(tolower(v) == "any" | v == "", ANY, as.integer(round(as.numeric(v) * 100)))
    return(vi)
  } else {
    if (is.integer(v)) return(v)
    return(as.integer(round(as.numeric(v) * 100)))
  }
}
sort_pair <- function(a1, a2, ANY) {
  a1 <- as.integer(a1); a2 <- as.integer(a2)
  swap <- ((a1 == ANY & a2 != ANY) | (a1 > a2 & a2 != ANY))
  list(a1 = ifelse(swap, a2, a1), a2 = ifelse(swap, a1, a2))
}

q <- NULL
if (grepl("\\.(rds|RDS)$", opt$query)) {
  qx <- readRDS(opt$query)
  if (!is.list(qx) || !all(c("q1","q2") %in% names(qx))) {
    stop("Query RDS must provide list(q1, q2)")
  }
  if (!(length(qx$q1)==L && length(qx$q2)==L)) stop("Query lengths must match L")
  qp <- sort_pair(as.integer(qx$q1), as.integer(qx$q2), ANY)
  q <- list(q1 = qp$a1, q2 = qp$a2)
} else {
  x <- read.csv(opt$query, stringsAsFactors = FALSE, check.names = FALSE)
  nms <- tolower(names(x)); names(x) <- nms
  lc <- which(nms %in% c("locus","marker"))[1]
  a1c <- which(nms %in% c("allele1","a1","q1"))[1]
  a2c <- which(nms %in% c("allele2","a2","q2"))[1]
  if (is.na(lc) || is.na(a1c) || is.na(a2c)) stop("CSV needs columns: Locus(+), Allele1/A1/Q1, Allele2/A2/Q2")
  mm <- match(db$locus_ids, x[[lc]])
  if (any(is.na(mm))) stop("CSV loci must cover all db$locus_ids")
  qp <- sort_pair(norm_scale100(x[[a1c]][mm]), norm_scale100(x[[a2c]][mm]), ANY)
  q <- list(q1 = qp$a1, q2 = qp$a2)
}
t1 <- proc.time()[["elapsed"]]
load_q_sec <- t1 - t0

# ---- idf mask bits (uint32_t; bit j=1 means enabled locus j) ----
make_all_ones <- function(L) {
  bits <- 0L
  for (j in 0:(L-1)) {
    if (j >= 31) break  # uint32_t safety; DMP uses 21 loci
    bits <- bitwOr(bits, bitwShiftL(1L, j))
  }
  bits
}
idf_mask_bits <- make_all_ones(L)
mask_name <- "normal"

if (!is.na(opt$idf_csv) && nzchar(opt$idf_csv) && file.exists(opt$idf_csv)) {
  mask_name <- "gf_idf_mask"
  dfm <- read.csv(opt$idf_csv, stringsAsFactors = FALSE, check.names = FALSE)
  nms <- tolower(names(dfm)); names(dfm) <- nms
  # tolerant columns: locus + (bit|mask|enabled) where 1=enabled,0=disabled
  lc <- which(nms %in% c("locus","marker","locus_id"))[1]
  bc <- which(nms %in% c("bit","mask","enabled","enable","use"))[1]
  if (is.na(lc) || is.na(bc)) {
    warning("idf_csv lacks recognizable columns; falling back to all-ones")
  } else {
    idf_mask_bits <- 0L
    mm <- match(db$locus_ids, dfm[[lc]])
    if (any(is.na(mm))) {
      warning("idf_csv does not cover all loci; missing ones are treated as enabled")
    }
    for (j in 0:(L-1)) {
      if (j >= 31) break  # uint32_t
      v <- 1L
      if (!is.na(mm[j+1])) {
        vv <- suppressWarnings(as.integer(dfm[[bc]][ mm[j+1] ]))
        if (!is.na(vv)) v <- ifelse(vv != 0L, 1L, 0L)
      }
      if (v == 1L) idf_mask_bits <- bitwOr(idf_mask_bits, bitwShiftL(1L, j))
    }
  }
}

if (as.integer(opt$force_all_loci) == 1L) {
  idf_mask_bits <- make_all_ones(L)
  mask_name <- "normal"
}

# ---- run core ----
rep_mode <- tolower(opt$report)
compute_scores <- (rep_mode %in% c("all","top"))
compute_hist   <- (rep_mode %in% c("all","hist"))

comp_sec <- 0
peak_mib <- NA_real_

# ensure out dir
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)

res <- data.frame(SampleID=character(0), Score=integer(0), stringsAsFactors = FALSE)
if (compute_scores) {
  if (requireNamespace("peakRAM", quietly=TRUE)) {
    pm <- peakRAM::peakRAM({
      opts_cpp <- list(
        any_code = as.integer(opt$any_code),
        score_min = as.integer(opt$score_min),
        n_cap = as.integer(opt$n_cap),
        force_all_loci = as.logical(opt$force_all_loci)
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
    comp_sec <- comp_sec + as.numeric(pm$Elapsed_Time_sec)
    peak_mib <- pm$Peak_RAM_Used_MiB
  } else {
    t0 <- proc.time()[["elapsed"]]
    opts_cpp <- list(
      any_code = as.integer(opt$any_code),
      score_min = as.integer(opt$score_min),
      n_cap = as.integer(opt$n_cap),
      force_all_loci = as.logical(opt$force_all_loci)
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
    comp_sec <- comp_sec + (t1 - t0)
    # peak_mib stays NA if peakRAM not available
  }
  # write scores (already sorted desc by C++)
  write.csv(res, opt$out, row.names = FALSE)
  if (as.integer(opt$debug) == 1L) cat(sprintf("[OK] wrote %s (%d rows)\n", opt$out, nrow(res)))
}

# histogram
hdf <- NULL
if (compute_hist) {
  hist_path <- file.path(dirname(opt$out), sub("^scores", "hist", basename(opt$out)))
  if (opt$score_min > 0L) {
    hist_path <- file.path(dirname(opt$out), sub("^hist", sprintf("hist_t%d", as.integer(opt$score_min)), basename(hist_path)))
  }
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
    comp_sec <- comp_sec + as.numeric(pmh$Elapsed_Time_sec)
    if (is.na(peak_mib)) peak_mib <- pmh$Peak_RAM_Used_MiB else peak_mib <- max(peak_mib, pmh$Peak_RAM_Used_MiB, na.rm=TRUE)
  } else {
    t0 <- proc.time()[["elapsed"]]
    hdf <- dmp_hist_cpp(
      A1 = A1m, A2 = A2m,
      q1 = q$q1, q2 = q$q2,
      score_table = SCORE_TABLE,
      idf_mask_bits = as.integer(idf_mask_bits),
      opts = list(any_code = as.integer(opt$any_code), force_all_loci = as.logical(opt$force_all_loci))
    )
    t1 <- proc.time()[["elapsed"]]
    comp_sec <- comp_sec + (t1 - t0)
  }
  write.csv(hdf, hist_path, row.names = FALSE)
  if (as.integer(opt$debug) == 1L) cat(sprintf("[OK] wrote %s (%d rows)\n", hist_path, nrow(hdf)))
}

# ---- opts_*.json ----
opts_list <- list(
  db_path = normalizePath(opt$db, winslash = "/"),
  query_path = normalizePath(opt$query, winslash = "/"),
  mask_name = mask_name,
  mask_ratio = NA_real_,    # Phase2以降で使用
  mask_seed = NA_integer_,  # Phase2以降で使用
  any_code = as.integer(opt$any_code),
  force_all_loci = as.logical(opt$force_all_loci),
  h2a_on = as.logical(opt$h2a_on)
)
json_path <- file.path(dirname(opt$out), sprintf("opts_%s.json", file_path_sans_ext(basename(opt$out))))
writeLines(jsonlite::toJSON(opts_list, auto_unbox = TRUE, pretty = TRUE), json_path)
if (as.integer(opt$debug) == 1L) cat(sprintf("[OK] wrote %s\n", json_path))

# ---- time.log ----
# name format: "{kit}/{mask}/{ratio}/{h2a}/{DBsize}"
kit <- "GF"
ratio_str <- if (identical(mask_name, "normal")) "0.0" else "1.0"  # Phase1: simple 0.0 or 1.0
h2a_str <- if (as.integer(opt$h2a_on) == 1L) "on" else "off"
dbsize <- as.integer(S)
core_name <- sprintf("%s/%s/%s/%s/%d", kit, mask_name, ratio_str, h2a_str, dbsize)

total_sec <- proc.time()[["elapsed"]] - t0_total
header <- "LOAD_DB_SEC,LOAD_Q_SEC,COMP_SEC,TOTAL_SEC,PEAK_MiB,name"
line <- sprintf("%.3f,%.3f,%.3f,%.3f,%.1f,%s",
                load_db_sec, load_q_sec, comp_sec, total_sec,
                ifelse(is.na(peak_mib), NA_real_, peak_mib), core_name)
logp <- file.path(dirname(opt$out), "time.log")
if (!file.exists(logp)) writeLines(header, logp)
cat(paste0(line, "\n"), file = logp, append = TRUE)
if (as.integer(opt$debug) == 1L) cat("[DBG] done\n")

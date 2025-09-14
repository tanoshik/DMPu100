# scripts/matcher_onepass.R
suppressPackageStartupMessages({
  library(optparse)
  library(Rcpp)
})

# Rcpp core
sourceCpp("src/dmp_match.cpp")

# SCORE_TABLE (DMP spec)
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
  make_option("--debug",   type="integer",   default=1L)
)
opt <- parse_args(OptionParser(option_list=option_list))
stopifnot(!is.null(opt$db), !is.null(opt$query))
if (!file.exists(opt$db))    stop(sprintf("db not found: %s", opt$db))
if (!file.exists(opt$query)) stop(sprintf("query not found: %s", opt$query))

# ---- DB load (前提：正しいRDS) ----
db <- readRDS(opt$db)
need <- c("A1","A2","sample_ids","locus_ids")
miss <- setdiff(need, names(db))
if (length(miss)) stop("db missing fields: ", paste(miss, collapse=","))
if (!is.integer(db$A1)) db$A1 <- apply(db$A1, 2, as.integer)
if (!is.integer(db$A2)) db$A2 <- apply(db$A2, 2, as.integer)
if (!"homo_mask" %in% names(db)) {
  S <- nrow(db$A1); L <- ncol(db$A1)
  hm <- integer(S)
  for (s in seq_len(S)) {
    bits <- 0L
    for (l in 0:(L-1)) if (db$A1[s,l+1] == db$A2[s,l+1]) bits <- bitwOr(bits, bitwShiftL(1L, l))
    hm[s] <- bits
  }
  db$homo_mask <- hm
}
L <- ncol(db$A1)

# ---- Query load: サニティの encode_one と同仕様 ----
load_query <- function(path, locus_ids, any_code=9999L) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    q <- readRDS(path)
    if (!is.list(q) || !all(c("q1","q2") %in% names(q))) stop("RDS query must have q1,q2")
    if (!is.integer(q$q1)) q$q1 <- as.integer(q$q1)
    if (!is.integer(q$q2)) q$q2 <- as.integer(q$q2)
    return(list(q1=q$q1, q2=q$q2))
  }
  if (ext != "csv") stop("Unsupported query extension: ", ext)
  
  df <- read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)
  req <- c("Locus","allele1","allele2")
  if (!all(req %in% names(df))) stop("CSV must have columns: Locus, allele1, allele2")
  
  encode_one <- function(xx, ANY=as.integer(any_code)) {
    if (is.na(xx)) return(ANY)
    if (is.character(xx)) {
      sx <- trimws(tolower(xx))
      if (sx == "any") return(ANY)
      vx <- suppressWarnings(as.numeric(sx))
      if (!is.na(vx)) return(as.integer(round(vx * 100)))
      stop("Unrecognized allele in query: ", xx)
    } else if (is.numeric(xx)) {
      if (abs(xx - round(xx)) > 1e-9) return(as.integer(round(xx * 100)))
      xi <- as.integer(round(xx))
      if (xi %% 100 == 0) return(xi) else return(as.integer(xi * 100))
    } else {
      stop("Unsupported allele type: ", typeof(xx))
    }
  }
  
  q1c <- vapply(df$allele1, encode_one, integer(1))
  q2c <- vapply(df$allele2, encode_one, integer(1))
  
  # locus順に並べ替え（完全一致前提）
  idx <- match(locus_ids, df$Locus)
  if (any(is.na(idx))) {
    miss <- locus_ids[is.na(idx)]
    stop("Query does not cover loci: ", paste(miss, collapse=", "))
  }
  list(q1=as.integer(q1c[idx]), q2=as.integer(q2c[idx]))
}

q <- load_query(opt$query, db$locus_ids, any_code=opt$any_code)

# ---- IDF mask bits ----
idf_mask_bits <- if (!is.na(opt$idf_csv) && nzchar(opt$idf_csv) && file.exists(opt$idf_csv)) {
  idf <- read.csv(opt$idf_csv, stringsAsFactors=FALSE, check.names=FALSE)
  if (!all(c("locus","enabled") %in% tolower(names(idf))))
    stop("idf_csv must have columns: locus, enabled")
  names(idf) <- tolower(names(idf))
  bits <- 0L
  for (l in 0:(L-1)) {
    loc <- tolower(as.character(db$locus_ids[l+1]))
    row <- idf[tolower(idf$locus)==loc, , drop=FALSE]
    en <- if (nrow(row)==1) as.integer(row$enabled[[1]]) else 1L
    if (is.na(en)) en <- 1L
    if (en != 0L) bits <- bitwOr(bits, bitwShiftL(1L, l))
  }
  bits
} else {
  bitwShiftL(1L, L) - 1L
}

if (opt$debug)
  cat(sprintf("[DBG] SxL=%dx%d idf=0x%X force_all=%d\n",
              nrow(db$A1), ncol(db$A1), as.integer(idf_mask_bits), as.integer(opt$force_all_loci)))

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

res <- dmp_match_cpp(
  A1 = db$A1, A2 = db$A2,
  q1 = q$q1, q2 = q$q2,
  score_table = SCORE_TABLE,
  homo_mask = db$homo_mask,
  idf_mask_bits = as.integer(idf_mask_bits),
  opts = opts_cpp,
  sample_ids = db$sample_ids
)

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.csv(res, opt$out, row.names = FALSE)
cat(sprintf("[OK] wrote %s (%d rows)\n", opt$out, nrow(res)))

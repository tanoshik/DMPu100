# scripts/cli/dmp_apply_idf_mask.R
# Pre-apply IDF% mask to a DB RDS and write a masked DB RDS.
# Input : list(sample_ids, locus_ids, A1, A2) integers (x100), ANY=9999
# Output: same schema as input (pre-masked A1/A2)
# Mask  : rows = top floor(ratio*S) by order(crc32(sample_id) ^ seed)
#         cols = loci with enabled == 0 in idf_csv
# Notes : ASCII-only; snake_case CLI; no external packages except optparse.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--in_db"),     type="character", help="input DB RDS path"),
  make_option(c("--out_db"),    type="character", help="output masked DB RDS path"),
  make_option(c("--idf_csv"),   type="character", help="mask CSV with columns locus,enabled"),
  make_option(c("--ratio"),     type="double",    default=0.0, help="ratio in [0,1]"),
  make_option(c("--seed"),      type="integer",   default=0L,  help="32-bit integer"),
  make_option(c("--any_code"),  type="integer",   default=9999L, help="ANY code"),
  make_option(c("--debug"),     type="integer",   default=1L)
)
opt <- parse_args(OptionParser(option_list=option_list))

# Guards
if (is.null(opt$in_db)  || !nzchar(opt$in_db))  stop("--in_db is required")
if (is.null(opt$out_db) || !nzchar(opt$out_db)) stop("--out_db is required")
if (is.null(opt$idf_csv)|| !nzchar(opt$idf_csv))stop("--idf_csv is required")
if (!file.exists(opt$in_db))  stop(sprintf("in_db not found: %s", opt$in_db))
if (!file.exists(opt$idf_csv)) stop(sprintf("idf_csv not found: %s", opt$idf_csv))
if (is.na(opt$ratio) || opt$ratio < 0 || opt$ratio > 1) stop("--ratio must be in [0,1]")

# Pure-R CRC32
crc32_table <- local({
  poly <- 0xEDB88320
  tab <- integer(256)
  for (i in 0:255) {
    crc <- i
    for (j in 1:8) {
      if (bitwAnd(crc,1L)!=0L) crc <- bitwXor(bitwShiftR(crc,1L), poly) else crc <- bitwShiftR(crc,1L)
    }
    tab[i+1] <- crc
  }
  tab
})
crc32_vec <- function(xs) {
  vapply(xs, function(s) {
    b <- charToRaw(as.character(s))
    crc <- 0xFFFFFFFF
    for (by in b) {
      idx <- bitwAnd(bitwXor(crc, as.integer(by)), 0xFF) + 1L
      crc <- bitwXor(bitwShiftR(crc,8L), crc32_table[idx])
    }
    bitwXor(crc, 0xFFFFFFFF)
  }, integer(1))
}

# Load DB
db <- readRDS(opt$in_db)
req <- c("sample_ids","locus_ids","A1","A2")
if (!is.list(db) || !all(req %in% names(db))) stop("DB RDS must be list(sample_ids, locus_ids, A1, A2)")

sample_ids <- as.character(db$sample_ids)
locus_ids  <- as.character(db$locus_ids)
A1 <- db$A1; A2 <- db$A2

S <- length(sample_ids); L <- length(locus_ids)
to_i_mat <- function(X){ X <- as.matrix(X); storage.mode(X) <- "integer"; X }
A1m <- to_i_mat(A1); A2m <- to_i_mat(A2)
if (nrow(A1m)!=S || ncol(A1m)!=L || nrow(A2m)!=S || ncol(A2m)!=L) stop("A1/A2 shape mismatch")

ANY <- as.integer(opt$any_code)

# Load IDF mask CSV
dfm <- read.csv(opt$idf_csv, stringsAsFactors=FALSE, check.names=FALSE)
nms <- tolower(names(dfm)); names(dfm) <- nms
lc <- which(nms %in% c("locus","marker","locus_id"))[1]
ec <- which(nms %in% c("enabled","enable","use","mask","bit"))[1]
if (is.na(lc) || is.na(ec)) stop("idf_csv must have columns: locus, enabled")

mm <- match(locus_ids, dfm[[lc]])
if (any(is.na(mm))) stop("idf_csv does not cover all loci")
enabled <- suppressWarnings(as.integer(dfm[[ec]][mm]))
enabled[is.na(enabled)] <- 1L
disabled_cols <- which(enabled == 0L)

if (length(disabled_cols) == 0L) {
  if (opt$debug==1L) cat("[INFO] no disabled loci; copying DB\n")
  saveRDS(db, opt$out_db)
  quit(status=0)
}

# Row selection by ratio
rows <- integer(0)
if (opt$ratio > 0) {
  h <- bitwXor(crc32_vec(sample_ids), as.integer(opt$seed))
  ord <- order(h, method="radix")
  k <- floor(opt$ratio * S)
  if (k > 0L) rows <- ord[seq_len(min(k, S))]
}

# Apply pre-mask (ANY)
if (length(rows) > 0L) {
  A1m[rows, disabled_cols] <- ANY
  A2m[rows, disabled_cols] <- ANY
  if (opt$debug==1L) cat(sprintf("[OK] masked rows=%d/%d, disabled_cols=%d\n", length(rows), S, length(disabled_cols)))
} else {
  if (opt$debug==1L) cat("[INFO] ratio==0 or k==0; no rows masked\n")
}

# Write masked DB (schema-preserving)
db_out <- list(sample_ids = sample_ids, locus_ids = locus_ids, A1 = A1m, A2 = A2m)
saveRDS(db_out, opt$out_db)
if (opt$debug==1L) cat(sprintf("[OK] wrote masked DB: %s\n", normalizePath(opt$out_db, winslash="/")))

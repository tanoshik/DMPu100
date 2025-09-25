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
  make_option(c("--ratio"),     type="double",    default=0.0,    help="ratio in [0,1]"),
  make_option(c("--seed"),      type="integer",   default=0L,     help="32-bit integer"),
  make_option(c("--any_code"),  type="integer",   default=9999L,  help="ANY code"),
  make_option(c("--debug"),     type="integer",   default=1L)
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---- guards ----
if (is.null(opt$in_db)  || !nzchar(opt$in_db))  stop("--in_db is required")
if (is.null(opt$out_db) || !nzchar(opt$out_db)) stop("--out_db is required")
if (is.null(opt$idf_csv)|| !nzchar(opt$idf_csv))stop("--idf_csv is required")
if (!file.exists(opt$in_db))  stop(sprintf("in_db not found: %s", opt$in_db))
if (!file.exists(opt$idf_csv)) stop(sprintf("idf_csv not found: %s", opt$idf_csv))
if (is.na(opt$ratio) || opt$ratio < 0 || opt$ratio > 1) stop("--ratio must be in [0,1]")

# ---- CRC32 (純R, 256-table) ----
# R の integer は 32bit 符号付きなので、unsigned 定数は 2^32 補数の符号付き値を使う
# 0xEDB88320 (unsigned) = -306674912 (signed 32bit)
# 0xFFFFFFFF (unsigned) = -1 (signed 32bit)
crc32_table <- local({
  poly <- -306674912L  # 0xEDB88320
  tab <- integer(256)
  for (i in 0:255) {
    crc <- as.integer(i)
    for (j in 1:8) {
      if (bitwAnd(crc, 1L) != 0L) crc <- bitwXor(bitwShiftR(crc, 1L), poly) else crc <- bitwShiftR(crc, 1L)
    }
    tab[i + 1L] <- crc
  }
  tab
})
crc32_vec <- function(xs) {
  vapply(xs, function(s) {
    b <- charToRaw(if (is.na(s)) "" else as.character(s))
    crc <- -1L  # 0xFFFFFFFF
    for (by in b) {
      idx <- bitwAnd(bitwXor(crc, as.integer(by)), 255L) + 1L  # 0xFF = 255
      crc <- bitwXor(bitwShiftR(crc, 8L), crc32_table[idx])
    }
    bitwXor(crc, -1L)  # 最終反転
  }, integer(1))
}

# ---- load DB ----
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

# ---- load IDF mask CSV ----
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

# ---- fast path: ratio==0 or no disabled cols -> passthrough ----
ratio <- as.numeric(opt$ratio)
if (is.na(ratio) || ratio <= 0 || length(disabled_cols) == 0L) {
  if (opt$debug==1L) cat(sprintf("[DBG] passthrough (ratio<=0 or no disabled cols); S=%d, L=%d\n", S, L))
  saveRDS(db, file = opt$out_db)
  if (opt$debug==1L) cat(sprintf("[OK] wrote masked DB (passthrough): %s\n", normalizePath(opt$out_db, winslash="/")))
  quit(save="no", status=0)
}

# ---- row selection by ratio (crc32 ^ seed, radix order) ----
seed <- as.integer(opt$seed)
h <- bitwXor(crc32_vec(sample_ids), seed)
if (any(is.na(h))) stop("crc32_vec produced NA (unexpected)")
ord <- order(h, method="radix")
k <- max(0L, min(S, floor(ratio * S)))
rows <- if (k > 0L) ord[seq_len(k)] else integer(0)

# ---- apply pre-mask (ANY) ----
if (length(rows) > 0L && length(disabled_cols) > 0L) {
  A1m[rows, disabled_cols] <- ANY
  A2m[rows, disabled_cols] <- ANY
  if (opt$debug==1L) cat(sprintf("[OK] masked rows=%d/%d, disabled_cols=%d\n", length(rows), S, length(disabled_cols)))
} else {
  if (opt$debug==1L) cat("[INFO] nothing masked (rows or disabled_cols empty)\n")
}

# ---- write masked DB (schema-preserving) ----
db_out <- list(sample_ids = sample_ids, locus_ids = locus_ids, A1 = A1m, A2 = A2m)
saveRDS(db_out, opt$out_db)
if (opt$debug==1L) cat(sprintf("[OK] wrote masked DB: %s\n", normalizePath(opt$out_db, winslash="/")))

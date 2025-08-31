# quickstart/test_data_integrity.R
# No multibyte chars. Minimal integrity checks for DB RDS and Query CSV.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  db   = "data/virtual_db_u100_S1000_seed123.rds",
  csv  = "data/query_profile_seed123.csv",
  out  = "output/sample/integrity_report.txt",
  any_code = 9999L,
  max_scan = 5e6
)
for (a in args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2) opt[[sub("^--?", "", kv[1])]] <- kv[2]
}
opt$any_code <- as.integer(opt$any_code)

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
logf <- file(opt$out, "wt"); on.exit(close(logf), add = TRUE)
say <- function(...) cat(..., "\n", file = logf, sep = "")

ok_all <- TRUE
fail <- function(msg) { say(paste0("[FAIL] ", msg)); ok_all <<- FALSE }

must_read_rds <- function(path) {
  if (!file.exists(path)) stop("RDS not found: ", path, call. = FALSE)
  readRDS(path)
}
must_read_csv <- function(path) {
  if (!file.exists(path)) stop("CSV not found: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE)
}

# Read locus order
if (!file.exists("data/locus_order.rds")) {
  fail("data/locus_order.rds not found"); loc_order <- character()
} else {
  loc_order <- readRDS("data/locus_order.rds")
  say("[OK] locus_order.rds loaded: ", length(loc_order), " loci")
}

# --- Check DB RDS ---
try({
  db <- must_read_rds(opt$db)
  need <- c("sample_ids","locus_ids","A1","A2")
  if (!all(need %in% names(db))) fail("db_index missing fields")
  if (!is.character(db$sample_ids) || !is.character(db$locus_ids)) fail("ids should be character")
  if (!is.matrix(db$A1) || !is.matrix(db$A2)) fail("A1/A2 should be matrices")
  if (!is.integer(db$A1) || !is.integer(db$A2)) fail("A1/A2 should be integer matrices")
  if (nrow(db$A1) != length(db$sample_ids) || nrow(db$A2) != length(db$sample_ids)) fail("A1/A2 row mismatch")
  if (ncol(db$A1) != length(db$locus_ids)  || ncol(db$A2) != length(db$locus_ids))  fail("A1/A2 col mismatch")
  if (length(loc_order) > 0 && !identical(db$locus_ids, as.character(loc_order))) {
    fail("db$locus_ids is not identical to locus_order.rds")
  }
  cells <- as.double(nrow(db$A1)) * as.double(ncol(db$A1))
  say("[OK] db_index shape: S=", nrow(db$A1), " L=", ncol(db$A1), " cells=", format(cells, big.mark=","))
  
  check_block <- function(A1, A2, any_code, max_scan) {
    total <- length(A1); k <- min(total, as.integer(max_scan))
    idx <- if (total <= k) seq_len(total) else sample.int(total, k)
    a1 <- A1[idx]; a2 <- A2[idx]
    bad_order <- sum(a1 > a2, na.rm = TRUE)
    bad_any   <- sum(a1 == any_code & a2 != any_code, na.rm = TRUE)
    c(bad_order = bad_order, bad_any = bad_any, scanned = length(idx))
  }
  chk <- check_block(db$A1, db$A2, opt$any_code, opt$max_scan)
  if (chk["bad_order"] > 0) fail(paste0("allele order violated: count=", chk["bad_order"]))
  if (chk["bad_any"]   > 0) fail(paste0("ANY should be right: count=",   chk["bad_any"]))
  say("[OK] allele ordering rules passed (scanned ", chk["scanned"], " cells)")
}, silent = TRUE)

# --- Check Query CSV ---
try({
  qraw <- must_read_csv(opt$csv)
  nms <- names(qraw)
  nms <- sub("(?i)^locus$","Locus", nms, perl=TRUE)
  nms <- sub("(?i)^allele1$","allele1", nms, perl=TRUE)
  nms <- sub("(?i)^allele2$","allele2", nms, perl=TRUE)
  names(qraw) <- nms
  need <- c("Locus","allele1","allele2")
  if (!all(need %in% names(qraw))) fail("CSV columns must include Locus, allele1, allele2")
  q <- qraw[,need,drop=FALSE]
  if (length(loc_order) > 0) {
    if (!setequal(q$Locus, loc_order)) fail("CSV loci set differs from locus_order.rds")
    o <- match(loc_order, q$Locus); q <- q[o,,drop=FALSE]
    if (!identical(q$Locus, loc_order)) fail("CSV order differs from locus_order.rds")
  }
  map_any <- function(x, any_code=opt$any_code) {
    x <- as.character(x); x[is.na(x) | x==""] <- "any"
    x[tolower(x)=="any"] <- as.character(any_code)
    suppressWarnings(as.integer(x))
  }
  q$allele1 <- map_any(q$allele1)
  q$allele2 <- map_any(q$allele2)
  bad_order <- sum(q$allele1 > q$allele2, na.rm=TRUE)
  bad_any   <- sum(q$allele1 == opt$any_code & q$allele2 != opt$any_code, na.rm=TRUE)
  if (bad_order > 0) fail(paste0("query allele order violated: count=", bad_order))
  if (bad_any   > 0) fail(paste0("query ANY should be right: count=", bad_any))
  say("[OK] query CSV structure passed: ", nrow(q), " loci")
}, silent = TRUE)

# --- final verdict (use braces!) ---
if (ok_all) {
  say("[PASS] integrity checks passed")
} else {
  say("[FAIL] integrity checks failed")
}

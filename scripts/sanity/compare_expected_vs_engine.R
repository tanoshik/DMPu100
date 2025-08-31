# scripts/sanity/compare_expected_vs_engine.R
# No multibyte chars. Compare two Top3 per-locus CSVs strictly, after canonical ordering.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  expected = "output/sample/match_detail_S1000_top3.csv",
  engine   = "output/sample/match_detail_S1000_top3_from_engine.csv",
  locus_order_rds = "data/locus_order.rds",
  out_diff = "output/sample/top3_expected_vs_engine_diff.csv"
)
for (a in args) { kv <- strsplit(a,"=",fixed=TRUE)[[1]]; if (length(kv)==2) opt[[sub("^--?","",kv[1])]] <- kv[2] }

stopifnot(file.exists(opt$expected), file.exists(opt$engine), file.exists(opt$locus_order_rds))
loc <- readRDS(opt$locus_order_rds)

canon_order <- function(df) {
  # normalize colnames
  n <- names(df)
  n <- sub("(?i)^locus$","Locus", n, perl=TRUE)
  n <- sub("(?i)^sampleid$","SampleID", n, perl=TRUE)
  n <- sub("(?i)^db_allele1$","r1", n, perl=TRUE)
  n <- sub("(?i)^db_allele2$","r2", n, perl=TRUE)
  names(df) <- n
  # keep expected columns if present
  keep <- intersect(c("Locus","q1","q2","r1","r2","bits","code","score","SampleID"), names(df))
  df <- df[, keep, drop=FALSE]
  # order: SampleID asc, then Locus by locus_order
  if (!"SampleID" %in% names(df)) stop("column 'SampleID' not found")
  if (!"Locus" %in% names(df)) stop("column 'Locus' not found")
  df$Locus <- factor(df$Locus, levels = loc, ordered = TRUE)
  df <- df[order(df$SampleID, df$Locus), , drop=FALSE]
  rownames(df) <- NULL
  df
}

A <- read.csv(opt$expected, stringsAsFactors = FALSE)
B <- read.csv(opt$engine, stringsAsFactors = FALSE)
A <- canon_order(A); B <- canon_order(B)

# align columns
cols <- union(names(A), names(B))
for (cc in cols) { if (!cc %in% names(A)) A[[cc]] <- NA; if (!cc %in% names(B)) B[[cc]] <- NA }
A <- A[, cols, drop=FALSE]; B <- B[, cols, drop=FALSE]

neq <- A != B
neq[is.na(neq)] <- FALSE
num_diff <- sum(neq)

if (nrow(A) != nrow(B) || ncol(A) != ncol(B)) {
  cat("[FAIL] shape differs: A(", nrow(A), "x", ncol(A), ") vs B(", nrow(B), "x", ncol(B), ")\n", sep="")
} else if (num_diff == 0) {
  cat("[PASS] expected vs engine: identical after canonical ordering\n")
} else {
  cat("[FAIL] expected vs engine: ", num_diff, " cells differ -> ", opt$out_diff, "\n", sep="")
  out <- cbind(src="expected", A)
  out2 <- cbind(src="engine", B)
  # annotate differences only rows where any col differs
  rows <- which(rowSums(neq) > 0)
  out_all <- rbind(out[rows, , drop=FALSE], out2[rows, , drop=FALSE])
  dir.create(dirname(opt$out_diff), recursive=TRUE, showWarnings=FALSE)
  write.csv(out_all, opt$out_diff, row.names=FALSE)
}

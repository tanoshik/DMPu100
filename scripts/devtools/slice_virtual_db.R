# scripts/devtools/slice_virtual_db.R
# No multibyte chars. Slice first N samples from an indexed DB RDS.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  src  = "output/data/virtual_db_u100_S2000000_seed123.rds",
  size = 1000000L,
  out  = "output/data/virtual_db_u100_S1000000_seed123.rds"
)
for (a in args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2) opt[[sub("^--?","",kv[1])]] <- kv[2]
}
opt$size <- as.integer(opt$size)

db <- readRDS(opt$src)
n  <- min(opt$size, length(db$sample_ids))
idx <- seq_len(n)

db2 <- list(
  sample_ids = db$sample_ids[idx],
  locus_ids  = db$locus_ids,
  A1 = db$A1[idx, , drop=FALSE],
  A2 = db$A2[idx, , drop=FALSE]
)
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(db2, file = opt$out)
cat("[OK] sliced: ", opt$out, " S=", n, " L=", ncol(db2$A1), "\n", sep = "")

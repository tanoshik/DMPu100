# scripts/devtools/slice_query_fromdb.R
# No multibyte chars. Extract one sample from indexed DB RDS into a query CSV.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  db       = "output/data/virtual_db_u100_S2000000_seed123.rds",
  sample   = 1L,   # 1-based index or SampleID string
  out_csv  = "output/sample/query_profile_seed123.csv",
  any_code = 9999L
)
for (a in args) {
  kv <- strsplit(a, "=", fixed=TRUE)[[1]]
  if (length(kv) == 2) opt[[sub("^--?","",kv[1])]] <- kv[2]
}
if (suppressWarnings(!is.na(as.integer(opt$sample)))) {
  opt$sample <- as.integer(opt$sample)
}

opt$any_code <- as.integer(opt$any_code)
dir.create(dirname(opt$out_csv), recursive = TRUE, showWarnings = FALSE)

# -- load DB index --
if (!file.exists(opt$db)) stop("DB not found: ", opt$db)
db <- readRDS(opt$db)
need <- c("sample_ids","locus_ids","A1","A2")
if (!all(need %in% names(db))) stop("invalid db_index structure")

S <- length(db$sample_ids)
if (is.integer(opt$sample)) {
  idx <- opt$sample
  if (idx < 1L || idx > S) stop("sample index out of range: ", idx)
} else {
  idx <- match(as.character(opt$sample), db$sample_ids)
  if (is.na(idx)) stop("SampleID not found: ", opt$sample)
}

# -- get one sample across loci --
a1 <- db$A1[idx, , drop=TRUE]
a2 <- db$A2[idx, , drop=TRUE]
loci <- as.character(db$locus_ids)

# optional: enforce allele1 <= allele2 and ANY on right
ANY <- opt$any_code
swap <- (a1 == ANY & a2 != ANY) | (a1 != ANY & a2 != ANY & a1 > a2)
if (any(swap, na.rm = TRUE)) {
  tmp <- a1[swap]; a1[swap] <- a2[swap]; a2[swap] <- tmp
}

# build df
df <- data.frame(
  Locus   = loci,
  allele1 = as.integer(a1),
  allele2 = as.integer(a2),
  stringsAsFactors = FALSE
)

# sanity: if locus_order.rds exists, reorder to match (should already match)
if (file.exists("data/locus_order.rds")) {
  locus_order <- readRDS("data/locus_order.rds")
  o <- match(locus_order, df$Locus)
  if (all(!is.na(o))) df <- df[o, , drop=FALSE]
}

write.csv(df, file = opt$out_csv, row.names = FALSE)
cat("[OK] wrote query CSV from DB: ", opt$out_csv,
    " (sample=", if (is.integer(opt$sample)) idx else opt$sample, ")\n", sep = "")

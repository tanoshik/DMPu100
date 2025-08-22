# scripts/devtools/make_virtual_db_small.R
# Generate a small virtual database (indexed for run_match_fast).
# Output: RDS (list(sample_ids, locus_ids, A1[int], A2[int])) ready for matcher.

suppressWarnings({suppressPackageStartupMessages({
  # microbenchmark は不要。依存なしで動くようにしています。
})})

ANY_CODE <- 9999L

# ---- helpers ----
.detect_freq_col <- function(df) {
  cand <- tolower(names(df))
  hit  <- which(cand %in% c("frequency","freq","p","prob","probability"))
  if (length(hit) == 0) stop("No frequency column found in freq_table.rds")
  names(df)[hit[1]]
}

.draw_genotype_one <- function(freq_df, locus) {
  # locusごとに Allele / freq を取り出して 2アレルを独立サンプル（HW前提の簡易）
  ln <- if ("Locus" %in% names(freq_df)) "Locus" else if ("locus" %in% names(freq_df)) "locus" else stop("freq table must have Locus")
  an <- if ("Allele" %in% names(freq_df)) "Allele" else if ("allele" %in% names(freq_df)) "allele" else stop("freq table must have Allele")
  fc <- .detect_freq_col(freq_df)
  
  sub <- freq_df[freq_df[[ln]] == locus, c(an, fc), drop = FALSE]
  if (!nrow(sub)) stop("No rows in freq table for locus: ", locus)
  
  alle <- as.character(sub[[an]])
  p    <- suppressWarnings(as.numeric(sub[[fc]]))
  if (!all(is.finite(p))) stop("Non-numeric freq at locus: ", locus)
  p <- p / sum(p)
  
  a1 <- sample(alle, size = 1L, prob = p, replace = TRUE)
  a2 <- sample(alle, size = 1L, prob = p, replace = TRUE)
  
  # normalize order: numeric ascending (no ANY in synthetic generation)
  n1 <- suppressWarnings(as.integer(a1))
  n2 <- suppressWarnings(as.integer(a2))
  if (!is.finite(n1) || !is.finite(n2)) stop("Non-integer allele in freq table at locus: ", locus)
  
  if (n1 <= n2) c(n1, n2) else c(n2, n1)
}

.build_index_small <- function(n = 1000, locus_order, freq_df, seed = 123) {
  set.seed(seed)
  S <- as.integer(n)
  L <- as.integer(length(locus_order))
  if (L == 0L) stop("locus_order is empty")
  
  A1 <- matrix(ANY_CODE, nrow = S, ncol = L)
  A2 <- matrix(ANY_CODE, nrow = S, ncol = L)
  
  for (j in seq_len(L)) {
    loc <- locus_order[j]
    for (i in seq_len(S)) {
      geno <- .draw_genotype_one(freq_df, loc)
      A1[i, j] <- geno[1]
      A2[i, j] <- geno[2]
    }
  }
  sample_ids <- sprintf("V%07d", seq_len(S))
  list(sample_ids = sample_ids, locus_ids = as.character(locus_order), A1 = A1, A2 = A2)
}

# ---- main ----
make_virtual_db_small <- function(n = 1000, seed = 123,
                                  out = sprintf("data/virtual_db_u100_S%d_seed%d.rds", n, seed),
                                  freq_rds = "data/freq_table.rds",
                                  locus_rds = "data/locus_order.rds") {
  # load deps
  if (!file.exists(freq_rds))  stop("freq_table.rds not found: ", freq_rds)
  if (!file.exists(locus_rds)) stop("locus_order.rds not found: ", locus_rds)
  
  freq_df    <- readRDS(freq_rds)
  locus_order <- readRDS(locus_rds)
  
  # build
  db_index <- .build_index_small(n = n, locus_order = locus_order, freq_df = freq_df, seed = seed)
  
  # save
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(db_index, out, compress = "xz")
  message(sprintf("[virtual-db small] wrote: %s (S=%d, L=%d)", out, length(db_index$sample_ids), length(db_index$locus_ids)))
  return(out)
}

# ---- auto-run on source (set to FALSE if you don't want) ----
AUTO_RUN <- TRUE
if (isTRUE(AUTO_RUN)) {
  try({
    p <- make_virtual_db_small()
    message("[virtual-db small] auto-run done: ", p)
  }, silent = FALSE)
}

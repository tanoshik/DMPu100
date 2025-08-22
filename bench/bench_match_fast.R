# bench/bench_match_fast.R
# Minimal bench harness for run_match_fast() with matrix-index DB.

suppressPackageStartupMessages({
  library(microbenchmark)
})

source("scripts/scoring_fast.R", local = TRUE)
source("scripts/matcher_fast.R", local = TRUE)

# tiny query builder (from locus_order, all-any except one locus example)
build_query_any_plus_one <- function(locus_order, picked, a1, a2, ANY_CODE = 9999L) {
  stopifnot(length(locus_order) > 0)
  df <- data.frame(
    Locus   = locus_order,
    allele1 = rep(ANY_CODE, length(locus_order)),
    allele2 = rep(ANY_CODE, length(locus_order)),
    stringsAsFactors = FALSE
  )
  df$allele1[df$Locus == picked] <- as.integer(a1)
  df$allele2[df$Locus == picked] <- as.integer(a2)
  df
}

run_one_bench <- function(db_rds, picked_locus = "D3S1358", a1 = 15, a2 = 17) {
  dbx <- readRDS(db_rds)
  stopifnot(is.list(dbx), is.matrix(dbx$A1), is.matrix(dbx$A2))
  lo  <- as.character(dbx$locus_ids)
  qdf <- build_query_any_plus_one(lo, picked_locus, a1, a2)
  
  mb <- microbenchmark(
    match = {
      rs <- run_match_fast(qdf, dbx, pre_add_for_any_any = 2L, include_bits_in_detail = FALSE)
    },
    times = 10L
  )
  print(mb)
  rs
}

# Example usage:
# rs <- run_one_bench("data/virtual_db_u100_S100000_seed20250822.rds")
# head(rs$summary)

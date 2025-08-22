# scripts/devtools/make_virtual_db_large.R
# Chunked generator for very large virtual DBs using freq_table.rds.
# Writes several chunk RDS and a manifest; later we can bind them or slice subsets.

suppressPackageStartupMessages({
  library(tools)
})

ANY_CODE <- 9999L

detect_freq_col <- function(df) {
  cn <- tolower(names(df))
  hits <- which(cn %in% c("frequency","freq","p","prob","probability"))
  if (length(hits)) names(df)[hits[1]] else stop("no frequency column in freq_table")
}

draw_geno_vec <- function(freq_df, locus, fcol, n) {
  sub <- freq_df[freq_df$Locus == locus, , drop=FALSE]
  p   <- suppressWarnings(as.numeric(sub[[fcol]]))
  a   <- as.character(sub$Allele); p <- p / sum(p)
  a1 <- sample(a, size = n, prob = p, replace = TRUE)
  a2 <- sample(a, size = n, prob = p, replace = TRUE)
  n1 <- suppressWarnings(as.numeric(a1))
  n2 <- suppressWarnings(as.numeric(a2))
  swap <- which(n1 > n2)
  if (length(swap)) {
    t <- n1[swap]; n1[swap] <- n2[swap]; n2[swap] <- t
  }
  cbind(as.integer(n1), as.integer(n2))
}

make_virtual_db_large <- function(n_total = 2e6, chunk_size = 100000,
                                  seed = 20250822,
                                  freq_rds = "data/freq_table.rds",
                                  locus_rds = "data/locus_order.rds",
                                  out_prefix = "data/virtual_db_u100_large") {
  stopifnot(file.exists(freq_rds), file.exists(locus_rds))
  set.seed(seed)
  ft <- readRDS(freq_rds)
  lo <- readRDS(locus_rds)
  fcol <- detect_freq_col(ft)
  
  L <- length(lo)
  nchunks <- ceiling(n_total / chunk_size)
  manifest <- data.frame(
    chunk = integer(0), file = character(0), rows = integer(0),
    stringsAsFactors = FALSE
  )
  dir.create("data", showWarnings = FALSE, recursive = TRUE)
  
  for (k in seq_len(nchunks)) {
    n <- if (k < nchunks) chunk_size else (n_total - chunk_size * (nchunks - 1L))
    SIDs <- sprintf("V%07d", (chunk_size*(k-1L) + 1L):(chunk_size*(k-1L) + n))
    A1 <- matrix(ANY_CODE, nrow = n, ncol = L)
    A2 <- matrix(ANY_CODE, nrow = n, ncol = L)
    
    for (j in seq_along(lo)) {
      loc <- lo[j]
      g   <- draw_geno_vec(ft, loc, fcol, n)
      A1[, j] <- g[, 1]; A2[, j] <- g[, 2]
    }
    
    dbx <- list(sample_ids = SIDs, locus_ids = lo, A1 = A1, A2 = A2)
    out_rds <- sprintf("%s_chunk%03d.rds", out_prefix, k)
    saveRDS(dbx, out_rds)
    manifest[nrow(manifest)+1, ] <- list(k, out_rds, n)
    message(sprintf("chunk %d/%d -> %s (n=%d)", k, nchunks, out_rds, n))
  }
  
  manfile <- sprintf("%s_manifest.csv", out_prefix)
  write.csv(manifest, manfile, row.names = FALSE)
  message("manifest: ", manfile)
  invisible(manfile)
}

if (sys.nframe() == 0) {
  make_virtual_db_large()
}

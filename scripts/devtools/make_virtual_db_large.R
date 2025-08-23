# scripts/devtools/make_virtual_db_large.R
# Generate a large virtual DB (u100 format) from allele frequency table.
# - Uses chunked sampling per locus
# - Stops early with clear messages if freq table or locus order missing
# - No multibyte chars in code/comments.

# ====== CONFIG (edit here) ======
N_TOTAL    <- 2000000L      # total #samples
CHUNK_SIZE <- 100000L       # rows per chunk
SEED       <- 123L          # base seed
ANY_CODE   <- 9999L
AUTO_RUN   <- TRUE          # source() -> auto run

# output path (gzip compressed RDS)
OUT_PATH <- file.path("data",
                      sprintf("virtual_db_u100_S%d_seed%d.rds", as.integer(N_TOTAL), as.integer(SEED))
)

# ====== helpers ======
.ts <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S] ")

.fail <- function(msg) {
  message(.ts(), " [virtual-db large] FAILED: ", msg)
  stop(msg, call. = FALSE)
}

.detect_freq_col <- function(df) {
  cn <- tolower(names(df))
  hit <- which(cn %in% c("frequency", "freq", "p", "prob", "probability"))
  if (length(hit) == 0L) return(NULL)
  names(df)[hit[1L]]
}

load_freq_df <- function() {
  # prefer RDS
  if (file.exists("data/freq_table.rds")) {
    ft <- tryCatch(readRDS("data/freq_table.rds"), error = function(e) NULL)
    if (is.data.frame(ft) && nrow(ft) > 0L) return(ft)
  }
  # fallback CSV (DMP legacy name)
  if (file.exists("data/Allele-count_GF_Japanese.csv")) {
    ft <- tryCatch(
      read.csv("data/Allele-count_GF_Japanese.csv", stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.data.frame(ft) && nrow(ft) > 0L) return(ft)
  }
  NULL
}

load_locus_order <- function(freq_df) {
  if (file.exists("data/locus_order.rds")) {
    lo <- tryCatch(readRDS("data/locus_order.rds"), error = function(e) NULL)
    if (is.null(lo)) return(unique(as.character(freq_df$Locus)))
    return(as.character(lo))
  }
  # fallback to freq order
  unique(as.character(freq_df$Locus))
}

build_per_locus_prob <- function(freq_df, locus_ids) {
  fcol <- .detect_freq_col(freq_df)
  if (is.null(fcol)) .fail("frequency column not found in freq table")
  if (!all(c("Locus", "Allele") %in% names(freq_df))) {
    .fail("freq table must contain columns: Locus, Allele, and a frequency column")
  }
  # per-locus list: list of list(alleles=int, probs=numeric)
  lapply(locus_ids, function(L) {
    rows <- which(as.character(freq_df$Locus) == as.character(L))
    if (length(rows) == 0L) .fail(paste0("no frequency rows for locus: ", L))
    alle <- suppressWarnings(as.integer(freq_df$Allele[rows]))
    pr   <- suppressWarnings(as.numeric(freq_df[[fcol]][rows]))
    # sanitize
    ok <- is.finite(alle) & is.finite(pr) & pr > 0
    alle <- alle[ok]; pr <- pr[ok]
    if (length(alle) == 0L) .fail(paste0("empty/invalid frequency entries at locus: ", L))
    # normalize probs
    s <- sum(pr)
    if (!is.finite(s) || s <= 0) .fail(paste0("non-positive prob sum at locus: ", L))
    pr <- pr / s
    list(alleles = alle, probs = pr)
  })
}

# sample 2 alleles (diploid) for m individuals at one locus
sample_two_alleles <- function(m, alleles, probs) {
  # draw 2*m alleles i.i.d., then split to A1/A2
  idx <- sample.int(length(alleles), size = 2L * m, replace = TRUE, prob = probs)
  a   <- alleles[idx]
  list(A1 = a[seq_len(m)], A2 = a[m + seq_len(m)])
}

# ====== main generator ======
make_virtual_db_large <- function(n_total = N_TOTAL, chunk_size = CHUNK_SIZE,
                                  seed = SEED, out_path = OUT_PATH,
                                  any_code = ANY_CODE) {
  
  if (!dir.exists("data")) dir.create("data", recursive = TRUE, showWarnings = FALSE)
  
  message(.ts(), " virtual-db large start: n=", n_total, ", chunk=", chunk_size, ", seed=", seed)
  
  freq_df <- load_freq_df()
  if (is.null(freq_df) || !is.data.frame(freq_df) || nrow(freq_df) == 0L) {
    .fail("freq table not found or empty (data/freq_table.rds or data/Allele-count_GF_Japanese.csv)")
  }
  
  locus_ids <- load_locus_order(freq_df)
  if (length(locus_ids) == 0L) .fail("locus order missing and could not infer from freq table")
  
  L <- length(locus_ids)
  per_locus <- build_per_locus_prob(freq_df, locus_ids)
  
  # Pre-allocate full matrices and fill by chunks
  S <- as.integer(n_total)
  A1 <- matrix(any_code, nrow = S, ncol = L,
               dimnames = list(NULL, locus_ids))
  A2 <- matrix(any_code, nrow = S, ncol = L,
               dimnames = list(NULL, locus_ids))
  
  # chunked fill
  n_chunks <- ceiling(S / chunk_size)
  t0 <- proc.time()[3L]
  
  for (ck in seq_len(n_chunks)) {
    i1 <- (ck - 1L) * chunk_size + 1L
    i2 <- min(S, ck * chunk_size)
    m  <- i2 - i1 + 1L
    set.seed(seed + ck) # deterministic per chunk
    
    # per locus sampling for this chunk
    for (j in seq_len(L)) {
      pl <- per_locus[[j]]
      samp <- sample_two_alleles(m, pl$alleles, pl$probs)
      A1[i1:i2, j] <- samp$A1
      A2[i1:i2, j] <- samp$A2
    }
    
    message(.ts(), " chunk ", ck, "/", n_chunks, " done (rows ", i1, "..", i2, ")")
  }
  
  # quick sanity (not all 9999)
  any_ratio <- (sum(A1 == any_code) + sum(A2 == any_code)) / (length(A1) + length(A2))
  if (any_ratio > 0.01) {
    message(.ts(), " [warn] ANY_CODE rate = ", sprintf("%.4f", any_ratio),
            " (expected near 0 with proper freq table)")
  }
  
  sample_ids <- sprintf("V%07d", seq_len(S))
  out <- list(sample_ids = sample_ids,
              locus_ids  = as.character(locus_ids),
              A1 = A1, A2 = A2)
  
  saveRDS(out, file = out_path, compress = "gzip")
  t1 <- proc.time()[3L]
  message(.ts(), " wrote: ", out_path,
          " (S=", S, ", L=", L, ", compress=gzip) ")
  message(.ts(), " elapsed: ", sprintf("%.2f", t1 - t0), " sec (",
          sprintf("%.2f", (t1 - t0)/60), " min) ")
  
  invisible(out_path)
}

# ====== autorun on source() ======
if (isTRUE(AUTO_RUN)) {
  tryCatch({
    path <- make_virtual_db_large()
    message(.ts(), " [virtual-db large] auto-run done: ", path)
  }, error = function(e) {
    message(.ts(), " [virtual-db large] FAILED: ", conditionMessage(e))
  })
}

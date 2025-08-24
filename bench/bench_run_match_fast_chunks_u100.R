# scripts/bench/bench_run_match_fast_chunks_u100.R
# Benchmark run_match_fast() with chunked DB and configurable workers.
# - No multibyte characters in code/comments.
# - Working directory assumed: project root (DMPu100).
# - Logging: append to output/bench/bench_chunks_log.csv (resumable).
# - Safe for small "smoke test" and scalable to 1M by config.
#
# Usage:
#   source("scripts/bench/bench_run_match_fast_chunks_u100.R")

suppressPackageStartupMessages({
  # base R only; optionally parallel if available
})

# -------------------------
# CONFIG (edit here)
# -------------------------
CONFIG <- list(
  # DB RDS path (relative to project root)
  db_rds_path   = "data/virtual_db_u100_S1000_seed123.rds",  # replace for 1M bench
  # Optional: query selection strategy: "first", "random", or "fixed_id"
  query_pick    = "first",
  fixed_query_id = NA_character_,
  # chunk sizes for test (records per chunk; for 1M bench use larger, e.g., 100000, 200000, 250000, 500000)
  chunk_sizes   = c(200L, 300L, 500L),
  # workers set; for light test keep 1; for 1M try c(1,2,4)
  workers_list  = c(1L),
  # repeat count per setting (stability)
  repeats       = 1L,
  # output log csv
  log_csv       = "output/bench/bench_chunks_log.csv",
  # optional detail out (disabled by default for speed)
  write_detail  = FALSE,
  # detail path pattern if enabled
  detail_dir    = "output/bench/details",
  # random seed for stable query selection if needed
  seed          = 123L,
  # enable resume: skip if (db_file, chunk_size, workers, run_idx) already logged
  enable_resume = TRUE
)

# -------------------------
# Helpers
# -------------------------

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_now_str <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

sha7_head <- function() {
  h <- tryCatch(system("git rev-parse --short=7 HEAD", intern = TRUE, ignore.stderr = TRUE),
                error = function(e) NA_character_)
  if (is.character(h) && length(h) == 1L && nchar(h) > 0) h else "unknown"
}

split_indices <- function(n, chunk_size) {
  if (n <= 0 || chunk_size <= 0) return(list())
  starts <- seq.int(1L, n, by = chunk_size)
  ends   <- pmin(starts + chunk_size - 1L, n)
  Map(seq.int, starts, ends)
}

read_rds_any <- function(path) {
  readRDS(path)
}

# Try to infer locus count from locus_order.rds
get_locus_count <- function() {
  p <- "data/locus_order.rds"
  if (!file.exists(p)) return(NA_integer_)
  loci <- tryCatch(readRDS(p), error = function(e) NULL)
  if (is.null(loci)) return(NA_integer_)
  as.integer(length(loci))
}

ensure_sources <- function() {
  # Load project functions needed by run_match_fast
  # Adjust if your file names differ.
  srcs <- c(
    "scripts/utils_profile.R",
    "scripts/utils_profile_modular.R",
    "scripts/scoring_fast.R",
    "scripts/matcher.R",
    "scripts/io_profiles.R",
    "scripts/utils_freq.R"
  )
  for (s in srcs) {
    if (file.exists(s)) {
      tryCatch(source(s), error = function(e) {
        message("[WARN] Failed to source: ", s, " : ", conditionMessage(e))
      })
    }
  }
}

pick_query_from_db <- function(db_df) {
  # Expect a long DF with columns: SampleID, Locus, allele1, allele2 (case-insensitive tolerated upstream)
  # We will pick one SampleID and build a query list/df compatible with run_match_fast().
  # Strategy controlled by CONFIG$query_pick.
  sid_col <- "SampleID"
  if (!sid_col %in% names(db_df)) {
    # try case-insensitive
    nms <- names(db_df)
    hit <- nms[tolower(nms) == "sampleid"]
    if (length(hit) == 1L) names(db_df)[names(db_df) == hit] <- "SampleID"
  }
  stopifnot("SampleID" %in% names(db_df))
  set.seed(CONFIG$seed)
  cand_ids <- unique(db_df$SampleID)
  
  qid <- switch(
    CONFIG$query_pick,
    "random"   = sample(cand_ids, 1L),
    "fixed_id" = {
      if (isTRUE(nzchar(CONFIG$fixed_query_id)) && CONFIG$fixed_query_id %in% cand_ids) CONFIG$fixed_query_id else cand_ids[1L]
    },
    "first"    = cand_ids[1L],
    cand_ids[1L]
  )
  
  q_df <- db_df[db_df$SampleID == qid, , drop = FALSE]
  list(query_sample_id = qid, query_df = q_df)
}

prepare_db_for_match <- function(db_df, qid) {
  # Remove query sample from DB to avoid self-match bias
  db_df[db_df$SampleID != qid, , drop = FALSE]
}

# Estimate "loci processed" as (#loci in locus_order) * (#rows processed per sample pair)
# For a quick throughput proxy we use loci_count * num_pairs_processed, where
# num_pairs_processed approximates n_chunk_samples * 1 query.
estimate_loci_processed <- function(n_records_chunk, loci_count) {
  if (is.na(loci_count)) return(NA_real_)
  as.numeric(loci_count) * as.numeric(n_records_chunk)
}

append_log <- function(csv_path, row_df) {
  safe_dir_create(dirname(csv_path))
  if (!file.exists(csv_path)) {
    utils::write.table(row_df, file = csv_path, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    suppressWarnings(utils::write.table(row_df, file = csv_path, sep = ",",
                                        row.names = FALSE, col.names = FALSE, append = TRUE))
  }
}

already_logged <- function(csv_path, key_cols, key_vals) {
  if (!file.exists(csv_path)) return(FALSE)
  df <- tryCatch(utils::read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df)) return(FALSE)
  if (!all(key_cols %in% names(df))) return(FALSE)
  # exact match on all key columns
  sel <- rep(TRUE, nrow(df))
  for (i in seq_along(key_cols)) {
    col <- key_cols[i]
    val <- key_vals[[i]]
    sel <- sel & (as.character(df[[col]]) == as.character(val))
  }
  any(sel)
}

# -------------------------
# Main
# -------------------------

ensure_sources()

db_path <- CONFIG$db_rds_path
if (!file.exists(db_path)) {
  stop(sprintf("DB RDS not found: %s (wd=%s)", db_path, getwd()))
}

message("[INFO] Loading DB: ", db_path)
db_obj <- read_rds_any(db_path)

# Accept either long DF directly or list/rds with element "db_df"
if (is.data.frame(db_obj)) {
  db_df <- db_obj
} else if (is.list(db_obj) && !is.null(db_obj$db_df) && is.data.frame(db_obj$db_df)) {
  db_df <- db_obj$db_df
} else {
  stop("Unsupported DB RDS structure. Expect data.frame or list with $db_df.")
}

# Normalize column names minimal (allele1/allele2 case)
nms <- names(db_df)
names(db_df)[tolower(nms) == "allele1"] <- "allele1"
names(db_df)[tolower(nms) == "allele2"] <- "allele2"
names(db_df)[tolower(nms) == "locus"]   <- "Locus"
names(db_df)[tolower(nms) == "sampleid"]<- "SampleID"

stopifnot(all(c("SampleID", "Locus", "allele1", "allele2") %in% names(db_df)))

# Pick query
pq <- pick_query_from_db(db_df)
qid <- pq$query_sample_id
q_df <- pq$query_df

# Prepare DB excluding query sample
db_df2 <- prepare_db_for_match(db_df, qid)

loci_count <- get_locus_count()
git_hash   <- sha7_head()

message("[INFO] Query SampleID: ", qid)
message("[INFO] Records (DB excl. query): ", nrow(db_df2))
message("[INFO] Loci count (from locus_order.rds): ", loci_count)
message("[INFO] Git short hash: ", git_hash)

# Index per SampleID to chunk by samples (prefer even split by samples)
sample_ids <- unique(db_df2$SampleID)
n_samples  <- length(sample_ids)

# Build combinations
combis <- expand.grid(
  chunk_size = as.integer(CONFIG$chunk_sizes),
  workers    = as.integer(CONFIG$workers_list),
  run_idx    = seq_len(as.integer(CONFIG$repeats)),
  stringsAsFactors = FALSE
)

safe_dir_create(dirname(CONFIG$log_csv))
if (isTRUE(CONFIG$write_detail)) safe_dir_create(CONFIG$detail_dir)

for (row in seq_len(nrow(combis))) {
  cs   <- combis$chunk_size[row]
  wks  <- combis$workers[row]
  ridx <- combis$run_idx[row]
  
  key_cols <- c("db_file", "chunk_size", "workers", "run_idx")
  key_vals <- list(basename(db_path), cs, wks, ridx)
  if (CONFIG$enable_resume && already_logged(CONFIG$log_csv, key_cols, key_vals)) {
    message(sprintf("[SKIP] Logged already: db=%s chunk=%d workers=%d run=%d",
                    basename(db_path), cs, wks, ridx))
    next
  }
  
  # Build chunk list by sample IDs
  # Approximate records per sample: use split to get counts
  # For chunking by samples, we split indices of sample_ids by chunk size in samples
  idx_chunks <- split_indices(n_samples, cs)
  n_chunks   <- length(idx_chunks)
  if (n_chunks == 0L) {
    message("[WARN] No chunks for chunk_size=", cs, " (n_samples=", n_samples, ")")
    next
  }
  
  # Prepare query structure for run_match_fast
  # If run_match_fast expects certain structure, adjust here as needed.
  query_profile <- q_df
  
  # Execution start
  t0 <- proc.time()[["elapsed"]]
  total_rows <- 0L
  total_loci <- 0
  detail_rows <- list()
  
  for (ci in seq_len(n_chunks)) {
    sidx <- idx_chunks[[ci]]
    sids <- sample_ids[sidx]
    sub_df <- db_df2[db_df2$SampleID %in% sids, , drop = FALSE]
    
    # Run matching (single-thread by default; multi-worker orchestration can be added as needed)
    # Here we assume run_match_fast(query_df, db_df) returns a data.frame with per-row scores
    # Replace with the actual call signature used in your project.
    t_chunk0 <- proc.time()[["elapsed"]]
    res <- NULL
    err <- NULL
    try({
      # If parallel workers >1 are desired for 1M bench, orchestrate here.
      # For the light test we run single call. You can parallelize by splitting sub_df further.
      res <- run_match_fast(query_profile = query_profile, db_profiles_df = sub_df)
    }, silent = TRUE)
    if (is.null(res)) {
      err <- "run_match_fast_failed"
      n_res <- 0L
    } else {
      n_res <- nrow(res)
    }
    t_chunk1 <- proc.time()[["elapsed"]]
    elapsed_chunk <- t_chunk1 - t_chunk0
    
    # Estimate loci processed (proxy)
    loci_proc <- estimate_loci_processed(n_records_chunk = nrow(sub_df), loci_count = loci_count)
    
    total_rows <- total_rows + n_res
    total_loci <- if (is.na(loci_proc)) total_loci else total_loci + loci_proc
    
    # Optional detail write (disabled by default)
    if (isTRUE(CONFIG$write_detail) && n_res > 0L) {
      fdetail <- file.path(CONFIG$detail_dir,
                           sprintf("detail_%s_cs%d_w%d_run%d_chunk%03d.csv",
                                   tools::file_path_sans_ext(basename(db_path)),
                                   cs, wks, ridx, ci))
      utils::write.csv(res, fdetail, row.names = FALSE)
    }
    
    # Per-chunk log append (progress visibility)
    per_chunk <- data.frame(
      timestamp    = safe_now_str(),
      git_hash     = git_hash,
      db_file      = basename(db_path),
      query_id     = as.character(qid),
      chunk_size   = as.integer(cs),
      workers      = as.integer(wks),
      run_idx      = as.integer(ridx),
      chunk_id     = as.integer(ci),
      chunks_total = as.integer(n_chunks),
      n_db_rows    = as.integer(nrow(sub_df)),
      n_results    = as.integer(n_res),
      sec_chunk    = as.numeric(elapsed_chunk),
      rows_per_sec = if (elapsed_chunk > 0) nrow(sub_df) / elapsed_chunk else NA_real_,
      loci_proc    = as.numeric(loci_proc),
      loci_per_sec = if (!is.na(loci_proc) && elapsed_chunk > 0) loci_proc / elapsed_chunk else NA_real_,
      status       = if (is.null(err)) "OK" else err,
      stringsAsFactors = FALSE
    )
    append_log(CONFIG$log_csv, per_chunk)
  }
  
  t1 <- proc.time()[["elapsed"]]
  elapsed_total <- t1 - t0
  
  # Final summary row for this combination
  final_row <- data.frame(
    timestamp    = safe_now_str(),
    git_hash     = git_hash,
    db_file      = basename(db_path),
    query_id     = as.character(qid),
    chunk_size   = as.integer(cs),
    workers      = as.integer(wks),
    run_idx      = as.integer(ridx),
    chunk_id     = as.integer(NA),
    chunks_total = as.integer(n_chunks),
    n_db_rows    = as.integer(NA),
    n_results    = as.integer(total_rows),
    sec_chunk    = as.numeric(elapsed_total),
    rows_per_sec = as.numeric(NA),
    loci_proc    = as.numeric(total_loci),
    loci_per_sec = if (!is.na(total_loci) && elapsed_total > 0) total_loci / elapsed_total else NA_real_,
    status       = "SUMMARY",
    stringsAsFactors = FALSE
  )
  append_log(CONFIG$log_csv, final_row)
  
  message(sprintf("[DONE] db=%s chunk=%d workers=%d run=%d elapsed=%.2fs",
                  basename(db_path), cs, wks, ridx, elapsed_total))
}

message("[OK] Benchmark finished. Log: ", CONFIG$log_csv)

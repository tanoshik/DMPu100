# scripts/bench/bench_run_match_fast_chunks_u100.R
# Chunked bench runner for run_match_fast() with resumable logging.
# No multibyte chars in code/comments.

suppressWarnings({
  options(stringsAsFactors = FALSE, scipen = 999)
})

# ========= User Config (edit here) =========
cfg <- list(
  # DB RDS path (relative to project root)
  db_rds_candidates = c(
    "output/virtual/virtual_db_u100_S1000_seed123.rds",  # light test
    "data/virtual_db_u100_S1000_seed123.rds"             # fallback
    # For 1M bench, replace with your real file below, e.g.:
    # "output/virtual/virtual_db_u100_S1000000_seed123.rds"
  ),
  
  # Query source: "first_sample_in_db" or path to a prepared query RDS/CSV
  query_source = "first_sample_in_db",   # or set to "data/query_profile.rds" etc.
  
  # Chunk sizes to try (light test will iterate these)
  chunk_sizes = c(200, 300, 500),
  
  # Workers (sequential only for now; keep 1)
  workers = 1,
  
  # Output dir for logs
  out_dir = "output/bench",
  
  # Benchmark label (used in log filename)
  bench_label = "u100_chunks",
  
  # Whether to write detail CSV for each chunk (heavy; keep FALSE for benches)
  write_detail = FALSE
)
# ==========================================

# ---- helpers ----
safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp_compact <- function() format(Sys.time(), "%Y%m%d%H%M%S")

find_first_existing <- function(paths) {
  for (p in paths) if (file.exists(p)) return(p)
  return(NA_character_)
}

require_script <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing script: %s", path))
  source(path, local = TRUE)
}

msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", sprintf(...), "\n", sep = "")

# ---- load core project scripts (relative) ----
require_script("scripts/scoring_fast.R")          # expects run_match_fast(...)
# If needed by your run_match_fast environment, also load modular utils:
# require_script("scripts/utils_profile_modular.R")

# ---- load DB ----
db_rds_path <- find_first_existing(cfg$db_rds_candidates)
if (is.na(db_rds_path)) {
  stop("DB RDS not found. Please update cfg$db_rds_candidates.")
}
msg("DB file: %s", db_rds_path)
db_obj <- readRDS(db_rds_path)

# Expect long-form data.frame: SampleID, Locus, allele1, allele2 (case-insensitive ok)
normalize_db_cols <- function(df) {
  nms <- names(df)
  map <- list(
    SampleID = c("^sampleid$","^SampleID$"),
    Locus    = c("^locus$","^Locus$"),
    allele1  = c("^allele1$","^Allele1$"),
    allele2  = c("^allele2$","^Allele2$")
  )
  for (std in names(map)) {
    if (!any(grepl(paste(map[[std]], collapse="|"), nms))) {
      stop(sprintf("DB missing required column %s (case-insensitive)", std))
    }
    hit <- which(grepl(paste(map[[std]], collapse="|"), nms))[1]
    names(df)[hit] <- std
  }
  df
}

if (is.data.frame(db_obj)) {
  db_df <- normalize_db_cols(db_obj)
} else if (is.list(db_obj)) {
  # If list, try to rbind
  df_list <- lapply(db_obj, function(x) {
    if (!is.data.frame(x)) return(NULL)
    normalize_db_cols(x)
  })
  db_df <- do.call(rbind, df_list)
  if (!is.data.frame(db_df)) stop("Unsupported DB structure (list not of data.frames).")
} else {
  stop("Unsupported DB object.")
}

# Collect sample ids
if (!"SampleID" %in% names(db_df)) stop("DB frame has no SampleID.")
sample_ids <- unique(db_df$SampleID)
n_samples  <- length(sample_ids)
msg("DB samples: %d", n_samples)

# ---- build query ----
build_query_from_first_sample <- function(db_df) {
  sid <- sample_ids[1]
  qdf <- db_df[db_df$SampleID == sid, c("Locus","allele1","allele2")]
  # Ensure correct types
  qdf$allele1 <- as.character(qdf$allele1)
  qdf$allele2 <- as.character(qdf$allele2)
  qdf
}

load_query <- function(source) {
  if (identical(source, "first_sample_in_db")) return(build_query_from_first_sample(db_df))
  if (!file.exists(source)) stop(sprintf("Query source not found: %s", source))
  obj <- readRDS(source)
  if (is.data.frame(obj)) {
    # Expect Locus, allele1, allele2
    need <- c("Locus","allele1","allele2")
    miss <- setdiff(need, names(obj))
    if (length(miss)) stop(sprintf("Query frame missing cols: %s", paste(miss, collapse=",")))
    return(obj[, need])
  }
  stop("Unsupported query source.")
}

query_df <- load_query(cfg$query_source)
msg("Query loci: %d", nrow(query_df))

# ---- bench core ----
safe_mkdir(cfg$out_dir)

make_log_path <- function(chunk_size) {
  sprintf("%s/bench_%s_cs%s_w%d_%s.csv",
          cfg$out_dir, cfg$bench_label, chunk_size, cfg$workers, timestamp_compact())
}

# Resume support: read last completed chunk index from existing log (if any).
# We treat a "series" as same filename; to resume, point to the same log file.
# For simple usage, each run creates a new log. For resume, re-use the same path.
get_last_completed_idx <- function(log_path) {
  if (!file.exists(log_path)) return(0L)
  x <- tryCatch(read.csv(log_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(0L)
  if (!"chunk_idx" %in% names(x)) return(0L)
  max(as.integer(x$chunk_idx), na.rm = TRUE)
}

append_log <- function(log_path, row) {
  write.table(row, file = log_path, sep = ",", row.names = FALSE, col.names = !file.exists(log_path), append = TRUE)
}

calc_kpis <- function(chunk_df, qdf, res_df, elapsed_sec) {
  # Estimate counts
  # samples in chunk:
  ns <- length(unique(chunk_df$SampleID))
  # loci processed ~ rows in qdf * ns
  nl <- nrow(qdf) * ns
  list(
    n_samples = ns,
    n_loci = nl,
    sps = if (elapsed_sec > 0) ns / elapsed_sec else NA_real_,
    lps = if (elapsed_sec > 0) nl / elapsed_sec else NA_real_
  )
}

subset_db_by_chunk <- function(db_df, ids, idx, chunk_size) {
  start <- (idx - 1L) * chunk_size + 1L
  end   <- min(idx * chunk_size, length(ids))
  if (start > end) return(NULL)
  sel_ids <- ids[start:end]
  db_df[db_df$SampleID %in% sel_ids, , drop = FALSE]
}

run_one_chunk <- function(idx, chunk_size, ids, db_df, qdf, detail = FALSE) {
  chunk_df <- subset_db_by_chunk(db_df, ids, idx, chunk_size)
  if (is.null(chunk_df) || nrow(chunk_df) == 0) return(NULL)
  
  t0 <- proc.time()[3]
  # run_match_fast is expected from scripts/scoring_fast.R
  # It should accept prepared query df and prepared db df.
  res <- run_match_fast(qdf, chunk_df)
  t1 <- proc.time()[3]
  elapsed <- as.numeric(t1 - t0)
  
  k <- calc_kpis(chunk_df, qdf, res, elapsed)
  
  # Optional detail output per chunk (heavy)
  detail_path <- NA_character_
  if (detail && is.data.frame(res) && nrow(res)) {
    safe_mkdir(sprintf("%s/detail", cfg$out_dir))
    detail_path <- sprintf("%s/detail/detail_cs%d_chunk%05d_%s.csv",
                           cfg$out_dir, chunk_size, idx, timestamp_compact())
    try(write.csv(res, detail_path, row.names = FALSE), silent = TRUE)
  }
  
  list(
    elapsed = elapsed,
    n_samples = k$n_samples,
    n_loci = k$n_loci,
    sps = k$sps,
    lps = k$lps,
    detail_path = detail_path
  )
}

# ---- main routine ----
do_series_for_chunk_size <- function(chunk_size) {
  log_path <- make_log_path(chunk_size)
  msg("Log: %s", log_path)
  
  total_chunks <- ceiling(n_samples / chunk_size)
  done_idx <- get_last_completed_idx(log_path)
  if (done_idx > 0) msg("Resuming from chunk %d (of %d)", done_idx + 1L, total_chunks)
  
  total_rows_done <- 0L
  if (file.exists(log_path)) {
    old <- tryCatch(read.csv(log_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(old) && nrow(old) && "rows_done_cum" %in% names(old)) {
      total_rows_done <- max(as.integer(old$rows_done_cum), na.rm = TRUE)
    }
  }
  
  for (idx in seq.int(done_idx + 1L, total_chunks)) {
    msg("Chunk %d / %d (size=%d)", idx, total_chunks, chunk_size)
    info <- run_one_chunk(idx, chunk_size, sample_ids, db_df, query_df, detail = isTRUE(cfg$write_detail))
    if (is.null(info)) break
    
    rows_in <- nrow(subset_db_by_chunk(db_df, sample_ids, idx, chunk_size))
    total_rows_done <- total_rows_done + rows_in
    
    row <- data.frame(
      ts = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      chunk_idx = idx,
      chunk_size = chunk_size,
      rows_in = rows_in,
      elapsed_sec = round(info$elapsed, 6),
      samples = info$n_samples,
      loci = info$n_loci,
      samples_per_sec = round(info$sps, 3),
      loci_per_sec = round(info$lps, 3),
      rows_done_cum = total_rows_done,
      detail_csv = ifelse(is.na(info$detail_path), "", info$detail_path),
      stringsAsFactors = FALSE
    )
    append_log(log_path, row)
  }
  
  msg("Done chunk_size=%d. Log at: %s", chunk_size, log_path)
}

# ---- run light test series (S1000) ----
for (cs in cfg$chunk_sizes) {
  do_series_for_chunk_size(cs)
}

msg("All series finished.")

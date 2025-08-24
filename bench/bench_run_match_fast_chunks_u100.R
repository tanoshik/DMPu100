# No multibyte characters in code/comments.

# DMPu100 bench: incremental CSV append + resume + KPI (+ optional parallel)
# NOTE: run_match_fast() signature is based on CURRENT RUN LOGS (assumption).
#       If your local function signature differs, adjust the call near "worker" section.

########################################
# ========== CONFIG (EDIT HERE) ==========
########################################

# --- data & output ---
db_rds_path   <- "data/virtual_db_u100_S1000_seed123.rds"   # <= 1K quick test
result_dir    <- "test/bench"
log_dir       <- "logs"

# --- bench grid (small defaults for quick test) ---
chunk_sizes   <- c(100, 200, 250, 500)  # <= adjust for bigger DB later
enable_parallel   <- FALSE               # set TRUE to test parallel
parallel_workers  <- max(1L, parallel::detectCores() - 1L)   # ex: CPU-1

# --- matcher behavior ---
ANY_CODE      <- 9999L
include_bits_in_detail <- FALSE         # fixed FALSE for speed
pre_add_for_any_any    <- 2L            # design doc default

# --- metadata for KPI (旧DMP互換) ---
seed          <- 1L
bench_type    <- "gf_like"
cores         <- if (enable_parallel) parallel_workers else 1L
include_io    <- FALSE                  # change to TRUE if measuring IO separately

# --- optional: reuse same CSV file name to resume across sessions ---
# leave empty "" to auto-generate time-stamped name each run.
manual_run_tag <- ""   # e.g., "resume_1M_try1" to append into same CSV

########################################
# ======== PROJECT SOURCES (REQUIRED) ========
########################################

# These scripts must exist in your repo. If not, please provide them before running.
# - scripts/scoring_fast.R
# - scripts/matcher_fast.R
source("scripts/scoring_fast.R",  local = TRUE)
source("scripts/matcher_fast.R",  local = TRUE)

########################################
# ========== OPTIONAL DEPENDENCIES ==========
########################################
if (enable_parallel) {
  suppressWarnings(suppressMessages({
    library(future)
    library(future.apply)
  }))
}

########################################
# ========== HELPERS ==========
########################################

safe_dir_create <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
ts_tag <- function() format(Sys.time(), "%Y%m%d%H%M%S")

# map "any" to ANY_CODE (integer), keep integers as-is
.map_any <- function(x, any_code = ANY_CODE) {
  if (is.null(x)) return(integer(0))
  if (is.factor(x)) x <- as.character(x)
  out <- suppressWarnings(as.integer(x))
  is_any <- !is.na(x) & tolower(as.character(x)) == "any"
  out[is_any] <- any_code
  out
}

# minimal query shape: Locus, allele1, allele2  (+ sorted by locus_order if available)
.prepare_query_min <- function(q_raw, locus_order = NULL, any_code = ANY_CODE) {
  n <- names(q_raw)
  names(q_raw) <- sub("(?i)^locus$","Locus", n, perl=TRUE)
  names(q_raw) <- sub("(?i)^allele1$","allele1", names(q_raw), perl=TRUE)
  names(q_raw) <- sub("(?i)^allele2$","allele2", names(q_raw), perl=TRUE)
  stopifnot(all(c("Locus","allele1","allele2") %in% names(q_raw)))
  q <- q_raw[, c("Locus","allele1","allele2"), drop = FALSE]
  q$Locus   <- as.character(q$Locus)
  q$allele1 <- .map_any(q$allele1, any_code)
  q$allele2 <- .map_any(q$allele2, any_code)
  if (!is.null(locus_order) && length(locus_order) > 0L) {
    ord <- match(q$Locus, locus_order)
    q <- q[order(ord), , drop = FALSE]
  }
  rownames(q) <- NULL
  q
}

# build db_index (list: sample_ids, locus_ids, A1, A2) from long DF
.build_db_index_v1 <- function(db_prep, locus_order, any_code = ANY_CODE) {
  req <- c("SampleID","Locus","Allele1","Allele2")
  stopifnot(is.data.frame(db_prep), all(req %in% names(db_prep)))
  SIDs <- unique(db_prep$SampleID)
  LIDs <- as.character(locus_order)
  S <- length(SIDs); L <- length(LIDs)
  A1 <- matrix(any_code, nrow = S, ncol = L, dimnames = list(NULL, LIDs))
  A2 <- matrix(any_code, nrow = S, ncol = L, dimnames = list(NULL, LIDs))
  Lmap <- setNames(seq_along(LIDs), LIDs)
  sid_map <- setNames(seq_along(SIDs), SIDs)
  for (k in seq_len(nrow(db_prep))) {
    sid <- db_prep$SampleID[k]
    loc <- as.character(db_prep$Locus[k])
    i <- sid_map[[sid]]; j <- Lmap[[loc]]
    if (is.na(i) || is.na(j)) next
    A1[i, j] <- .map_any(db_prep$Allele1[k], any_code)
    A2[i, j] <- .map_any(db_prep$Allele2[k], any_code)
  }
  list(sample_ids = SIDs, locus_ids = LIDs, A1 = A1, A2 = A2)
}

load_locus_order <- function() {
  if (file.exists("data/locus_order.rds")) {
    lo <- tryCatch(readRDS("data/locus_order.rds"), error = function(e) NULL)
    if (is.character(lo)) return(lo)
  }
  NULL
}

# Load either index RDS or long DF RDS; also returns a query builder
load_db_index_or_df <- function(path) {
  obj <- readRDS(path)
  lo <- load_locus_order()
  if (is.list(obj) && all(c("sample_ids","locus_ids","A1","A2") %in% names(obj))) {
    make_query_from_index <- function(idx) {
      sidx <- 1L
      LIDs <- idx$locus_ids
      q <- data.frame(
        Locus   = LIDs,
        allele1 = idx$A1[sidx, ],
        allele2 = idx$A2[sidx, ],
        stringsAsFactors = FALSE
      )
      .prepare_query_min(q, LIDs, ANY_CODE)
    }
    return(list(db_index = obj, make_query = function() make_query_from_index(obj), locus_order = obj$locus_ids))
  }
  if (is.data.frame(obj)) {
    n <- names(obj)
    names(obj) <- sub("(?i)^sampleid$","SampleID", n, perl=TRUE)
    names(obj) <- sub("(?i)^locus$","Locus", names(obj), perl=TRUE)
    names(obj) <- sub("(?i)^allele1$","Allele1", names(obj), perl=TRUE)
    names(obj) <- sub("(?i)^allele2$","Allele2", names(obj), perl=TRUE)
    LIDs <- lo; if (is.null(LIDs)) LIDs <- unique(as.character(obj$Locus))
    idx <- .build_db_index_v1(obj, LIDs, ANY_CODE)
    make_query_from_df <- function(df) {
      sid <- df$SampleID[1]
      q <- subset(df, SampleID == sid, c("Locus","Allele1","Allele2"))
      names(q) <- c("Locus","allele1","allele2")
      .prepare_query_min(q, LIDs, ANY_CODE)
    }
    return(list(db_index = idx, make_query = function() make_query_from_df(obj), locus_order = LIDs))
  }
  stop("Unsupported RDS format. Expect index list or long data.frame.")
}

get_result_nrow <- function(res) {
  if (is.data.frame(res)) return(nrow(res))
  if (is.list(res)) {
    if (!is.null(res$match_scores) && is.data.frame(res$match_scores)) return(nrow(res$match_scores))
    if (!is.null(res$detail) && is.data.frame(res$detail)) return(nrow(res$detail))
    if (!is.null(res$scores) && is.atomic(res$scores)) return(length(res$scores))
    return(NA_integer_)
  }
  if (is.atomic(res)) return(length(res))
  NA_integer_
}

append_row <- function(df_row, csv_path) {
  exists <- file.exists(csv_path)
  utils::write.table(
    df_row, file = csv_path, sep = ",",
    row.names = FALSE, col.names = !exists, append = exists, qmethod = "double"
  )
  flush.console()
}

# chunk worker (serial callable; also used inside parallel)
run_one_chunk <- function(cs, k, db_index, q, tag, db_rds_path) {
  t0 <- proc.time()[3L]
  
  # slice
  s_from <- (k-1L)*cs + 1L
  s_to   <- min(length(db_index$sample_ids), k*cs)
  sl <- s_from:s_to
  idx_slice <- list(
    sample_ids = db_index$sample_ids[sl],
    locus_ids  = db_index$locus_ids,
    A1 = db_index$A1[sl, , drop = FALSE],
    A2 = db_index$A2[sl, , drop = FALSE]
  )
  
  # prep timing (slice + query reuse, so prep is mostly slicing)
  t_prep1 <- proc.time()[3L]
  ms_prep <- (t_prep1 - t0) * 1000
  
  ms_io <- 0
  ms_sched <- 0
  
  # ---- MATCHER CALL (adjust signature if needed) ----
  t_work0 <- proc.time()[3L]
  res <- run_match_fast(
    q,
    idx_slice,
    pre_add_for_any_any = pre_add_for_any_any,
    any_code = ANY_CODE,
    include_bits_in_detail = include_bits_in_detail
  )
  t_work1 <- proc.time()[3L]
  gc_before <- gc(reset = TRUE)
  invisible(gc())
  gc_after  <- gc()
  had_gc <- any(gc_after[, "used"] < gc_before[, "used"])  # 簡易判定
  row$gc_triggered <- had_gc
  
  # ---------------------------------------------------
  
  ms_worker <- (t_work1 - t_work0) * 1000
  ms_reduce <- 0
  ms_total  <- (t_work1 - t0) * 1000
  elapsed   <- as.numeric(t_work1 - t0)
  
  samples <- length(sl)
  loci_total <- samples * length(idx_slice$locus_ids)
  thr_samples_sec <- if (elapsed > 0) samples / elapsed else NA_real_
  thr_loci_sec    <- if (elapsed > 0) loci_total / elapsed else NA_real_
  res_rows <- get_result_nrow(res)
  
  row <- data.frame(
    ts = format(Sys.time(), "%Y%m%d%H%M%S"),
    N = samples,
    seed = seed,
    type = bench_type,
    cores = if (enable_parallel) parallel_workers else 1L,
    chunk_size = cs,
    include_io = include_io,
    ms_prep = round(ms_prep, 3),
    ms_io = round(ms_io, 3),
    ms_sched = round(ms_sched, 3),
    ms_worker = round(ms_worker, 3),
    ms_reduce = round(ms_reduce, 3),
    ms_total = round(ms_total, 3),
    rows = res_rows,
    note = "bench_u100",
    # extensions
    run_id = tag,
    dataset_id = basename(db_rds_path),
    chunk_index = k,
    throughput_samples_per_sec = thr_samples_sec,
    throughput_loci_per_sec    = thr_loci_sec,
    workers = if (enable_parallel) parallel_workers else 1L,
    stringsAsFactors = FALSE
  )
  
  logline <- sprintf(
    "  [done] size=%6d idx=%4d  N=%6d  worker=%.2fs total=%.2fs  thr=%.0f samp/s, %.0f loci/s  rows=%d\n",
    cs, k, samples, row$ms_worker/1000, row$ms_total/1000,
    thr_samples_sec, thr_loci_sec, res_rows
  )
  
  list(row = row, logline = logline)
}

########################################
# ========== MAIN ==========
########################################

safe_dir_create(result_dir); safe_dir_create(log_dir)

info <- load_db_index_or_df(db_rds_path)
db_index   <- info$db_index
locus_order <- info$locus_order
make_query <- info$make_query

S <- length(db_index$sample_ids)
L <- length(db_index$locus_ids)
cat("[BENCH] DB index S=", S, " L=", L, "\n", sep="")

# Query externalization (build ONCE)
q <- make_query()

# CSV target (support manual tag to truly resume across sessions)
tag <- if (nzchar(manual_run_tag)) manual_run_tag else ts_tag()
out_csv <- file.path(result_dir, sprintf("bench_run_match_fast_%s.csv", tag))
cat("[BENCH] output CSV: ", out_csv, "\n", sep="")

# resume map for this run (if appending to same CSV, previous rows are respected)
completed <- list()
if (file.exists(out_csv)) {
  done <- tryCatch(read.csv(out_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(done) && nrow(done) > 0) {
    keys <- paste(done$chunk_size, done$chunk_index, sep=":")
    completed <- as.list(stats::setNames(rep(TRUE, length(keys)), keys))
  }
}
is_done <- function(cs, k) {
  key <- paste(cs, k, sep=":")
  isTRUE(completed[[key]])
}

for (cs in chunk_sizes) {
  n_chunks <- max(1L, ceiling(S / cs))
  cat("[BENCH] chunk_size=", cs, " chunks=", n_chunks, "\n", sep="")
  
  if (!enable_parallel) {
    for (k in seq_len(n_chunks)) {
      key <- paste(cs, k, sep=":")
      if (is_done(cs, k)) { cat("  - skip chunk (already in CSV): size=", cs, " idx=", k, "\n", sep=""); next }
      out <- run_one_chunk(cs, k, db_index, q, tag, db_rds_path)
      append_row(out$row, out_csv); cat(out$logline); completed[[key]] <- TRUE; gc()
    }
  } else {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(multisession, workers = parallel_workers)
    
    todo <- which(!vapply(seq_len(n_chunks), function(k) is_done(cs, k), logical(1)))
    if (length(todo) > 0) {
      res_list <- future_lapply(todo, function(k)
        run_one_chunk(cs, k, db_index, q, tag, db_rds_path))
      for (i in seq_along(todo)) {
        k <- todo[i]; key <- paste(cs, k, sep=":")
        append_row(res_list[[i]]$row, out_csv); cat(res_list[[i]]$logline); completed[[key]] <- TRUE
      }
      gc()
    } else {
      cat("  - all chunks already in CSV\n")
    }
  }
}

cat("[BENCH DONE] summary => ", out_csv, "\n", sep="")

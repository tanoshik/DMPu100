# bench/bench_run_match_fast_chunks_u100.R
# Chunked matching bench with explicit run plan (no full cartesian).
# Root working dir = project root (DMPu100).
# No multibyte characters in code/comments.

# ==================== CONFIG (EDIT HERE) ====================
bench_type   <- "vdb_gf"   # gf_like -> vdb_gf (virtual DB built from real frequencies)
db_rds_path  <- "data/virtual_db_u100_S1000_seed123.rds"  # 1M RDS (relative path)

# Output root and tag
out_dir      <- "output/bench/1M_vdb_gf"
run_tag      <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ---- IMPORTANT: explicit run plan (list only what you want to run) ----
# Each row is one run: a (chunk_size, workers) pair. Only these pairs will execute.
# Example below: narrower than a 4x4 grid.
run_plan <- data.frame(
  chunk_size = c(100000L, 200000L, 250000L, 500000L),
  workers    = c(8L,      4L,      4L,      2L),
  stringsAsFactors = FALSE
)

# You can override run_plan via env or CLI:
#   Env: RUN_PLAN="100000x8,200000x4,500000x2"
#   CLI: Rscript bench/bench_run_match_fast_chunks_u100.R --plan=100000x8,200000x4
#
# Resume/log settings
resume_enabled     <- TRUE
log_flush_every    <- 1L       # per chunk
state_save_every   <- 1L       # per chunk
gc_every_n_chunks  <- 2L       # do not over-GC

# Output weight: Summary-only for 1M
write_detail_csv   <- FALSE

# One-time strict schema check for vdb_gf
strict_schema_check_once <- TRUE
# ==================== END CONFIG ============================

# ==================== BOOTSTRAP =============================
`%||%` <- function(x, y) if (is.null(x)) y else x

safe_show <- function(path) {
  p <- tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE),
                error = function(e) path)
  if (is.na(p) || p == "") path else p
}

ensure_dir <- function(d) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Allow overriding run_plan via ENV or CLI
parse_plan_string <- function(s) {
  # format: "100000x8,200000x4,500000x2"
  if (is.null(s) || is.na(s) || !nzchar(s)) return(NULL)
  parts <- strsplit(s, ",", fixed = TRUE)[[1]]
  cs <- integer(); wk <- integer()
  for (p in parts) {
    kv <- strsplit(trimws(p), "x", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Bad RUN_PLAN item: ", p)
    cs <- c(cs, as.integer(kv[1]))
    wk <- c(wk, as.integer(kv[2]))
  }
  data.frame(chunk_size = cs, workers = wk, stringsAsFactors = FALSE)
}

# CLI override
parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(plan = NULL)
  if (length(args) == 0) return(out)
  for (a in args) {
    if (startsWith(a, "--plan=")) {
      out$plan <- sub("^--plan=", "", a)
    } else if (startsWith(a, "--db=")) {
      db_rds_path <<- sub("^--db=", "", a)
    } else if (startsWith(a, "--out=")) {
      out_dir <<- sub("^--out=", "", a)
    } else if (startsWith(a, "--tag=")) {
      run_tag <<- sub("^--tag=", "", a)
    }
  }
  out
}

# Source project scripts if available
maybe_source <- function(p) if (file.exists(p)) try(suppressWarnings(source(p)), silent = TRUE)

maybe_source("scripts/utils_profile.R")
maybe_source("scripts/utils_profile_modular.R")
maybe_source("scripts/scoring_fast.R")
maybe_source("scripts/matcher.R")

match_fun_name <- if (exists("run_match_fast")) "run_match_fast" else if (exists("matcher_fast")) "matcher_fast" else NA
if (is.na(match_fun_name)) {
  stop("No run_match_fast() or matcher_fast() found. Please source scripts/scoring_fast.R or scripts/matcher.R.")
}
match_fun <- get(match_fun_name)

# ==================== IO PREP =================================
run_root <- file.path(out_dir, run_tag)
ensure_dir(run_root)

log_path  <- file.path(run_root, "bench.log")
state_dir <- file.path(run_root, "state")
ensure_dir(state_dir)

summary_path_for <- function(cs, wk) file.path(run_root, sprintf("summary_chunk%06d_w%02d.csv", as.integer(cs), as.integer(wk)))
detail_path_for  <- function(cs, wk) file.path(run_root, sprintf("detail_chunk%06d_w%02d.csv",  as.integer(cs), as.integer(wk)))
state_path_for   <- function(cs, wk) file.path(state_dir, sprintf("state_chunk%06d_w%02d.rds",  as.integer(cs), as.integer(wk)))

log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_path, append = TRUE)
}

# ==================== DATA LOADER ============================
# Expect list(sample_ids, locus_ids, A1, A2) for vdb_gf
load_vdb_gf <- function(rds_path) {
  obj <- readRDS(rds_path)
  if (!is.list(obj) || !all(c("sample_ids","locus_ids","A1","A2") %in% names(obj))) {
    stop("RDS must be a list with names: sample_ids, locus_ids, A1, A2")
  }
  sample_ids <- obj$sample_ids
  locus_ids  <- obj$locus_ids
  A1 <- obj$A1; A2 <- obj$A2
  
  to_mat <- function(x) {
    if (is.matrix(x)) return(x)
    if (is.data.frame(x)) return(as.matrix(x))
    if (is.atomic(x) && length(x) == length(sample_ids) * length(locus_ids)) {
      return(matrix(x, nrow = length(sample_ids), ncol = length(locus_ids), byrow = TRUE))
    }
    stop("A1/A2 must be matrix/data.frame [N x L] or atomic vector of length N*L.")
  }
  A1m <- to_mat(A1)
  A2m <- to_mat(A2)
  
  if (nrow(A1m) != length(sample_ids) || ncol(A1m) != length(locus_ids)) stop("A1 shape mismatch")
  if (nrow(A2m) != length(sample_ids) || ncol(A2m) != length(locus_ids)) stop("A2 shape mismatch")
  
  list(sample_ids = sample_ids, locus_ids = locus_ids, A1 = A1m, A2 = A2m,
       N = length(sample_ids), L = length(locus_ids))
}

# Make long DF for sample row range [i_start:i_end]
db_chunk_df <- function(db, i_start, i_end) {
  i_start <- max(1L, i_start); i_end <- min(db$N, i_end)
  idx <- i_start:i_end
  sid <- db$sample_ids[idx]
  L <- db$L
  SampleID <- rep(sid, each = L)
  Locus    <- rep(db$locus_ids, times = length(idx))
  A1v <- as.vector(t(db$A1[idx, , drop = FALSE]))
  A2v <- as.vector(t(db$A2[idx, , drop = FALSE]))
  df <- data.frame(SampleID = SampleID, Locus = Locus,
                   allele1 = as.character(A1v),
                   allele2 = as.character(A2v),
                   stringsAsFactors = FALSE)
  df
}

query_from_sample <- function(db, sample_index = 1L) {
  if (sample_index < 1L || sample_index > db$N) stop("query sample_index out of range")
  df <- db_chunk_df(db, sample_index, sample_index)
  df[, c("Locus","allele1","allele2")]
}

# ==================== RESULT NORMALIZER ======================
normalize_match_result <- function(res) {
  if (is.list(res) && "summary" %in% names(res) && is.data.frame(res$summary)) {
    summary <- res$summary
  } else if (is.data.frame(res) && all(c("SampleID","Score") %in% names(res))) {
    summary <- res
  } else if (is.list(res) && "scores" %in% names(res) && is.data.frame(res$scores)) {
    summary <- res$scores
  } else {
    stop("Cannot find summary scores in result.")
  }
  detail <- NULL
  cand_names <- c("detail","log","match_log","rows")
  for (nm in cand_names) {
    if (is.list(res) && nm %in% names(res) && is.data.frame(res[[nm]])) { detail <- res[[nm]]; break }
  }
  list(summary = summary, detail = detail)
}

# ==================== CORE RUNNERS ===========================
run_for_pair <- function(db, pair_chunk, pair_workers, q_df) {
  cs <- as.integer(pair_chunk)
  wk <- as.integer(pair_workers)
  state_path   <- state_path_for(cs, wk)
  summary_path <- summary_path_for(cs, wk)
  detail_path  <- detail_path_for(cs, wk)
  
  state <- if (file.exists(state_path) && resume_enabled) readRDS(state_path) else list(next_start = 1L, chunk_count = 0L)
  next_start <- as.integer(state$next_start %||% 1L)
  chunk_count <- as.integer(state$chunk_count %||% 0L)
  
  total <- db$N
  wrote_summary_header <- file.exists(summary_path) && (file.info(summary_path)$size > 0)
  wrote_detail_header  <- file.exists(detail_path)  && (file.info(detail_path)$size  > 0)
  
  log_line(sprintf("[pair cs=%d, w=%d] resume from %d/%d", cs, wk, next_start, total))
  
  repeat {
    if (next_start > total) {
      log_line(sprintf("[pair cs=%d, w=%d] completed all samples.", cs, wk))
      break
    }
    end_i <- min(total, next_start + cs - 1L)
    db_df <- db_chunk_df(db, next_start, end_i)
    
    # call matcher; pass workers if supported
    t0 <- proc.time()[3]
    res <- try({
      fml <- try(formals(match_fun), silent = TRUE)
      if (!inherits(fml, "try-error") && "workers" %in% names(fml)) {
        match_fun(query_df = q_df, db_df = db_df, workers = wk)
      } else {
        match_fun(query_df = q_df, db_df = db_df)
      }
    }, silent = TRUE)
    t1 <- proc.time()[3]
    
    if (inherits(res, "try-error")) stop("Matcher call failed: ", as.character(res))
    
    norm <- normalize_match_result(res)
    s <- norm$summary
    s$ChunkSize  <- cs
    s$Workers    <- wk
    s$IdxStart   <- next_start
    s$IdxEnd     <- end_i
    s$ElapsedSec <- round(as.numeric(t1 - t0), 3)
    
    write.table(s, file = summary_path, sep = ",", row.names = FALSE, col.names = !wrote_summary_header, append = wrote_summary_header)
    wrote_summary_header <- TRUE
    
    if (isTRUE(write_detail_csv) && !is.null(norm$detail)) {
      d <- norm$detail
      d$ChunkSize <- cs
      d$Workers   <- wk
      d$IdxStart  <- next_start
      d$IdxEnd    <- end_i
      write.table(d, file = detail_path, sep = ",", row.names = FALSE, col.names = !wrote_detail_header, append = wrote_detail_header)
      wrote_detail_header <- TRUE
    }
    
    if ((end_i - next_start + 1L) >= log_flush_every) {
      log_line(sprintf("[pair cs=%d, w=%d] processed [%d..%d] in %.3f sec (rows=%d)",
                       cs, wk, next_start, end_i, (t1 - t0), nrow(db_df)))
    }
    
    # update state for resume
    next_start <- end_i + 1L
    chunk_count <- chunk_count + 1L
    if (chunk_count %% state_save_every == 0L) {
      saveRDS(list(next_start = next_start, chunk_count = chunk_count), state_path)
    }
    if (chunk_count %% gc_every_n_chunks == 0L) {
      gc(FALSE)
    }
    flush.console()
  }
}

# ==================== ENTRYPOINT =============================
main <- function() {
  # CLI overrides
  cli <- parse_cli()
  env_plan <- parse_plan_string(Sys.getenv("RUN_PLAN", ""))
  if (!is.null(cli$plan)) env_plan <- parse_plan_string(cli$plan)
  if (!is.null(env_plan)) {
    run_plan <<- env_plan
  }
  
  # Basic logging
  log_line("=== bench start ===")
  log_line("bench_type: ", bench_type)
  log_line("DB: ", db_rds_path)
  log_line("Output: ", safe_show(run_root))
  log_line("Matcher: ", match_fun_name)
  
  # Load DB once
  if (bench_type == "vdb_gf") {
    db <- load_vdb_gf(db_rds_path)
  } else {
    stop("Unsupported bench_type: ", bench_type)
  }
  log_line(sprintf("Loaded DB: N=%d samples, L=%d loci", db$N, db$L))
  log_line("First SampleID: ", db$sample_ids[1])
  
  # Optional one-time strict schema checks
  if (strict_schema_check_once) {
    if (!is.character(db$sample_ids)) stop("sample_ids must be character")
    if (!is.character(db$locus_ids))  stop("locus_ids must be character")
  }
  
  # Build query from the first sample
  q_df <- query_from_sample(db, 1L)
  
  # Validate run_plan
  if (!all(c("chunk_size","workers") %in% names(run_plan))) {
    stop("run_plan must have columns: chunk_size, workers")
  }
  if (nrow(run_plan) == 0L) {
    stop("run_plan is empty. Provide pairs or set RUN_PLAN/--plan.")
  }
  
  # Execute pairs sequentially (explicit plan; no cartesian)
  for (i in seq_len(nrow(run_plan))) {
    cs <- as.integer(run_plan$chunk_size[i])
    wk <- as.integer(run_plan$workers[i])
    log_line(sprintf("--- run (%d/%d): chunk=%d, workers=%d ---", i, nrow(run_plan), cs, wk))
    run_for_pair(db, cs, wk, q_df)
  }
  
  log_line("=== bench done ===")
}

# Prepare dirs
ensure_dir(out_dir)
ensure_dir(run_root)

# Auto-run when sourced or called via Rscript
if (sys.nframe() == 0L) main()

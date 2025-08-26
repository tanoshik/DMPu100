# No multibyte characters.

########################################
# ========== CONFIG ==========
########################################
db_rds_path   <- "data/virtual_db_u100_S1000_seed123.rds"
result_dir    <- "test/bench"
log_dir       <- "logs"

ANY_CODE      <- 9999L
pre_add_for_any_any <- 2L
include_bits_in_detail <- FALSE

bench_type    <- "vdb_gf"
seed          <- 1L

# default plan (used when RUN_PLAN/--plan not provided)
chunk_sizes   <- c(100000L, 200000L, 250000L, 500000L)
parallel_workers <- max(1L, parallel::detectCores() - 1L)

manual_run_tag <- ""       # empty => timestamp
chunk_timeout_sec <- 1800L # can be overridden by --timeout=

########################################
# ========== DEPS / SOURCES ==========
########################################
suppressWarnings(suppressMessages({
  library(future)
  library(future.apply)
  library(ps)   # <- added
}))
options(future.globals.maxSize = 8 * 1024^3)
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

source("scripts/scoring_fast.R",  local = TRUE)
source("scripts/matcher_fast.R",  local = TRUE)

########################################
# ========== UTILS ==========
########################################
`%||%` <- function(x,y) if (is.null(x)) y else x
.ensure_dir <- function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)
ts_tag <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

log_file <- function() { .ensure_dir(log_dir); file.path(log_dir,"bench.log") }
log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(..., collapse=""))
  cat(msg,"\n"); cat(msg,"\n", file=log_file(), append=TRUE)
}

append_row <- function(row, csv_path, retries=25L, sleep_sec=0.2) {
  .ensure_dir(dirname(csv_path))
  has_file <- file.exists(csv_path) && (file.info(csv_path)$size %||% 0) > 0
  tmp <- tempfile("bench_row_", tmpdir=dirname(csv_path), fileext=".csv")
  on.exit(try(unlink(tmp), silent=TRUE), add=TRUE)
  write.table(row, file=tmp, sep=",", row.names=FALSE, col.names=!has_file, append=FALSE)
  for (i in seq_len(retries)) {
    ok <- try({
      if (!has_file) file.rename(tmp, csv_path) else { file.append(csv_path, tmp); TRUE }
    }, silent=TRUE)
    if (!inherits(ok,"try-error") && isTRUE(ok)) { if (file.exists(tmp)) try(unlink(tmp), silent=TRUE); return(invisible(TRUE)) }
    Sys.sleep(sleep_sec)
  }
  stop("Failed to append: ", csv_path)
}

.map_any <- function(x, any_code=ANY_CODE) {
  if (is.null(x)) return(integer(0))
  if (is.factor(x)) x <- as.character(x)
  out <- suppressWarnings(as.integer(x))
  is_any <- !is.na(x) & tolower(as.character(x))=="any"
  out[is_any] <- any_code
  out
}

.prepare_query_min <- function(q_raw, locus_order=NULL, any_code=ANY_CODE) {
  n <- names(q_raw)
  names(q_raw) <- sub("(?i)^locus$","Locus", n, perl=TRUE)
  names(q_raw) <- sub("(?i)^allele1$","allele1", names(q_raw), perl=TRUE)
  names(q_raw) <- sub("(?i)^allele2$","allele2", names(q_raw), perl=TRUE)
  stopifnot(all(c("Locus","allele1","allele2") %in% names(q_raw)))
  q <- q_raw[,c("Locus","allele1","allele2"),drop=FALSE]
  q$Locus <- as.character(q$Locus)
  q$allele1 <- .map_any(q$allele1, any_code)
  q$allele2 <- .map_any(q$allele2, any_code)
  if (!is.null(locus_order) && length(locus_order)>0L) {
    ord <- match(q$Locus, locus_order); q <- q[order(ord), , drop=FALSE]
  }
  rownames(q) <- NULL; q
}

.build_db_index_v1 <- function(db_prep, locus_order, any_code=ANY_CODE) {
  req <- c("SampleID","Locus","Allele1","Allele2")
  stopifnot(is.data.frame(db_prep), all(req %in% names(db_prep)))
  SIDs <- unique(db_prep$SampleID); LIDs <- as.character(locus_order)
  S <- length(SIDs); L <- length(LIDs)
  A1 <- matrix(any_code, nrow=S, ncol=L, dimnames=list(NULL,LIDs))
  A2 <- matrix(any_code, nrow=S, ncol=L, dimnames=list(NULL,LIDs))
  Lmap <- setNames(seq_along(LIDs), LIDs); sid_map <- setNames(seq_along(SIDs), SIDs)
  for (k in seq_len(nrow(db_prep))) {
    sid <- db_prep$SampleID[k]; loc <- as.character(db_prep$Locus[k])
    i <- sid_map[[sid]]; j <- Lmap[[loc]]; if (is.na(i) || is.na(j)) next
    A1[i,j] <- .map_any(db_prep$Allele1[k], any_code)
    A2[i,j] <- .map_any(db_prep$Allele2[k], any_code)
  }
  list(sample_ids=SIDs, locus_ids=LIDs, A1=A1, A2=A2)
}

load_locus_order <- function() {
  if (file.exists("data/locus_order.rds")) {
    lo <- tryCatch(readRDS("data/locus_order.rds"), error=function(e) NULL)
    if (is.character(lo)) return(lo)
  }
  NULL
}

load_rds_as_index <- function(path) {
  obj <- readRDS(path); lo <- load_locus_order()
  if (is.list(obj) && all(c("sample_ids","locus_ids","A1","A2") %in% names(obj))) {
    return(list(db_index=obj, make_query=function() {
      data.frame(Locus=obj$locus_ids, allele1=obj$A1[1,,drop=TRUE], allele2=obj$A2[1,,drop=TRUE], stringsAsFactors=FALSE)
    }))
  }
  if (is.data.frame(obj)) {
    n <- names(obj)
    names(obj) <- sub("(?i)^sampleid$","SampleID", n, perl=TRUE)
    names(obj) <- sub("(?i)^locus$","Locus", names(obj), perl=TRUE)
    names(obj) <- sub("(?i)^allele1$","Allele1", names(obj), perl=TRUE)
    names(obj) <- sub("(?i)^allele2$","Allele2", names(obj), perl=TRUE)
    LIDs <- lo; if (is.null(LIDs)) LIDs <- unique(as.character(obj$Locus))
    idx <- .build_db_index_v1(obj, LIDs, ANY_CODE)
    return(list(db_index=idx, make_query=function() {
      sid <- obj$SampleID[1]
      q <- subset(obj, SampleID==sid, c("Locus","Allele1","Allele2"))
      names(q) <- c("Locus","allele1","allele2"); q
    }))
  }
  stop("Unsupported RDS")
}

make_slice <- function(db_index, cs, k) {
  s_from <- (k-1L)*cs + 1L; s_to <- min(length(db_index$sample_ids), k*cs)
  sl_id <- s_from:s_to
  list(
    sample_ids = db_index$sample_ids[sl_id],
    locus_ids  = db_index$locus_ids,
    A1 = db_index$A1[sl_id, , drop=FALSE],
    A2 = db_index$A2[sl_id, , drop=FALSE]
  )
}

.make_common_cols <- function() {
  data.frame(
    ts = character(), N = integer(), seed = integer(), type = character(),
    cores = integer(), chunk_size = integer(), include_io = logical(),
    ms_prep = numeric(), ms_io = numeric(), ms_sched = numeric(),
    ms_worker = numeric(), ms_reduce = numeric(), ms_total = numeric(),
    rows = integer(), note = character(), run_id = character(),
    dataset_id = character(), chunk_index = integer(),
    throughput_samples_per_sec = numeric(), throughput_loci_per_sec = numeric(),
    workers = integer(),
    # existing metrics
    mem_rss_mb_before = numeric(), mem_rss_mb_after = numeric(), mem_peak_mb = numeric(),
    cpu_user_sec = numeric(), cpu_system_sec = numeric(), cpu_total_sec = numeric(),
    cpu_util_pct_of_one_core = numeric(), cpu_util_pct_per_core_est = numeric(),
    io_read_bytes = numeric(), io_write_bytes = numeric(),
    io_read_MBps = numeric(), io_write_MBps = numeric(),
    pid = integer(), host = character(),
    # NEW:
    elapsed_sec = numeric(),           # 実測経過（worker内での t_work1 - t0）
    iowait_est_pct = numeric(),        # ≈ max(0, (elapsed - cpu_total)/elapsed*100)
    stringsAsFactors = FALSE
  )
}

make_start_row <- function(cs, k, N, tag, dataset_id, wk) {
  # Create a single-row DF with all columns (including new metrics) initialized
  df <- .make_common_cols()
  df[1, ] <- NA
  df$ts[1] <- format(Sys.time(), "%Y%m%d%H%M%S")
  df$N[1] <- N; df$seed[1] <- seed; df$type[1] <- paste0(bench_type, "_START")
  df$cores[1] <- wk; df$chunk_size[1] <- cs; df$include_io[1] <- FALSE
  df$note[1] <- "START"; df$run_id[1] <- tag; df$dataset_id[1] <- dataset_id
  df$chunk_index[1] <- k; df$workers[1] <- wk
  df
}

get_rows <- function(res) {
  if (is.data.frame(res)) return(nrow(res))
  if (is.list(res) && !is.null(res$summary)) return(nrow(res$summary))
  NA_integer_
}

########################################
# ========== CLI / PLAN ==========
########################################
parse_plan_string <- function(s) {
  if (is.null(s) || !nzchar(s)) return(NULL)
  parts <- strsplit(s, ",", fixed=TRUE)[[1]]
  cs <- integer(); wk <- integer()
  for (p in parts) {
    kv <- strsplit(trimws(p), "x", fixed=TRUE)[[1]]
    if (length(kv)!=2) stop("Bad --plan item: ", p)
    cs <- c(cs, as.integer(kv[1])); wk <- c(wk, as.integer(kv[2]))
  }
  data.frame(chunk_size=cs, workers=wk, stringsAsFactors=FALSE)
}

apply_overrides <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  plan_cli <- NULL
  for (a in args) {
    if (startsWith(a,"--plan=")) plan_cli <- sub("^--plan=","",a)
    else if (startsWith(a,"--db="))  db_rds_path <<- sub("^--db=","",a)
    else if (startsWith(a,"--tag=")) manual_run_tag <<- sub("^--tag=","",a)
    else if (startsWith(a,"--timeout=")) chunk_timeout_sec <<- as.integer(sub("^--timeout=","",a))
    # NOTE: --out= は無視して test/bench に固定
  }
  plan_env <- Sys.getenv("RUN_PLAN","")
  plan <- parse_plan_string(plan_cli %||% plan_env)
  if (!is.null(plan)) return(plan)
  data.frame(chunk_size=chunk_sizes, workers=parallel_workers, stringsAsFactors=FALSE)
}

########################################
# ========== PROCESS METRICS (FIX) ==========
########################################
.get_proc_metrics <- function() {
  h <- ps::ps_handle(Sys.getpid())
  
  # ps_* は OS により「名前付きベクトル（atomic）」として返ることがある
  mem <- ps::ps_memory_info(h)
  cpu <- ps::ps_cpu_times(h)
  io  <- tryCatch(ps::ps_io_counters(h), error=function(e) NULL)
  
  get_num <- function(x, key) {
    if (is.null(x)) return(NA_real_)
    v <- tryCatch(x[[key]], error=function(e) NA)
    if (is.null(v)) return(NA_real_)
    suppressWarnings(as.numeric(v))
  }
  
  list(
    rss = get_num(mem, "rss"),                 # bytes
    user = get_num(cpu, "user"),               # sec (cumulative)
    system = get_num(cpu, "system"),           # sec (cumulative)
    read_bytes  = if (!is.null(io)) get_num(io, "read_bytes")  else NA_real_,
    write_bytes = if (!is.null(io)) get_num(io, "write_bytes") else NA_real_,
    pid  = ps::ps_pid(h),
    host = (Sys.info()[["nodename"]] %||% Sys.getenv("COMPUTERNAME", "unknown-host"))
  )
}

########################################
# ========== WORKER ==========
########################################
run_slice_worker <- function(slice, q, tag, dataset_id, cs, wk,
                             timeout_sec, project_root, seed_param) {
  setwd(project_root)
  assign("ANY_CODE", ANY_CODE, envir=globalenv())
  source("scripts/scoring_fast.R", local=FALSE)
  source("scripts/matcher_fast.R", local=FALSE)
  
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)
  setTimeLimit(cpu = Inf, elapsed = timeout_sec, transient = TRUE)
  
  k <- if (!is.null(attr(slice, "chunk_index"))) attr(slice, "chunk_index") else NA_integer_
  
  # capture process metrics before
  gc()
  m0 <- .get_proc_metrics()
  t0_wall <- proc.time()[3L]
  
  ok <- TRUE; err <- NA_character_; res_row <- NULL
  
  tryCatch({
    t_work0 <- proc.time()[3L]
    res <- run_match_fast(
      q,
      slice,
      pre_add_for_any_any = pre_add_for_any_any,
      any_code = ANY_CODE,
      include_bits_in_detail = include_bits_in_detail
    )
    t_work1 <- proc.time()[3L]
    
    # capture after metrics
    gc()
    m1 <- .get_proc_metrics()
    
    ms_total <- (t_work1 - t_work0) * 1000
    
    N <- length(slice$sample_ids)
    L <- length(slice$locus_ids)
    elapsed <- (t_work1 - t0_wall)
    
    # throughput
    thr_s <- if (elapsed>0) N/elapsed else NA_real_
    thr_L <- if (elapsed>0) (N*L)/elapsed else NA_real_
    rows <- get_rows(res)
    
    # deltas
    cpu_user   <- if (is.finite(m1$user) && is.finite(m0$user))   m1$user   - m0$user   else NA_real_
    cpu_system <- if (is.finite(m1$system) && is.finite(m0$system)) m1$system - m0$system else NA_real_
    cpu_total  <- if (is.finite(cpu_user) && is.finite(cpu_system)) cpu_user + cpu_system else NA_real_
    
    rbytes <- if (is.finite(m1$read_bytes) && is.finite(m0$read_bytes)) m1$read_bytes - m0$read_bytes else NA_real_
    wbytes <- if (is.finite(m1$write_bytes) && is.finite(m0$write_bytes)) m1$write_bytes - m0$write_bytes else NA_real_
    
    rss_before_mb <- if (is.finite(m0$rss)) m0$rss / (1024^2) else NA_real_
    rss_after_mb  <- if (is.finite(m1$rss)) m1$rss / (1024^2) else NA_real_
    mem_peak_mb   <- max(rss_before_mb, rss_after_mb, na.rm = TRUE)
    
    cpu_util_pct_of_one_core <- if (elapsed > 0 && is.finite(cpu_total)) (cpu_total / elapsed) * 100 else NA_real_
    cpu_util_pct_per_core_est <- if (elapsed > 0 && is.finite(cpu_total) && wk > 0) (cpu_total / elapsed / wk) * 100 else NA_real_
    
    io_read_MBps  <- if (elapsed > 0 && is.finite(rbytes)) (rbytes / (1024^2)) / elapsed else NA_real_
    io_write_MBps <- if (elapsed > 0 && is.finite(wbytes)) (wbytes / (1024^2)) / elapsed else NA_real_
    iowait_est_pct <- if (elapsed > 0 && is.finite(cpu_total))
      max(0, (elapsed - cpu_total) / elapsed * 100) else NA_real_
    
    res_row <- data.frame(
      ts = format(Sys.time(), "%Y%m%d%H%M%S"),
      N = N, seed = seed_param, type = bench_type, cores = wk,
      chunk_size = cs, include_io = FALSE,
      ms_prep = NA_real_, ms_io = NA_real_, ms_sched = NA_real_,
      ms_worker = NA_real_, ms_reduce = NA_real_, ms_total = round(ms_total,3),
      rows = rows, note = "OK",
      run_id = tag, dataset_id = dataset_id, chunk_index = k,
      throughput_samples_per_sec = thr_s, throughput_loci_per_sec = thr_L,
      workers = wk,
      mem_rss_mb_before = rss_before_mb, mem_rss_mb_after = rss_after_mb, mem_peak_mb = mem_peak_mb,
      cpu_user_sec = cpu_user, cpu_system_sec = cpu_system, cpu_total_sec = cpu_total,
      cpu_util_pct_of_one_core = cpu_util_pct_of_one_core,
      cpu_util_pct_per_core_est = cpu_util_pct_per_core_est,
      io_read_bytes = rbytes, io_write_bytes = wbytes,
      io_read_MBps = io_read_MBps, io_write_MBps = io_write_MBps,
      pid = as.integer(m1$pid), host = as.character(m1$host),
      elapsed_sec = elapsed,                 # NEW
      iowait_est_pct = iowait_est_pct,       # NEW
      stringsAsFactors = FALSE
    )
  }, error=function(e){ ok <<- FALSE; err <<- conditionMessage(e) })
  
  if (!ok || is.null(res_row)) {
    N <- length(slice$sample_ids)
    res_row <- .make_common_cols()
    res_row[1, ] <- NA
    res_row$ts[1] <- format(Sys.time(), "%Y%m%d%H%M%S")
    res_row$N[1] <- N; res_row$seed[1] <- seed_param
    res_row$type[1] <- paste0(bench_type,"_ERR")
    res_row$cores[1] <- wk; res_row$chunk_size[1] <- cs
    res_row$include_io[1] <- FALSE; res_row$rows[1] <- NA_integer_
    res_row$note[1] <- paste0("ERROR: ", err %||% "unknown")
    res_row$run_id[1] <- tag; res_row$dataset_id[1] <- dataset_id; res_row$chunk_index[1] <- k
    res_row$workers[1] <- wk
  }
  res_row
}

########################################
# ========== MAIN ==========
########################################
main <- function() {
  .ensure_dir(result_dir); .ensure_dir(log_dir)
  plan <- apply_overrides()
  
  project_root <- normalizePath(".", winslash="/", mustWork=FALSE)
  dataset_id <- basename(db_rds_path)
  
  info <- load_rds_as_index(db_rds_path)
  db_index <- info$db_index
  q <- .prepare_query_min(info$make_query(), db_index$locus_ids, ANY_CODE)
  
  S <- length(db_index$sample_ids); L <- length(db_index$locus_ids)
  cat("[BENCH] DB index S=",S," L=",L,"\n", sep="")
  tag <- if (nzchar(manual_run_tag)) manual_run_tag else ts_tag()
  out_csv <- file.path(result_dir, sprintf("bench_run_match_fast_%s.csv", tag))
  cat("[BENCH] output CSV: ", out_csv, "\n", sep="")
  
  # resume table: consider OK rows as completed
  completed <- list()
  if (file.exists(out_csv)) {
    done <- tryCatch(read.csv(out_csv, stringsAsFactors=FALSE), error=function(e) NULL)
    if (!is.null(done) && nrow(done)>0) {
      ok_rows <- subset(done, type==bench_type & (is.na(note) | note=="OK"))
      if (nrow(ok_rows)>0) {
        ok_rows$workers <- ok_rows$workers %||% NA_integer_
        keys <- paste(ok_rows$chunk_size, ok_rows$workers, ok_rows$chunk_index, sep=":")
        completed <- as.list(stats::setNames(rep(TRUE,length(keys)), keys))
      }
    }
  }
  is_done <- function(cs,wk,k) isTRUE(completed[[paste(cs,wk,k,sep=":")]])
  
  for (i in seq_len(nrow(plan))) {
    cs <- as.integer(plan$chunk_size[i]); wk <- as.integer(plan$workers[i])
    n_chunks <- max(1L, ceiling(S/cs))
    log_line(sprintf("[PAIR] chunk=%d workers=%d chunks=%d", cs, wk, n_chunks))
    
    todo <- which(!vapply(seq_len(n_chunks), function(k) is_done(cs,wk,k), logical(1)))
    if (length(todo)==0L) { cat("  - all chunks already in CSV for this pair\n"); next }
    
    # START rows (visible in CSV)
    for (k in todo) {
      N <- length(((k-1L)*cs+1L):min(S, k*cs))
      append_row(make_start_row(cs,k,N,tag,dataset_id,wk), out_csv)
      cat(sprintf("[CHUNK][START] cs=%d k=%d N=%d\n", cs, k, N))
    }
    
    if (wk <= 1L) {
      for (k in todo) {
        sl <- make_slice(db_index, cs, k); attr(sl,"chunk_index") <- k
        row <- run_slice_worker(sl, q, tag, dataset_id, cs, wk,
                                chunk_timeout_sec, project_root, seed_param = seed)
        append_row(row, out_csv)
        if (identical(row$note, "OK")) completed[[paste(cs,wk,k,sep=":")]] <- TRUE
        cat(if (identical(row$note,"OK"))
          sprintf("[CHUNK] cs=%d k=%d N=%d rows=%s ms_total=%s thr_s=%s thr_L=%s\n",
                  cs,k,row$N, ifelse(is.na(row$rows),"NA",as.character(row$rows)),
                  as.character(row$ms_total),
                  as.character(row$throughput_samples_per_sec),
                  as.character(row$throughput_loci_per_sec))
          else sprintf("[CHUNK][FAIL] cs=%d k=%d msg=%s\n", cs,k,row$note))
        gc()
      }
    } else {
      old <- future::plan(); on.exit(future::plan(old), add=TRUE)
      future::plan(multisession, workers=wk)
      
      remaining_k <- todo
      in_flight   <- list()
      last_ping   <- Sys.time()
      
      submit_next <- function() {
        while (length(in_flight) < wk && length(remaining_k) > 0) {
          k <- remaining_k[1]; remaining_k <<- remaining_k[-1]
          sl <- make_slice(db_index, cs, k); attr(sl,"chunk_index") <- k
          fut <- future::future(
            run_slice_worker(sl, q, tag, dataset_id, cs, wk,
                             chunk_timeout_sec, project_root, seed_param = seed),
            globals = TRUE
          )
          in_flight[[as.character(k)]] <<- fut
        }
      }
      
      submit_next()
      log_line(sprintf("[SUBMIT] batch=%d (cs=%d, workers=%d)", length(in_flight), cs, wk))
      
      while (length(in_flight) > 0) {
        ready_k <- names(in_flight)[ vapply(in_flight, future::resolved, logical(1)) ]
        if (length(ready_k) == 0) {
          if (as.numeric(difftime(Sys.time(), last_ping, units="secs")) >= 30) {
            log_line(sprintf("[WAIT] cs=%d workers=%d inflight=%d remain=%d",
                             cs, wk, length(in_flight), length(remaining_k)))
            last_ping <- Sys.time()
          }
          Sys.sleep(1)
          next
        }
        
        for (kk in ready_k) {
          k <- as.integer(kk)
          fut <- in_flight[[kk]]
          row <- tryCatch(
            future::value(fut),
            error = function(e) {
              df <- .make_common_cols(); df[1, ] <- NA
              df$ts[1] <- format(Sys.time(),"%Y%m%d%H%M%S"); df$seed[1] <- seed
              df$type[1] <- paste0(bench_type,"_ERR"); df$cores[1] <- wk
              df$chunk_size[1] <- cs; df$include_io[1] <- FALSE
              df$note[1] <- paste0("ERROR: ", conditionMessage(e))
              df$run_id[1] <- tag; df$dataset_id[1] <- dataset_id
              df$chunk_index[1] <- k; df$workers[1] <- wk
              df
            }
          )
          append_row(row, out_csv)
          if (identical(row$note, "OK"))
            completed[[paste(cs,wk,k,sep=":")]] <- TRUE
          
          cat(if (identical(row$note,"OK"))
            sprintf("[CHUNK] cs=%d k=%d N=%d rows=%s ms_total=%s thr_s=%s thr_L=%s\n",
                    cs,k,row$N, ifelse(is.na(row$rows),"NA",as.character(row$rows)),
                    as.character(row$ms_total),
                    as.character(row$throughput_samples_per_sec),
                    as.character(row$throughput_loci_per_sec))
            else sprintf("[CHUNK][FAIL] cs=%d k=%d msg=%s\n", cs,k,row$note))
          
          in_flight[[kk]] <- NULL
          submit_next()
        }
        gc()
      }
      log_line(sprintf("[BATCH DONE] cs=%d workers=%d", cs, wk))
    }
  }
  cat("[BENCH DONE] summary => ", out_csv, "\n", sep="")
}

if (sys.nframe()==0L) main()

# bench/summarize_bench_kpi.R
# No multibyte characters in code/comments.
# in R-Terminal, use like: Rscript bench/summarize_bench_kpi.R --src=test/bench/bench_run_match_fast_stress_2129.csv

suppressWarnings(suppressMessages({
  # base only
  # tools is base R; no extra packages required
}))

BENCH_DIR <- "test/bench"

list_bench_csv <- function(dir = BENCH_DIR) {
  if (!dir.exists(dir)) stop("bench dir not found: ", dir)
  fs <- list.files(dir, pattern = "^bench_run_match_fast_\\d{8}_\\d{6}\\.csv$", full.names = TRUE)
  fs[order(file.info(fs)$mtime, decreasing = TRUE)]
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  src <- NULL
  for (a in args) {
    if (startsWith(a, "--src=")) src <- sub("^--src=", "", a)
  }
  list(src = src)
}

read_target_bench <- function(src_opt) {
  if (!is.null(src_opt) && nzchar(src_opt)) {
    if (!file.exists(src_opt)) stop("src csv not found: ", src_opt)
    message("[KPI] using (explicit): ", src_opt)
    return(list(df = read.csv(src_opt, stringsAsFactors = FALSE), path = src_opt))
  }
  fs <- list_bench_csv()
  if (length(fs) == 0) stop("no bench csv in ", BENCH_DIR)
  latest <- fs[1]
  message("[KPI] using (latest): ", latest)
  list(df = read.csv(latest, stringsAsFactors = FALSE), path = latest)
}

num <- function(x) suppressWarnings(as.numeric(x))

aggregate_kpi <- function(df) {
  ok <- subset(df, (is.na(note) | note == "OK") & type %in% c("vdb_gf"))
  if (nrow(ok) == 0) stop("no OK rows to summarize")
  
  # normalize expected columns
  cols <- c("chunk_size","workers","chunk_index","N",
            "throughput_samples_per_sec","throughput_loci_per_sec",
            "ms_total",
            "mem_peak_mb","mem_rss_mb_before","mem_rss_mb_after",
            "cpu_total_sec","cpu_user_sec","cpu_system_sec",
            "io_read_bytes","io_write_bytes",
            "io_read_MBps","io_write_MBps")
  for (c in cols) if (!c %in% names(ok)) ok[[c]] <- NA
  
  ok$chunk_size  <- num(ok$chunk_size)
  ok$workers     <- num(ok$workers)
  ok$chunk_index <- num(ok$chunk_index)
  ok$N           <- num(ok$N)
  ok$thr_s       <- num(ok$throughput_samples_per_sec)
  ok$thr_L       <- num(ok$throughput_loci_per_sec)
  ok$ms_total    <- num(ok$ms_total)
  
  ok$mem_peak_mb <- num(ok$mem_peak_mb)
  ok$mem_bef_mb  <- num(ok$mem_rss_mb_before)
  ok$mem_aft_mb  <- num(ok$mem_rss_mb_after)
  
  ok$cpu_total   <- num(ok$cpu_total_sec)
  ok$cpu_user    <- num(ok$cpu_user_sec)
  ok$cpu_sys     <- num(ok$cpu_system_sec)
  
  ok$io_rb       <- num(ok$io_read_bytes)
  ok$io_wb       <- num(ok$io_write_bytes)
  ok$io_rMBps    <- num(ok$io_read_MBps)
  ok$io_wMBps    <- num(ok$io_write_MBps)
  
  ok$key <- paste(ok$chunk_size, ok$workers, sep="x")
  
  agg <- do.call(rbind, lapply(split(ok, ok$key), function(g) {
    data.frame(
      plan = unique(g$key)[1],
      chunk_size = unique(g$chunk_size)[1],
      workers    = unique(g$workers)[1],
      chunks     = length(unique(g$chunk_index)),
      N_total    = sum(g$N, na.rm=TRUE),
      
      thr_samples_mean = mean(g$thr_s, na.rm=TRUE),
      thr_samples_p95  = suppressWarnings(quantile(g$thr_s, 0.95, na.rm=TRUE)),
      thr_samples_min  = min(g$thr_s, na.rm=TRUE),
      
      thr_loci_mean    = mean(g$thr_L, na.rm=TRUE),
      
      ms_total_mean    = mean(g$ms_total, na.rm=TRUE),
      
      mem_peak_mb_mean   = mean(g$mem_peak_mb, na.rm=TRUE),
      mem_peak_mb_p95    = suppressWarnings(quantile(g$mem_peak_mb, 0.95, na.rm=TRUE)),
      mem_before_mb_mean = mean(g$mem_bef_mb, na.rm=TRUE),
      mem_after_mb_mean  = mean(g$mem_aft_mb, na.rm=TRUE),
      
      cpu_total_sec_sum = sum(g$cpu_total, na.rm=TRUE),
      cpu_user_sec_sum  = sum(g$cpu_user,  na.rm=TRUE),
      cpu_sys_sec_sum   = sum(g$cpu_sys,   na.rm=TRUE),
      
      io_read_MBps_mean  = mean(g$io_rMBps, na.rm=TRUE),
      io_write_MBps_mean = mean(g$io_wMBps, na.rm=TRUE),
      
      stringsAsFactors = FALSE
    )
  }))
  rownames(agg) <- NULL
  agg
}

write_summary <- function(agg, src_csv) {
  # derive tag safely: strip prefix + extension once
  base <- basename(src_csv)
  base_noext <- tools::file_path_sans_ext(base)      # remove .csv
  tag <- sub("^bench_run_match_fast_", "", base_noext) # remove prefix
  out <- file.path(BENCH_DIR, sprintf("bench_kpi_%s.csv", tag))
  write.csv(agg, out, row.names = FALSE)
  message("[KPI] wrote: ", out)
  out
}

main <- function() {
  a <- parse_args()
  obj <- read_target_bench(a$src)
  agg <- aggregate_kpi(obj$df)
  print(agg)
  write_summary(agg, obj$path)
}

if (sys.nframe()==0L) main()

# scripts/devtools/summarize_bench_times.R
# No multibyte chars. Requires: data.table

suppressPackageStartupMessages({
  library(data.table)
})

markers_dir <- "output/bench_markers"
results_dir <- "output/bench_runs"
out_csv     <- "output/bench_runs/timings_summary.csv"

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

# ---- helpers ----

# return list(type, S, N, impl, run)
parse_meta <- function(fname) {
  b <- basename(fname)
  
  # Normalize: remove extension
  stem <- sub("\\.[^.]+$", "", b)
  
  # Common fields init
  type <- NA_character_
  S    <- NA_integer_
  N    <- NA_integer_   # TopN; for top3 we'll set N=3
  impl <- NA_character_ # "R" | "CPP" | NA
  run  <- NA_integer_
  
  # Heuristics by prefix
  if (startsWith(stem, "start_")) {
    type <- if (grepl("^start_score_", stem)) "score" else "detail"
    rest <- sub("^start_(score|detail)_", "", stem)
  } else {
    # result
    if (startsWith(stem, "score_")) {
      type <- "score"
      rest <- sub("^score_", "", stem)
    } else if (startsWith(stem, "detail_")) {
      type <- "detail"
      rest <- sub("^detail_", "", stem)
    } else {
      return(list(type=NA,S=NA,N=NA,impl=NA,run=NA))
    }
  }
  
  # Expect S{digits}
  m <- regexpr("S([0-9]+)", rest, perl = TRUE)
  if (m > 0) {
    S <- as.integer(sub("S", "", regmatches(rest, m)))
  }
  
  # Top3 / TopN / TopN{N}
  if (grepl("top3", rest)) {
    N <- 3L
  } else {
    mN <- regexpr("topN([0-9]+)", rest, perl = TRUE)
    if (mN > 0) {
      N <- as.integer(sub("topN", "", regmatches(rest, mN)))
    } else if (grepl("topN($|[^0-9])", rest)) {
      # "topN" with no numeric → N=NA (後で --top の既定値等に依存)
      N <- NA_integer_
    }
  }
  
  # impl
  if (grepl("_CPP(_|$)", rest)) impl <- "CPP"
  if (grepl("_R(_|$)",   rest)) impl <- ifelse(is.na(impl), "R", impl)
  
  # run
  mrun <- regexpr("run([0-9]+)", rest, perl = TRUE)
  if (mrun > 0) {
    run <- as.integer(sub("run", "", regmatches(rest, mrun)))
  }
  
  list(type=type, S=S, N=N, impl=impl, run=run)
}

scan_dir <- function(d) {
  if (!dir.exists(d)) return(data.table())
  f <- list.files(d, full.names = TRUE)
  if (!length(f)) return(data.table())
  
  meta <- lapply(f, parse_meta)
  dt <- rbindlist(lapply(seq_along(f), function(i) {
    data.table(file=f[i],
               kind=if (startsWith(basename(f[i]), "start_")) "start" else "result",
               type = meta[[i]]$type,
               S    = meta[[i]]$S,
               N    = meta[[i]]$N,
               impl = meta[[i]]$impl,
               run  = meta[[i]]$run)
  }), fill=TRUE)
  
  # attach times
  fi <- file.info(dt$file)
  dt[, mtime := fi$mtime]
  dt
}

# ---- main ----

starts  <- scan_dir(markers_dir)
results <- scan_dir(results_dir)

# keep only recognized
starts  <- starts[!is.na(type)]
results <- results[!is.na(type)]

# Keys to join: (type,S,N,impl,run)
# Note: impl/run can be NA; join on those too.
setkeyv(starts,  c("type","S","N","impl","run","mtime"))
setkeyv(results, c("type","S","N","impl","run","mtime"))

# Best-effort join: exact on all keys first.
merged <- results[starts, on = .(type, S, N, impl, run), nomatch = 0L,
                  .(type, S, N, impl, run,
                    start_file = i.file, start_time = i.mtime,
                    result_file= x.file, result_time= x.mtime)]

# Fallback 1: if not matched due to impl/ run differing in naming,
# try ignoring impl & run (take nearest result after the start).
if (nrow(merged) == 0L) {
  # Coarse grouping
  coarse <- results[starts, on = .(type, S, N), allow.cartesian = TRUE,
                    .(type, S, N, impl_start = i.impl, run_start = i.run,
                      start_file = i.file, start_time = i.mtime,
                      result_file= x.file, result_time= x.mtime,
                      impl = x.impl, run = x.run)]
  # Keep only result_time >= start_time, and pick the earliest
  coarse <- coarse[result_time >= start_time]
  setorder(coarse, type, S, N, start_time, result_time)
  # Deduplicate per start_file (first match)
  coarse <- coarse[!duplicated(start_file)]
  # Harmonize cols
  merged <- coarse[, .(type, S, N, impl, run, start_file, start_time, result_file, result_time)]
}

if (!nrow(merged)) {
  message("[WARN] No matched start/result pairs found. Check file names.")
  fwrite(data.table(), out_csv)
  quit(status = 0)
}

merged[, dt_sec := as.numeric(difftime(result_time, start_time, units = "secs"))]
merged[, size_label := sprintf("S%g", S)]
merged[, top_label  := fifelse(is.na(N), "NA", sprintf("TopN%d", N))]
merged[, impl := ifelse(is.na(impl), "NA", impl)]

setorder(merged, type, S, N, impl, start_time)

# Write detailed table
fwrite(merged[, .(type, S, N, impl, run,
                  start_time, result_time, dt_sec,
                  start_file, result_file)], out_csv)

# Also print a small summary in console
cat("\n== Timings Summary (seconds, median by group) ==\n")
summ <- merged[, .(n=.N, sec_median=median(dt_sec), sec_mean=mean(dt_sec)),
               by=.(type, S, N, impl)]
setorder(summ, type, S, N, impl)
print(summ)

invisible(TRUE)

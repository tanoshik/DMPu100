# bench/bench_smoke_cpp_toggle_1M.R
# Smoke benchmark for run_match_fast() with use_cpp=FALSE vs TRUE on 1M DB.
# Writes to its own CSV (does NOT touch the legacy KPI file).
# No multibyte characters in code/comments.

suppressPackageStartupMessages({})

ANY_CODE <- 9999L

# ---- Load core scripts
source("scripts/scoring_fast.R")   # R impl (score_2x2 etc.)
source("scripts/matcher_fast.R")   # run_match_fast()
source("scripts/scoring_cpp.R")    # ensures compute_scores_uint16() is compiled/loaded

# ---- Helpers
.now_str <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
.ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

read_query_min <- function(csv = "data/query_profile.csv", locus_order = NULL, any_code = ANY_CODE) {
  if (!file.exists(csv)) stop(sprintf("Missing file: %s", csv))
  q <- read.csv(csv, stringsAsFactors = FALSE)
  names(q) <- sub("(?i)^locus$",    "Locus",    names(q), perl = TRUE)
  names(q) <- sub("(?i)^allele1$",  "Allele1",  names(q), perl = TRUE)
  names(q) <- sub("(?i)^allele2$",  "Allele2",  names(q), perl = TRUE)
  req <- c("Locus","Allele1","Allele2")
  if (!all(req %in% names(q))) stop("query CSV must have columns: Locus, Allele1, Allele2")
  
  map_any <- function(x) {
    x <- as.character(x); x[is.na(x) | x==""] <- "any"
    x[tolower(x)=="any"] <- as.character(ANY_CODE)
    suppressWarnings(as.integer(x))
  }
  q$allele1 <- map_any(q$Allele1)
  q$allele2 <- map_any(q$Allele2)
  
  if (!is.null(locus_order) && length(locus_order)) {
    ord <- match(q$Locus, locus_order)
    q <- q[order(ord), , drop = FALSE]
  }
  q[, c("Locus","allele1","allele2")]
}

# db_index must be a list(sample_ids, locus_ids, A1, A2)
read_db_indexed <- function(rds_path) {
  if (!file.exists(rds_path)) stop(sprintf("Missing file: %s", rds_path))
  db <- readRDS(rds_path)
  exp <- c("sample_ids","locus_ids","A1","A2")
  if (!all(exp %in% names(db))) stop("DB RDS must be an indexed list with: sample_ids, locus_ids, A1, A2")
  stopifnot(is.character(db$sample_ids), is.character(db$locus_ids))
  stopifnot(is.matrix(db$A1), is.matrix(db$A2))
  db
}

bench_once <- function(db_rds,
                       include_bits_in_detail = FALSE,
                       use_cpp = FALSE,
                       pre_add_for_any_any = 2L,
                       any_code = ANY_CODE) {
  db_index <- read_db_indexed(db_rds)
  q_prep   <- read_query_min("data/query_profile.csv", db_index$locus_ids, any_code)
  
  is_any_any <- (q_prep$allele1 == any_code & q_prep$allele2 == any_code)
  n_eval_loci <- sum(!is_any_any)
  n_samples <- length(db_index$sample_ids)
  est_evals <- as.double(n_eval_loci) * as.double(n_samples)
  
  t0 <- proc.time()[["elapsed"]]
  res <- run_match_fast(q_prep, db_index,
                        pre_add_for_any_any = pre_add_for_any_any,
                        any_code = any_code,
                        include_bits_in_detail = include_bits_in_detail,
                        use_cpp = use_cpp)
  t1 <- proc.time()[["elapsed"]]
  elapsed <- t1 - t0
  
  list(
    elapsed_sec = elapsed,
    est_evals   = est_evals,
    throughput_evals_per_sec = if (elapsed > 0) est_evals / elapsed else NA_real_,
    n_samples = n_samples,
    n_loci = length(db_index$locus_ids),
    result = res
  )
}

write_kpi <- function(kpi_row, out_csv) {
  hdr <- !file.exists(out_csv)
  con <- file(out_csv, open = if (hdr) "w" else "a", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  if (hdr) writeLines(paste(names(kpi_row), collapse = ","), con)
  writeLines(paste(kpi_row, collapse = ","), con)
}

# ---- Params (1M DB; do NOT touch legacy KPI) ----
db_rds  <- "data/virtual_db_u100_S1000000_seed123.rds"   # ← 1Mに固定
include_bits <- FALSE
out_csv <- "test/bench/bench_kpi_cpp_toggle_smoke_1M.csv"  # ← 新規ファイル（既存KPIに追記しない）

.ensure_dir(dirname(out_csv))
cat(sprintf("[SMOKE-1M] DB=%s  include_bits=%s\n", db_rds, include_bits))

# R path
cat("[SMOKE-1M] use_cpp=FALSE ...\n")
r <- bench_once(db_rds, include_bits_in_detail = include_bits, use_cpp = FALSE)
cat(sprintf("  elapsed=%.3f sec, est_evals=%.0f, thr=%.1f eval/s\n",
            r$elapsed_sec, r$est_evals, r$throughput_evals_per_sec))

row_r <- c(
  timestamp = .now_str(),
  db_rds = db_rds,
  db_samples = r$n_samples,
  n_loci = r$n_loci,
  use_cpp = "FALSE",
  include_bits = if (include_bits) "TRUE" else "FALSE",
  elapsed_sec = sprintf("%.6f", r$elapsed_sec),
  est_evals = sprintf("%.0f", r$est_evals),
  throughput_eval_per_sec = sprintf("%.3f", r$throughput_evals_per_sec)
)
write_kpi(row_r, out_csv)

# CPP path
cat("[SMOKE-1M] use_cpp=TRUE ...\n")
c <- bench_once(db_rds, include_bits_in_detail = include_bits, use_cpp = TRUE)
cat(sprintf("  elapsed=%.3f sec, est_evals=%.0f, thr=%.1f eval/s\n",
            c$elapsed_sec, c$est_evals, c$throughput_evals_per_sec))

row_c <- c(
  timestamp = .now_str(),
  db_rds = db_rds,
  db_samples = c$n_samples,
  n_loci = c$n_loci,
  use_cpp = "TRUE",
  include_bits = if (include_bits) "TRUE" else "FALSE",
  elapsed_sec = sprintf("%.6f", c$elapsed_sec),
  est_evals = sprintf("%.0f", c$est_evals),
  throughput_eval_per_sec = sprintf("%.3f", c$throughput_evals_per_sec)
)
write_kpi(row_c, out_csv)

cat("\n[SMOKE-1M] Done. KPI appended to: ", out_csv, "\n", sep = "")

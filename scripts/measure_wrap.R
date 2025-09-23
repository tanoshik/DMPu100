# scripts/measure_wrap.R
# Lightweight external measurer for time.log (LOAD_DB_SEC, LOAD_Q_SEC, TOTAL_SEC, PEAK_MiB)
# Usage (either):
#   Rscript scripts/measure_wrap.R --db <rds> --query <rds> --out_dir <dir>
#   Rscript scripts/measure_wrap.R --out_dir <dir>   # infer db/query from opts_scores.json

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--db",      type="character", default=NA),
  make_option("--query",   type="character", default=NA),
  make_option("--out_dir", type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$out_dir)) {
  stop("usage: Rscript scripts/measure_wrap.R --db <rds> --query <rds> --out_dir <dir>  OR  --out_dir <dir> (infer)")
}
out_dir <- opt$out_dir
if (!dir.exists(out_dir)) stop(sprintf("out_dir not found: %s", out_dir))

# infer db/query from opts if missing
db_path <- opt$db
q_path  <- opt$query
opts_json <- file.path(out_dir, "opts_scores.json")
if ((is.na(db_path) || is.na(q_path)) && file.exists(opts_json)) {
  js <- tryCatch(jsonlite::fromJSON(opts_json, simplifyVector = TRUE), error = function(e) NULL)
  if (!is.null(js)) {
    if (is.na(db_path) && !is.null(js$db_path)) db_path <- js$db_path
    if (is.na(q_path)  && !is.null(js$query_path)) q_path <- js$query_path
  }
}
if (is.na(db_path) || is.na(q_path)) {
  stop("cannot determine db/query: pass --db/--query or put opts_scores.json in out_dir")
}
if (!file.exists(db_path)) stop(sprintf("db not found: %s", db_path))
if (!file.exists(q_path))  stop(sprintf("query not found: %s", q_path))

# measure
peak_mib <- NA_real_
load_db  <- NA_real_
load_q   <- NA_real_
total    <- NA_real_

if (requireNamespace("peakRAM", quietly=TRUE)) {
  pm <- peakRAM::peakRAM({
    tdb <- system.time({ db <- readRDS(db_path) })
    tqr <- system.time({ q  <- readRDS(q_path) })
    load_db <<- as.numeric(tdb["elapsed"])
    load_q  <<- as.numeric(tqr["elapsed"])
  })
  peak_mib <- as.numeric(pm$Peak_RAM_Used_MiB)
  total    <- as.numeric(pm$Elapsed_Time_sec)
} else {
  t0  <- proc.time()[["elapsed"]]
  tdb <- system.time({ db <- readRDS(db_path) })
  tqr <- system.time({ q  <- readRDS(q_path) })
  load_db <- as.numeric(tdb["elapsed"])
  load_q  <- as.numeric(tqr["elapsed"])
  total   <- proc.time()[["elapsed"]] - t0
}

# name: reuse last line if present
logp <- file.path(out_dir, "time.log")
name <- "GF/normal/0.0/off/NA"
if (file.exists(logp)) {
  lines <- readLines(logp, warn = FALSE)
  if (length(lines) >= 2L) {
    last <- lines[length(lines)]
    toks <- strsplit(last, ",")[[1]]
    if (length(toks) >= 6L) name <- toks[6]
  }
}

header <- "LOAD_DB_SEC,LOAD_Q_SEC,COMP_SEC,TOTAL_SEC,PEAK_MiB,name"
line   <- sprintf("%.3f,%.3f,%s,%.3f,%.1f,%s",
                  load_db, load_q, "NA", total,
                  ifelse(is.na(peak_mib), NA_real_, peak_mib), name)

if (!file.exists(logp)) writeLines(header, logp)
cat(paste0(line, "\n"), file = logp, append = TRUE)

cat(sprintf("[OK] measure updated: LOAD_DB=%.3f, LOAD_Q=%.3f, TOTAL=%.3f, PEAK=%s\n",
            load_db, load_q, total, ifelse(is.na(peak_mib), "NA", sprintf("%.1f", peak_mib))))

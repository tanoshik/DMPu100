# scripts/cli/run_consistency_all.R
# Minimal launcher for nightly consistency battery.
# No multibyte chars in code/comments.

suppressPackageStartupMessages({
  library(optparse)
})

opt <- OptionParser() |>
  add_option(c("-m","--mode"), type="character", default="quick",
             help="consistency mode: quick or full") |>
  add_option(c("-q","--query"), type="character", default="data/query_profile_seed123.csv",
             help="query CSV path") |>
  add_option(c("-o","--outdir"), type="character", default="test/consistency",
             help="output dir for logs") |>
  parse_args()

dbs <- c(
  "data/virtual_db_u100_S1000_seed123.rds",
  "data/virtual_db_u100_S100000_seed123.rds",
  "data/virtual_db_u100_S1000000_seed123.rds",
  "data/virtual_db_u100_S2000000_seed123.rds"
)

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

src <- "scripts/cli/dmp_consistency_battery.R"   # exists in dist :contentReference[oaicite:4]{index=4}
source(src)

for (db in dbs) {
  if (!file.exists(db)) { message("[SKIP] ", db); next }
  message("[RUN] ", db, " mode=", opt$mode)
  try({
    dmp_consistency_battery(
      db_path   = db,
      query_csv = opt$query,
      mode      = opt$mode,  # quick: sample subset / full: all
      out_dir   = opt$outdir
    )
  }, silent = TRUE)
}
message("[DONE] consistency all")

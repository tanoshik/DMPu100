# scripts/cli/dmp_make_virtual_db.R
suppressPackageStartupMessages({ library(optparse) })

opt <- OptionParser() |>
  add_option("--size", type="integer", default=100000,
             help="db size (profiles)") |>
  add_option("--seed", type="integer", default=123,
             help="random seed") |>
  add_option("--out", type="character", default=NULL,
             help="output RDS path (optional: copy of default indexed RDS)") |>
  parse_args()

source("scripts/devtools/make_virtual_db_small.R")  # included in dist :contentReference[oaicite:1]{index=1}

set.seed(opt$seed)
invisible(make_virtual_db_small(n = opt$size))  # 既定RDSを自動生成する

# 既定RDS（＝インデックス済み）のパス
default_rds <- sprintf("data/virtual_db_u100_S%d_seed%d.rds", opt$size, opt$seed)
if (!file.exists(default_rds)) stop("default RDS not found: ", default_rds)

# --out が指定されていればコピーして提供
if (!is.null(opt$out) && normalizePath(opt$out, winslash="/", mustWork=FALSE) !=
    normalizePath(default_rds, winslash="/", mustWork=FALSE)) {
  dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
  file.copy(default_rds, opt$out, overwrite = TRUE)
  cat("[OK] wrote:", opt$out, "\n")
} else {
  cat("[OK] default indexed RDS:", default_rds, "\n")
}

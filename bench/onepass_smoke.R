# bench/onepass_smoke.R
suppressPackageStartupMessages({ library(tools) })

# Example:
# Rscript bench/onepass_smoke.R data/virtual_db_u100_S1000_seed101.rds data/query_001.rds output/scores.csv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: Rscript bench/onepass_smoke.R <db.rds> <query.rds> <out.csv>")
}
db <- args[[1]]; qu <- args[[2]]; out <- args[[3]]

cmd <- sprintf(
  'Rscript scripts/matcher_onepass.R --db %s --query %s --out %s --mode topn --top_n 200 --display_limit 200',
  shQuote(db), shQuote(qu), shQuote(out)
)
cat("[RUN]", cmd, "\n")
st <- system(cmd)
if (st != 0) stop("matcher_onepass failed")

if (!file.exists(out)) stop("scores.csv not found after run")
cat("[OK] smoke success; file:", out, "size:", file.info(out)$size, "bytes\n")

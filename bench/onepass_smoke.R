# bench/onepass_smoke.R
suppressPackageStartupMessages({ library(tools) })

# Example:
# Rscript bench/onepass_smoke.R data/virtual_db_u100_S1000_seed101.rds data/query_profile_seed101.csv output/scores.csv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: Rscript bench/onepass_smoke.R <db.rds> <query.csv|rds> <out.csv>")
}
db <- args[[1]]; qu <- args[[2]]; out <- args[[3]]

cmd <- sprintf('Rscript scripts/matcher_onepass.R --db "%s" --query "%s" --out "%s" --report all --n_cap 200',
               db, qu, out)
cat(cmd, "\n")
system(cmd)

# test/smoke_cli_entry.R
# Quick smoke test for new CLI scripts
# No multibyte chars in code/comments.

cat("=== Smoke test: CLI entrypoints ===\n")

# 1. Make a tiny virtual DB (this writes the indexed default RDS)
cat("[STEP] make_virtual_db (S=1000)\n")
system("Rscript scripts/cli/dmp_make_virtual_db.R --size=1000 --seed=123")

# ここがポイント：インデックス済みの既定RDSを使う
db_path <- "data/virtual_db_u100_S1000_seed123.rds"

stopifnot(file.exists(db_path))

# 2. Run match once (engine=r)
cat("[STEP] run_match (engine=r)\n")

# build args as '--key=value' tokens and use system2
args <- c(
  "scripts/cli/dmp_run_match.R",
  sprintf("--db=%s",   normalizePath(db_path, winslash = "/", mustWork = TRUE)),
  sprintf("--query=%s", "data/query_profile_seed123.csv"),
  sprintf("--out=%s",   "output/smoke_summary.csv")
)
rc <- system2("Rscript", args)
if (rc != 0) stop("dmp_run_match.R failed with code=", rc)

stopifnot(file.exists("output/smoke_summary.csv"))

# 3. Consistency battery quick-mode on 1k (skip large DBs)
cat("[STEP] consistency battery (quick, 1k only)\n")
system("Rscript scripts/cli/run_consistency_all.R --mode=quick --query=data/query_profile_seed123.csv --outdir=test/consistency_smoke")

cat("[OK] smoke_cli_entry finished\n")

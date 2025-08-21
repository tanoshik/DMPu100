# test/bench_run_match_u100.R
# Minimal bench for u100 run_match_fast.

source("scripts/scoring_fast.R", local = TRUE)
source("scripts/matcher_fast.R", local = TRUE)

# Load or mock query/db
q <- read.csv("data/query_profile.csv", stringsAsFactors = FALSE)
db <- read.csv("data/database_profile.csv", stringsAsFactors = FALSE)

# Optional: locus_order
locus_order <- if (file.exists("data/locus_order.rds")) readRDS("data/locus_order.rds") else unique(q$Locus)

cat("Samples in DB:", length(unique(db$SampleID)), "\n")
t0 <- proc.time()
res <- run_match_fast(q, db, locus_order = locus_order)
dt <- (proc.time() - t0)[["elapsed"]]
cat("Elapsed:", dt, "sec\n")
print(head(res$summary, 10))

# Save outputs compat with GUI
dir.create("output", showWarnings = FALSE)
write.csv(res$detail,  "output/match_log.csv", row.names = FALSE)
write.csv(res$summary, "output/match_scores.csv", row.names = FALSE)

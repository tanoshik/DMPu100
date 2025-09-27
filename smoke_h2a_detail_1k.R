#!/usr/bin/env Rscript
# Minimal one-shot to produce raw_detail.csv for 1K DB with h2a=on.

suppressPackageStartupMessages({ library(jsonlite) })

stop_on_fail <- function(ok, msg) if (!ok) stop(msg)

root <- getwd()
db  <- file.path("data", "virtual_db_u100_S1000_seed101.rds")
q   <- file.path("data", "query_profile_seed101.csv")
any <- 9999L
od  <- file.path("output", "h2a_detail_1k_on")
od_detail <- file.path(od, "detail")
dir.create(od_detail, recursive=TRUE, showWarnings=FALSE)

stop_on_fail(file.exists(db),  sprintf("DB not found: %s", db))
stop_on_fail(file.exists(q),   sprintf("Query not found: %s", q))

# 1) scores (n_cap=1000) with h2a=on
scores_csv <- file.path(od, "scores.csv")
cmd1 <- c("scripts/matcher_onepass.R",
          "--db", db, "--query", q,
          "--out", scores_csv,
          "--report", "all",
          "--n_cap", "1000",
          "--any_code", as.character(any),
          "--h2a_on", "1",
          "--debug", "0")
message("[RUN] Rscript ", paste(cmd1, collapse=" "))
st1 <- system2("Rscript", cmd1, stdout=TRUE, stderr=TRUE)
writeLines(st1, file.path(od, "runner_scores.log"))

# 2) raw_detail for every SampleID present in scores.csv
opt_json <- file.path(od, "opts_scores.json")
stop_on_fail(file.exists(opt_json), sprintf("Missing %s (step1 did not finish?)", opt_json))

cmd2 <- c("scripts/emit/emit_detail.R",
          "--opt_path", opt_json,
          "--scores_path", scores_csv,
          "--out_dir", od_detail,
          "--any_code", as.character(any))
message("[RUN] Rscript ", paste(cmd2, collapse=" "))
st2 <- system2("Rscript", cmd2, stdout=TRUE, stderr=TRUE)
writeLines(st2, file.path(od, "runner_detail.log"))

# 3) sanity: check that per-sample sum(detail$score) == scores$Score
scores <- read.csv(scores_csv, stringsAsFactors=FALSE, check.names=FALSE)
raw    <- read.csv(file.path(od_detail, "raw_detail.csv"), stringsAsFactors=FALSE, check.names=FALSE)
agg <- aggregate(score ~ SampleID, raw, sum)
chk <- merge(scores[c("SampleID","Score")], agg, by="SampleID", all.x=TRUE)
ok  <- all.equal(chk$Score, chk$score)
if (!isTRUE(ok)) {
  writeLines(c("[WARN] score mismatch detected", capture.output(print(chk))), file.path(od, "sanity_mismatch.txt"))
  message("[WARN] some Score mismatched; see sanity_mismatch.txt")
} else {
  write.csv(chk, file.path(od, "sanity_check.csv"), row.names=FALSE)
  message("[OK] sanity_check.csv written")
}

message(sprintf("[DONE] raw_detail.csv -> %s", file.path(od_detail, "raw_detail.csv")))

#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(jsonlite) })

pick_db <- function() {
  cands <- c("data/virtual_db_u100_S2000000_seed101.rds",
             "data/virtual_db_u100_S1000_seed101.rds")
  for (p in cands) if (file.exists(p)) return(p)
  stop("No DB RDS found: need 2M or 1k sample RDS")
}

db  <- pick_db()
q   <- "data/query_profile_seed101.csv"
any <- 9999L
out <- "output/smoke_h2a_2M"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

run_one <- function(h2a_on) {
  tag <- if (h2a_on==1L) "on" else "off"
  od  <- file.path(out, tag)
  dir.create(od, recursive=TRUE, showWarnings=FALSE)
  cmd <- c("scripts/matcher_onepass.R",
           "--db", db, "--query", q,
           "--out", file.path(od, "scores.csv"),
           "--report", "all",
           "--n_cap", "1000",
           "--any_code", as.character(any),
           "--h2a_on", as.character(h2a_on),
           "--debug", "0")
  t0 <- Sys.time()
  st <- tryCatch(system2("Rscript", cmd, stdout=TRUE, stderr=TRUE),
                 error=function(e) sprintf("ERROR: %s", conditionMessage(e)))
  t1 <- Sys.time()
  writeLines(st, file.path(od, "runner.log"))
  data.frame(
    mode        = tag,
    db_path     = normalizePath(db, winslash="/", mustWork=FALSE),
    elapsed_sec = as.numeric(difftime(t1, t0, units="secs")),
    scores_sz   = if (file.exists(file.path(od,"scores.csv"))) file.info(file.path(od,"scores.csv"))$size else NA_real_,
    hist_sz     = if (file.exists(file.path(od,"hist.csv")))   file.info(file.path(od,"hist.csv"))$size   else NA_real_
  )
}

res <- rbind(run_one(0L), run_one(1L))
write.csv(res, file.path(out, "summary.csv"), row.names=FALSE)
cat(sprintf("[2M-h2a] wrote %s\n", file.path(out, "summary.csv")))

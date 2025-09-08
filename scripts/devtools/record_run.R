# scripts/devtools/record_run.R
record_run <- function(seed, size, db_path, query_path, sum_sec, detail_sec, match_flag) {
  dir.create("output/bench_ledger", recursive = TRUE, showWarnings = FALSE)
  out <- "output/bench_ledger/bench_runs.csv"
  if (!file.exists(out)) {
    writeLines(
      "seed,size,db_path,query_path,sum_sec,detail_sec,top1_match,ts,commit,host,os,cpu,ram_gb,gpu",
      out
    )
  }
  row <- sprintf(
    "%d,%d,%s,%s,%g,%g,%s,%s,%s,%s,%s,%s,%s,%s",
    as.integer(seed), as.integer(size),
    normalizePath(db_path, winslash="/", mustWork = FALSE),
    normalizePath(query_path, winslash="/", mustWork = FALSE),
    as.numeric(sum_sec), as.numeric(detail_sec), as.character(match_flag),
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    tryCatch(system('git rev-parse --short HEAD', intern=TRUE), error=function(e) "unknown"),
    Sys.info()[["nodename"]],
    sprintf("%s %s", Sys.info()[["sysname"]], Sys.info()[["release"]]),
    "Ryzen 7 5700G","32","RTX 4060 Ti"
  )
  write(row, file = out, append = TRUE)
}

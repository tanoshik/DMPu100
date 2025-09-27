# smoke_r0to1_2M.R
# Run r0 / r0.5 / r1 on the 2M DB with report=all.
# ASCII only. No extra packages needed.

normalize <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

root <- normalize(".")
setwd(root)

# --- paths (edit if needed) ---
DB_2M  <- "data/virtual_db_u100_S2000000_seed101.rds"
QCSV   <- "data/query_profile_seed101.csv"
IDFCSV <- "data/gf_idf_mask.csv"
OUTDIR <- "output/smoke_r0to1_2M"
ANY    <- 9999L
SEED   <- 0L

# Rscript resolver (PATH優先、未検出ならWindowsの典型を順に)
rs <- Sys.which("Rscript")
if (!nzchar(rs)) {
  cand <- c("C:/Program Files/R/R-4.5.1/bin/Rscript.exe",
            "C:/Program Files/R/R-4.5.0/bin/Rscript.exe")
  rs <- cand[file.exists(cand)][1]
}
stopifnot(nzchar(rs) && file.exists(rs))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

label_ratio <- function(x){
  if (x == 0) "r0" else if (x == 1) "r1" else "r0p5"
}
masked_name <- function(x){
  if (x == 0) {
    "data/_masked_r0_2M_seed0.rds"
  } else if (x == 1) {
    "data/_masked_r1_2M_seed0.rds"
  } else {
    "data/_masked_r0p5_2M_seed0.rds"
  }
}

run_r <- function(args, log_file = NULL) {
  # args: character vector; first element should be the script path or '-e'
  if (is.null(log_file)) {
    system2(rs, c("--vanilla", args), stdout = "", stderr = "")
  } else {
    con <- file(log_file, open = "wt")
    on.exit(close(con), add = TRUE)
    system2(rs, c("--vanilla", args), stdout = con, stderr = con)
  }
}

summ_rows <- list()

for (ratio in c(0, 0.5, 1)) {
  tag <- label_ratio(ratio)
  out_sub <- file.path(OUTDIR, paste0("ratio", if (ratio == 0.5) "0p5" else ratio))
  dir.create(out_sub, showWarnings = FALSE, recursive = TRUE)
  
  masked <- masked_name(ratio)
  
  message(sprintf("[2M] STEP1 mask ratio=%s -> %s", ratio, masked))
  log_mask <- file.path(out_sub, "_log_mask.txt")
  run_r(c("scripts/cli/dmp_apply_idf_mask.R",
          "--in_db",  shQuote(DB_2M),
          "--out_db", shQuote(masked),
          "--idf_csv", shQuote(IDFCSV),
          "--ratio", as.character(ratio),
          "--seed",  as.character(SEED),
          "--any_code", as.character(ANY),
          "--debug", "0"), log_file = log_mask)
  if (!file.exists(masked)) stop("masked DB not created: ", masked)
  
  message(sprintf("[2M] STEP2 match ratio=%s (report=all)", ratio))
  out_scores <- file.path(out_sub, "scores.csv")
  run_r(c("scripts/matcher_onepass.R",
          "--db",    shQuote(masked),
          "--query", shQuote(QCSV),
          "--out",   shQuote(out_scores),
          "--report","all",
          "--n_cap","2000000",
          "--any_code", as.character(ANY),
          "--h2a_on","0",
          "--debug","0"))
  
  # quick summary from scores
  if (file.exists(out_scores)) {
    df <- tryCatch(read.csv(out_scores, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) {
      sc_col <- if ("Score" %in% names(df)) "Score" else names(df)[2]
      sc <- suppressWarnings(as.integer(df[[sc_col]]))
      sc <- sc[!is.na(sc)]
      if (length(sc) > 0) {
        srow <- data.frame(
          run    = paste0("ratio=", ratio),
          n      = length(sc),
          min    = min(sc),
          p25    = as.numeric(quantile(sc, 0.25, names = FALSE, type = 7)),
          median = as.numeric(quantile(sc, 0.50, names = FALSE, type = 7)),
          mean   = mean(sc),
          p75    = as.numeric(quantile(sc, 0.75, names = FALSE, type = 7)),
          max    = max(sc),
          stringsAsFactors = FALSE
        )
        summ_rows[[length(summ_rows)+1]] <- srow
      }
    }
  }
}

# write summary CSV
if (length(summ_rows) > 0) {
  summ <- do.call(rbind, summ_rows)
  write.csv(summ, file.path(OUTDIR, "summary_stats.csv"), row.names = FALSE)
  message("[2M] wrote summary: ", normalize(file.path(OUTDIR, "summary_stats.csv")))
}

message("[2M] done.")

# scripts/devtools/merge_bench_runs.R
suppressWarnings(suppressMessages({
  if (!requireNamespace("readr", quietly=TRUE)) install.packages("readr", repos="https://cran.rstudio.com")
  if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr", repos="https://cran.rstudio.com")
}))

library(readr); library(dplyr)

csv_out   <- "output/bench_ledger/bench_runs.csv"
log_f     <- "output/bench_logs/bench_sanity_all.log"    # メイン
sum_f_opt <- "output/bench_logs/sanity_summary.csv"      # あれば補助で読む

# 1) ログから拾う（第一候補）
parse_from_log <- function(path) {
  if (!file.exists(path)) return(tibble())
  txt <- readr::read_lines(path)
  # 例: [START] S=100000000 ...
  #     [ELAPSED] sum=1871s  detail=1742s
  #     [DONE] S=100000000 ...
  s_pat  <- "^\\[START\\] S=([0-9]+)\\b"
  e_pat  <- "^\\[ELAPSED\\] sum=([0-9]+)s\\s+detail=([0-9]+)s"
  d_pat  <- "^\\[DONE\\] S=([0-9]+)\\b"
  
  sizes <- as.integer(gsub(s_pat, "\\1", grep(s_pat, txt, value=TRUE)))
  # ELAPSED はサイズが無いので、直近の START のサイズを引き継ぐ
  rows <- tibble()
  curr_size <- NA_integer_
  for (ln in txt) {
    if (grepl(s_pat, ln)) {
      curr_size <- as.integer(sub(s_pat, "\\1", ln))
    } else if (grepl(e_pat, ln) && !is.na(curr_size)) {
      m <- regexec(e_pat, ln); mm <- regmatches(ln, m)[[1]]
      sum_sec    <- as.numeric(mm[2])
      detail_sec <- as.numeric(mm[3])
      rows <- bind_rows(rows, tibble(size=curr_size, sum_sec=sum_sec, detail_sec=detail_sec))
    }
  }
  # 重複は最後のもの優先でユニーク化
  rows %>% group_by(size) %>% summarise(sum_sec=last(sum_sec), detail_sec=last(detail_sec), .groups="drop")
}

# 2) sanity_summary.csv がある場合は補助（列名がケースバイケースなので素直に try）
parse_from_summary <- function(path) {
  if (!file.exists(path)) return(tibble())
  df <- tryCatch(readr::read_csv(path, show_col_types=FALSE), error=function(e) tibble())
  if (!nrow(df)) return(tibble())
  
  # size 列が無い時の救済（"S=..." 文字列があれば拾う）
  if (!("size" %in% names(df))) {
    cand <- names(df)[grepl("^S$", names(df), ignore.case=TRUE)]
    if (length(cand) == 1) df <- df %>% mutate(size = !!sym(cand))
  }
  # sum, detail 列の名前も場合によって違うかも
  s_col <- names(df)[grepl("sum(_sec)?$", names(df), ignore.case=TRUE)][1]
  d_col <- names(df)[grepl("detail(_sec)?$", names(df), ignore.case=TRUE)][1]
  if (is.na(s_col) || is.na(d_col) || !("size" %in% names(df))) return(tibble())
  
  df %>% transmute(size = as.integer(.data$size),
                   sum_sec = as.numeric(.data[[s_col]]),
                   detail_sec = as.numeric(.data[[d_col]])) %>%
    distinct(size, .keep_all = TRUE)
}

log_rows <- parse_from_log(log_f)
sum_rows <- parse_from_summary(sum_f_opt)

src <- bind_rows(log_rows, sum_rows) %>%
  group_by(size) %>% summarise(sum_sec=last(sum_sec), detail_sec=last(detail_sec), .groups="drop")

if (!nrow(src)) stop("[FATAL] no parsable rows in log or csv")

# 追加レコードを構築（共通メタは空でもOK：あとで埋めてOKとのこと）
seed  <- 123L
add <- src %>%
  mutate(
    seed       = seed,
    db_path    = sprintf("data/virtual_db_u100_S%s_seed%s.rds", size, seed),
    query_path = sprintf("data/query_profile_seed%s_S1000.csv", seed),
    top1_match = NA,   # 不明なら後で埋める
    ts         = format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'),
    commit     = tryCatch(system("git rev-parse --short HEAD", intern=TRUE), error=function(e) "unknown"),
    host       = NA_character_, os = NA_character_,
    cpu        = NA_character_, ram_gb = NA_character_, gpu = NA_character_
  ) %>%
  select(seed,size,db_path,query_path,sum_sec,detail_sec,top1_match,ts,commit,host,os,cpu,ram_gb,gpu)

dir.create(dirname(csv_out), recursive = TRUE, showWarnings = FALSE)
if (file.exists(csv_out)) {
  runs <- suppressWarnings(readr::read_csv(csv_out, show_col_types = FALSE))
  # 列揃え
  cols <- c("seed","size","db_path","query_path","sum_sec","detail_sec","top1_match",
            "ts","commit","host","os","cpu","ram_gb","gpu")
  miss <- setdiff(cols, names(runs))
  if (length(miss)) for (k in miss) runs[[k]] <- NA
  runs <- runs[, cols]
  add  <- add[,  cols]
  
  # 既存 (seed,size) は追記しない
  key_new <- paste(add$seed, add$size)
  key_old <- paste(runs$seed, runs$size)
  add <- add[!key_new %in% key_old, , drop = FALSE]
  
  out <- bind_rows(runs, add)
} else {
  out <- add
}

readr::write_csv(out, csv_out)
message("[OK] merged: ", nrow(add), " new rows into ", csv_out)

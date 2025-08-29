# scripts/cli/dmp_run_match.R
# Single-shot matcher CLI (engine=r first). No multibyte chars.

suppressPackageStartupMessages({ library(optparse) })

opt <- OptionParser() |>
  add_option("--db",     type="character", help="DB path (.rds or .csv)") |>
  add_option("--query",  type="character", default="data/query_profile_seed123.csv",
             help="Query CSV path") |>
  add_option("--engine", type="character", default="r",
             help="r or cpp (future)") |>
  add_option("--out",    type="character", default="output/summary.csv",
             help="output summary CSV") |>
  parse_args()

if (is.null(opt$db) || !file.exists(opt$db)) stop("DB path missing or not found")
if (!file.exists(opt$query)) stop("Query CSV missing")

# Load project functions
source("scripts/scoring_fast.R")
source("scripts/matcher_fast.R")

# --- read DB as db_index (RDS 推奨) ---
db_obj <- if (grepl("\\.rds$", opt$db, ignore.case = TRUE)) {
  readRDS(opt$db)
} else {
  # CSV 経路の場合は別途 indexer が必要（現状は RDS 前提推奨）
  stop("CSV DB not supported in this CLI yet: please use indexed RDS")
}

# （デバッグに役立つ軽ログ）
if (exists(".assert_db_index")) {
  tryCatch({ .assert_db_index(db_obj) },
           error = function(e) stop("Invalid db_index: ", conditionMessage(e)))
}

# --- read Query CSV -> data.frame ---
q_df <- read.csv(opt$query, stringsAsFactors = FALSE, check.names = FALSE)

# 1) ヘッダ正規化（堅牢）
raw_names <- names(q_df)
canon_key <- tolower(trimws(gsub("\\uFEFF", "", raw_names)))
map <- setNames(raw_names, canon_key)
alias <- list(
  locus    = c("locus","loc","marker","markers","loci"),
  allele1  = c("allele1","allele_1","a1","allele 1","allele-1"),
  allele2  = c("allele2","allele_2","a2","allele 2","allele-2"),
  sampleid = c("sampleid","sample_id","id","sample")
)
resolve <- function(keys){ hit <- intersect(keys, canon_key); if(length(hit)) map[hit[1]] else NA_character_ }
col_Locus   <- resolve(alias$locus)
col_allele1 <- resolve(alias$allele1)
col_allele2 <- resolve(alias$allele2)

if (any(is.na(c(col_Locus, col_allele1, col_allele2)))) {
  # locus_order を使って推測（先頭3列ヒューリスティック）
  if (ncol(q_df) >= 3) {
    locs <- try(readRDS("data/locus_order.rds"), silent = TRUE)
    if (!inherits(locs,"try-error") && is.character(locs)) {
      match_score <- sapply(seq_len(min(5, ncol(q_df))), function(i){
        sum(tolower(trimws(q_df[[i]])) %in% tolower(locs), na.rm=TRUE)
      })
      guess_loc_col <- if (any(match_score>0)) which.max(match_score) else 1
    } else guess_loc_col <- 1
    idx <- setdiff(1:3, guess_loc_col)
    col_Locus   <- names(q_df)[guess_loc_col]
    col_allele1 <- names(q_df)[idx[1]]
    col_allele2 <- names(q_df)[idx[2]]
  }
}

need <- c("Locus","allele1","allele2")
still_missing <- c(if (is.na(col_Locus)) "Locus",
                   if (is.na(col_allele1)) "allele1",
                   if (is.na(col_allele2)) "allele2")
if (length(still_missing)) {
  cat("[ERR] Query header names: ", paste(raw_names, collapse=", "), "\n")
  stop("Query CSV missing required columns (after normalization): ",
       paste(still_missing, collapse=", "))
}

# 2) 取り出し（文字型）
tmp <- data.frame(
  Locus   = as.character(q_df[[col_Locus]]),
  allele1 = as.character(q_df[[col_allele1]]),
  allele2 = as.character(q_df[[col_allele2]]),
  stringsAsFactors = FALSE, check.names = FALSE
)

# 3) 文字→数値（ANYは9999）。空文字/NAはANY扱い。
ANY_CODE <- 9999L
to_allele_num <- function(x) {
  x0 <- trimws(tolower(x))
  x0[x0=="" | is.na(x0)] <- "any"
  is_any <- (x0 == "any")
  suppressWarnings(num <- as.numeric(x0))
  # 非数値で "any" 以外が混ざっていないか検査
  bad <- which(is.na(num) & !is_any)
  if (length(bad)) {
    uniq <- unique(x[bad])
    stop("Non-numeric, non-any allele(s) detected: ",
         paste(head(uniq, 10), collapse=", "),
         if (length(uniq) > 10) "...", call. = FALSE)
  }
  num[is_any] <- ANY_CODE
  # 整数化（小数も想定されるため numeric のままでOK：is.numeric() を満たす）
  num
}

tmp$allele1 <- to_allele_num(tmp$allele1)
tmp$allele2 <- to_allele_num(tmp$allele2)

# 4) 最終 q_df として確定
q_df <- tmp

cat("[DBG] Query columns normalized as: ", paste(names(q_df), collapse=", "), "\n")
cat("[DBG] Allele1/Allele2 class: ", class(q_df$allele1), "/", class(q_df$allele2), "\n")
# --- call run_match_fast with NAMED arguments (順序取り違え防止) ---
res <- run_match_fast(
  q_df                   = q_df,
  db_index               = db_obj,
  pre_add_for_any_any    = TRUE,
  include_bits_in_detail = FALSE,
  use_cpp                = (tolower(opt$engine) == "cpp"),
  block_size             = 16000
)

# --- accept list(summary=...) or data.frame ---
summary_df <-
  if (is.list(res) && !is.null(res$summary)) res$summary else
    if (is.data.frame(res)) res else
      stop("Unexpected result shape; need data.frame or list(summary=...)")

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
write.csv(summary_df, opt$out, row.names = FALSE)
cat("[OK] wrote:", opt$out, "\n")

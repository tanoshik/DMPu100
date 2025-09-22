# scripts/fix_query.R
# Lenient ingestion for query CSV (補完＋標準化＋監査ログ)
# 出力: onepass がそのまま読める RDS (list(q1, q2)) と監査 CSV/JSON

suppressPackageStartupMessages({
  library(jsonlite)
})

fix_query <- function(in_csv, out_rds="query_fixed.rds",
                      audit_csv="query_audit.csv", audit_json="query_audit.json",
                      any_code=9999L, freq_table_path="data/freq_table.rds",
                      do_audit=TRUE) {
  if (!file.exists(in_csv)) stop(sprintf("input CSV not found: %s", in_csv))
  if (!file.exists(freq_table_path)) stop(sprintf("freq_table not found: %s", freq_table_path))
  freq_tbl <- readRDS(freq_table_path)
  
  x <- read.csv(in_csv, stringsAsFactors=FALSE, check.names=FALSE)
  nms <- tolower(names(x)); names(x) <- nms
  lc  <- which(nms %in% c("locus","marker"))[1]
  a1c <- which(nms %in% c("allele1","a1","q1"))[1]
  a2c <- which(nms %in% c("allele2","a2","q2"))[1]
  if (is.na(lc) || is.na(a1c) || is.na(a2c)) stop("CSV needs columns: Locus, Allele1, Allele2")
  
  # 出力ベクトル
  loci <- x[[lc]]
  q1 <- integer(length(loci))
  q2 <- integer(length(loci))
  audit <- if (isTRUE(do_audit)) {
    data.frame(Locus=loci, Raw1=x[[a1c]], Raw2=x[[a2c]],
               Fix1=NA_character_, Fix2=NA_character_, stringsAsFactors=FALSE)
  } else NULL
  
  # 小数→ユニット文字列（.00 は落とす）
  to_units_str <- function(v, any_code) {
    if (v == any_code) return("any")
    s <- sprintf("%.2f", as.numeric(v) / 100)
    while (substr(s, nchar(s), nchar(s)) == "0") s <- substr(s, 1, nchar(s)-1)
    if (substr(s, nchar(s), nchar(s)) == ".") s <- substr(s, 1, nchar(s)-1)
    s
  }
  
  # freq_table 照合：存在しない (locus, allele[unit]) は ANY
  if (is.data.frame(freq_tbl)) {
    nmsf <- tolower(names(freq_tbl)); names(freq_tbl) <- nmsf
    lc_f <- which(nmsf %in% c("locus","marker","locus_id"))[1]
    ac_f <- which(nmsf %in% c("allele","a","value","val"))[1]
    if (is.na(lc_f) || is.na(ac_f)) stop("freq_table must have columns: locus(+), allele")
    names(freq_tbl)[lc_f] <- "locus"
    names(freq_tbl)[ac_f] <- "allele"
  } else {
    stop("freq_table.rds must be a data.frame for this script")
  }
  
  to_int <- function(val, locus) {
    if (is.na(val) || val=="" || tolower(val)=="any") return(any_code)
    num <- suppressWarnings(as.numeric(val))
    if (is.na(num)) return(any_code)
    ok <- any(freq_tbl$locus == locus & suppressWarnings(as.numeric(freq_tbl$allele)) == num)
    if (!ok) return(any_code)
    as.integer(round(num*100))
  }
  
  for (i in seq_along(loci)) {
    a1 <- to_int(x[[a1c]][i], loci[i])
    a2 <- to_int(x[[a2c]][i], loci[i])
    # any右寄せ＋昇順
    if (a1==any_code & a2!=any_code) { tmp<-a1; a1<-a2; a2<-tmp }
    if (a1>a2 & a2!=any_code) { tmp<-a1; a1<-a2; a2<-tmp }
    q1[i] <- a1; q2[i] <- a2
    if (isTRUE(do_audit)) {
      audit$Fix1[i] <- to_units_str(a1, any_code)
      audit$Fix2[i] <- to_units_str(a2, any_code)
    }
  }
  
  # 出力先ディレクトリ作成
  ensure_parent <- function(p) { d <- dirname(p); if (!dir.exists(d)) dir.create(d, recursive=TRUE, showWarnings=FALSE) }
  ensure_parent(out_rds)
  saveRDS(list(q1=q1,q2=q2), out_rds)
  
  if (isTRUE(do_audit)) {
    ensure_parent(audit_csv); ensure_parent(audit_json)
    write.csv(audit, audit_csv, row.names=FALSE)
    writeLines(jsonlite::toJSON(audit, auto_unbox=TRUE, pretty=TRUE), audit_json)
    message(sprintf("[OK] wrote %s, %s, %s", out_rds, audit_csv, audit_json))
  } else {
    message(sprintf("[OK] wrote %s", out_rds))
  }
  invisible(list(q1=q1,q2=q2))
}

# ---- CLI entry ----
if (sys.nframe() == 0) {
  # 超簡易引数パーサ（--key value）
  parse_args <- function() {
    a <- commandArgs(trailingOnly = TRUE)
    kv <- list(); i <- 1
    while (i <= length(a)) {
      k <- a[i]
      if (substr(k,1,2) == "--") {
        key <- substring(k,3)
        if (i+1 <= length(a) && substr(a[i+1],1,2) != "--") { kv[[key]] <- a[i+1]; i <- i + 2 }
        else { kv[[key]] <- "1"; i <- i + 1 }
      } else i <- i + 1
    }
    kv
  }
  args <- parse_args()
  
  # kebab-case → snake_case のエイリアス解決（A.4 準拠・後方互換）
  alias <- function(from, to) { if (!is.null(args[[from]])) args[[to]] <<- args[[from]] }
  alias("out-rds",       "out_rds")
  alias("no-audit",      "no_audit")
  alias("any-code",      "any_code")
  alias("freq-table",    "freq_table")
  alias("out-csv",       "out_csv")
  alias("summary-json",  "summary_json")
  
  # 必須
  in_csv  <- args[["query"]]
  out_rds <- args[["out_rds"]]
  if (is.null(in_csv) || is.null(out_rds)) {
    stop(paste0(
      "usage: Rscript scripts/fix_query.R ",
      "--query <csv> --out-rds <rds> ",
      "[--any-code 9999] [--freq-table data/freq_table.rds] [--no-audit 1]\n",
      "  (aliases accepted: --out_rds, --any_code, --freq_table, --no_audit)"
    ))
  }
  
  any_code <- if (!is.null(args[["any_code"]])) as.integer(args[["any_code"]]) else 9999L
  freq_table_path <- if (!is.null(args[["freq_table"]])) args[["freq_table"]] else "data/freq_table.rds"
  do_audit <- is.null(args[["no_audit"]]) || as.integer(args[["no_audit"]]) == 0L
  
  # 監査ファイル名（明示なしなら out_rds と同じ場所）
  base_dir <- dirname(out_rds)
  audit_csv  <- if (!is.null(args[["out_csv"]])) args[["out_csv"]] else file.path(base_dir, "query_audit.csv")
  audit_json <- if (!is.null(args[["summary_json"]])) args[["summary_json"]] else file.path(base_dir, "query_audit.json")
  
  fix_query(in_csv = in_csv,
            out_rds = out_rds,
            audit_csv = audit_csv,
            audit_json = audit_json,
            any_code = any_code,
            freq_table_path = freq_table_path,
            do_audit = do_audit)
}

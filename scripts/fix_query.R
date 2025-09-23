#!/usr/bin/env Rscript
# scripts/fix_query.R
# Lenient ingestion for query CSV: 補完(ANY) + 標準化 + 監査(任意)
# 出力: onepass がそのまま読める RDS(list(q1,q2))
suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

# --- CLI ---
option_list <- list(
  make_option("--query",      type="character", help="path to query CSV"),
  make_option("--out_rds",    type="character", help="path to output RDS(list(q1,q2))"),
  make_option("--db",         type="character", default=NA, help="(optional) DB RDS for locus order"),
  make_option("--freq_table", type="character", default="data/freq_table.rds", help="freq_table.rds"),
  make_option("--any_code",   type="integer",   default=9999L, help="ANY code (int)"),
  make_option("--no_audit",   type="integer",   default=0L,    help="1=監査CSV/JSONを出力しない"),
  make_option("--audit_csv",  type="character", default=NA,    help="(opt) audit CSV path"),
  make_option("--audit_json", type="character", default=NA,    help="(opt) audit JSON path")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$query) || is.null(opt$out_rds)) {
  stop("usage: Rscript scripts/fix_query.R --query <csv> --out_rds <rds> [--db <rds>] [--any_code 9999] [--freq_table data/freq_table.rds] [--no_audit 1]")
}
if (!file.exists(opt$query)) stop(sprintf("input CSV not found: %s", opt$query))
if (!file.exists(opt$freq_table)) stop(sprintf("freq_table not found: %s", opt$freq_table))

# --- load freq_table + locus universe ---
freq_tbl <- readRDS(opt$freq_table) # 期待: data.frame(Locus, Allele) など
if (!all(c("Locus","Allele") %in% names(freq_tbl))) {
  stop("freq_table.rds must contain columns: Locus, Allele")
}
# locus universe: DB優先、なければfreq_table順
locus_univ <- NULL
L <- NULL
if (!is.na(opt$db) && nzchar(opt$db)) {
  if (!file.exists(opt$db)) stop(sprintf("db not found: %s", opt$db))
  db <- readRDS(opt$db)
  if (!all(c("locus_ids") %in% names(db))) stop("DB RDS must contain locus_ids")
  locus_univ <- as.character(db$locus_ids)
  L <- length(locus_univ)
} else {
  locus_univ <- unique(as.character(freq_tbl$Locus))
  L <- length(locus_univ)
}

# --- read query CSV ---
x <- read.csv(opt$query, stringsAsFactors=FALSE, check.names=FALSE)
nms <- tolower(names(x)); names(x) <- nms
lc  <- which(nms %in% c("locus","marker"))[1]
a1c <- which(nms %in% c("allele1","a1","q1"))[1]
a2c <- which(nms %in% c("allele2","a2","q2"))[1]
if (is.na(lc) || is.na(a1c) || is.na(a2c)) stop("CSV needs columns: Locus, Allele1, Allele2")

ANY <- as.integer(opt$any_code)
q1 <- rep.int(ANY, L)
q2 <- rep.int(ANY, L)

# 文字列→数値(×100 int)、freq_tableに無い(ローカス×等位)はANY化
to_i100_checked <- function(val, locus) {
  if (is.na(val) || val=="" || tolower(val)=="any") return(ANY)
  num <- suppressWarnings(as.numeric(val))
  if (is.na(num)) return(ANY)
  # freq_table に存在しなければ ANY
  ok <- any(freq_tbl$Locus == locus & as.numeric(freq_tbl$Allele) == num)
  if (!ok) return(ANY)
  as.integer(round(num * 100))
}
# ANY右寄せ + 昇順
ord_pair <- function(a1, a2, ANY) {
  swap <- ((a1 == ANY & a2 != ANY) | (a1 > a2 & a2 != ANY))
  a1n <- ifelse(swap, a2, a1)
  a2n <- ifelse(swap, a1, a2)
  c(a1n, a2n)
}

# マッピング（不足ローカスは ANY のまま）
loc_in <- as.character(x[[lc]])
for (i in seq_len(nrow(x))) {
  locus <- loc_in[i]
  pos <- match(locus, locus_univ)
  if (is.na(pos)) next  # DBに存在しない locus は黙ってスキップ（strict入口じゃないので警告は出さない）
  a1 <- to_i100_checked(x[[a1c]][i], locus)
  a2 <- to_i100_checked(x[[a2c]][i], locus)
  pr <- ord_pair(a1, a2, ANY)
  q1[pos] <- pr[1]; q2[pos] <- pr[2]
}

# 監査
emit_audit <- (as.integer(opt$no_audit) != 1L)
if (emit_audit) {
  audit <- data.frame(
    Locus = locus_univ,
    Fix1  = ifelse(q1==ANY, "any", sub("\\.00$","",format(q1/100, nsmall=2))),
    Fix2  = ifelse(q2==ANY, "any", sub("\\.00$","",format(q2/100, nsmall=2))),
    stringsAsFactors = FALSE
  )
  aud_csv <- ifelse(is.na(opt$audit_csv) || !nzchar(opt$audit_csv), "query_audit.csv", opt$audit_csv)
  aud_json <- ifelse(is.na(opt$audit_json) || !nzchar(opt$audit_json), "query_audit.json", opt$audit_json)
  write.csv(audit, aud_csv, row.names=FALSE)
  writeLines(jsonlite::toJSON(audit, auto_unbox=TRUE, pretty=TRUE), aud_json)
}

saveRDS(list(q1=as.integer(q1), q2=as.integer(q2)), opt$out_rds)
cat(sprintf("[OK] wrote %s\n", opt$out_rds))

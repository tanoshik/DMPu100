# scripts/fix_query.R
# Lenient ingestion for query CSV (補完＋標準化＋監査ログ)
# 出力: onepass がそのまま読める RDS (list(q1, q2)) と監査 CSV/JSON

suppressPackageStartupMessages({
  library(jsonlite)
})

fix_query <- function(in_csv, out_rds="query_fixed.rds",
                      audit_csv="query_audit.csv", audit_json="query_audit.json",
                      any_code=9999L, freq_table_path="data/freq_table.rds") {
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
  audit <- data.frame(Locus=loci, Raw1=x[[a1c]], Raw2=x[[a2c]],
                      Fix1=NA, Fix2=NA, stringsAsFactors=FALSE)

  to_int <- function(val, locus) {
    if (is.na(val) || val=="" || tolower(val)=="any") return(any_code)
    num <- suppressWarnings(as.numeric(val))
    if (is.na(num)) return(any_code)
    # freq_tableに存在しない場合はANY
    if (!any(freq_tbl$Locus==locus & freq_tbl$Allele==num)) return(any_code)
    as.integer(round(num*100))
  }

  for (i in seq_along(loci)) {
    a1 <- to_int(x[[a1c]][i], loci[i])
    a2 <- to_int(x[[a2c]][i], loci[i])
    if (a1==any_code & a2!=any_code) { tmp<-a1; a1<-a2; a2<-tmp } # anyは右寄せ
    if (a1>a2 & a2!=any_code) { tmp<-a1; a1<-a2; a2<-tmp }        # 昇順
    q1[i] <- a1; q2[i] <- a2
    audit$Fix1[i] <- ifelse(a1==any_code,"any",sub("\.00$","",format(a1/100,nsmall=2)))
    audit$Fix2[i] <- ifelse(a2==any_code,"any",sub("\.00$","",format(a2/100,nsmall=2)))
  }

  saveRDS(list(q1=q1,q2=q2), out_rds)
  write.csv(audit, audit_csv, row.names=FALSE)
  writeLines(toJSON(audit, auto_unbox=TRUE, pretty=TRUE), audit_json)
  message(sprintf("[OK] wrote %s, %s, %s", out_rds, audit_csv, audit_json))
  invisible(list(q1=q1,q2=q2))
}

# scripts/fix_db.R
# DB RDS の標準化ユーティリティ (A1<=A2, ANY右寄せ, x100未満自動補正)

fix_db <- function(in_rds, out_rds="db_fixed.rds", any_code=9999L) {
  if (!file.exists(in_rds)) stop(sprintf("db RDS not found: %s", in_rds))
  db <- readRDS(in_rds)
  if (!is.list(db) || !all(c("sample_ids","locus_ids","A1","A2") %in% names(db))) {
    stop("DB RDS must contain list(sample_ids, locus_ids, A1, A2)")
  }
  A1 <- as.matrix(db$A1)
  A2 <- as.matrix(db$A2)

  # x100未適用を自動補正 (<100を検出したら×100)
  if (any(A1[1,] < 100 | A2[1,] < 100, na.rm=TRUE)) {
    message("[WARN] values <100 detected in first row; applying x100 scaling")
    A1 <- A1*100L; A2 <- A2*100L
  }

  ANY <- as.integer(any_code)
  for (i in seq_len(nrow(A1))) {
    for (j in seq_len(ncol(A1))) {
      a1 <- A1[i,j]; a2 <- A2[i,j]
      if (a1==ANY & a2!=ANY) { tmp<-a1; a1<-a2; a2<-tmp }
      if (a1>a2 & a2!=ANY) { tmp<-a1; a1<-a2; a2<-tmp }
      A1[i,j]<-a1; A2[i,j]<-a2
    }
  }
  db$A1 <- A1; db$A2 <- A2
  saveRDS(db, out_rds)
  message(sprintf("[OK] wrote %s", out_rds))
  invisible(db)
}

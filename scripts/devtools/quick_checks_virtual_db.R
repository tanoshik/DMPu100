# devtools/quick_checks_virtual_db.R

source("scripts/scoring_fast.R")
source("scripts/matcher_fast.R")

quick_check_virtual_db <- function(rds_path, sample_n = 10000L) {
  cat("[quick-check] loading", rds_path, "\n")
  db <- readRDS(rds_path)
  
  stopifnot(all(dim(db$A1) == dim(db$A2)))
  stopifnot(length(db$sample_ids) == nrow(db$A1))
  stopifnot(length(db$locus_ids)  == ncol(db$A1))
  
  any1 <- sum(db$A1 == 9999L)
  any2 <- sum(db$A2 == 9999L)
  cat("ANY codes: A1=", any1, " A2=", any2, "\n")
  
  # ヘテロ率 vs freq_table（上位数座位）
  ft <- readRDS("data/freq_table.rds")
  ln <- if ("Locus" %in% names(ft)) "Locus" else "locus"
  an <- if ("Allele" %in% names(ft)) "Allele" else "allele"
  fc <- tolower(names(ft))
  fcol <- names(ft)[match(TRUE, fc %in% c("frequency","freq","p","prob","probability"))]
  
  het_exp <- tapply(ft[[fcol]], ft[[ln]], function(p) 1 - sum(p^2))
  het_obs <- sapply(db$locus_ids, function(L) {
    a1 <- db$A1[, L == db$locus_ids]
    a2 <- db$A2[, L == db$locus_ids]
    mean(a1 != a2)
  })
  delta <- het_obs - het_exp[names(het_obs)]
  print(round(delta[1:5], 4))
  
  # 部分サンプルで run_match_fast
  take <- seq_len(min(sample_n, nrow(db$A1)))
  db_sub <- list(
    sample_ids = db$sample_ids[take],
    locus_ids  = db$locus_ids,
    A1 = db$A1[take, ],
    A2 = db$A2[take, ]
  )
  q <- data.frame(Locus=db$locus_ids, allele1="any", allele2="any", stringsAsFactors=FALSE)
  q$allele1[q$Locus=="D3S1358"] <- 15L
  q$allele2[q$Locus=="D3S1358"] <- 17L
  
  cat("running run_match_fast on", length(take), "samples...\n")
  print(system.time({
    res <- run_match_fast(q, db_sub,
                          pre_add_for_any_any = 2L,
                          include_bits_in_detail = FALSE)
  }))
  cat("[quick-check] done\n")
}

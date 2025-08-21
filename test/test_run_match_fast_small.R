# test/test_run_match_fast_small.R
# Small end-to-end test for run_match_fast() with new API (db_index matrices + int query)

# --- load core ---
source("scripts/scoring_fast.R")
source("scripts/matcher_fast.R")

# --- tiny helpers (test-local) ---
ANY_CODE <- 9999L

map_any <- function(x, any_code=ANY_CODE) {
  x <- as.character(x); x[is.na(x) | x==""] <- "any"
  x[tolower(x)=="any"] <- as.character(any_code)
  suppressWarnings(as.integer(x))
}

prepare_query_min <- function(q_raw, locus_order=NULL, any_code=ANY_CODE) {
  n <- names(q_raw)
  names(q_raw) <- sub("(?i)^locus$","Locus",n,perl=TRUE)
  names(q_raw) <- sub("(?i)^allele1$","allele1",names(q_raw),perl=TRUE)
  names(q_raw) <- sub("(?i)^allele2$","allele2",names(q_raw),perl=TRUE)
  q <- q_raw[, c("Locus","allele1","allele2"), drop=FALSE]
  q$Locus   <- as.character(q$Locus)
  q$allele1 <- map_any(q$allele1, any_code)
  q$allele2 <- map_any(q$allele2, any_code)
  if (!is.null(locus_order) && length(locus_order)>0L) {
    ord <- match(q$Locus, locus_order)
    q <- q[order(ord), , drop=FALSE]
  }
  rownames(q) <- NULL
  q
}

build_db_index_v1 <- function(db_prep, locus_order, any_code = ANY_CODE) {
  # db_prep: data.frame(SampleID,Locus,Allele1,Allele2) with integers; any->9999
  SIDs <- unique(db_prep$SampleID)
  LIDs <- as.character(locus_order)
  S <- length(SIDs); L <- length(LIDs)
  A1 <- matrix(any_code, nrow = S, ncol = L)
  A2 <- matrix(any_code, nrow = S, ncol = L)
  sid2i <- setNames(seq_len(S), SIDs)
  loc2j <- setNames(seq_len(L), LIDs)
  for (i in seq_len(nrow(db_prep))) {
    si <- sid2i[[ db_prep$SampleID[i] ]]
    lj <- loc2j[[ db_prep$Locus[i] ]]
    if (!is.na(si) && !is.na(lj)) {
      A1[si, lj] <- as.integer(db_prep$Allele1[i])
      A2[si, lj] <- as.integer(db_prep$Allele2[i])
    }
  }
  list(sample_ids = SIDs, locus_ids = LIDs, A1 = A1, A2 = A2)
}

# --- build a tiny DB (3 samples x 2 loci) in "prep" format ---
locus_order <- c("D8S1179", "D21S11")

db_prep <- data.frame(
  SampleID = c("Sample1","Sample1","Sample2","Sample2","Sample3","Sample3"),
  Locus    = c("D8S1179","D21S11","D8S1179","D21S11","D8S1179","D21S11"),
  Allele1  = c(12, 28, 14, 29, 13, 30),
  Allele2  = c(13, 30, 14, 31, 14, 30),
  stringsAsFactors = FALSE
)
# any/NA を整数 9999 に統一（今回は全て数値なので不要だが一応同じ経路を通す）
db_prep$Allele1 <- map_any(db_prep$Allele1)
db_prep$Allele2 <- map_any(db_prep$Allele2)

db_index <- build_db_index_v1(db_prep, locus_order)

# --- build a query (as data.frame) ---
# 例：D8S1179=(12,14), D21S11=(28,30)
q_raw <- data.frame(
  Locus   = locus_order,
  allele1 = c("12","28"),
  allele2 = c("14","30"),
  stringsAsFactors = FALSE
)
q_prep <- prepare_query_min(q_raw, locus_order)

# --- run & inspect ---
res <- run_match_fast(q_prep, db_index, pre_add_for_any_any = 0L)  # この小テストは事前加点なしで
cat("=== summary ===\n"); print(res$summary)
cat("=== detail ===\n");  print(res$detail)

# quick sanity checks
stopifnot(is.data.frame(res$summary), is.data.frame(res$detail))
stopifnot(all(c("SampleID","Score") %in% names(res$summary)))
stopifnot(all(c("SampleID","Locus","DB_Allele1","DB_Allele2","Score") %in% names(res$detail)))
cat("✅ small fast matcher test passed\n")

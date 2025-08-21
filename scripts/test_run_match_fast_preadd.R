# test_run_match_fast_preadd.R
# Real-DB test: 21 loci, 20 are (any,any) -> pre-add = +2 * 20 = +40.
# Verifies that total score = per-sample score on the single tested locus + pre-add.

# --- core ---
source("scripts/scoring_fast.R")
source("scripts/matcher_fast.R")

# --- helpers (local) ---
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

# --- load real DB & locus order ---
# Expect: data/locus_order.rds, data/database_profile.csv
if (file.exists("data/locus_order.rds")) {
  locus_order <- readRDS("data/locus_order.rds")
} else {
  stop("locus_order.rds not found at data/locus_order.rds")
}
db_prep <- read.csv("data/database_profile.csv", stringsAsFactors = FALSE)

# normalize/limit columns
cn <- names(db_prep)
names(db_prep) <- sub("(?i)^sample[_ ]?id$","SampleID",cn,perl=TRUE)
names(db_prep) <- sub("(?i)^locus$","Locus",names(db_prep),perl=TRUE)
names(db_prep) <- sub("(?i)^allele1$","Allele1",names(db_prep),perl=TRUE)
names(db_prep) <- sub("(?i)^allele2$","Allele2",names(db_prep),perl=TRUE)
db_prep <- db_prep[, c("SampleID","Locus","Allele1","Allele2"), drop=FALSE]
db_prep$Allele1 <- map_any(db_prep$Allele1)
db_prep$Allele2 <- map_any(db_prep$Allele2)
# keep only loci present in locus_order
db_prep <- db_prep[ db_prep$Locus %in% locus_order, , drop = FALSE ]

# build index
db_index <- build_db_index_v1(db_prep, locus_order)

# --- construct query: 1 locus specified, others any/any ---
# choose a target locus that definitely exists in db_prep
present_loci <- intersect(locus_order, unique(db_prep$Locus))
if (length(present_loci) == 0L) stop("No overlapping loci between DB and locus_order")
target_locus <- present_loci[1]

# pick the first non-any genotype in DB for that locus as query genotype
row_t <- which(db_prep$Locus == target_locus)[1]
q_a1 <- db_prep$Allele1[row_t]
q_a2 <- db_prep$Allele2[row_t]
if (is.na(q_a1)) q_a1 <- ANY_CODE
if (is.na(q_a2)) q_a2 <- ANY_CODE

q_raw <- data.frame(
  Locus   = locus_order,
  allele1 = rep("any", length(locus_order)),
  allele2 = rep("any", length(locus_order)),
  stringsAsFactors = FALSE
)
q_raw$allele1[ q_raw$Locus == target_locus ] <- as.character(q_a1)
q_raw$allele2[ q_raw$Locus == target_locus ] <- as.character(q_a2)
q_prep <- prepare_query_min(q_raw, locus_order)

# --- run ---
res <- run_match_fast(q_prep, db_index, pre_add_for_any_any = 2L)
cat("=== Target locus & genotype ===\n")
cat(sprintf("Locus=%s, Query=(%s,%s)\n", target_locus, as.character(q_a1), as.character(q_a2)))
pre_add <- 2L * (length(locus_order) - 1L)
cat(sprintf("Pre-add total = %d (2 x %d loci)\n\n", pre_add, length(locus_order) - 1L))

cat("=== summary (top 10) ===\n"); print(utils::head(res$summary, 10))
cat("=== detail (should have exactly one row per sample, locus=", target_locus, ") ===\n", sep="")
print(utils::head(res$detail, 10))

# --- validate: summary == aggregated detail score + pre_add ---
# only the target locus should appear in detail
stopifnot(all(res$detail$Locus == target_locus))
# aggregate per-sample detail score
agg <- aggregate(Score ~ SampleID, data = res$detail, FUN = sum, drop = FALSE)
# merge with summary
chk <- merge(res$summary, agg, by = "SampleID", all.x = TRUE, sort = FALSE)
# NA (should not happen) -> 0
chk$Score.y[is.na(chk$Score.y)] <- 0L
# expect equality
ok_vec <- (chk$Score.x == (chk$Score.y + pre_add))
if (!all(ok_vec)) {
  bad <- chk[!ok_vec, ]
  print(bad)
  stop("Total score != (detail score + pre-add) for some samples")
}

# rownames reset already handled in matcher; ensure clean
row.names(res$summary) <- NULL
row.names(res$detail)  <- NULL

cat("\n✅ Passed: total score equals (single-locus detail + pre-add=", pre_add, ") for all samples.\n", sep="")

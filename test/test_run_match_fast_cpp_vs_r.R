# test/test_run_match_fast_cpp_vs_r.R
# Purpose:
#   Compare run_match_fast(use_cpp=FALSE) vs run_match_fast(use_cpp=TRUE)
#   using an already-indexed virtual DB:
#     data/virtual_db_u100_S1000_seed123.rds
#   Structure: list(sample_ids, locus_ids, A1, A2)
#
# Prereqs:
#   - scripts/scoring_fast.R      (R implementation: score_2x2 etc.)
#   - scripts/matcher_fast.R      (run_match_fast)
#   - scripts/scoring_cpp.R       (Rcpp loader: compute_scores_uint16)
#
# No multibyte characters in code/comments.

# --- load core scripts ---
source("scripts/scoring_fast.R")
source("scripts/matcher_fast.R")
source("scripts/scoring_cpp.R")  # will compile src/matcher_fast.cpp on first call

ANY_CODE <- 9999L

# --- tiny helpers (self-contained) ---
.map_any <- function(x, any_code = ANY_CODE) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "any"
  x[tolower(x) == "any"] <- as.character(any_code)
  suppressWarnings(as.integer(x))
}

.prepare_query_df_min <- function(csv = "data/query_profile.csv",
                                  locus_order = NULL,
                                  any_code = ANY_CODE) {
  if (!file.exists(csv)) {
    stop(sprintf("Missing query CSV: %s", csv))
  }
  q <- read.csv(csv, stringsAsFactors = FALSE)
  # normalize minimal expected columns (case-insensitive)
  names(q) <- sub("(?i)^locus$","Locus",    names(q), perl=TRUE)
  names(q) <- sub("(?i)^allele1$","Allele1",names(q), perl=TRUE)
  names(q) <- sub("(?i)^allele2$","Allele2",names(q), perl=TRUE)
  req <- c("Locus","Allele1","Allele2")
  if (!all(req %in% names(q))) {
    stop("query CSV must have columns: Locus, Allele1, Allele2")
  }
  # encode to integers (ANY -> any_code)
  q$allele1 <- .map_any(q$Allele1, any_code)
  q$allele2 <- .map_any(q$Allele2, any_code)
  if (!is.null(locus_order) && length(locus_order)) {
    ord <- match(q$Locus, locus_order)
    q <- q[order(ord), , drop = FALSE]
  }
  rownames(q) <- NULL
  q[, c("Locus","allele1","allele2")]
}

# --- load inputs ---
db_index <- readRDS("data/virtual_db_u100_S1000_seed123.rds")
# expected names: sample_ids, locus_ids, A1, A2 (already indexed)
exp <- c("sample_ids","locus_ids","A1","A2")
if (!all(exp %in% names(db_index))) {
  stop("virtual DB must be a list with: sample_ids, locus_ids, A1, A2")
}
# type checks (light)
stopifnot(is.character(db_index$sample_ids))
stopifnot(is.character(db_index$locus_ids))
stopifnot(is.matrix(db_index$A1))
stopifnot(is.matrix(db_index$A2))

locus_order <- db_index$locus_ids
q_prep <- .prepare_query_df_min("data/query_profile.csv", locus_order, any_code = ANY_CODE)

# --- run R & Cpp paths identically ---
res_r <- run_match_fast(q_prep, db_index,
                        pre_add_for_any_any = 2L,
                        any_code = ANY_CODE,
                        include_bits_in_detail = TRUE,
                        use_cpp = FALSE)

res_c <- run_match_fast(q_prep, db_index,
                        pre_add_for_any_any = 2L,
                        any_code = ANY_CODE,
                        include_bits_in_detail = TRUE,
                        use_cpp = TRUE)

# --- compare ---
# 1) summary scores (order-insensitive by SampleID)
r_sum <- res_r$summary[, c("SampleID","Score")]
c_sum <- res_c$summary[, c("SampleID","Score")]
r_sum <- r_sum[order(r_sum$SampleID), , drop = FALSE]
c_sum <- c_sum[order(c_sum$SampleID), , drop = FALSE]
rownames(r_sum) <- NULL; rownames(c_sum) <- NULL

if (!identical(r_sum, c_sum)) {
  # show a small diff and stop
  cat("\n--- Summary diff head (R) ---\n"); print(head(r_sum, 10))
  cat("\n--- Summary diff head (CPP) ---\n"); print(head(c_sum, 10))
  stop("Mismatch in summary scores between R and Cpp paths.")
}

# 2) detail rows (order-insensitive by SampleID,Locus)
keycols <- intersect(c("SampleID","Locus","DB_Allele1","DB_Allele2","Score","Code","Bits"),
                     intersect(names(res_r$detail), names(res_c$detail)))
r_key <- res_r$detail[order(res_r$detail$SampleID, res_r$detail$Locus), keycols, drop=FALSE]
c_key <- res_c$detail[order(res_c$detail$SampleID, res_c$detail$Locus), keycols, drop=FALSE]
rownames(r_key) <- NULL; rownames(c_key) <- NULL

if (!identical(r_key, c_key)) {
  # show rows that differ (limited print)
  diff_idx <- which(rowSums(r_key != c_key) > 0)
  cat(sprintf("\nDetail mismatch rows: %d\n", length(diff_idx)))
  print(utils::head(r_key[diff_idx, , drop=FALSE], 20))
  print(utils::head(c_key[diff_idx, , drop=FALSE], 20))
  stop("Mismatch in detail rows between R and Cpp paths.")
}

cat("\n✅ Passed: run_match_fast(use_cpp=FALSE) == run_match_fast(use_cpp=TRUE) on virtual_db_u100_S1000_seed123.rds\n")

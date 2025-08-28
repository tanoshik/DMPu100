# scripts/matcher_fast.R
# DMPu100 matcher core with blockized C++ option. No multibyte chars.

# Requires scoring core (R) and optional Rcpp wrapper
source("scripts/scoring_fast.R",  local = FALSE)
if (file.exists("scripts/scoring_cpp.R")) {
  source("scripts/scoring_cpp.R", local = FALSE)
}

# Fallbacks
if (!exists("ANY_CODE", mode = "numeric")) ANY_CODE <- 9999L
if (!exists("SCORE_TABLE", mode = "numeric")) {
  # index by code 0..15
  SCORE_TABLE <- as.integer(c(0,1,1,1, 1,1,2,2, 1,2,1,2, 1,2,2,2))
}

# Expect db_index to be:
#   list(
#     sample_ids = <chr length S>,
#     locus_ids  = <chr length L>,
#     A1 = <int matrix [S x L]>,  # any=9999
#     A2 = <int matrix [S x L]>   # any=9999
#   )
.assert_db_index <- function(db_index) {
  ok <- is.list(db_index) &&
    all(c("sample_ids","locus_ids","A1","A2") %in% names(db_index)) &&
    is.character(db_index$sample_ids) &&
    is.character(db_index$locus_ids) &&
    is.matrix(db_index$A1) && is.integer(db_index$A1) &&
    is.matrix(db_index$A2) && is.integer(db_index$A2) &&
    nrow(db_index$A1) == length(db_index$sample_ids) &&
    nrow(db_index$A2) == length(db_index$sample_ids) &&
    ncol(db_index$A1) == length(db_index$locus_ids) &&
    ncol(db_index$A2) == length(db_index$locus_ids)
  if (!ok) stop("Invalid db_index structure")
  invisible(TRUE)
}

# Vectorized R scorer for a block (fallback when use_cpp=FALSE)
.score_block_R <- function(q1, q2, r1, r2, any_code = ANY_CODE) {
  m0 <- (r1 == q1) | (r1 == any_code) | (q1 == any_code)
  m1 <- (r2 == q1) | (r2 == any_code) | (q1 == any_code)
  m2 <- (r1 == q2) | (r1 == any_code) | (q2 == any_code)
  m3 <- (r2 == q2) | (r2 == any_code) | (q2 == any_code)
  code <- as.integer(m0) + bitwShiftL(as.integer(m1), 1L) +
    bitwShiftL(as.integer(m2), 2L) + bitwShiftL(as.integer(m3), 3L)
  SCORE_TABLE[code + 1L]
}

# Main runner
# q_df: data.frame with columns Locus, allele1, allele2 (integers; ANY_CODE allowed)
# db_index: as described above (may be a slice)
# pre_add_for_any_any: add this score and skip comparison for (any,any) locus
# include_bits_in_detail: if TRUE, include per-locus rows in detail
# use_cpp: if TRUE, use Rcpp kernel; else use R fallback
# block_size: number of samples per block for C++/R vector calls
run_match_fast <- function(q_df, db_index,
                           pre_add_for_any_any = 2L,
                           include_bits_in_detail = TRUE,
                           use_cpp = FALSE,
                           block_size = 16000L) {
  .assert_db_index(db_index)
  stopifnot(is.data.frame(q_df), all(c("Locus","allele1","allele2") %in% names(q_df)))
  stopifnot(is.integer(q_df$allele1) || is.numeric(q_df$allele1))
  stopifnot(is.integer(q_df$allele2) || is.numeric(q_df$allele2))
  if (!is.integer(q_df$allele1)) q_df$allele1 <- as.integer(q_df$allele1)
  if (!is.integer(q_df$allele2)) q_df$allele2 <- as.integer(q_df$allele2)
  
  # Align query order to db_index$locus_ids
  LIDs <- db_index$locus_ids
  o <- match(LIDs, as.character(q_df$Locus))
  if (any(is.na(o))) stop("Query does not cover all loci in db_index$locus_ids")
  q1_vec <- q_df$allele1[o]; q2_vec <- q_df$allele2[o]
  
  S <- length(db_index$sample_ids)
  L <- length(LIDs)
  any_code <- if (exists("ANY_CODE")) ANY_CODE else 9999L
  
  # Prepare accumulators
  total_scores <- integer(S); total_scores[] <- 0L
  detail_list <- if (isTRUE(include_bits_in_detail)) vector("list", L) else NULL
  
  # Optionally ensure Rcpp compiled
  if (isTRUE(use_cpp)) {
    if (!exists("score_block_cpp", mode = "function")) {
      if (file.exists("scripts/scoring_cpp.R")) source("scripts/scoring_cpp.R", local = FALSE)
    }
    if (!exists("score_block_cpp", mode = "function")) {
      stop("score_block_cpp() not available (scripts/scoring_cpp.R)")
    }
    ensure_rcpp_compiled(rebuild = FALSE)
  }
  
  # Iterate loci
  for (j in seq_len(L)) {
    q1 <- q1_vec[j]; q2 <- q2_vec[j]
    # Pre-add and skip heavy work for (any,any)
    if (q1 == any_code && q2 == any_code) {
      if (pre_add_for_any_any != 0L) total_scores <- total_scores + as.integer(pre_add_for_any_any)
      # Detail rows for (any,any) are skipped by default to reduce output size.
      next
    }
    
    # Block scan over samples
    i <- 1L
    while (i <= S) {
      b_to <- min(S, i + as.integer(block_size) - 1L)
      r1_blk <- db_index$A1[i:b_to, j]
      r2_blk <- db_index$A2[i:b_to, j]
      
      if (isTRUE(use_cpp)) {
        sc <- score_block_cpp(q1, q2, r1_blk, r2_blk, any_code = any_code)
      } else {
        sc <- .score_block_R(q1, q2, r1_blk, r2_blk, any_code = any_code)
      }
      # accumulate
      total_scores[i:b_to] <- total_scores[i:b_to] + as.integer(sc)
      
      if (isTRUE(include_bits_in_detail)) {
        blk_detail <- data.frame(
          SampleID = db_index$sample_ids[i:b_to],
          Locus    = LIDs[j],
          DB_Allele1 = as.integer(r1_blk),
          DB_Allele2 = as.integer(r2_blk),
          Score = as.integer(sc),
          stringsAsFactors = FALSE
        )
        if (is.null(detail_list[[j]])) {
          detail_list[[j]] <- blk_detail
        } else {
          detail_list[[j]] <- rbind(detail_list[[j]], blk_detail)
        }
      }
      
      i <- b_to + 1L
    }
  }
  
  # Build outputs
  summary <- data.frame(
    SampleID = db_index$sample_ids,
    Score    = as.integer(total_scores),
    stringsAsFactors = FALSE
  )
  # Order: Score desc, SampleID asc
  if (nrow(summary) > 0L) {
    o <- order(-summary$Score, summary$SampleID)
    summary <- summary[o, , drop = FALSE]
    rownames(summary) <- NULL
  }
  
  detail <- if (isTRUE(include_bits_in_detail)) {
    do.call(rbind, detail_list)
  } else {
    data.frame()
  }
  
  used_flag <- if (isTRUE(use_cpp)) "fast_cpp_blk" else "fast_blk"
  list(detail = detail, summary = summary, used = used_flag)
}

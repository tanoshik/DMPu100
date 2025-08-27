# scripts/scoring_cpp.R
# Thin R wrapper to compile/load the Rcpp kernel and provide a simple API.
# No multibyte characters in code/comments.

ensure_rcpp_compiled <- function() {
  # Compile on the fly if not loaded.
  # You can call this in app startup or before running matches.
  src_file <- file.path("src", "matcher_fast.cpp")
  if (!file.exists(src_file)) {
    stop(sprintf("Missing source file: %s", src_file))
  }
  # If the function is already available, skip compilation
  if (!exists("compute_scores_uint16", mode = "function")) {
    message("[Rcpp] Compiling src/matcher_fast.cpp ...")
    Rcpp::sourceCpp(src_file, verbose = FALSE, rebuild = TRUE)
    message("[Rcpp] Done.")
  }
}

# Public helper: score a prepared data.frame by columns.
# Expect integer columns named q1, q2, r1, r2 (uint16-like).
# any_code defaults to 9999 (DMP ANY_CODE).
score_block_cpp <- function(df, any_code = 9999L) {
  stopifnot(is.data.frame(df))
  req <- c("q1", "q2", "r1", "r2")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) {
    stop(sprintf("score_block_cpp(): missing columns: %s", paste(miss, collapse = ", ")))
  }
  # Ensure integer type
  q1 <- as.integer(df$q1); q2 <- as.integer(df$q2)
  r1 <- as.integer(df$r1); r2 <- as.integer(df$r2)
  
  ensure_rcpp_compiled()
  compute_scores_uint16(q1, q2, r1, r2, as.integer(any_code))
}

# Example integration:
# Replace the inner per-row scoring of run_match_fast() with this call
# after you have constructed an intermediate data.frame "blk" that has
# columns q1, q2, r1, r2 (integer-encoded alleles).
#
# blk$Score <- score_block_cpp(blk, any_code = 9999L)
#
# You can then proceed to aggregate or output as before.

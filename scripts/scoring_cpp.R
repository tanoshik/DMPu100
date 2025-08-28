# scripts/scoring_cpp.R
# No multibyte characters.

ensure_rcpp_compiled <- function(rebuild = FALSE) {
  src_file <- file.path("src", "matcher_fast.cpp")
  if (!file.exists(src_file)) stop(sprintf("Missing source file: %s", src_file))
  
  need <- isTRUE(rebuild) || !exists("compute_scores_uint16", mode = "function")
  
  # sanity probe: even if the name exists, the DLL may be unloaded
  if (!need) {
    probe_ok <- FALSE
    try({
      invisible(compute_scores_uint16(
        q1 = 9999L, q2 = 9999L,
        r1 = as.integer(9999L),
        r2 = as.integer(9999L),
        any_code = 9999L
      ))
      probe_ok <- TRUE
    }, silent = TRUE)
    if (!probe_ok) need <- TRUE
  }
  
  if (need) {
    message("[Rcpp] Compiling src/matcher_fast.cpp ...")
    Rcpp::sourceCpp(src_file, verbose = FALSE, rebuild = TRUE, cacheDir = "src/.rcpp_cache")
    message("[Rcpp] Done.")
  }
}

score_block_cpp <- function(q1, q2, r1, r2, any_code = 9999L) {
  q1 <- as.integer(q1); q2 <- as.integer(q2)
  r1 <- as.integer(r1); r2 <- as.integer(r2)
  ensure_rcpp_compiled(rebuild = FALSE)
  compute_scores_uint16(q1, q2, r1, r2, as.integer(any_code))
}

# scripts/scoring_cpp.R
# No multibyte characters.

ensure_rcpp_compiled <- function(rebuild = FALSE) {
  src_file <- file.path("src", "matcher_fast.cpp")
  if (!file.exists(src_file)) stop(sprintf("Missing source file: %s", src_file))
  
  need <- isTRUE(rebuild) || !exists("compute_scores_uint16", mode = "function")
  
  # すでに関数が見えていないなら、まずは静かにキャッシュロードを試す（rebuild=FALSE）
  if (need && !isTRUE(rebuild)) {
    try({
      Rcpp::sourceCpp(src_file, verbose = FALSE, rebuild = FALSE, cacheDir = "src/.rcpp_cache")
    }, silent = TRUE)
    # ここでロードできたら need を再評価
    need <- !exists("compute_scores_uint16", mode = "function")
  }
  
  # まだ無ければ初回コンパイル（このときのみメッセージを出す）
  if (need) {
    message("[Rcpp] Compiling src/matcher_fast.cpp ...")
    Rcpp::sourceCpp(src_file, verbose = FALSE, rebuild = TRUE, cacheDir = "src/.rcpp_cache")
    message("[Rcpp] Done.")
  }
  
  # プローブ（DLLがUnloadされていた等の保険）
  if (!exists("compute_scores_uint16", mode = "function")) {
    stop("Failed to load/compile compute_scores_uint16")
  }
  invisible(TRUE)
}

score_block_cpp <- function(q1, q2, r1, r2, any_code = 9999L) {
  q1 <- as.integer(q1); q2 <- as.integer(q2)
  r1 <- as.integer(r1); r2 <- as.integer(r2)
  ensure_rcpp_compiled(rebuild = FALSE)
  compute_scores_uint16(q1, q2, r1, r2, as.integer(any_code))
}

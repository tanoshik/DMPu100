# No multibyte chars.
# Make 1M subset from the 2M index RDS (list: sample_ids, locus_ids, A1, A2)
make_subset_from_2M <- function(in_path = "data/virtual_db_u100_S2000000_seed123.rds",
                                out_path = "data/virtual_db_u100_S1000000_seed123.rds",
                                S_out = 1000000L,
                                seed = 123L) {
  stopifnot(file.exists(in_path))
  obj <- readRDS(in_path)
  stopifnot(is.list(obj), all(c("sample_ids","locus_ids","A1","A2") %in% names(obj)))
  S <- length(obj$sample_ids)
  if (S_out > S) stop("S_out > S in source DB")
  
  set.seed(seed)
  idx <- sort(sample.int(S, S_out))
  sub <- list(
    sample_ids = obj$sample_ids[idx],
    locus_ids  = obj$locus_ids,
    A1 = obj$A1[idx, , drop = FALSE],
    A2 = obj$A2[idx, , drop = FALSE]
  )
  saveRDS(sub, out_path, compress = "gzip")
  message("[SUBSET] written: ", out_path)
  invisible(out_path)
}

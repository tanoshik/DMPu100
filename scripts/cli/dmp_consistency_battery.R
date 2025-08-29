# scripts/cli/dmp_consistency_battery.R
# No multibyte characters.

dmp_consistency_battery <- function(
    db_rds,
    query_csv = "data/query_profile_seed123.csv",
    rows      = 1000L,        # quick size
    seed      = 123L,
    block_size = 16000L,
    out_dir   = "test/consistency",
    emit_index_rds = FALSE
) {
  ANY_CODE <- 9999L
  `%||%` <- function(x,y) if (is.null(x)) y else x
  .ensure_dir <- function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .ensure_dir(out_dir)
  
  message("[CONSIST] load deps")
  source("scripts/scoring_fast.R",  local = TRUE)
  source("scripts/matcher_fast.R",  local = TRUE)
  
  # load DB index or build from CSV-like RDS
  message("[CONSIST] load RDS: ", db_rds)
  load_locus_order <- function() {
    if (file.exists("data/locus_order.rds")) {
      lo <- tryCatch(readRDS("data/locus_order.rds"), error=function(e) NULL)
      if (is.character(lo)) return(lo)
    }
    NULL
  }
  .map_any <- function(x, any_code=ANY_CODE) {
    if (is.null(x)) return(integer(0))
    if (is.factor(x)) x <- as.character(x)
    out <- suppressWarnings(as.integer(x))
    is_any <- !is.na(x) & tolower(as.character(x))=="any"
    out[is_any] <- any_code
    out
  }
  .build_db_index_v1 <- function(db_prep, locus_order, any_code=ANY_CODE) {
    req <- c("SampleID","Locus","Allele1","Allele2")
    stopifnot(is.data.frame(db_prep), all(req %in% names(db_prep)))
    SIDs <- unique(db_prep$SampleID); LIDs <- as.character(locus_order)
    S <- length(SIDs); L <- length(LIDs)
    A1 <- matrix(any_code, nrow=S, ncol=L, dimnames=list(NULL,LIDs))
    A2 <- matrix(any_code, nrow=S, ncol=L, dimnames=list(NULL,LIDs))
    Lmap <- setNames(seq_along(LIDs), LIDs); sid_map <- setNames(seq_along(SIDs), SIDs)
    for (k in seq_len(nrow(db_prep))) {
      sid <- db_prep$SampleID[k]; loc <- as.character(db_prep$Locus[k])
      i <- sid_map[[sid]]; j <- Lmap[[loc]]; if (is.na(i) || is.na(j)) next
      A1[i,j] <- .map_any(db_prep$Allele1[k], any_code)
      A2[i,j] <- .map_any(db_prep$Allele2[k], any_code)
    }
    list(sample_ids=SIDs, locus_ids=LIDs, A1=A1, A2=A2)
  }
  load_rds_as_index <- function(path) {
    obj <- readRDS(path); lo <- load_locus_order()
    if (is.list(obj) && all(c("sample_ids","locus_ids","A1","A2") %in% names(obj))) {
      return(list(
        db_index=obj,
        make_query=function() {
          data.frame(Locus=obj$locus_ids, allele1=obj$A1[1,,drop=TRUE], allele2=obj$A2[1,,drop=TRUE], stringsAsFactors=FALSE)
        }
      ))
    }
    if (is.data.frame(obj)) {
      n <- names(obj)
      names(obj) <- sub("(?i)^sampleid$","SampleID", n, perl=TRUE)
      names(obj) <- sub("(?i)^locus$","Locus", names(obj), perl=TRUE)
      names(obj) <- sub("(?i)^allele1$","Allele1", names(obj), perl=TRUE)
      names(obj) <- sub("(?i)^allele2$","Allele2", names(obj), perl=TRUE)
      LIDs <- lo; if (is.null(LIDs)) LIDs <- unique(as.character(obj$Locus))
      idx <- .build_db_index_v1(obj, LIDs, ANY_CODE)
      return(list(
        db_index=idx,
        make_query=function() {
          sid <- obj$SampleID[1]
          q <- subset(obj, SampleID==sid, c("Locus","Allele1","Allele2"))
          names(q) <- c("Locus","allele1","allele2"); q
        }
      ))
    }
    stop("Unsupported RDS")
  }
  
  info <- load_rds_as_index(db_rds)
  db_index <- info$db_index
  S <- length(db_index$sample_ids); L <- length(db_index$locus_ids)
  
  # subset rows for quick mode
  set.seed(seed)
  pick <- sort(sample.int(S, size = min(rows, S)))
  slice <- list(
    sample_ids = db_index$sample_ids[pick],
    locus_ids  = db_index$locus_ids,
    A1 = db_index$A1[pick, , drop=FALSE],
    A2 = db_index$A2[pick, , drop=FALSE]
  )
  
  # make query
  .prepare_query_min <- function(q_raw, locus_order=NULL, any_code=ANY_CODE) {
    n <- names(q_raw)
    names(q_raw) <- sub("(?i)^locus$","Locus", n, perl=TRUE)
    names(q_raw) <- sub("(?i)^allele1$","allele1", names(q_raw), perl=TRUE)
    names(q_raw) <- sub("(?i)^allele2$","allele2", names(q_raw), perl=TRUE)
    stopifnot(all(c("Locus","allele1","allele2") %in% names(q_raw)))
    q <- q_raw[,c("Locus","allele1","allele2"),drop=FALSE]
    q$Locus <- as.character(q$Locus)
    q$allele1 <- .map_any(q$allele1, ANY_CODE)
    q$allele2 <- .map_any(q$allele2, ANY_CODE)
    if (!is.null(locus_order) && length(locus_order)>0L) {
      ord <- match(q$Locus, locus_order); q <- q[order(ord), , drop=FALSE]
    }
    rownames(q) <- NULL; q
  }
  if (file.exists(query_csv)) {
    q_raw <- read.csv(query_csv, stringsAsFactors = FALSE)
    q <- .prepare_query_min(q_raw, db_index$locus_ids, ANY_CODE)
  } else {
    q <- info$make_query()
  }
  
  # run R
  message("[CONSIST] run R engine ...")
  rs_R <- run_match_fast(q, slice,
                         include_bits_in_detail = FALSE,
                         use_cpp = FALSE,
                         block_size = as.integer(block_size))
  # run Rcpp
  message("[CONSIST] run Rcpp engine ...")
  suppressWarnings(suppressMessages({
    if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp not installed")
  }))
  # compile once in master
  Rcpp::sourceCpp("src/matcher_fast.cpp", rebuild = TRUE, cacheDir = "src/.rcpp_cache")
  rs_C <- run_match_fast(q, slice,
                         include_bits_in_detail = FALSE,
                         use_cpp = TRUE,
                         block_size = as.integer(block_size))
  
  # compare summary scores
  sm_R <- rs_R$summary; sm_C <- rs_C$summary
  names(sm_R)[names(sm_R)=="Score"] <- "Score_R"
  names(sm_C)[names(sm_C)=="Score"] <- "Score_C"
  key <- intersect(names(sm_R), names(sm_C)); key <- key[key %in% c("SampleID","Locus")]
  if (!"SampleID" %in% key) key <- "SampleID"
  
  merged <- merge(sm_R[, c("SampleID","Score_R")], sm_C[, c("SampleID","Score_C")], by="SampleID", all=TRUE)
  merged$equal <- merged$Score_R == merged$Score_C
  ok <- all(merged$equal %||% TRUE)
  
  out_tag <- paste0("CONSIST_", format(Sys.time(), "%Y%m%d%H%M%S"))
  csv_R <- file.path(out_dir, paste0("bench_run_match_fast__CONSIST_R_", out_tag, ".csv"))
  csv_C <- file.path(out_dir, paste0("bench_run_match_fast__CONSIST_CPP_", out_tag, ".csv"))
  csv_cmp <- file.path(out_dir, paste0("bench_run_match_fast__CONSIST_DIFF_", out_tag, ".csv"))
  
  write.csv(sm_R, csv_R, row.names = FALSE)
  write.csv(sm_C, csv_C, row.names = FALSE)
  write.csv(merged, csv_cmp, row.names = FALSE)
  
  message("[CONSIST] equal = ", ok, "  rows=", nrow(merged))
  if (emit_index_rds) {
    saveRDS(slice, file.path(out_dir, paste0("db_index_slice_", out_tag, ".rds")))
  }
  invisible(list(equal = ok, diff = merged, R=rs_R, C=rs_C,
                 files = list(R=csv_R, C=csv_C, cmp=csv_cmp)))
}

# CLI entry
if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly = TRUE)
  kv <- as.list(setNames(rep(NA_character_, length(args)), rep("", length(args))))
  for (a in args) {
    if (grepl("^--", a)) {
      s <- sub("^--", "", a)
      k <- sub("=.*$", "", s); v <- sub("^[^=]+=", "", s)
      kv[[k]] <- v
    }
  }
  get <- function(k, def=NULL) {
    v <- kv[[k]]
    if (is.na(v) || is.null(v)) def else v
  }
  db_rds <- get("db", "data/virtual_db_u100_S100000_seed123.rds")
  rows   <- as.integer(get("rows", "1000"))
  seed   <- as.integer(get("seed", "123"))
  blk    <- as.integer(get("block_size", "16000"))
  outd   <- get("out", "test/consistency")
  emit   <- tolower(get("emit_index_rds","false")) %in% c("1","true","yes","on")
  
  res <- dmp_consistency_battery(
    db_rds = db_rds, rows = rows, seed = seed,
    block_size = blk, out_dir = outd, emit_index_rds = emit
  )
  cat("[CONSIST][DONE] equal=", res$equal, "  out_dir=", outd, "\n", sep="")
}

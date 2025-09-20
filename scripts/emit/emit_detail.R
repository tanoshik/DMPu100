# scripts/emit/emit_detail.R
# ASCII-only; JP comments allowed for dev.

suppressPackageStartupMessages({
  library(jsonlite)
})

# ---- public I/F ----
# emit_detail(opt_path, scores_path, out_dir, mode=c("raw","detail"), any_code=9999L)
emit_detail <- function(opt_path, scores_path, out_dir, mode=c("raw","detail"), any_code=9999L) {
  mode <- match.arg(mode)
  if (!file.exists(opt_path))   stop(sprintf("opt_path not found: %s", opt_path))
  if (!file.exists(scores_path)) stop(sprintf("scores_path not found: %s", scores_path))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- load opts ----
  opts <- jsonlite::fromJSON(opt_path, simplifyVector = TRUE)
  if (is.null(opts$db_path) || is.null(opts$query_path)) {
    stop("opts json must contain db_path and query_path")
  }
  db_path <- opts$db_path
  q_path  <- opts$query_path
  if (!file.exists(db_path)) stop(sprintf("db_path not found: %s", db_path))
  if (!file.exists(q_path))  stop(sprintf("query_path not found: %s", q_path))
  
  # ---- load DB ----
  db <- readRDS(db_path) # list(sample_ids, locus_ids, A1, A2) with integers (×100), ANY=9999
  sample_ids <- as.character(db$sample_ids)
  locus_ids  <- as.character(db$locus_ids)
  S <- length(sample_ids); L <- length(locus_ids)
  A1 <- as.matrix(db$A1); storage.mode(A1) <- "integer"
  A2 <- as.matrix(db$A2); storage.mode(A2) <- "integer"
  if (nrow(A1)!=S || ncol(A1)!=L || nrow(A2)!=S || ncol(A2)!=L) stop("DB A1/A2 shape mismatch")
  
  # ---- load query (RDS or CSV), align to locus_ids, standardize (ascending + ANY right) ----
  ANY <- as.integer(any_code)
  to_pair <- function(a1, a2) {
    # place ANY to right; otherwise ascending
    swap <- ((a1 == ANY & a2 != ANY) | (a1 > a2 & a2 != ANY))
    a1n <- ifelse(swap, a2, a1)
    a2n <- ifelse(swap, a1, a2)
    list(a1=a1n, a2=a2n)
  }
  
  if (grepl("\\.(rds|RDS)$", q_path)) {
    qx <- readRDS(q_path)
    q1 <- as.integer(qx$q1); q2 <- as.integer(qx$q2)
    if (length(q1)!=L || length(q2)!=L) stop("query length must match L")
    qp <- to_pair(q1, q2); q1 <- qp$a1; q2 <- qp$a2
  } else {
    qcsv <- read.csv(q_path, stringsAsFactors = FALSE, check.names = FALSE)
    nms <- tolower(names(qcsv)); names(qcsv) <- nms
    lc <- which(nms %in% c("locus","marker"))[1]
    a1c <- which(nms %in% c("allele1","a1","q1"))[1]
    a2c <- which(nms %in% c("allele2","a2","q2"))[1]
    if (is.na(lc) || is.na(a1c) || is.na(a2c)) stop("query CSV needs Locus & A1 & A2")
    
    mm <- match(locus_ids, qcsv[[lc]])
    if (any(is.na(mm))) stop("query CSV does not cover all loci")
    
    # strict: numeric-only; NO lenient mapping here
    as_i100 <- function(v) {
      num <- suppressWarnings(as.numeric(v))
      if (any(is.na(num))) stop("query CSV contains non-numeric allele(s) under strict mode")
      as.integer(round(num*100))
    }
    qp <- to_pair(as_i100(qcsv[[a1c]][mm]), as_i100(qcsv[[a2c]][mm]))
    q1 <- qp$a1; q2 <- qp$a2
  }
  
  # ---- select target SampleID(s) from scores_path ----
  scids <- read.csv(scores_path, stringsAsFactors = FALSE, check.names = FALSE)
  # 「候補が複数あれば左端」を踏襲
  id_col <- which(tolower(names(scids)) %in% c("sampleid","id","sid"))[1]
  if (is.na(id_col)) stop("scores_path needs a column containing SampleID(s)")
  ids <- as.character(scids[[id_col]])
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) < 1L) stop("no SampleID in scores_path")
  mmS <- match(ids, sample_ids)
  if (any(is.na(mmS))) stop("some SampleID in scores_path not found in DB")
  
  # ---- SCORE_TABLE（0832固定）----
  SCORE_TABLE <- as.integer(c(
    0,1,1,1,
    1,1,2,2,
    1,2,1,2,
    1,2,2,2
  ))
  
  # ---- ANY match fn ----
  m <- function(x, y) { (x == ANY) | (y == ANY) | (x == y) }
  
  # ---- emit ----
  if (mode == "raw") {
    outp <- file.path(out_dir, "raw_detail.csv")
    # preallocate vectors (append per locus)
    out_list <- vector("list", length(ids))
    for (k in seq_along(ids)) {
      sidx <- mmS[k]
      r1 <- A1[sidx, ]
      r2 <- A2[sidx, ]
      # per-locus bits / code / score
      b0 <- as.integer(m(q1, r1))
      b1 <- as.integer(m(q1, r2))
      b2 <- as.integer(m(q2, r1))
      b3 <- as.integer(m(q2, r2))
      code <- b0 + 2L*b1 + 4L*b2 + 8L*b3
      # 4-bit binary string, MSB=b3 .. LSB=b0  -> "b3b2b1b0"
      bits <- paste0(b3, b2, b1, b0)
      sc   <- SCORE_TABLE[code + 1L]
      
      df <- data.frame(
        Locus = locus_ids,
        q1 = q1,
        q2 = q2,
        r1 = r1,
        r2 = r2,
        bits = bits,
        code = code,
        score = sc,                # ← sanity_topN に合わせて lower-case
        SampleID = rep(ids[k], L),
        stringsAsFactors = FALSE
      )
      out_list[[k]] <- df
    }
    out_df <- do.call(rbind, out_list)
    # 出力はローカス順・SampleID順のまま（必要ならここで並べ替え可）
    write.csv(out_df, outp, row.names = FALSE)
    return(invisible(outp))
  } else {
    stop("mode='detail' is not implemented in Phase1")
  }
}

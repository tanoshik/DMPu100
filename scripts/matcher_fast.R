# scripts/matcher_fast.R
# DMPu100 matcher core. No multibyte chars.

# Requires scoring core
source("scripts/scoring_fast.R", local = FALSE)

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
    is.matrix(db_index$A1) && storage.mode(db_index$A1) %in% c("integer","double") &&
    is.matrix(db_index$A2) && storage.mode(db_index$A2) %in% c("integer","double") &&
    nrow(db_index$A1) == nrow(db_index$A2) &&
    ncol(db_index$A1) == ncol(db_index$A2) &&
    length(db_index$sample_ids) == nrow(db_index$A1) &&
    length(db_index$locus_ids)  == ncol(db_index$A1)
  if (!ok) stop("db_index must be list(sample_ids, locus_ids, A1[int], A2[int]) with matching dims")
  invisible(TRUE)
}

# Query must be data.frame(Locus, allele1, allele2), integer alleles (any=9999)
.assert_query_df <- function(q) {
  if (!is.data.frame(q) || !all(c("Locus","allele1","allele2") %in% names(q)))
    stop("q_prep must be data.frame(Locus, allele1, allele2)")
  invisible(TRUE)
}

# locus name -> db_index column
.match_locus_to_col <- function(locus_ids) {
  idx <- seq_along(locus_ids)
  names(idx) <- locus_ids
  idx
}

# Main matcher (fast).
# Args:
#   q_prep: data.frame(Locus, allele1, allele2) normalized, integer alleles, any=9999
#   db_index: list(sample_ids, locus_ids, A1, A2) as asserted above
#   pre_add_for_any_any: integer, default 2 (u100 spec)
#   any_code: integer, default ANY_CODE (9999)
#   include_bits_in_detail: logical, include Bits/Code columns in detail
# Returns:
#   list(detail=<df>, summary=<df>, used="fast")
run_match_fast <- function(q_prep,
                           db_index,
                           pre_add_for_any_any = 2L,
                           any_code = ANY_CODE,
                           include_bits_in_detail = TRUE) {
  .assert_db_index(db_index)
  .assert_query_df(q_prep)
  
  # Query-side skip (+2) for (any,any)
  q <- q_prep
  is_any_any <- (q$allele1 == any_code & q$allele2 == any_code)
  q$skip <- is_any_any
  q$preadd <- ifelse(is_any_any, as.integer(pre_add_for_any_any), 0L)
  pre_add_total <- sum(q$preadd)
  
  # Map query loci to db_index columns; keep only existing loci
  Lmap <- .match_locus_to_col(db_index$locus_ids)
  keep <- q$Locus %in% names(Lmap)
  if (!all(keep)) {
    # drop unknown loci on DB side
    q <- q[keep, , drop = FALSE]
    row.names(q) <- NULL
  }
  if (nrow(q) == 0L) {
    return(list(
      detail  = data.frame(SampleID=character(0), Locus=character(0),
                           DB_Allele1=integer(0), DB_Allele2=integer(0), Score=integer(0),
                           Bits=character(0), Code=integer(0), stringsAsFactors = FALSE),
      summary = data.frame(SampleID=character(0), Score=integer(0), stringsAsFactors = FALSE),
      used = "fast"
    ))
  }
  
  # Pre-allocate holders (list of data.frame rows) —— ここが重要：行ごとに data.frame を積む
  detail_rows  <- vector("list", length = 0L)
  summary_rows <- vector("list", length = length(db_index$sample_ids))
  
  # Iterate DB samples
  S <- length(db_index$sample_ids)
  for (sidx in seq_len(S)) {
    tot <- 0L
    
    # Iterate query loci
    for (k in seq_len(nrow(q))) {
      if (q$skip[k]) {
        next  # (any,any) locus: skipped and already pre-added
      }
      j <- Lmap[[ q$Locus[k] ]]
      if (is.null(j) || is.na(j)) next
      
      q1 <- as.integer(q$allele1[k])
      q2 <- as.integer(q$allele2[k])
      r1 <- as.integer(db_index$A1[sidx, j])
      r2 <- as.integer(db_index$A2[sidx, j])
      
      # score 2x2
      sp <- score_2x2(q1, q2, r1, r2)
      sc <- as.integer(sp$score)
      tot <- tot + sc
      
      # append row as **data.frame** (atomic) to avoid list-columns
      if (isTRUE(include_bits_in_detail)) {
        detail_rows[[length(detail_rows)+1L]] <- data.frame(
          SampleID   = as.character(db_index$sample_ids[sidx]),
          Locus      = as.character(q$Locus[k]),
          DB_Allele1 = r1,
          DB_Allele2 = r2,
          Score      = sc,
          Bits       = as.character(sp$bits0123),
          Code       = as.integer(sp$code),
          stringsAsFactors = FALSE
        )
      } else {
        detail_rows[[length(detail_rows)+1L]] <- data.frame(
          SampleID   = as.character(db_index$sample_ids[sidx]),
          Locus      = as.character(q$Locus[k]),
          DB_Allele1 = r1,
          DB_Allele2 = r2,
          Score      = sc,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # summary row as data.frame（こちらも atomic で積む）
    summary_rows[[sidx]] <- data.frame(
      SampleID = as.character(db_index$sample_ids[sidx]),
      Score    = as.integer(tot + pre_add_total),
      stringsAsFactors = FALSE
    )
  }
  
  # Bind to data.frames（list of data.frame → data.frame, list列にならない）
  detail <- if (length(detail_rows)) {
    do.call(rbind, detail_rows)
  } else {
    data.frame(SampleID=character(0), Locus=character(0),
               DB_Allele1=integer(0), DB_Allele2=integer(0), Score=integer(0),
               Bits=character(0), Code=integer(0), stringsAsFactors = FALSE)
  }
  summary <- do.call(rbind, summary_rows)
  
  # Types and ordering
  detail$SampleID <- as.character(detail$SampleID)
  detail$Locus    <- as.character(detail$Locus)
  suppressWarnings(detail$DB_Allele1 <- as.integer(detail$DB_Allele1))
  suppressWarnings(detail$DB_Allele2 <- as.integer(detail$DB_Allele2))
  suppressWarnings(detail$Score      <- as.integer(detail$Score))
  if ("Code" %in% names(detail)) suppressWarnings(detail$Code <- as.integer(detail$Code))
  summary$SampleID <- as.character(summary$SampleID)
  suppressWarnings(summary$Score <- as.integer(summary$Score))
  
  if (nrow(summary) > 0L) {
    summary <- summary[order(-summary$Score, summary$SampleID), , drop = FALSE]
  }
  row.names(summary) <- NULL
  row.names(detail)  <- NULL
  
  list(detail = detail, summary = summary, used = "fast")
}

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
    is.character(db_index$locus_ids)  &&
    is.matrix(db_index$A1) && is.integer(db_index$A1) &&
    is.matrix(db_index$A2) && is.integer(db_index$A2) &&
    identical(dim(db_index$A1), dim(db_index$A2))
  if (!ok) stop("db_index must be list(sample_ids, locus_ids, A1[int], A2[int]) with matching dims")
}

# Expect q_prep to be a data.frame:
#   columns: Locus (chr), allele1 (int), allele2 (int)
#   values: ANY_CODE(9999) used for any, no NAs for alleles
.assert_query_df <- function(q_prep) {
  req <- c("Locus","allele1","allele2")
  if (!is.data.frame(q_prep) || !all(req %in% names(q_prep)))
    stop("q_prep must be data.frame with Locus, allele1, allele2")
  if (!is.integer(q_prep$allele1) || !is.integer(q_prep$allele2))
    stop("q_prep alleles must be integer (use prepare_query_* before)")
}

# Build a simple integer lookup for locus -> column index
.build_loc_col_map <- function(locus_ids) {
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
  loc_map <- .build_loc_col_map(db_index$locus_ids)
  q$col_j <- as.integer(unname(loc_map[q$Locus]))
  keep <- !is.na(q$col_j)
  if (!all(keep)) q <- q[keep, , drop = FALSE]
  if (nrow(q) == 0L) {
    # No comparable loci -> everyone gets just pre_add_total
    summary <- data.frame(SampleID = db_index$sample_ids,
                          Score = rep.int(pre_add_total, length(db_index$sample_ids)),
                          stringsAsFactors = FALSE)
    return(list(detail = data.frame(SampleID=character(0), Locus=character(0),
                                    DB_Allele1=integer(0), DB_Allele2=integer(0),
                                    Score=integer(0), Bits=character(0), Code=integer(0),
                                    stringsAsFactors = FALSE),
                summary = summary[order(-summary$Score, summary$SampleID), , drop = FALSE],
                used = "fast"))
  }
  
  S <- length(db_index$sample_ids)
  Lsel <- which(!q$skip)        # indices in q (rows) that we actually compare
  detail_rows  <- vector("list", 0L)
  summary_rows <- vector("list", S)
  
  A1 <- db_index$A1
  A2 <- db_index$A2
  
  # Loop over samples
  for (sidx in seq_len(S)) {
    tot <- 0L
    
    # Compare only non-skip loci
    if (length(Lsel) > 0L) {
      for (k in Lsel) {
        j <- q$col_j[k]
        q1 <- q$allele1[k]; q2 <- q$allele2[k]
        r1 <- A1[sidx, j];  r2 <- A2[sidx, j]
        if (is.na(r1)) r1 <- any_code
        if (is.na(r2)) r2 <- any_code
        
        sp <- score_2x2(q1, q2, r1, r2)
        sc <- sp$score
        tot <- tot + sc
        
        if (include_bits_in_detail) {
          detail_rows[[length(detail_rows)+1L]] <- list(
            SampleID = db_index$sample_ids[sidx],
            Locus    = q$Locus[k],
            DB_Allele1 = r1,
            DB_Allele2 = r2,
            Score      = as.integer(sc),
            Bits       = sp$bits0123,
            Code       = as.integer(sp$code)
          )
        } else {
          detail_rows[[length(detail_rows)+1L]] <- list(
            SampleID = db_index$sample_ids[sidx],
            Locus    = q$Locus[k],
            DB_Allele1 = r1,
            DB_Allele2 = r2,
            Score      = as.integer(sc)
          )
        }
      }
    }
    
    summary_rows[[sidx]] <- list(
      SampleID = db_index$sample_ids[sidx],
      Score    = as.integer(tot + pre_add_total)
    )
  }
  
  # Assemble data.frames
  detail <- if (length(detail_rows)) {
    as.data.frame(do.call(rbind, detail_rows), stringsAsFactors = FALSE)
  } else {
    data.frame(SampleID=character(0), Locus=character(0),
               DB_Allele1=integer(0), DB_Allele2=integer(0), Score=integer(0),
               Bits=character(0), Code=integer(0), stringsAsFactors = FALSE)
  }
  summary <- as.data.frame(do.call(rbind, summary_rows), stringsAsFactors = FALSE)
  
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

# scripts/matcher_detail.R
# No multibyte characters.

SCORE_TABLE <- c(
  0, # 0000 (code 0)
  1, # 0001 (1)
  1, # 0010 (2)
  1, # 0011 (3)
  1, # 0100 (4)
  1, # 0101 (5)
  2, # 0110 (6)
  2, # 0111 (7)
  1, # 1000 (8)
  2, # 1001 (9)
  1, # 1010 (10)
  2, # 1011 (11)
  1, # 1100 (12)
  2, # 1101 (13)
  2, # 1110 (14)
  2  # 1111 (15)
)

# decode integer allele to human-readable (e.g., 1320 -> "13.2", 1500 -> "15")
decode_allele <- function(a, ANY) {
  if (is.na(a)) return("")
  if (as.character(a) == as.character(ANY)) return("any")
  v <- as.numeric(a) / 100
  if (abs(v - round(v)) < 1e-9) as.character(as.integer(round(v))) else format(v, trim = TRUE, nsmall = 1)
}

# fetch sample's r1/r2 across loci as numeric vector (length L)
.fetch_sample_locus_vecs <- function(A, s) {
  if (is.matrix(A)) {
    A[s, ]
  } else if (is.list(A)) {
    vapply(A, `[`, numeric(1), s)
  } else {
    stop("DB A1/A2 must be matrix or list")
  }
}

# one-locus scoring bits with ANY wildcard
.bits_for_locus <- function(q1, q2, r1, r2, ANY) {
  # ANY matches everything
  m <- function(x, y) (x == ANY) || (y == ANY) || (x == y)
  b0 <- as.integer(m(q1, r1))
  b1 <- as.integer(m(q1, r2))
  b2 <- as.integer(m(q2, r1))
  b3 <- as.integer(m(q2, r2))
  code <- b0 + 2L * b1 + 4L * b2 + 8L * b3
  list(
    b0 = b0, b1 = b1, b2 = b2, b3 = b3,
    code = code,
    score = SCORE_TABLE[code + 1L]
  )
}

# TopN detail builder (O(TopN * L))
build_detail_for_ids <- function(q_df, db, sample_ids,
                                 any_code = 9999L,
                                 decode = TRUE,
                                 include_code = TRUE,
                                 include_bits = TRUE) {
  if (is.null(db$A1) || is.null(db$A2)) stop("db index must have A1/A2")
  if (is.null(db$locus_ids)) {
    # fall back: infer L from columns or length
    L <- if (is.matrix(db$A1)) ncol(db$A1) else length(db$A1)
    locus_ids <- q_df$Locus
  } else {
    locus_ids <- db$locus_ids
    L <- length(locus_ids)
  }
  if (length(sample_ids) < 1L) return(data.frame())
  
  # map SampleID -> row index
  if (is.null(db$sample_ids)) stop("db index missing sample_ids")
  ridx <- match(sample_ids, db$sample_ids)
  if (anyNA(ridx)) {
    miss <- sample_ids[is.na(ridx)]
    stop(sprintf("unknown SampleID(s) in DB: %s", paste(miss, collapse = ",")))
  }
  
  out_list <- vector("list", length(sample_ids))
  for (k in seq_along(sample_ids)) {
    s <- ridx[k]
    r1v <- .fetch_sample_locus_vecs(db$A1, s)
    r2v <- .fetch_sample_locus_vecs(db$A2, s)
    
    # ensure query is aligned to locus_ids
    if (!is.null(locus_ids)) {
      ord <- match(locus_ids, q_df$Locus)
      if (anyNA(ord)) stop("Query loci do not cover DB locus_ids")
      q1v <- as.integer(q_df$allele1[ord])
      q2v <- as.integer(q_df$allele2[ord])
      locv <- locus_ids
    } else {
      q1v <- as.integer(q_df$allele1)
      q2v <- as.integer(q_df$allele2)
      locv <- q_df$Locus
    }
    
    # per-locus bits/score
    b0 <- integer(L); b1 <- integer(L); b2 <- integer(L); b3 <- integer(L)
    code <- integer(L); sc <- integer(L)
    for (i in seq_len(L)) {
      z <- .bits_for_locus(q1v[i], q2v[i], r1v[i], r2v[i], any_code)
      b0[i] <- z$b0; b1[i] <- z$b1; b2[i] <- z$b2; b3[i] <- z$b3
      code[i] <- z$code; sc[i] <- z$score
    }
    
    # assemble dataframe for this sample
    df <- data.frame(
      Locus = locv,
      q1 = q1v, q2 = q2v, r1 = r1v, r2 = r2v,
      stringsAsFactors = FALSE
    )
    if (include_bits) {
      df$bits <- paste0(b3, b2, b1, b0) # [b3 b2 b1 b0]
    }
    if (include_code) {
      df$code <- code
    }
    df$score <- sc
    df$SampleID <- sample_ids[k]
    
    if (decode) {
      df$q1 <- vapply(df$q1, decode_allele, "", ANY = any_code)
      df$q2 <- vapply(df$q2, decode_allele, "", ANY = any_code)
      df$r1 <- vapply(df$r1, decode_allele, "", ANY = any_code)
      df$r2 <- vapply(df$r2, decode_allele, "", ANY = any_code)
    }
    
    out_list[[k]] <- df
  }
  
  do.call(rbind, out_list)
}

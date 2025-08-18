# scripts/utils_profile.R
# Utilities for preparing STR profiles for matching.
# - No multibyte characters in code/comments.
# - ANY_CODE = 9999 for "any" wildcard.
# - Allele pair rule: numeric ascending, ANY_CODE always on the right.
# - Locus ordering is enforced by locus_order.
# - Homozygous -> any conversion is optional and independent per caller.

# ---------- constants ----------

ANY_CODE <- 9999L

# ---------- helpers: type/codec ----------

is_any_code <- function(x) {
  # TRUE if value equals ANY_CODE
  !is.na(x) & as.integer(x) == ANY_CODE
}

to_int_or_any <- function(x, lenient = FALSE) {
  # Convert allele tokens to integer; map "any"/"ANY" to ANY_CODE.
  # If lenient = TRUE and token cannot be parsed as numeric, map to ANY_CODE.
  if (is.null(x)) return(ANY_CODE)
  if (length(x) == 0) return(ANY_CODE)
  if (is.na(x)) return(ANY_CODE)
  
  if (is.character(x)) {
    tx <- trimws(x)
    if (tx == "" || tolower(tx) == "any") return(ANY_CODE)
    suppressWarnings(nx <- as.integer(round(as.numeric(tx))))
    if (is.na(nx)) {
      if (lenient) return(ANY_CODE) else stop(sprintf("Non-numeric allele: '%s'", x))
    }
    return(nx)
  }
  
  if (is.numeric(x)) {
    nx <- as.integer(round(x))
    if (is.na(nx)) {
      if (lenient) return(ANY_CODE) else stop("Non-numeric allele (NA)")
    }
    return(nx)
  }
  
  if (is.logical(x)) {
    if (lenient) return(ANY_CODE) else stop("Unsupported allele type (logical)")
  }
  
  # fallback
  if (lenient) return(ANY_CODE)
  stop(sprintf("Unsupported allele type: %s", typeof(x)))
}

encode_any <- function(x, lenient = FALSE) {
  # Vectorized mapping to integer; "any" -> ANY_CODE
  vapply(x, to_int_or_any, integer(1), lenient = lenient)
}

decode_any <- function(x) {
  # Map ANY_CODE to "any"; keep others as character of integer
  out <- as.character(as.integer(x))
  out[is_any_code(x)] <- "any"
  out
}

# ---------- allele pair normalization ----------

order_pair_any_last <- function(a1, a2) {
  # Return ordered pair (left <= right), with ANY_CODE always on the right.
  # For non-any values, ascending numeric order.
  a1 <- as.integer(a1); a2 <- as.integer(a2)
  
  # if both ANY
  if (is_any_code(a1) && is_any_code(a2)) return(c(ANY_CODE, ANY_CODE))
  
  # if one ANY -> put ANY on the right
  if (is_any_code(a1) && !is_any_code(a2)) return(c(a2, ANY_CODE))
  if (!is_any_code(a1) && is_any_code(a2)) return(c(a1, ANY_CODE))
  
  # both non-any -> ascending
  if (a1 <= a2) c(a1, a2) else c(a2, a1)
}

std_allele_pair <- function(a1, a2, lenient = FALSE, homo_to_any = FALSE) {
  # Standardize one allele pair:
  # - encode tokens to integers with ANY_CODE support
  # - reorder so ANY_CODE is on the right, otherwise ascending
  # - optionally convert homozygous into (allele, ANY_CODE)
  x1 <- encode_any(a1, lenient = lenient)
  x2 <- encode_any(a2, lenient = lenient)
  
  if (homo_to_any && !is_any_code(x1) && !is_any_code(x2) && x1 == x2) {
    return(c(x1, ANY_CODE))
  }
  
  order_pair_any_last(x1, x2)
}

# ---------- locus order utilities ----------

load_locus_order <- function(path = file.path("data", "locus_order.rds")) {
  # Expecting a character vector (unique locus names)
  if (!file.exists(path)) stop(sprintf("locus_order RDS not found: %s", path))
  lo <- readRDS(path)
  if (!is.character(lo)) stop("locus_order.rds must be a character vector")
  unique(lo)
}

ensure_all_loci <- function(df, locus_order, fill_any = TRUE) {
  # Ensure that all loci in locus_order exist in df.
  # If missing and fill_any = TRUE, append rows with any,any.
  stopifnot(is.data.frame(df))
  if (!all(c("Locus","Allele1","Allele2") %in% names(df))) {
    stop("df must have columns: Locus, Allele1, Allele2")
  }
  
  present <- unique(df$Locus)
  missing <- setdiff(locus_order, present)
  
  if (length(missing) == 0L) return(df)
  
  if (!fill_any) {
    # Keep as is but warn
    warning(sprintf("Missing loci not filled: %s", paste(missing, collapse = ", ")))
    return(df)
  }
  
  add <- data.frame(
    Locus   = missing,
    Allele1 = rep(ANY_CODE, length(missing)),
    Allele2 = rep(ANY_CODE, length(missing)),
    stringsAsFactors = FALSE
  )
  rbind(df, add)
}

order_by_locus <- function(df, locus_order) {
  # Order rows by locus_order; unknown loci go to the end in original order.
  idx <- match(df$Locus, locus_order)
  o <- order(is.na(idx), idx, seq_len(nrow(df)))
  df[o, , drop = FALSE]
}

# ---------- profile preparation (data.frame) ----------

prepare_profile_df <- function(
    profile_df,
    locus_order = load_locus_order(),
    homo_to_any = FALSE,
    lenient = FALSE,
    fill_missing_any = TRUE
) {
  # Normalize a profile given as data.frame with at least:
  #   Locus, Allele1, Allele2
  # Optional columns are kept but not used here.
  if (!is.data.frame(profile_df)) stop("profile_df must be a data.frame")
  req <- c("Locus","Allele1","Allele2")
  if (!all(req %in% names(profile_df))) {
    stop("profile_df must contain columns: Locus, Allele1, Allele2")
  }
  
  # Copy minimal view
  df <- data.frame(
    Locus   = as.character(profile_df$Locus),
    Allele1 = profile_df$Allele1,
    Allele2 = profile_df$Allele2,
    stringsAsFactors = FALSE
  )
  
  # Standardize pairs row-wise
  if (nrow(df) > 0L) {
    ap <- mapply(
      function(a1, a2) std_allele_pair(a1, a2, lenient = lenient, homo_to_any = homo_to_any),
      df$Allele1, df$Allele2, SIMPLIFY = TRUE
    )
    # mapply returns 2 x N matrix
    df$Allele1 <- as.integer(ap[1L, ])
    df$Allele2 <- as.integer(ap[2L, ])
  }
  
  # Ensure all loci, then order by locus_order
  df2 <- ensure_all_loci(df, locus_order, fill_any = fill_missing_any)
  df3 <- order_by_locus(df2, locus_order)
  
  rownames(df3) <- NULL
  df3
}

# ---------- profile preparation (list) ----------

# Accepts either:
# - named list: names are loci, each element length 2 vector (a1,a2)
# - data.frame: forwarded to prepare_profile_df
prepare_profile <- function(
    profile,
    locus_order = load_locus_order(),
    homo_to_any = FALSE,
    lenient = FALSE,
    fill_missing_any = TRUE
) {
  if (is.data.frame(profile)) {
    return(prepare_profile_df(
      profile_df = profile,
      locus_order = locus_order,
      homo_to_any = homo_to_any,
      lenient = lenient,
      fill_missing_any = fill_missing_any
    ))
  }
  
  if (is.list(profile)) {
    if (is.null(names(profile)) || any(names(profile) == "")) {
      stop("list profile must be a named list with locus names")
    }
    # Build DF from list
    loci <- names(profile)
    vals <- lapply(profile, function(v) {
      v <- as.list(v)
      if (length(v) < 2L) v <- c(v, list(ANY_CODE))
      c(v[[1L]], v[[2L]])
    })
    mat <- do.call(rbind, vals)
    df <- data.frame(
      Locus   = loci,
      Allele1 = mat[,1],
      Allele2 = mat[,2],
      stringsAsFactors = FALSE
    )
    return(prepare_profile_df(
      profile_df = df,
      locus_order = locus_order,
      homo_to_any = homo_to_any,
      lenient = lenient,
      fill_missing_any = fill_missing_any
    ))
  }
  
  stop("profile must be a data.frame or a named list")
}

# ---------- pretty-print helpers (optional) ----------

format_profile_df <- function(df) {
  # Return a copy with Allele1/2 rendered as character, ANY_CODE -> "any"
  if (!all(c("Locus","Allele1","Allele2") %in% names(df))) return(df)
  out <- df
  out$Allele1 <- decode_any(out$Allele1)
  out$Allele2 <- decode_any(out$Allele2)
  out
}

# Return the number of loci, safely.
# locus_order: character vector of loci (required)
get_locus_count <- function(locus_order) {
  if (is.null(locus_order)) return(0L)
  as.integer(length(locus_order))
}

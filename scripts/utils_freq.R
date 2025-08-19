# scripts/utils_freq.R
# Frequency utilities for GlobalFiler 21-locus profiles.
# No multibyte characters in code/comments.

# Load frequency table:
# - Prefer RDS at data/freq_table.rds
# - Fallback to CSV at data/Allele-count_GF_Japanese.csv (columns: Locus, Allele, Count)
load_freq_table <- function(
    rds_path  = "data/freq_table.rds",
    csv_path  = "data/Allele-count_GF_Japanese.csv"
) {
  ft <- NULL
  if (file.exists(rds_path)) {
    ft <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  }
  if (is.null(ft) && file.exists(csv_path)) {
    ft <- tryCatch(utils::read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
  }
  if (is.null(ft)) return(NULL)
  nm <- names(ft)
  nm <- sub("^marker$", "Locus", nm, ignore.case = TRUE)
  nm <- sub("^allele$", "Allele", nm, ignore.case = TRUE)
  nm <- sub("^count$",  "Count",  nm, ignore.case = TRUE)
  names(ft) <- nm
  keep <- intersect(c("Locus", "Allele", "Count"), names(ft))
  ft <- ft[, keep, drop = FALSE]
  ft$Locus  <- as.character(ft$Locus)
  ft$Allele <- as.character(ft$Allele)
  suppressWarnings({ ft$Count <- as.numeric(ft$Count) })
  ft$Count[is.na(ft$Count)] <- 0
  ft
}

# Per-locus genotype frequency with "any" support (neutral handling).
# Rules:
# - if both alleles are "any" -> return 1 (neutral)
# - if exactly one allele is "any" -> return p(other)
# - if allele not found in table (p==0) -> return 1 (neutral, lenient)
# - homozygote -> p^2
# - heterozygote -> 2*p1*p2
calc_freq_locus <- function(locus, allele1, allele2, freq_table) {
  if (is.null(freq_table)) return(NA_real_)
  sub <- freq_table[freq_table$Locus == locus, , drop = FALSE]
  tot <- sum(sub$Count, na.rm = TRUE)
  if (tot <= 0) return(1)
  get_p <- function(a) {
    if (is.null(a) || is.na(a) || !nzchar(a)) return(0)
    cnt <- sub$Count[sub$Allele == a]
    if (length(cnt) == 0 || is.na(cnt)) return(0)
    as.numeric(cnt) / tot
  }
  a1 <- as.character(allele1); a2 <- as.character(allele2)
  if (identical(a1, "any") && identical(a2, "any")) return(1)
  if (identical(a1, "any") && !identical(a2, "any")) return(get_p(a2))
  if (!identical(a1, "any") && identical(a2, "any")) return(get_p(a1))
  p1 <- get_p(a1); p2 <- get_p(a2)
  if (p1 == 0 || p2 == 0) return(1)
  if (identical(a1, a2)) return(p1^2)
  2 * p1 * p2
}

# Append per-locus Freq and return product (TotalFreq) as attribute.
calc_freq_loci_df <- function(profile_df, freq_table) {
  if (is.null(profile_df) || nrow(profile_df) == 0) return(profile_df)
  need <- c("Locus", "Allele1", "Allele2")
  if (!all(need %in% names(profile_df))) return(profile_df)
  profile_df$Freq <- mapply(
    function(lc, a1, a2) calc_freq_locus(lc, a1, a2, freq_table),
    profile_df$Locus, profile_df$Allele1, profile_df$Allele2
  )
  total <- prod(profile_df$Freq, na.rm = TRUE)
  attr(profile_df, "TotalFreq") <- total
  profile_df
}

# Compute only total frequency from a prepared query profile (long format).
calc_total_frequency <- function(profile_df, freq_table) {
  if (is.null(profile_df) || nrow(profile_df) == 0) return(NA_real_)
  need <- c("Locus", "Allele1", "Allele2")
  if (!all(need %in% names(profile_df))) return(NA_real_)
  vals <- mapply(
    function(lc, a1, a2) calc_freq_locus(lc, a1, a2, freq_table),
    profile_df$Locus, profile_df$Allele1, profile_df$Allele2
  )
  prod(vals, na.rm = TRUE)
}

# Format as scientific notation string.
format_total_frequency <- function(x, digits = 3) {
  if (is.null(x) || is.na(x)) return("NA")
  format(x, scientific = TRUE, digits = digits)
}

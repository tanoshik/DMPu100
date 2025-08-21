# scripts/scoring_fast.R
# DMPu100 scoring core. No multibyte chars.

# Any / wildcard code
ANY_CODE <- 9999L

# Bit layout (internal):
#   b0 = (q1 vs r1)
#   b1 = (q1 vs r2)
#   b2 = (q2 vs r1)
#   b3 = (q2 vs r2)
# Code:
#   code = b0 + 2*b1 + 4*b2 + 8*b3
# Display for logs:
#   "b0b1b2b3" (0123 order)
#
# Score table (code 0..15) per u100 spec:
# 0000:0
# 0001:1  0010:1  0011:1
# 0100:1  0101:1  0110:2  0111:2
# 1000:1  1001:2  1010:1  1011:2
# 1100:1  1101:2  1110:2  1111:2
SCORE_TABLE <- as.integer(c(
  0, # 0000
  1, # 0001
  1, # 0010
  1, # 0011
  1, # 0100
  1, # 0101
  2, # 0110
  2, # 0111
  1, # 1000
  2, # 1001
  1, # 1010
  2, # 1011
  1, # 1100
  2, # 1101
  2, # 1110
  2  # 1111
))

# Single allele match with ANY support.
# Returns TRUE if a==b or either is ANY_CODE. NA never matches.
allele_match <- function(a, b) {
  if (is.na(a) || is.na(b)) return(FALSE)
  if (identical(a, ANY_CODE) || identical(b, ANY_CODE)) return(TRUE)
  identical(a, b)
}

# Locus scoring: 2x2 allele grid
# Args:
#   q1,q2: integers (ANY_CODE allowed)
#   r1,r2: integers (ANY_CODE allowed)
# Returns:
#   list(code=<int 0..15>, score=<int 0/1/2>, bits0123=<chr "b0b1b2b3">)
score_2x2 <- function(q1, q2, r1, r2) {
  b0 <- as.integer(allele_match(q1, r1))
  b1 <- as.integer(allele_match(q1, r2))
  b2 <- as.integer(allele_match(q2, r1))
  b3 <- as.integer(allele_match(q2, r2))
  code <- b0 + bitwShiftL(b1, 1L) + bitwShiftL(b2, 2L) + bitwShiftL(b3, 3L)
  score <- SCORE_TABLE[code + 1L]
  bits0123 <- paste0(b0, b1, b2, b3)
  list(code = as.integer(code), score = as.integer(score), bits0123 = bits0123)
}

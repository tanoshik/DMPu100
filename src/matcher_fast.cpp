// src/matcher_fast.cpp
// Minimal Rcpp kernel for DnaMatchProject scoring.
// No multibyte characters in code/comments.

#include <Rcpp.h>
using namespace Rcpp;

/*
 DMP bit definition (exactly matches the project spec):
 b0 = (q1 vs r1)
 b1 = (q1 vs r2)
 b2 = (q2 vs r1)
 b3 = (q2 vs r2)
 Internal code: code = b0*1 + b1*2 + b2*4 + b3*8  (b0 is LSB)
 Display order for 4-bit string is [b3 b2 b1 b0] (not used here).
 ANY_CODE (e.g., 9999) always matches.
 
 Score table (SCORE_TABLE[code]):
 code(Display) : score
 0 (0000) -> 0
 1 (0001) -> 1
 2 (0010) -> 1
 3 (0011) -> 1
 4 (0100) -> 1
 5 (0101) -> 1
 6 (0110) -> 2
 7 (0111) -> 2
 8 (1000) -> 1
 9 (1001) -> 2
 10 (1010) -> 1
 11 (1011) -> 2
 12 (1100) -> 1
 13 (1101) -> 2
 14 (1110) -> 2
 15 (1111) -> 2
 
 Note:
 - Inputs are integer-encoded alleles (uint16-like in R integer), with ANY_CODE used as wildcard.
 - Allele ordering convention ("any" on the right, ascending) is assumed to be applied beforehand on R side.
 */

// Compute a single bit: whether allele a matches allele b considering ANY_CODE.
// Returns 1 if match, 0 otherwise.
inline int bit_match_int(int a, int b, int any_code) {
  if (a == any_code || b == any_code) return 1;
  return (a == b) ? 1 : 0;
}

// [[Rcpp::export]]
IntegerVector compute_scores_uint16(IntegerVector q1,
                                    IntegerVector q2,
                                    IntegerVector r1,
                                    IntegerVector r2,
                                    int any_code = 9999) {
  R_xlen_t n = q1.size();
  if (q2.size() != n || r1.size() != n || r2.size() != n) {
    stop("All input vectors must have the same length.");
  }
  
  // DMP score table for codes 0..15 (matches project spec exactly)
  static int SCORE_TABLE[16] = {
    0, // 0000
    1, // 0001
    1, // 0010
    1, // 0011
    1, // 0100
    1, // 0101
    2, // 0110
    2, // 0111
    1, // 1000
    2, // 1001
    1, // 1010
    2, // 1011
    1, // 1100
    2, // 1101
    2, // 1110
    2  // 1111
  };
  
  IntegerVector out(n);
  
  for (R_xlen_t i = 0; i < n; ++i) {
    // Load alleles
    int a_q1 = q1[i];
    int a_q2 = q2[i];
    int a_r1 = r1[i];
    int a_r2 = r2[i];
    
    // Bits (b0,b1,b2,b3), then code = b0 + 2*b1 + 4*b2 + 8*b3
    int b0 = bit_match_int(a_q1, a_r1, any_code);
    int b1 = bit_match_int(a_q1, a_r2, any_code);
    int b2 = bit_match_int(a_q2, a_r1, any_code);
    int b3 = bit_match_int(a_q2, a_r2, any_code);
    
    int code = b0 + (b1 << 1) + (b2 << 2) + (b3 << 3);
    out[i] = SCORE_TABLE[code];
  }
  
  return out;
}

/*** R
# (Optional) quick sanity check demo from R if you sourceCpp this file standalone:
# q1 <- c(10L, 11L, 9999L, 14L)
# q2 <- c(12L, 11L, 9999L, 15L)
# r1 <- c(10L, 9L,   14L, 15L)
# r2 <- c(12L, 11L,  15L, 14L)
# compute_scores_uint16(q1, q2, r1, r2, 9999L)
*/

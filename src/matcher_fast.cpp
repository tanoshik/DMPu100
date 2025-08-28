// src/matcher_fast.cpp
// Rcpp kernel for DMPu100. No multibyte characters.

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <array>
using namespace Rcpp;

/*
 Scoring spec:
 b0=(q1 vs r1), b1=(q1 vs r2), b2=(q2 vs r1), b3=(q2 vs r2)
 code = b0 + 2*b1 + 4*b2 + 8*b3
 Score table (index = code):
 0:0, 1:1, 2:1, 3:1,
 4:1, 5:1, 6:2, 7:2,
 8:1, 9:2,10:1,11:2,
 12:1,13:2,14:2,15:2
 ANY_CODE matches everything.
 */

// [[Rcpp::export]]
Rcpp::IntegerVector compute_scores_uint16(const int q1,
                                          const int q2,
                                          const Rcpp::IntegerVector& r1,
                                          const Rcpp::IntegerVector& r2,
                                          const int any_code) {
  const R_xlen_t n = r1.size();
  if (r2.size() != n) {
    Rcpp::stop("r1 and r2 must have the same length");
  }
  
  Rcpp::IntegerVector out(n);
  
  const bool q1_any = (q1 == any_code);
  const bool q2_any = (q2 == any_code);
  
  static const std::array<int,16> score_table = {
    0,1,1,1, 1,1,2,2, 1,2,1,2, 1,2,2,2
  };
  
  for (R_xlen_t i = 0; i < n; ++i) {
    const int a1 = r1[i];
    const int a2 = r2[i];
    
    const bool a1_any = (a1 == any_code);
    const bool a2_any = (a2 == any_code);
    
    const bool b0 = q1_any || a1_any || (q1 == a1);
    const bool b1 = q1_any || a2_any || (q1 == a2);
    const bool b2 = q2_any || a1_any || (q2 == a1);
    const bool b3 = q2_any || a2_any || (q2 == a2);
    
    const int code =
      static_cast<int>(b0) +
      (static_cast<int>(b1) << 1) +
      (static_cast<int>(b2) << 2) +
      (static_cast<int>(b3) << 3);
    
    out[i] = score_table[code];
  }
  return out;
}

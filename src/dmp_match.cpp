// src/dmp_match.cpp
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <cstdint>

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame dmp_match_cpp(
    IntegerMatrix A1,           // S x L
    IntegerMatrix A2,           // S x L
    IntegerVector q1,           // L
    IntegerVector q2,           // L
    IntegerVector score_table,  // length 16
    IntegerVector homo_mask,    // length S (uint32 packed in R)
    uint32_t idf_mask_bits,     // 21-bit mask (1=enabled)
    List opts,                  // any_code, h2a_Q, h2a_R, mode, top_n, threshold, display_limit
    CharacterVector sample_ids  // length S
) {
  const int S = A1.nrow();
  const int L = A1.ncol();
  const int any_code = opts.containsElementNamed("any_code") ? as<int>(opts["any_code"]) : 9999;
  const bool force_all_loci = opts.containsElementNamed("force_all_loci") ? as<bool>(opts["force_all_loci"]) : false;
  // const bool h2a_R = opts.containsElementNamed("h2a_R") ? as<bool>(opts["h2a_R"]) : false;
  const std::string mode = opts.containsElementNamed("mode") ? as<std::string>(opts["mode"]) : std::string("topn");
  const int top_n = opts.containsElementNamed("top_n") ? as<int>(opts["top_n"]) : 200;
  const int threshold = opts.containsElementNamed("threshold") ? as<int>(opts["threshold"]) : 0;
  const int display_limit = opts.containsElementNamed("display_limit") ? as<int>(opts["display_limit"]) : 200;
  
  // score buffer
  std::vector<int> scores(S, 0);
  
  // main loop: SoA style: for l then for s
  for (int l = 0; l < L; ++l) {
    const int q1l = q1[l];
    const int q2l = q2[l];
    const uint32_t use_locus = (idf_mask_bits >> l) & 1u;
    
    for (int s = 0; s < S; ++s) {
      // NOTE:
      // - h2a_R handling is expected to be done in R side preprocessing for A1/A2 or via homo_mask logic there.
      // - Here we keep branchless eq: any matches anything, or equal value.
      const int r1 = A1(s, l);
      const int r2 = A2(s, l);
      
      auto eq = [&](int a, int b) -> uint32_t {
        const uint32_t a_any = (a == any_code);
        const uint32_t b_any = (b == any_code);
        const uint32_t same  = (a == b);
        return (a_any | b_any | same);
      };
      
      const uint32_t b0 = eq(q1l, r1);
      const uint32_t b1 = eq(q1l, r2);
      const uint32_t b2 = eq(q2l, r1);
      const uint32_t b3 = eq(q2l, r2);
      const uint32_t code = (b0) | (b1 << 1) | (b2 << 2) | (b3 << 3);
      
      if (use_locus) {
        scores[s] += score_table[code];
      }
    }
  }
  
  // selection
  std::vector<int> idx;
  idx.reserve((mode == "topn") ? std::min(top_n, S) : std::min((display_limit > 0 ? display_limit : S), S));
  
  if (mode == "topn") {
    std::vector<int> order(S);
    std::iota(order.begin(), order.end(), 0);
    const int K = std::min(top_n, S);
    std::partial_sort(order.begin(), order.begin() + K, order.end(),
                      [&](int a, int b) {
                        if (scores[a] != scores[b]) return scores[a] > scores[b];
                        return std::string(sample_ids[a]) < std::string(sample_ids[b]);
                      });
    idx.assign(order.begin(), order.begin() + K);
  } else { // "threshold"
    for (int s = 0; s < S; ++s) {
      if (scores[s] >= threshold) {
        idx.push_back(s);
        if (display_limit > 0 && (int)idx.size() >= display_limit) break;
      }
    }
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
      if (scores[a] != scores[b]) return scores[a] > scores[b];
      return std::string(sample_ids[a]) < std::string(sample_ids[b]);
    });
  }
  
  const int K = (int)idx.size();
  CharacterVector out_id(K);
  IntegerVector out_sc(K);
  for (int i = 0; i < K; ++i) {
    const int s = idx[i];
    out_id[i] = sample_ids[s];
    out_sc[i] = scores[s];
  }
  
  return DataFrame::create(_["SampleID"] = out_id,
                           _["Score"]    = out_sc,
                           _["stringsAsFactors"] = false);
}

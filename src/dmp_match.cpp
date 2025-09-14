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
    uint32_t idf_mask_bits,     // bit=1: enabled locus
    List opts,                  // any_code, force_all_loci, mode, top_n, threshold, display_limit
    CharacterVector sample_ids  // length S
) {
  const int S = A1.nrow();
  const int L = A1.ncol();
  const int any_code = opts.containsElementNamed("any_code") ? as<int>(opts["any_code"]) : 9999;
  const bool force_all = opts.containsElementNamed("force_all_loci") ? as<bool>(opts["force_all_loci"]) : false;
  const std::string mode = opts.containsElementNamed("mode") ? as<std::string>(opts["mode"]) : std::string("topn");
  const int top_n = opts.containsElementNamed("top_n") ? as<int>(opts["top_n"]) : 200;
  const int threshold = opts.containsElementNamed("threshold") ? as<int>(opts["threshold"]) : 0;
  const int display_limit = opts.containsElementNamed("display_limit") ? as<int>(opts["display_limit"]) : top_n;
  
  // enabled loci
  std::vector<int> loci;
  loci.reserve(L);
  if (force_all) {
    for (int l = 0; l < L; ++l) loci.push_back(l);
  } else {
    for (int l = 0; l < L; ++l) if (((idf_mask_bits >> l) & 1U) != 0U) loci.push_back(l);
  }
  const int LE = (int)loci.size();
  
  // score LUT
  int scoreLUT[16];
  for (int i = 0; i < 16; ++i) scoreLUT[i] = score_table[i];
  
  // raw pointers to query
  const int* __restrict Q1 = &q1[0];
  const int* __restrict Q2 = &q2[0];
  
  std::vector<int> scores(S);
  
  // main loop (cache-friendly): outer=locus, inner=sample
  std::fill(scores.begin(), scores.end(), 0);
  for (int ei = 0; ei < LE; ++ei) {
    const int l  = loci[ei];
    const int x1 = Q1[l];
    const int x2 = Q2[l];
    const int* __restrict colA1 = &A1(0, l); // column pointer (contiguous)
    const int* __restrict colA2 = &A2(0, l);
    for (int s = 0; s < S; ++s) {
     const int r1 = colA1[s];
     const int r2 = colA2[s];
     const int m0 = (r1 == x1) || (r1 == any_code) || (x1 == any_code);
     const int m1 = (r2 == x1) || (r2 == any_code) || (x1 == any_code);
     const int m2 = (r1 == x2) || (r1 == any_code) || (x2 == any_code);
     const int m3 = (r2 == x2) || (r2 == any_code) || (x2 == any_code);
     const int code = m0 + (m1 << 1) + (m2 << 2) + (m3 << 3);
     scores[s] += scoreLUT[code];
   }
  }
  
  // selection
  std::vector<int> idx;
  idx.reserve(S);
  if (mode == "topn") {
    const int K = std::min(std::max(top_n, 0), S);
    std::vector<int> order(S);
    std::iota(order.begin(), order.end(), 0);
    
    // cache sample_ids to avoid string construction in comparator
    std::vector<std::string> sid(S);
    for (int i = 0; i < S; ++i) sid[i] = std::string(sample_ids[i]);
    
    std::stable_sort(order.begin(), order.end(), [&](int a, int b){
      if (scores[a] != scores[b]) return scores[a] > scores[b];
      return sid[a] < sid[b];
    });
    idx.assign(order.begin(), order.begin() + K);
  } else { // threshold
    for (int s = 0; s < S; ++s) {
      if (scores[s] >= threshold) {
        idx.push_back(s);
        if (display_limit > 0 && (int)idx.size() >= display_limit) break;
      }
    }
    std::stable_sort(idx.begin(), idx.end(), [&](int a, int b){
      if (scores[a] != scores[b]) return scores[a] > scores[b];
      return std::string(sample_ids[a]) < std::string(sample_ids[b]);
    });
  }
  
  const int K = (int)idx.size();
  CharacterVector out_id(K);
  IntegerVector  out_sc(K);
  for (int i = 0; i < K; ++i) {
    const int s = idx[i];
    out_id[i] = sample_ids[s];
    out_sc[i] = scores[s];
  }
  
  return DataFrame::create(
    _["SampleID"] = out_id,
    _["Score"]    = out_sc,
    _["stringsAsFactors"] = false
  );
}

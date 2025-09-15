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
  
  // ---- selection (bucketized; O(S) + intra-bucket stable) ----
  std::vector<int> idx;
  idx.reserve(S);
  
  // cache SampleID once (same as original)
  std::vector<std::string> sid(S);
  for (int i = 0; i < S; ++i) sid[i] = std::string(sample_ids[i]);
  
  // derive MAX_SC dynamically
  int max_per_locus = scoreLUT[0];
  for (int t = 1; t < 16; ++t) if (scoreLUT[t] > max_per_locus) max_per_locus = scoreLUT[t];
  const int MAX_SC = std::max(0, max_per_locus * LE);
  std::vector<std::vector<int>> buckets(MAX_SC + 1);
  for (int s = 0; s < S; ++s) {
    int sc = scores[s];
    if (sc < 0) sc = 0;
    if (sc > MAX_SC) sc = MAX_SC;
    buckets[sc].push_back(s);
  }
  
  if (mode == "topn" || mode == "top") {
    const int K = std::min(std::max(top_n, 0), S);
    for (int sc = MAX_SC; sc >= 0 && (int)idx.size() < K; --sc) {
      auto &b = buckets[sc];
      if (b.empty()) continue;
      std::stable_sort(b.begin(), b.end(), [&](int a, int b_){ return sid[a] < sid[b_]; });
      for (int v : b) {
        idx.push_back(v);
        if ((int)idx.size() >= K) break;
      }
    }
  } else { // threshold
    const int T = std::max(0, threshold);
    for (int sc = MAX_SC; sc >= T && (display_limit <= 0 || (int)idx.size() < display_limit); --sc) {
      auto &b = buckets[sc];
      if (b.empty()) continue;
      std::stable_sort(b.begin(), b.end(), [&](int a, int b_){ return sid[a] < sid[b_]; });
      if (display_limit > 0) {
        for (int v : b) {
          idx.push_back(v);
          if ((int)idx.size() >= display_limit) break;
        }
      } else {
        idx.insert(idx.end(), b.begin(), b.end());
      }
    }
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

// -------- exported: full-database histogram (Score,Count; always 0..MAX_SC) --------
// 0件のスコアも含め、必ず 1+MAX_SC 行を出す
// [[Rcpp::export]]
DataFrame dmp_hist_cpp(
    IntegerMatrix A1,           // S x L
    IntegerMatrix A2,           // S x L
    IntegerVector q1,           // L
    IntegerVector q2,           // L
    IntegerVector score_table,  // length 16
    uint32_t idf_mask_bits,     // bit=1: enabled locus
    List opts                   // any_code, force_all_loci
) {
  const int S = A1.nrow();
  const int L = A1.ncol();
  if (A2.nrow() != S || A2.ncol() != L) stop("A2 shape mismatch");
  if (q1.size() != L || q2.size() != L) stop("q size mismatch");
  if (score_table.size() != 16)         stop("score_table length must be 16");
  
  const int any_code = opts.containsElementNamed("any_code") ? as<int>(opts["any_code"]) : 9999;
  const bool force_all = opts.containsElementNamed("force_all_loci") ? as<bool>(opts["force_all_loci"]) : false;
  
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
  
  // compute MAX_SC
  int max_per_locus = scoreLUT[0];
  for (int i = 1; i < 16; ++i) if (scoreLUT[i] > max_per_locus) max_per_locus = scoreLUT[i];
  const int MAX_SC = std::max(0, max_per_locus * LE);
  
  // accumulate scores then histogram
  std::vector<int> scores(S, 0);
  for (int ei = 0; ei < LE; ++ei) {
    const int l  = loci[ei];
    const int x1 = Q1[l];
    const int x2 = Q2[l];
    const int* __restrict colA1 = &A1(0, l);
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
  
  // histogram over 0..MAX_SC (descending)
  std::vector<long long> hist(MAX_SC + 1, 0);
  for (int s = 0; s < S; ++s) {
    int sc = scores[s];
    if (sc < 0) sc = 0;
    if (sc > MAX_SC) sc = MAX_SC;
    hist[sc] += 1;
  }
  const int N = MAX_SC + 1;
  IntegerVector Score(N);
  NumericVector Count(N);
  for (int i = 0; i < N; ++i) {
    const int sc = MAX_SC - i;
    Score[i] = sc;
    Count[i] = (double)hist[sc];
  }
  return DataFrame::create(
    _["Score"] = Score,
    _["Count"] = Count,
    _["stringsAsFactors"] = false
  );
}

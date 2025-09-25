// src/dmp_match.cpp
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <cstdint>
#include <climits>  // INT_MIN
using namespace Rcpp;

// Helper: compute score code (2bit x 2 alleles -> 4bit index) is assumed pre-encoded via SCORE_TABLE on R side.

// [[Rcpp::export]]
DataFrame dmp_match_cpp(
    IntegerMatrix A1,           // S x L
    IntegerMatrix A2,           // S x L
    IntegerVector q1,           // L (scaled x100, ANY allowed)
    IntegerVector q2,           // L
    IntegerVector score_table,  // length 16
    uint32_t idf_mask_bits,     // bit=1: enabled locus
    List opts,                  // any_code, score_min(NA=INT_MIN), n_cap, force_all_loci
    CharacterVector sample_ids  // S
) {
  const int S = A1.nrow();
  const int L = A1.ncol();
  if (A2.nrow()!=S || A2.ncol()!=L) stop("A2 size mismatch");
  if (q1.size()!=L || q2.size()!=L) stop("q size mismatch");
  if ((int)sample_ids.size()!=S) stop("sample_ids size mismatch");
  
  // loci enabled by idf_mask_bits or force_all_loci
  std::vector<int> loci; loci.reserve(L);
  if (opts.containsElementNamed("force_all_loci") && as<bool>(opts["force_all_loci"])) {
    for (int j=0;j<L;++j) loci.push_back(j);
  } else {
    for (int j=0;j<L;++j) if ( (idf_mask_bits >> j) & 1u ) loci.push_back(j);
    if (loci.empty()) { for (int j=0;j<L;++j) loci.push_back(j); } // fallback: all
  }
  const int LE = (int)loci.size();
  
  // any code
  const int ANY = opts.containsElementNamed("any_code") ? as<int>(opts["any_code"]) : 9999;
  
  // score LUT
  int scoreLUT[16];
  for (int i=0; i<16; ++i) scoreLUT[i] = score_table[i];
  
  // raw pointers
  const int* Q1 = &q1[0];
  const int* Q2 = &q2[0];
  const int* A1p = &A1[0];
  const int* A2p = &A2[0];
  
  std::vector<int> scores(S,0);
  
  // main loop: outer=locus (enabled), inner=sample
  for (int ei=0; ei<LE; ++ei) {
    const int l = loci[ei];
    const int q1v = Q1[l];
    const int q2v = Q2[l];
    // q_any removed; ANY handled in per-bit comparisons
    for (int s=0; s<S; ++s) {
      const int a1 = A1p[s + l*S];
      const int a2 = A2p[s + l*S];
      // ANY on either side acts as wildcard (q or db)
      int b0 = ((q1v==ANY) || (a1==ANY) || (q1v==a1)) ? 1:0;
      int b1 = ((q1v==ANY) || (a2==ANY) || (q1v==a2)) ? 1:0;
      int b2 = ((q2v==ANY) || (a1==ANY) || (q2v==a1)) ? 1:0;
      int b3 = ((q2v==ANY) || (a2==ANY) || (q2v==a2)) ? 1:0;
      int code = (b0) | (b1<<1) | (b2<<2) | (b3<<3);
      scores[s] += scoreLUT[code];
    }
  }
  
  // unified selection (no R-side sorting): C++ で順序を確定する
  int n_cap = 200;
  if (opts.containsElementNamed("n_cap")) {
    Rcpp::RObject o = opts["n_cap"];
    if (!o.isNULL()) {
      int v = as<int>(o);
      if (v != NA_INTEGER) n_cap = v;
    }
  }
  int score_min = INT_MIN;
  if (opts.containsElementNamed("score_min")) {
    Rcpp::RObject o = opts["score_min"];
    if (!o.isNULL()) {
      int v = as<int>(o);
      if (v != NA_INTEGER) score_min = v;
    }
  }
  const bool use_threshold = (score_min != INT_MIN);
  
  struct Node { int sc; int idx; };
  struct Cmp { bool operator()(const Node& a, const Node& b) const { return a.sc > b.sc; } }; // min-heap
  
  // 事前に ID を C++ 側へコピー（NA は空文字に）
  std::vector<std::string> sid(S);
  for (int i = 0; i < S; ++i) {
    SEXP si = sample_ids[i];  // Rcpp::String 経由ではなく SEXP として扱う
    sid[i] = (si == NA_STRING) ? std::string()
      : std::string(Rf_translateCharUTF8(si)); // UTF-8 に正規化して取り出す
  }
  
  // 最終並びを確定する comparator（score desc, id asc）
  // ここから先は R に触れない
  auto less_final = [&](const Node& a, const Node& b) {
    if (a.sc != b.sc) return a.sc > b.sc;
    const int ia = (a.idx >= 0 && a.idx < S) ? a.idx : 0;
    const int ib = (b.idx >= 0 && b.idx < S) ? b.idx : 0;
    return sid[ia] < sid[ib];
  };
  
    std::vector<Node> pick;  // ここに最終候補を集約
    pick.reserve(n_cap>0 ? n_cap : S);

    if (n_cap <= 0) {
    // キャップなし：全件を閾値でふるい、あとで最終整列
    for (int s=0; s<S; ++s) {
      const int sc = scores[s] + (2 * (L - LE));  // mask_bit=0 loci add +2 each
      if (use_threshold && sc < score_min) continue;
      pick.push_back({sc, s});
     }
     std::stable_sort(pick.begin(), pick.end(), less_final);
    } else {
    // TopN 抽出：min-heapで上位n_cap候補を抽出 → 最終整列
    std::vector<Node> heap; heap.reserve(n_cap + 1);
    for (int s=0; s<S; ++s) {
      const int sc = scores[s] + (2 * (L - LE));
      if (use_threshold && sc < score_min) continue;
      if ((int)heap.size() < n_cap) { heap.push_back({sc,s}); std::push_heap(heap.begin(), heap.end(), Cmp()); }
      else if (sc > heap.front().sc) { std::pop_heap(heap.begin(), heap.end(), Cmp()); heap.back()={sc,s}; std::push_heap(heap.begin(), heap.end(), Cmp()); }
    }
    // heap から降順列を構成し、tie-break を ID 昇順で確定
    std::sort_heap(heap.begin(), heap.end(), Cmp()); // 昇順
    std::reverse(heap.begin(), heap.end());          // 降順に並び替え
    pick.assign(heap.begin(), heap.end());
    std::stable_sort(pick.begin(), pick.end(), less_final);
  }

  const int K = (int)pick.size();
  IntegerVector outScore(K);
  CharacterVector outId(K);
  for (int i=0;i<K;++i){
    outScore[i] = pick[i].sc;
    outId[i]    = sample_ids[ pick[i].idx ];
  }
  return DataFrame::create(_["SampleID"]=outId, _["Score"]=outScore, _["stringsAsFactors"]=false);
}

// [[Rcpp::export]]
DataFrame dmp_hist_cpp(
    IntegerMatrix A1,           // S x L
    IntegerMatrix A2,           // S x L
    IntegerVector q1,           // L
    IntegerVector q2,           // L
    IntegerVector score_table,  // length 16
    uint32_t idf_mask_bits,     // bit=1: enabled locus
    List opts                   // any_code, force_all_loci
){
  const int S = A1.nrow();
  const int L = A1.ncol();
  
  // loci enabled
  std::vector<int> loci; loci.reserve(L);
  if (opts.containsElementNamed("force_all_loci") && as<bool>(opts["force_all_loci"])) {
    for (int j=0;j<L;++j) loci.push_back(j);
  } else {
    for (int j=0;j<L;++j) if ( (idf_mask_bits >> j) & 1u ) loci.push_back(j);
    if (loci.empty()) { for (int j=0;j<L;++j) loci.push_back(j); }
  }
  const int LE = (int)loci.size();
  
  const int ANY = opts.containsElementNamed("any_code") ? as<int>(opts["any_code"]) : 9999;
  
  int scoreLUT[16];
  for (int i=0;i<16;++i) scoreLUT[i] = score_table[i];
  
  // compute max possible score for hist size
  // Assume per-locus max is max(score_table). Total MAX_SC = max_per_locus * L.
  int max_per = 0; for (int i=0;i<16;++i) if (scoreLUT[i] > max_per) max_per = scoreLUT[i];
  const int MAX_SC = max_per * L;
  
  std::vector<int> hist(MAX_SC + 1, 0);
  
  const int* Q1 = &q1[0];
  const int* Q2 = &q2[0];
  const int* A1p = &A1[0];
  const int* A2p = &A2[0];
  
  for (int s=0; s<S; ++s) {
    int sc = 0;
    for (int ei=0; ei<LE; ++ei) {
      const int l = loci[ei];
      const int q1v = Q1[l];
      const int q2v = Q2[l];
      const int a1 = A1p[s + l*S];
      const int a2 = A2p[s + l*S];
      // ANY on either side acts as wildcard (q or db)
      int b0 = ((q1v==ANY) || (a1==ANY) || (q1v==a1)) ? 1:0;
      int b1 = ((q1v==ANY) || (a2==ANY) || (q1v==a2)) ? 1:0;
      int b2 = ((q2v==ANY) || (a1==ANY) || (q2v==a1)) ? 1:0;
      int b3 = ((q2v==ANY) || (a2==ANY) || (q2v==a2)) ? 1:0;
      int code = (b0) | (b1<<1) | (b2<<2) | (b3<<3);
      sc += scoreLUT[code];
    }
    sc += (2 * (L - LE)); // mask_bit=0 loci add +2 each
    if (sc<0) sc=0; if (sc>MAX_SC) sc=MAX_SC;
    hist[sc] += 1;
  }
  
  // output Score (desc) and Count
  const int N = MAX_SC + 1;
  IntegerVector Score(N);
  NumericVector Count(N);
  for (int i=0;i<N;++i){
    const int sc = MAX_SC - i;
    Score[i] = sc;
    Count[i] = (double)hist[sc];
  }
  return DataFrame::create(_["Score"]=Score, _["Count"]=Count, _["stringsAsFactors"]=false);
}


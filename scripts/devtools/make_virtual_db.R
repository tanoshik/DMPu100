# --- make_virtual_db.R : GlobalFiler 21 loci Virtual DB generator (u100) ---
# hard requirements:
#   - data/freq_table.rds : 必須。無ければ停止
#   - data/locus_order.rds : 必須。無ければ停止
# features:
#   - 小数アレル（例: 12.2）に対応
#   - SCALE = 1 or 100 を自動決定（小数が含まれたら 100）
#   - 2倍体サンプリング（i.i.d.）をローカス毎に実施
#   - メタ情報を同梱（再現性・検証用）
#   - チャンク生成によりメモリ負荷を抑制

suppressWarnings({
  options(stringsAsFactors = FALSE)
})

library(tools)  # md5sum 用

# ====== config ======
ANY_CODE   <- 9999L
SEED       <- 123L
N_TOTAL    <- 2_000_000L     # デフォルトは 2M
CHUNK_SIZE <- 100_000L
AUTO_RUN   <- FALSE          # 手動呼び出しを基本に

# 出力パス
out_path_for <- function(n_total, seed) {
  file.path("data", sprintf("virtual_db_u100_S%d_seed%d.rds", as.integer(n_total), as.integer(seed)))
}

# ====== utils ======
ts_ <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")

fail <- function(msg) { ts_("[FAIL] ", msg); stop(msg, call. = FALSE) }

must_read_rds <- function(path, label) {
  if (!file.exists(path)) fail(sprintf("%s not found: %s", label, path))
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj)) fail(sprintf("failed to read %s: %s", label, path))
  obj
}

detect_freq_col <- function(df) {
  cn <- tolower(names(df))
  hit <- which(cn %in% c("frequency","freq","p","prob","probability"))
  if (length(hit) == 0L) return(NULL)
  names(df)[hit[1L]]
}

# 文字/数値混在 Allele を numeric に正規化（"12.2" -> 12.2）
norm_allele_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  suppressWarnings(as.numeric(as.character(x)))
}

# スケール決定（小数が1つでもあれば 100）
decide_scale <- function(freq_df, locus_ids, fcol) {
  rows <- freq_df$Locus %in% locus_ids
  alle <- norm_allele_num(freq_df$Allele[rows])
  if (any(!is.na(alle) & (alle != floor(alle)))) 100L else 1L
}

# 1ローカスの確率表を作成
build_per_locus_prob <- function(freq_df, locus_ids, fcol) {
  lapply(locus_ids, function(L) {
    rows <- which(as.character(freq_df$Locus) == as.character(L))
    if (length(rows) == 0L) fail(paste0("no frequency rows for locus: ", L))
    alle <- norm_allele_num(freq_df$Allele[rows])
    pr   <- suppressWarnings(as.numeric(freq_df[[fcol]][rows]))
    ok <- is.finite(alle) & is.finite(pr) & pr > 0
    alle <- alle[ok]; pr <- pr[ok]
    if (length(alle) == 0L) fail(paste0("empty/invalid frequency entries at locus: ", L))
    s <- sum(pr); if (!is.finite(s) || s <= 0) fail(paste0("non-positive prob sum at locus: ", L))
    pr <- pr / s
    list(alleles = alle, probs = pr)
  })
}

# 2m個のアレルを i.i.d. に引いて A1/A2 へ
sample_two_alleles <- function(m, alleles, probs) {
  idx <- sample.int(length(alleles), size = 2L*m, replace = TRUE, prob = probs)
  a   <- alleles[idx]
  list(A1 = a[seq_len(m)], A2 = a[m + seq_len(m)])
}

# 符号化（整数へ）。SCALE=1 なら整数化、SCALE=100 なら *100
encode_int <- function(x, SCALE) {
  y <- as.numeric(x)
  as.integer(round(y * SCALE))
}

# ====== main ======
make_virtual_db <- function(n_total = N_TOTAL, seed = SEED,
                            chunk_size = CHUNK_SIZE,
                            any_code = ANY_CODE,
                            out_path = out_path_for(n_total, seed)) {
  
  # 必須ファイル読み込み
  freq_path  <- "data/freq_table.rds"
  locus_path <- "data/locus_order.rds"
  if (!file.exists(freq_path))  fail("data/freq_table.rds is mandatory and missing")
  if (!file.exists(locus_path)) fail("data/locus_order.rds is mandatory and missing")
  
  freq_df <- must_read_rds(freq_path, "freq_table.rds")
  locus_ids <- must_read_rds(locus_path, "locus_order.rds")
  if (!is.data.frame(freq_df) || nrow(freq_df) == 0L) fail("freq_table.rds is empty/invalid")
  if (!all(c("Locus","Allele") %in% names(freq_df))) fail("freq_table must have Locus, Allele, and a frequency column")
  if (!is.character(locus_ids) || length(locus_ids) < 1L) fail("locus_order.rds must be a character vector")
  
  fcol <- detect_freq_col(freq_df)
  if (is.null(fcol)) fail("frequency column not found in freq_table")
  
  # ローカス順に合わせて確率表を構築
  per_locus <- build_per_locus_prob(freq_df, locus_ids, fcol)
  
  # SCALE 決定（小数が含まれたら 100）
  SCALE <- decide_scale(freq_df, locus_ids, fcol)
  ts_(sprintf("[info] SCALE decided = %d (1:integer only, 100:decimal present)", SCALE))
  
  # 出力ディレクトリ
  if (!dir.exists("data")) dir.create("data", recursive = TRUE, showWarnings = FALSE)
  
  set.seed(seed)
  S <- as.integer(n_total); L <- length(locus_ids)
  ts_(sprintf("start make_virtual_db: S=%d, L=%d, chunk=%d, seed=%d", S, L, chunk_size, seed))
  
  # 直書き（行=sample, 列=locus）で整数化を保持
  A1 <- matrix(any_code, nrow = S, ncol = L, dimnames = list(NULL, locus_ids))
  A2 <- matrix(any_code, nrow = S, ncol = L, dimnames = list(NULL, locus_ids))
  
  n_chunks <- ceiling(S / chunk_size)
  t0 <- proc.time()[3L]
  
  for (ck in seq_len(n_chunks)) {
    i1 <- (ck - 1L) * chunk_size + 1L
    i2 <- min(S, ck * chunk_size)
    m  <- i2 - i1 + 1L
    set.seed(seed + ck)
    
    for (j in seq_len(L)) {
      pl <- per_locus[[j]]
      samp <- sample_two_alleles(m, pl$alleles, pl$probs)  # numeric
      A1[i1:i2, j] <- encode_int(samp$A1, SCALE)
      A2[i1:i2, j] <- encode_int(samp$A2, SCALE)
    }
    ts_(sprintf("chunk %d/%d done (rows %d..%d)", ck, n_chunks, i1, i2))
  }
  
  # ANY_CODE 率の軽いチェック
  any_ratio <- (sum(A1 == any_code) + sum(A2 == any_code)) / ((length(A1) + length(A2)))
  if (any_ratio > 0) ts_(sprintf("[warn] ANY_CODE rate = %.6f (should be 0 for pure sampling)", any_ratio))
  
  # メタ
  meta <- list(
    engine      = "make_virtual_db.R",
    version     = "u100",
    scale       = SCALE,
    seed        = seed,
    n_total     = S,
    n_loci      = L,
    locus_order = locus_ids,
    freq_hash   = unname(md5sum(freq_path)),
    locus_hash  = unname(md5sum(locus_path)),
    timestamp   = as.character(Sys.time())
  )
  
  out <- list(
    sample_ids = sprintf("V%07d", seq_len(S)),
    locus_ids  = locus_ids,
    A1 = A1, A2 = A2,
    meta = meta
  )
  
  saveRDS(out, file = out_path, compress = "gzip")
  t1 <- proc.time()[3L]
  ts_(sprintf("wrote: %s (S=%d, L=%d, compress=gzip) in %.2f sec (%.2f min)",
              out_path, S, L, t1 - t0, (t1 - t0)/60))
  invisible(out_path)
}

# -- helper: 1k smoke（頻度の粗チェック） --------------------------------------
smoke_1k_check <- function(seed = SEED) {
  path <- out_path_for(1000L, seed)
  make_virtual_db(n_total = 1000L, seed = seed)
  db <- readRDS(path)
  A1 <- db$A1; A2 <- db$A2; L <- colnames(A1)
  
  # 観測頻度（各ローカスごとに *SCALE から逆変換）
  SCALE <- db$meta$scale
  obs <- lapply(seq_along(L), function(j){
    x <- c(A1[,j], A2[,j])
    x <- x[x != ANY_CODE]
    tab <- sort(table(round(x / SCALE, 2)), decreasing = TRUE)
    data.frame(Locus = L[j], Allele = as.numeric(names(tab)), Count = as.integer(tab), stringsAsFactors = FALSE)
  })
  obs_df <- do.call(rbind, obs)
  
  cat("[SMOKE] 1k generated at:", path, "\n")
  # 代表3ローカスだけ上位を出力
  pick <- head(unique(obs_df$Locus), 3)
  for (loc in pick) {
    cat("  Locus:", loc, "\n")
    print(head(subset(obs_df, Locus == loc), 10), row.names = FALSE)
  }
  invisible(path)
}

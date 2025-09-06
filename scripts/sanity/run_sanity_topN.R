# scripts/sanity/run_sanity_topN.R
# No multibyte characters. TopN-only detail extractor (O(TopN * L)).

# deps: use existing fast matcher and new detail helper
source("scripts/matcher_fast.R")   # provides run_match_fast()
source("scripts/matcher_detail.R") # provides build_detail_for_ids()

# --- default options (house style) ---
opt <- list(
  db       = "data/virtual_db_u100_S1000_seed123.rds",
  query    = "data/query_profile_seed123.csv",
  out      = "output/sample/match_detail_S1000_topN.csv",
  top      = 3L,
  use_cpp  = FALSE,
  any_code = 9999L
)

# --- tiny arg parser: supports --key=value or --key value ---
argv <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(argv)) {
  a <- argv[i]
  if (startsWith(a, "--")) {
    kv <- sub("^--", "", a)
    if (grepl("=", kv, fixed = TRUE)) {
      key <- sub("=.*$", "", kv)
      val <- sub("^[^=]*=", "", kv)
    } else {
      key <- kv
      if ((i + 1L) <= length(argv) && !startsWith(argv[i + 1L], "--")) {
        val <- argv[i + 1L]; i <- i + 1L
      } else {
        val <- TRUE
      }
    }
    if (key %in% c("db", "query", "out")) {
      opt[[key]] <- as.character(val)
    } else if (key %in% c("top", "any_code")) {
      opt[[key]] <- as.integer(val)
    } else if (key %in% c("use_cpp")) {
      opt[[key]] <- as.logical(val)
    } else {
      warning(sprintf("Unknown option ignored: --%s", key))
    }
  }
  i <- i + 1L
}

# --- validate args ---
req <- c("db","query","out")
miss <- req[!vapply(req, function(k) !is.null(opt[[k]]) && nzchar(as.character(opt[[k]])), FALSE)]
if (length(miss)) {
  stop(paste0(
    "Usage: Rscript scripts/sanity/run_sanity_topN.R ",
    "--db=... --query=... --out=... [--top=3] [--use_cpp=TRUE|FALSE] [--any_code=9999]\n"
  ))
}

# --- load DB index ---
if (!file.exists(opt$db)) stop(sprintf("DB not found: %s", opt$db))
db <- readRDS(opt$db)
if (is.null(db$A1) || is.null(db$A2)) stop("Invalid DB index RDS (missing A1/A2)")

# --- load query (human-readable or x100 both OK) ---
q_raw <- read.csv(opt$query, stringsAsFactors = FALSE)
need_cols <- c("Locus","allele1","allele2")
if (!all(need_cols %in% names(q_raw))) {
  stop("Query CSV must have columns: Locus, allele1, allele2")
}
enc_one <- function(x, ANY=opt$any_code) {
  if (is.na(x)) return(NA_integer_)
  if (is.character(x) && tolower(x) == "any") return(as.integer(ANY))
  as.integer(round(as.numeric(x) * 100))
}
q_df <- data.frame(
  Locus   = q_raw$Locus,
  allele1 = vapply(q_raw$allele1, enc_one, 0L),
  allele2 = vapply(q_raw$allele2, enc_one, 0L),
  stringsAsFactors = FALSE
)

# reorder query to DB locus order if provided
if (!is.null(db$locus_ids)) {
  ord <- match(db$locus_ids, q_df$Locus)
  if (anyNA(ord)) stop("Query loci do not cover DB locus_ids")
  q_df <- q_df[ord, , drop=FALSE]
}

# --- Stage 1: score-only (no detail) to get TopN ---
res <- run_match_fast(
  q_df = q_df,
  db_index = db,
  use_cpp = isTRUE(opt$use_cpp),
  include_bits_in_detail = FALSE
)
sum_tab <- res$summary
if (nrow(sum_tab) < opt$top) {
  warning(sprintf("Requested top=%d but only %d candidates exist", opt$top, nrow(sum_tab)))
}
top_ids <- head(sum_tab$SampleID, opt$top)

# --- Stage 2: build detail only for TopN (O(TopN * L)) ---
det <- build_detail_for_ids(
  q_df = q_df, db = db, sample_ids = top_ids,
  any_code = opt$any_code,
  decode = TRUE, include_code = TRUE, include_bits = TRUE
)

# --- order: keep TopN score order, then locus order ---
# sum_tab は run_match_fast() が返した合計スコアの降順テーブル
# top_ids は head(sum_tab$SampleID, opt$top) なので、これが「Score降順」の正解順
score_order <- top_ids  # already in score-desc / tie-broken by SampleID
det$SampleID <- factor(as.character(det$SampleID),
                       levels = score_order, ordered = TRUE)

if (!is.null(db$locus_ids)) {
  det <- det[order(det$SampleID, match(det$Locus, db$locus_ids)), , drop = FALSE]
} else {
  det <- det[order(det$SampleID, det$Locus), , drop = FALSE]
}

# write
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.csv(det, opt$out, row.names = FALSE)
cat("[OK] wrote (TopN detail):", opt$out, "rows=", nrow(det), "\n")

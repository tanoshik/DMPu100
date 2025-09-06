# scripts/sanity/topn_sum_csv.R
# No multibyte chars. Write Top-N total scores as CSV (lightweight).

# ----- defaults (house style) -----
defaults <- list(
  db    = "data/virtual_db_u100_S1000_seed123.rds",
  query = "data/query_profile_seed123.csv",
  out   = "output/sample/match_scores_S1000_top3.csv",
  top   = 3L,
  any_code = 9999L,
  allow_defaults = TRUE
)

opt <- defaults
keys_set <- character(0)

# ----- tiny arg parser: supports --k=v  and  --k v -----
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
        val <- "TRUE"
      }
    }
    # normalize value types by key
    if (key %in% c("db","query","out")) {
      opt[[key]] <- as.character(val)
      keys_set <- union(keys_set, key)
    } else if (key %in% c("top","any_code")) {
      opt[[key]] <- as.integer(val)
      keys_set <- union(keys_set, key)
    } else if (key %in% c("allow_defaults")) {
      lv <- tolower(as.character(val))
      opt[[key]] <- lv %in% c("1","true","t","yes","y")
      keys_set <- union(keys_set, key)
    } else {
      warning(sprintf("Unknown option ignored: --%s", key))
    }
    end_if <- TRUE
  }
  i <- i + 1L
}

# ----- validate requireds (with default-guard) -----
req <- c("db","query","out")
if (!isTRUE(opt$allow_defaults)) {
  # user must have *explicitly* set all required keys
  missing_explicit <- setdiff(req, keys_set)
  if (length(missing_explicit)) {
    stop(sprintf("[FATAL] The following options are still at defaults: %s\nPass explicit values (e.g. --db=path --query=path --out=path) or add --allow_defaults=TRUE if you really intend to use defaults.",
                 paste(missing_explicit, collapse=", ")))
  }
}

# ----- load resources -----
if (!file.exists(opt$db)) stop(sprintf("DB not found: %s", opt$db))
loc <- readRDS("data/locus_order.rds")
x <- readRDS(opt$db)
q <- read.csv(opt$query, stringsAsFactors = FALSE)

# ----- normalize query: allow human-readable (decimals/any) or legacy x100 -----
encode_one <- function(xx, ANY=as.integer(opt$any_code)) {
  if (is.na(xx)) return(ANY)             # NAはany寄せ（必要なら運用で変更）
  if (is.character(xx)) {
    sx <- trimws(tolower(xx))
    if (sx == "any") return(ANY)
    vx <- suppressWarnings(as.numeric(sx))
    if (!is.na(vx)) return(as.integer(round(vx * 100)))
    stop("Unrecognized allele format in query: ", xx)
  } else if (is.numeric(xx)) {
    # 100で割り切れる整数（例: 1500）は x100 とみなし、そのまま返す
    if (abs(xx - round(xx)) > 1e-9) return(as.integer(round(xx * 100)))
    xi <- as.integer(round(xx))
    if (xi %% 100 == 0) return(xi) else return(as.integer(xi * 100))
  } else {
    stop("Unsupported allele type in query: ", typeof(xx))
  }
}

need_cols <- c("Locus","allele1","allele2")
if (!all(need_cols %in% names(q))) {
  stop("query CSV must have columns: Locus, allele1, allele2")
}
q$allele1 <- vapply(q$allele1, encode_one, integer(1))
q$allele2 <- vapply(q$allele2, encode_one, integer(1))

# align query to locus order
names(q) <- c("Locus","allele1","allele2")
idx <- match(loc, q$Locus)
if (anyNA(idx)) {
  missing <- loc[is.na(idx)]
  stop("Query loci do not cover the DB locus set: ", paste(missing, collapse=", "))
}
q <- q[idx, , drop=FALSE]

# ----- scoring (bit pattern → SCORE_TABLE) -----
SCORE <- c(0,1,1,1,1,1,2,2,1,2,1,2,1,2,2,2)
ANY <- as.integer(opt$any_code)
S <- nrow(x$A1)
tot <- integer(S)

for (j in seq_along(loc)) {
  q1 <- as.integer(q$allele1[j]); q2 <- as.integer(q$allele2[j])
  r1 <- x$A1[, j]; r2 <- x$A2[, j]
  b0 <- (r1==q1)|(r1==ANY)|(q1==ANY)
  b1 <- (r2==q1)|(r2==ANY)|(q1==ANY)
  b2 <- (r1==q2)|(r1==ANY)|(q2==ANY)
  b3 <- (r2==q2)|(r2==ANY)|(q2==ANY)
  code <- (b0*1L) + (b1*2L) + (b2*4L) + (b3*8L)
  tot <- tot + SCORE[code + 1L]
}

# ----- select TopN & write -----
o <- order(-tot, x$sample_ids)
k <- min(as.integer(opt$top), length(o))
res <- data.frame(
  SampleID = x$sample_ids[o][seq_len(k)],
  Score    = tot[o][seq_len(k)],
  stringsAsFactors = FALSE
)

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.csv(res, opt$out, row.names = FALSE)
cat("[OK] wrote: ", opt$out, "\n", sep = "")

# scripts/sanity/run_sanity_top3.R
# No multibyte chars. Use run_match_fast() to emit per-locus Top-N detail (canonical sort: SampleID -> Locus).

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  db       = "data/virtual_db_u100_S1000_seed123.rds",
  query    = "data/query_profile_seed123.csv",
  out      = "output/sample/match_detail_S1000_top3_from_engine.csv",
  top      = 3L,
  use_cpp  = FALSE,
  any_code = 9999L
)
for (a in args) {
  kv <- strsplit(a, "=", fixed=TRUE)[[1]]
  if (length(kv)==2) opt[[sub("^--?","",kv[1])]] <- kv[2]
}
opt$top      <- as.integer(opt$top)
opt$use_cpp  <- as.logical(opt$use_cpp)
opt$any_code <- as.integer(opt$any_code)

source("scripts/scoring_fast.R")
source("scripts/matcher_fast.R")

if (!file.exists("data/locus_order.rds")) stop("missing data/locus_order.rds")
locus_order <- readRDS("data/locus_order.rds")
db <- readRDS(opt$db)
if (!identical(as.character(db$locus_ids), as.character(locus_order))) {
  stop("db$locus_ids must equal locus_order")
}

q <- read.csv(opt$query, stringsAsFactors = FALSE)
nn <- names(q)
nn <- sub("(?i)^locus$","Locus", nn, perl=TRUE)
nn <- sub("(?i)^allele1$","allele1", nn, perl=TRUE)
nn <- sub("(?i)^allele2$","allele2", nn, perl=TRUE)
names(q) <- nn
q <- q[, c("Locus","allele1","allele2")]
map_any <- function(x, ANY){
  x <- as.character(x); x[x==""] <- "any"
  x[tolower(x)=="any"] <- as.character(ANY)
  suppressWarnings(as.integer(x))
}
q$allele1 <- map_any(q$allele1, opt$any_code)
q$allele2 <- map_any(q$allele2, opt$any_code)
o <- match(locus_order, q$Locus)
if (any(is.na(o))) stop("query loci mismatch against locus_order")
q <- q[o,,drop=FALSE]; rownames(q) <- NULL

res <- run_match_fast(q_df = q, db_index = db,
                      include_bits_in_detail = TRUE,
                      use_cpp = opt$use_cpp)

sum_tab <- res$summary
sum_tab <- sum_tab[order(-sum_tab$Score, sum_tab$SampleID), ]
top_ids <- head(sum_tab$SampleID, opt$top)

det <- res$detail
names(det) <- sub("(?i)^locus$","Locus", names(det), perl=TRUE)
names(det) <- sub("(?i)^sampleid$","SampleID", names(det), perl=TRUE)
names(det) <- sub("(?i)^db_allele1$","DB_Allele1", names(det), perl=TRUE)
names(det) <- sub("(?i)^db_allele2$","DB_Allele2", names(det), perl=TRUE)

qmap <- data.frame(Locus = locus_order, q1=q$allele1, q2=q$allele2, stringsAsFactors = FALSE)
idx <- match(det$Locus, qmap$Locus)
det$q1 <- qmap$q1[idx]; det$q2 <- qmap$q2[idx]

need_cols <- c("code","bits","score")
missing_cols <- setdiff(need_cols, names(det))
if (length(missing_cols) > 0) {
  tmp <- mapply(function(q1,q2,r1,r2){
    s <- score_2x2(q1,q2,r1,r2)
    c(code=as.integer(s$code), bits=paste0(rev(strsplit(s$bits0123,"")[[1]]), collapse=""), score=as.integer(s$score))
  }, det$q1, det$q2, det$DB_Allele1, det$DB_Allele2)
  tmp <- as.data.frame(t(tmp), stringsAsFactors=FALSE)
  tmp$code  <- as.integer(tmp$code)
  tmp$score <- as.integer(tmp$score)
  det$code  <- tmp$code
  det$bits  <- tmp$bits
  det$score <- tmp$score
} else {
  det$bits <- vapply(det$bits, function(x){
    xs <- as.character(x)
    if (nchar(xs)==4) paste0(substr(xs,4,4),substr(xs,3,3),substr(xs,2,2),substr(xs,1,1)) else xs
  }, character(1))
}

det <- det[det$SampleID %in% top_ids, , drop=FALSE]

decode <- function(a, ANY=opt$any_code) {
  if (is.na(a)) return("")
  if (as.character(a) == as.character(ANY)) return("any")
  v <- as.numeric(a) / 100
  if (abs(v - round(v)) < 1e-9) as.character(as.integer(round(v))) else format(v, trim=TRUE, nsmall=1)
}

out <- data.frame(
  Locus = det$Locus,
  q1 = vapply(det$q1, decode, character(1)),
  q2 = vapply(det$q2, decode, character(1)),
  r1 = vapply(det$DB_Allele1, decode, character(1)),
  r2 = vapply(det$DB_Allele2, decode, character(1)),
  bits = det$bits,
  code = det$code,
  score = det$score,
  SampleID = det$SampleID,
  stringsAsFactors = FALSE
)

# --- canonical sort: SampleID -> Locus (by locus_order)
out$Locus <- factor(out$Locus, levels = locus_order, ordered = TRUE)
out <- out[order(out$SampleID, out$Locus), , drop=FALSE]
rownames(out) <- NULL

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file = opt$out, row.names = FALSE)
cat("[OK] wrote: ", opt$out, " rows=", nrow(out), "\n", sep="")

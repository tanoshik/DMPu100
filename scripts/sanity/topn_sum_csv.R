# scripts/sanity/topn_sum_csv.R
# No multibyte chars. Write Top-N total scores as CSV (lightweight).

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  db    = "data/virtual_db_u100_S1000_seed123.rds",
  query = "data/query_profile_seed123.csv",
  out   = "output/sample/match_scores_S1000_top3.csv",
  top   = 3L, any_code = 9999L
)
for (a in args){ kv <- strsplit(a,"=",TRUE)[[1]]; if(length(kv)==2) opt[[sub("^--?","",kv[1])]] <- kv[2] }
opt$top <- as.integer(opt$top); ANY <- as.integer(opt$any_code)

loc <- readRDS("data/locus_order.rds")
x <- readRDS(opt$db)
q <- read.csv(opt$query, stringsAsFactors=FALSE)
names(q) <- c("Locus","allele1","allele2")
q <- q[match(loc,q$Locus),]

SCORE <- c(0,1,1,1,1,1,2,2,1,2,1,2,1,2,2,2)
tot <- integer(nrow(x$A1))
for (j in seq_along(loc)) {
  q1 <- as.integer(q$allele1[j]); q2 <- as.integer(q$allele2[j])
  r1 <- x$A1[,j]; r2 <- x$A2[,j]
  b0 <- (r1==q1)|(r1==ANY)|(q1==ANY); b1 <- (r2==q1)|(r2==ANY)|(q1==ANY)
  b2 <- (r1==q2)|(r1==ANY)|(q2==ANY); b3 <- (r2==q2)|(r2==ANY)|(q2==ANY)
  code <- (b0*1)+(b1*2)+(b2*4)+(b3*8)
  tot <- tot + SCORE[code+1]
}
o <- order(-tot, x$sample_ids)
k <- min(opt$top, length(o))
res <- data.frame(SampleID = x$sample_ids[o][seq_len(k)],
                  Score    = tot[o][seq_len(k)],
                  stringsAsFactors=FALSE)
dir.create(dirname(opt$out), recursive=TRUE, showWarnings=FALSE)
write.csv(res, opt$out, row.names=FALSE)
cat("[OK] wrote: ", opt$out, "\n", sep="")

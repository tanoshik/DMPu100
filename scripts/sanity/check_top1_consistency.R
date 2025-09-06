# scripts/sanity/check_top1_consistency.R
# No multibyte chars. Compare Top1 between summary and detail.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/sanity/check_top1_consistency.R <summary_csv> <detail_csv>")
}
summary_csv <- args[1]
detail_csv  <- args[2]

s <- read.csv(summary_csv, stringsAsFactors = FALSE, check.names = FALSE)
d <- read.csv(detail_csv,  stringsAsFactors = FALSE, check.names = FALSE)

# summary Top1
if (!all(c("SampleID","Score") %in% names(s))) {
  stop("summary csv must have columns: SampleID, Score")
}
top1s <- as.character(s$SampleID[order(-s$Score, s$SampleID)][1])

# detail Top1 (aggregate score by SampleID)
if (!all(c("SampleID","score") %in% names(d))) {
  stop("detail csv must have columns: SampleID, score")
}
agg <- aggregate(score ~ SampleID, data = d, FUN = sum)
top1d <- as.character(agg$SampleID[order(-agg$score, agg$SampleID)][1])

cat("[CHECK] summary:", top1s, " detail:", top1d, " MATCH=", top1s == top1d, "\n", sep = "")
if (!identical(top1s, top1d)) quit(status = 2)

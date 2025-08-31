# scripts/sanity/compare_routes_top3.R
# No multibyte chars. Route parity check R vs CPP (canonical sort).

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  db    = "data/virtual_db_u100_S1000_seed123.rds",
  query = "data/query_profile_seed123.csv",
  top   = 3L
)
for (a in args) { kv <- strsplit(a,"=",fixed=TRUE)[[1]]; if (length(kv)==2) opt[[sub("^--?","",kv[1])]] <- kv[2] }
opt$top <- as.integer(opt$top)

loc <- readRDS("data/locus_order.rds")
canon <- function(df){
  n <- names(df)
  n <- sub("(?i)^locus$","Locus", n, perl=TRUE)
  n <- sub("(?i)^sampleid$","SampleID", n, perl=TRUE)
  names(df) <- n
  df$Locus <- factor(df$Locus, levels=loc, ordered=TRUE)
  df <- df[order(df$SampleID, df$Locus), , drop=FALSE]
  rownames(df) <- NULL
  df
}

dir.create("output/sample", recursive=TRUE, showWarnings=FALSE)
cmdR <- function(use_cpp, out){
  system2(command = file.path(R.home("bin"), "Rscript"),
          args = c("scripts/sanity/run_sanity_top3.R",
                   paste0("--db=", opt$db),
                   paste0("--query=", opt$query),
                   paste0("--out=", out),
                   paste0("--top=", opt$top),
                   paste0("--use_cpp=", if (use_cpp) "TRUE" else "FALSE")),
          stdout = TRUE, stderr = TRUE)
}

out_r  <- "output/sample/match_detail_S1000_top3_routeR.csv"
out_cpp<- "output/sample/match_detail_S1000_top3_routeCPP.csv"
invisible(cmdR(FALSE, out_r))
invisible(cmdR(TRUE , out_cpp))

xr <- canon(read.csv(out_r , stringsAsFactors=FALSE))
xc <- canon(read.csv(out_cpp, stringsAsFactors=FALSE))

# allow only exact equality after canonical sort
same <- identical(xr, xc)
if (same) {
  cat("[PASS] R vs CPP parity: identical on Top", opt$top, "\n", sep="")
} else {
  m <- merge(xr, xc, by=c("SampleID","Locus"), suffixes=c(".R",".CPP"), all=TRUE, sort=FALSE)
  write.csv(m, "output/sample/top3_route_diff.csv", row.names=FALSE)
  cat("[FAIL] R vs CPP parity mismatches -> output/sample/top3_route_diff.csv\n")
}

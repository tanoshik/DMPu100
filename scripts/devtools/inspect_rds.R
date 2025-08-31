# scripts/devtools/inspect_rds.R
# No multibyte chars. Minimal RDS inspector.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  path = "data/freq_table.rds",
  out  = "output/dev/inspect_freq_table.txt"
)
for (a in args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2) opt[[sub("^--?","",kv[1])]] <- kv[2]
}

if (!file.exists(opt$path)) {
  stop("file not found: ", opt$path)
}

if (!is.null(opt$out) && nzchar(opt$out)) {
  dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
  sink(opt$out); on.exit(sink(), add = TRUE)
}

cat("[inspect_rds] file: ", opt$path, "\n", sep = "")
x <- readRDS(opt$path)
cat("class: ", paste(class(x), collapse = ", "), "\n", sep = "")
cat("type : ", typeof(x), "\n", sep = "")

if (is.data.frame(x)) {
  cat("data.frame: rows=", nrow(x), " cols=", ncol(x), "\n", sep = "")
  cat("names: ", paste(names(x), collapse = ", "), "\n", sep = "")
  print(utils::head(x, 5))
} else if (is.list(x)) {
  cat("list: length=", length(x), "\n", sep = "")
  cat("names: ", paste(names(x), collapse = ", "), "\n", sep = "")
  str(x, max.level = 1)
} else if (is.matrix(x)) {
  d <- dim(x); cat("matrix: dim=", d[1], "x", d[2], " type=", typeof(x), "\n", sep = "")
  print(x[1:min(5, nrow(x)), 1:min(5, ncol(x)), drop = FALSE])
} else {
  str(x, max.level = 1)
}

if (!is.null(opt$out) && nzchar(opt$out)) {
  cat("[OK] wrote: ", opt$out, "\n", sep = "")
}

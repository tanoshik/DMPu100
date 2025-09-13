# inspect_rds.R
# Purpose: dump lightweight structure of an RDS for debugging.
# ASCII comments only.

args <- commandArgs(trailingOnly = TRUE)

parse_opt <- function(args) {
  out <- list(path = NULL, out = NULL)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("--path", "-p")) {
      i <- i + 1L; out$path <- args[[i]]
    } else if (a %in% c("--out", "-o")) {
      i <- i + 1L; out$out <- args[[i]]
    }
    i <- i + 1L
  }
  if (is.null(out$path)) stop("--path is required")
  if (is.null(out$out))  stop("--out is required")
  out
}

describe_obj <- function(obj) {
  cat("class:", paste(class(obj), collapse = ","), "\n")
  if (is.data.frame(obj)) {
    cat("type: data.frame\n")
    cat("nrow:", nrow(obj), " ncol:", ncol(obj), "\n")
    cat("cols:\n")
    for (nm in names(obj)) {
      cat(" -", nm, ":", class(obj[[nm]])[1], "\n")
    }
  } else if (is.list(obj)) {
    cat("type: list\n")
    cat("names:", paste(names(obj), collapse = ","), "\n")
    if (all(c("sample_ids","locus_ids","A1","A2") %in% names(obj))) {
      cat("detected: indexed_db (sample_ids,locus_ids,A1,A2)\n")
      si <- obj$sample_ids; li <- obj$locus_ids
      a1 <- obj$A1; a2 <- obj$A2
      cat(" sample_ids:", length(si), "\n")
      cat(" locus_ids :", length(li), "\n")
      if (is.matrix(a1)) cat(" A1 dim     :", paste(dim(a1), collapse = " x "), "\n")
      if (is.matrix(a2)) cat(" A2 dim     :", paste(dim(a2), collapse = " x "), "\n")
    }
  } else if (is.vector(obj)) {
    cat("type: vector\n")
    cat("length:", length(obj), "\n")
  } else if (is.matrix(obj)) {
    cat("type: matrix\n")
    cat("dim:", paste(dim(obj), collapse = " x "), "\n")
  } else {
    cat("type: unknown\n")
  }
}

main <- function() {
  opt <- parse_opt(args)
  if (!file.exists(opt$path)) stop("file not found: ", opt$path)
  obj <- readRDS(opt$path)
  
  dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
  con <- file(opt$out, open = "wt")
  sink(con, type = "output")
  on.exit({
    sink(type = "output")
    close(con)
  }, add = TRUE)
  
  cat("RDS path:", opt$path, "\n")
  describe_obj(obj)
  cat("OK\n")
}

if (sys.nframe() == 0) main()

# scripts/devtools/inspect_rds.R
# ASCII-only. Inspect an RDS and emit meta.json (+ optional sample.csv).
# Default output dir is unified to output/exports/<family>/ when path matches known families.
# CLI:
#   Rscript scripts/devtools/inspect_rds.R --path data/locus_order.rds
#   Rscript scripts/devtools/inspect_rds.R data/locus_layout.rds
#   (override) --out output/exports/custom_dir

suppressWarnings(suppressMessages({
  library(jsonlite)
  library(utils)
}))

`%||%` <- function(a,b){ if (is.null(a) || !length(a)) b else a }

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  res <- list(path = NULL, out = NULL, sample_n = 50L)
  if (!length(a)) return(res)
  i <- 1L
  while (i <= length(a)) {
    if (a[i] %in% c("--path","-p")) { res$path <- a[i+1]; i <- i + 2L; next }
    if (a[i] %in% c("--out","-o"))  { res$out  <- a[i+1]; i <- i + 2L; next }
    if (a[i] %in% c("--n"))         { res$sample_n <- as.integer(a[i+1]); i <- i + 2L; next }
    if (is.null(res$path)) { res$path <- a[i]; i <- i + 1L; next }
    if (is.null(res$out))  { res$out  <- a[i]; i <- i + 1L; next }
    i <- i + 1L
  }
  res
}

infer_family_from_path <- function(path) {
  bn <- basename(path)
  if (grepl("^locus_order\\.rds$", bn, ignore.case = TRUE))  return("locus_order")
  if (grepl("^locus_layout\\.rds$", bn, ignore.case = TRUE)) return("locus_layout")
  if (grepl("^freq_table\\.rds$", bn, ignore.case = TRUE))   return("freq_table")
  if (grepl("^virtual_db_u100_", bn, ignore.case = TRUE))    return("virtual_db_u100")
  return(NULL)
}

inspect_rds <- function(path, out_dir = NULL, sample_n = 50L) {
  if (is.null(path) || !nzchar(path)) stop("--path is required")
  if (!file.exists(path)) stop("File not found: ", path)
  
  fam <- infer_family_from_path(path)
  stem <- if (is.null(fam)) tools::file_path_sans_ext(basename(path)) else fam
  default_out <- file.path("output","exports", stem)
  out_dir <- out_dir %||% default_out
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  obj <- readRDS(path)
  meta <- list(
    file   = normalizePath(path, winslash="/", mustWork=FALSE),
    class  = class(obj),
    type   = typeof(obj),
    length = tryCatch(length(obj), error=function(e) NA_integer_),
    str    = capture.output(utils::str(obj, max.level = 1))
  )
  
  if (is.data.frame(obj)) {
    meta$kind <- "data.frame"; meta$dim <- dim(obj); meta$cols <- names(obj)
    n <- min(nrow(obj), sample_n)
    utils::write.csv(obj[seq_len(n), , drop = FALSE],
                     file = file.path(out_dir, paste0(stem, "_inspect.sample.csv")),
                     row.names = FALSE)
  } else if (is.vector(obj) && is.character(obj)) {
    meta$kind <- "character_vector"
    meta$preview <- head(obj, n = min(length(obj), 30L))
    writeLines(obj, con = file.path(out_dir, paste0(stem, "_inspect.vector.txt")), useBytes = TRUE)
  } else if (is.list(obj)) {
    meta$kind <- "list"; meta$names <- names(obj)
  } else {
    meta$kind <- "unknown"
  }
  
  jsonlite::write_json(meta, file.path(out_dir, paste0(stem, "_inspect.meta.json")), pretty=TRUE, auto_unbox=TRUE)
  message("[INSPECT] Wrote: ", file.path(out_dir, paste0(stem, "_inspect.meta.json")))
  invisible(meta)
}

if (sys.nframe() == 0L) {
  opt <- parse_args()
  inspect_rds(opt$path, opt$out, opt$sample_n)
}

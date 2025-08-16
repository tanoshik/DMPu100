# scripts/devtools/make_dist_zip.R
# Build a distribution zip one level above the project root.
# Also writes "meta" folder (sibling to the zip) containing:
#  - EXCLUDED_LIST.txt : list of excluded paths
#  - excluded_heads/<path>.head.txt : up to head_lines lines for selected excluded files
#
# No multibyte characters in code/comments.

# Do not call normalizePath for printing; just show raw absolute paths.
safe_show <- function(path) {
  # best-effort absolute without mustWork
  p <- tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE),
                error = function(e) path)
  if (is.na(p) || p == "") path else p
}

make_dist_zip <- function(project_root = getwd(),
                          zip_prefix   = "DMP_dist",
                          parent_dir   = dirname(getwd()),  # one level above project root
                          head_lines   = 50,                # number of lines for head capture
                          max_head_kb  = 128,               # per-file cap for head capture
                          dry_run      = FALSE,
                          verbose      = TRUE) {
  olwd <- getwd()
  on.exit(setwd(olwd), add = TRUE)
  setwd(project_root)
  
  # ---- helpers ----
  tryQuiet <- function(expr) suppressWarnings(suppressMessages(try(expr, silent = TRUE)))
  
  ts_jst <- function() {
    old_tz <- Sys.getenv("TZ", unset = "")
    on.exit(Sys.setenv(TZ = old_tz), add = TRUE)
    Sys.setenv(TZ = "Asia/Tokyo")
    format(Sys.time(), "%Y%m%d%H%M")
  }
  
  short_hash <- function() {
    out <- tryQuiet(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = TRUE))
    if (inherits(out, "try-error") || length(out) == 0 || grepl("fatal", paste(out, collapse = " "))) {
      paste(sample(c(0:9, letters[1:6]), 7, replace = TRUE), collapse = "")
    } else trimws(out[[1]])
  }
  
  is_under <- function(path, dir_prefix) startsWith(path, paste0(dir_prefix, "/"))
  
  # text-like extension check for head capture
  is_text_like <- function(path) {
    ext <- tolower(tools::file_ext(path))
    ext %in% c(
      "r","rmd","qmd","md","txt","csv","tsv","json","yaml","yml",
      "log","gitignore","renvignore","rprofile","rproj"
    )
  }
  
  # ---- inventory ----
  all <- list.files(".", recursive = TRUE, all.files = FALSE, full.names = FALSE)
  # drop directories
  all <- all[file.info(all)$isdir == FALSE]
  
  # whitelist: include scripts/ and data/ (but filter), plus top-level docs
  must_include <- character(0)
  if (file.exists("README.md"))  must_include <- c(must_include, "README.md")
  if (file.exists("LICENSE"))    must_include <- c(must_include, "LICENSE")
  if (file.exists("DESCRIPTION")) must_include <- c(must_include, "DESCRIPTION")
  
  include_dirs <- c("scripts", "data")
  keep_seed <- Reduce(`|`, lapply(include_dirs, function(d) is_under(all, d))) |
    (basename(all) %in% basename(must_include))
  
  # data/ extension filter (no .rds)
  allow_ext_data <- function(path) {
    if (!is_under(path, "data")) return(TRUE)
    ext <- tolower(tools::file_ext(path))
    nchar(ext) > 0 && ext %in% c("csv","txt","tsv","md","rda","rdata","json","yaml","yml")
  }
  
  # exclude dirs anywhere in path
  exclude_dirs <- c(
    "output","test",".git",".github",".Rproj.user",".idea",".vscode",
    ".renv","renv","packrat","docs","dist","release","tmp","bench","notes","meta"
  )
  
  is_excluded_dir <- function(path) any(startsWith(path, paste0(exclude_dirs, "/")))
  
  # exclude patterns
  exclude_patterns <- c(
    "\\.rds$","\\.RDS$","\\.zip$","~$","^\\.",     # common
    "^debug_.*\\.R$","^scratch_.*\\.R$",
    "^.*[/\\\\]debug_.*\\.R$","^.*[/\\\\]scratch_.*\\.R$"
  )
  is_excluded_pat <- function(path) any(vapply(exclude_patterns, function(rx) grepl(rx, path, perl = TRUE), logical(1)))
  
  # build keep list
  candidates <- all[keep_seed]
  candidates <- candidates[vapply(candidates, allow_ext_data, logical(1))]
  keep <- candidates[!vapply(candidates, is_excluded_dir,  logical(1))]
  keep <- keep[!vapply(keep,      is_excluded_pat, logical(1))]
  
  # excluded = all - keep
  excluded <- setdiff(all, keep)
  
  if (length(keep) == 0) stop("No files to include. Check filters.", call. = FALSE)
  
  # ---- output paths (one level above) ----
  stamp <- ts_jst()
  shash <- short_hash()
  zip_name  <- sprintf("%s_%s_%s.zip", zip_prefix, shash, stamp)
  zip_path  <- file.path(parent_dir, zip_name)
  
  meta_dir  <- file.path(parent_dir, sprintf("meta_%s_%s", shash, stamp))
  heads_dir <- file.path(meta_dir, "excluded_heads")
  if (!dry_run) {
    dir.create(meta_dir,  recursive = TRUE, showWarnings = FALSE)
    dir.create(heads_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # ---- meta: write excluded list and heads ----
  write_meta <- function() {
    # EXCLUDED_LIST.txt
    list_path <- file.path(meta_dir, "EXCLUDED_LIST.txt")
    cat(sprintf("# excluded files: %d\n", length(excluded)), file = list_path)
    if (length(excluded)) {
      cat(paste(excluded, collapse = "\n"), file = list_path, append = TRUE)
      cat("\n", file = list_path, append = TRUE)
    }
    
    # heads for text-like smallish files
    for (p in excluded) {
      if (!file.exists(p)) next
      if (!is_text_like(p)) next
      info <- file.info(p)
      if (is.na(info$size) || info$size > (max_head_kb * 1024)) next
      
      rel_out <- paste0(p, ".head.txt")
      rel_out <- gsub("[/\\\\]+", "_", rel_out)  # sanitize
      out_path <- file.path(heads_dir, rel_out)
      
      con <- file(p, open = "rt", encoding = "UTF-8")
      on.exit(tryQuiet(close(con)), add = TRUE)
      lines <- readLines(con, n = head_lines, warn = FALSE)
      tryQuiet(close(con))  # close immediately, regardless of loop
      cat(sprintf("# HEAD of excluded file: %s\n", p), file = out_path)
      if (length(lines)) cat(paste(lines, collapse = "\n"), file = out_path, append = TRUE)
    }
  }
  
  # PREVIEW
  if (verbose) {
    cat("== DIST PREVIEW ==\n")
    cat(sprintf("* project_root: %s\n", safe_show(project_root)))
    cat(sprintf("* parent_dir  : %s\n", safe_show(parent_dir)))
    cat(sprintf("* files kept  : %d\n", length(keep)))
    cat(sprintf("* files excl. : %d\n", length(excluded)))
    cat(sprintf("* zip target  : %s\n", safe_show(zip_path)), "\n")
  }
  
  if (dry_run) {
    cat("\n-- dry_run=TRUE; no zip/meta created. Preview lists follow --\n")
    cat(">> KEEP (head):\n"); print(utils::head(keep, 10))
    cat(">> EXCLUDED (head):\n"); print(utils::head(excluded, 10))
    return(invisible(list(keep = keep, excluded = excluded, zip = zip_path)))
  }
  
  # write meta before zipping
  write_meta()
  
  # ---- zip build ----
  zip_bin <- Sys.which("zip")
  if (nzchar(zip_bin)) {
    zres <- tryQuiet(utils::zip(zipfile = zip_path, files = keep, flags = "-q"))
    if (inherits(zres, "try-error")) stop("Failed to create zip via system zip.", call. = FALSE)
  } else if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(zipfile = zip_path, files = keep, include_directories = TRUE)
  } else {
    stop("No system 'zip' found and package {zip} not installed.", call. = FALSE)
  }
  
  if (!file.exists(zip_path)) stop("Zip not created; unknown error.", call. = FALSE)
  
  if (verbose) {
    sz_kb <- tryCatch(file.info(zip_path)$size / 1024, error = function(e) NA_real_)
    cat(sprintf(
      "== DONE ==\n* ZIP : %s (%.1f KB)\n* META: %s\n",
      safe_show(zip_path),
      sz_kb,
      safe_show(meta_dir)
    ))
  }
  
  invisible(list(zip = zip_path, meta = meta_dir, keep = keep, excluded = excluded))
}

# Run when invoked directly (only when sourced as a script)
if (sys.nframe() == 0) {
  make_dist_zip()
}

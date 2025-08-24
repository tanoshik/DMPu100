# scripts/devtools/make_dist_zip.R
# Build a distribution zip one level above the project root.
# Also writes "meta" folder (sibling to the zip) containing:
#  - EXCLUDED_LIST.txt : list of excluded paths
#  - excluded_heads/<path>.head.txt : up to head_lines lines for selected excluded files
#
# No multibyte characters in code/comments.

# --- tiny helpers -------------------------------------------------------------
safe_show <- function(path) {
  p <- tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE),
                error = function(e) path)
  if (is.na(p) || p == "") path else p
}

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

# normalize path separators to '/' before prefix judgment
.normalize_slashes <- function(x) gsub("[\\\\]+", "/", x)

is_under <- function(path, dir_prefix) {
  path <- .normalize_slashes(path)
  dir_prefix <- .normalize_slashes(dir_prefix)
  startsWith(path, paste0(dir_prefix, "/"))
}

is_text_like <- function(path) {
  ext <- tolower(tools::file_ext(path))
  ext %in% c(
    "r","rmd","qmd","md","txt","csv","tsv","json","yaml","yml","rdata","rda",
    "html","css","js","xml","ini","cfg","conf","toml","sh","bat","ps1","py"
  )
}

# --- main ---------------------------------------------------------------------
make_dist_zip <- function(
    project_root      = getwd(),
    zip_prefix        = "DMP_dist",
    parent_dir        = dirname(getwd()),
    head_lines        = 50,
    max_head_kb       = 128,
    include_dirs      = c("scripts", "data", "test", "output"),  # ★ test/output を既定で含める
    include_data_rds  = FALSE,
    extra_include     = character(0),
    extra_exclude_dirs    = character(0),
    extra_exclude_patterns= character(0),
    data_max_bytes        = 1024 * 1024,   # ★ data/ の閾値（1MB）
    stage_large_data_as_empty = TRUE,      # ★ 閾値超は中身なしファイルで出力
    dry_run           = FALSE,
    verbose           = TRUE
){
  olwd <- getwd()
  on.exit(setwd(olwd), add = TRUE)
  setwd(project_root)
  
  # ---- inventory ----
  all <- list.files(".", recursive = TRUE, all.files = FALSE, full.names = FALSE)
  finfo <- file.info(all)
  all   <- all[!finfo$isdir]                 # files only
  all   <- .normalize_slashes(all)
  
  # top-level docs
  must_include <- character(0)
  for (f in c("README.md","LICENSE","DESCRIPTION")) {
    if (file.exists(f)) must_include <- c(must_include, f)
  }
  
  # 初期シード（scripts, data, test, output + トップドキュメント）
  seed_dirs <- unique(.normalize_slashes(include_dirs))
  keep_seed <- rep(FALSE, length(all))
  for (d in seed_dirs) {
    keep_seed <- keep_seed | startsWith(all, paste0(.normalize_slashes(d), "/"))
  }
  if (length(must_include)) {
    keep_seed <- keep_seed | (all %in% .normalize_slashes(must_include))
  }
  if (length(extra_include)) {
    # files or directories
    # if directory: include all under it
    add <- unique(unlist(lapply(extra_include, function(p) {
      p <- .normalize_slashes(p)
      if (dir.exists(p)) {
        all[startsWith(all, paste0(p, "/"))]
      } else {
        p
      }
    }), use.names = FALSE))
    keep_seed <- keep_seed | (all %in% .normalize_slashes(add))
  }
  
  # data/ の拡張子フィルタ（RDSは既定では除外）
  allow_ext_data <- function(path) {
    if (!is_under(path, "data")) return(TRUE)
    ext <- tolower(tools::file_ext(path))
    if (include_data_rds) {
      nchar(ext) > 0 && ext %in% c("csv","txt","tsv","md","rda","rdata","json","yaml","yml","rds","rds.gz")
    } else {
      nchar(ext) > 0 && ext %in% c("csv","txt","tsv","md","rda","rdata","json","yaml","yml")
    }
  }
  
  # 除外ディレクトリ（bench/test/output を除外しないように修正）
  exclude_dirs <- c(
    ".git",".github",".Rproj.user",".idea",".vscode",
    ".renv","renv","packrat","docs","dist","release","tmp","notes","meta"
  )
  if (length(extra_exclude_dirs)) exclude_dirs <- unique(c(exclude_dirs, extra_exclude_dirs))
  
  is_excluded_dir <- function(path) {
    any(vapply(exclude_dirs, function(d) startsWith(path, paste0(d, "/")), logical(1)))
  }
  
  # 除外パターン
  exclude_patterns <- c(
    "\\.rds$","\\.RDS$","\\.zip$","~$","^\\.",     # common
    "^debug_.*\\.R$","^scratch_.*\\.R$",
    "^.*[/\\\\]debug_.*\\.R$","^.*[/\\\\]scratch_.*\\.R$"
  )
  if (include_data_rds) {
    exclude_patterns <- setdiff(exclude_patterns, c("\\.rds$","\\.RDS$"))
    exclude_patterns <- c(exclude_patterns, "^(?!data/).+\\.rds$")  # non-data rds excluded
  }
  if (length(extra_exclude_patterns)) exclude_patterns <- unique(c(exclude_patterns, extra_exclude_patterns))
  
  is_excluded_pat <- function(path) any(vapply(exclude_patterns, function(rx) grepl(rx, path, perl = TRUE), logical(1)))
  
  # build keep list
  candidates <- all[keep_seed]
  candidates <- candidates[vapply(candidates, allow_ext_data, logical(1))]
  keep <- candidates[!vapply(candidates, is_excluded_dir,  logical(1))]
  keep <- keep[!vapply(keep,      is_excluded_pat, logical(1))]
  keep <- sort(unique(keep))  # stable
  
  # excluded = all - keep
  excluded <- setdiff(all, keep)
  excluded <- sort(unique(excluded))
  
  # meta directory paths
  out_hash <- short_hash()
  out_ts   <- ts_jst()
  zip_name <- sprintf("%s_%s_%s.zip", zip_prefix, out_hash, out_ts)
  zip_path <- file.path(parent_dir, zip_name)
  meta_dir <- file.path(parent_dir, sprintf("meta_%s_%s", out_hash, out_ts))
  heads_dir<- file.path(meta_dir, "excluded_heads")
  
  # meta writer
  write_meta <- function() {
    dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
    # excluded list
    cat(paste(excluded, collapse = "\n"), file = file.path(meta_dir, "EXCLUDED_LIST.txt"))
    
    # heads for text-like smallish files
    for (p in excluded) {
      if (!file.exists(p)) next
      if (!is_text_like(p)) next
      info <- file.info(p)
      if (is.na(info$size) || info$size > (max_head_kb * 1024)) next
      
      rel_out  <- paste0(p, ".head.txt")
      rel_out  <- gsub("[/\\\\]+", "_", rel_out)  # sanitize
      out_path <- file.path(heads_dir, rel_out)
      
      con <- file(p, open = "rt", encoding = "UTF-8")
      on.exit(tryQuiet(close(con)), add = TRUE)
      lines <- tryCatch(readLines(con, n = head_lines, warn = FALSE), error = function(e) character(0))
      tryQuiet(close(con))
      on.exit(NULL, add = FALSE)
      
      dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
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
  }
  
  if (dry_run) {
    write_meta()
    if (verbose) {
      cat("== DRY-RUN ONLY ==\n")
      cat(sprintf("ZIP would be written to: %s\nMETA would be written to: %s\n", safe_show(zip_path), safe_show(meta_dir)))
    }
    return(invisible(list(keep = keep, excluded = excluded, zip = zip_path)))
  }
  
  # write meta before building staging/zip
  write_meta()
  
  # ---- staging (to handle large data files as empty) ----
  staging_dir <- file.path(tempdir(), sprintf("dmp_dist_stage_%s_%s", out_hash, out_ts))
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  
  # copy with special handling for data/*
  for (p in keep) {
    src  <- p
    dest <- file.path(staging_dir, p)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    
    if (is_under(p, "data") && stage_large_data_as_empty) {
      info <- file.info(src)
      if (!is.na(info$size) && info$size > data_max_bytes) {
        # create empty file with same name (ファイル名だけ)
        file.create(dest, showWarnings = FALSE)
        next
      }
    }
    file.copy(src, dest, overwrite = TRUE, copy.date = TRUE, copy.mode = TRUE)
  }
  
  # ---- zip build (from staging) ----
  owd2 <- getwd()
  on.exit(setwd(owd2), add = TRUE)
  setwd(staging_dir)
  
  zip_bin <- Sys.which("zip")
  if (nzchar(zip_bin)) {
    zres <- tryQuiet(utils::zip(zipfile = zip_path, files = keep, flags = "-q"))
    if (inherits(zres, "try-error")) stop("Failed to create zip via system zip.", call. = FALSE)
  } else if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(zipfile = zip_path, files = keep, include_directories = TRUE)
  } else {
    stop("No system 'zip' found and package {zip} not installed.", call. = FALSE)
  }
  
  # size
  if (file.exists(zip_path)) {
    sz_kb <- round(file.info(zip_path)$size / 1024, 1)
    if (verbose) cat(sprintf("== DONE ==\n* ZIP : %s (%.1f KB)\n* META: %s\n", safe_show(zip_path), sz_kb, safe_show(meta_dir)))
  } else {
    warning("Zip not found after build.")
  }
  
  invisible(list(zip = zip_path, meta = meta_dir, keep = keep, excluded = excluded))
}

# Run when invoked directly (only when sourced as a script)
if (sys.nframe() == 0) {
  make_dist_zip()
}

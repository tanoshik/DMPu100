# scripts/devtools/make_dist_zip.R
# Build a distribution zip one level above the project root.
# Includes data/, scripts/bench/, bench/, test/ (incl. test/bench), output/
# For data/: files larger than 1MB are replaced by empty stub files (same name).
# Writes meta folder with EXCLUDED_LIST.txt and LARGE_FILES.txt.

# ---- helpers ---------------------------------------------------------------

bytes_1mb <- 1024L * 1024L

safe_show <- function(path) {
  p <- tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE),
                error = function(e) path)
  if (is.na(p) || p == "") path else p
}

is_subpath <- function(p, root) {
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  startsWith(paste0(p, "/"), paste0(root, "/")) || identical(p, root)
}

list_files_recursive <- function(dir) {
  if (!dir.exists(dir)) return(character())
  files <- list.files(dir, recursive = TRUE, all.files = TRUE, full.names = TRUE)
  files[file.info(files)$isdir %in% FALSE]
}

# copy with parent dirs, create dirs if needed
copy_with_dirs <- function(src, dst_root, rel) {
  dst <- file.path(dst_root, rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dst, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  if (!ok) stop("Failed to copy: ", src, " -> ", dst)
  dst
}

# create empty stub file at dst_root/rel
create_empty_stub <- function(dst_root, rel) {
  dst <- file.path(dst_root, rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  con <- file(dst, open = "wb")
  close(con)
  dst
}

# ---- main ------------------------------------------------------------------

make_dist_zip <- function(project_root = getwd(),
                          zip_prefix   = "DMP_dist",
                          parent_dir   = dirname(getwd()),
                          head_lines   = 50,     # kept for compatibility (not used here)
                          data_large_threshold = bytes_1mb) {
  
  project_root <- safe_show(project_root)
  parent_dir   <- safe_show(parent_dir)
  
  # Staging dir for exact zip contents
  hash_stub <- substr(as.character(as.integer(as.numeric(Sys.time()))), 1, 8)
  ts <- format(Sys.time(), "%Y%m%d%H%M")
  zip_base <- sprintf("%s_%s_%s", zip_prefix, hash_stub, ts)
  staging <- file.path(tempdir(), paste0("staging_", zip_base))
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  
  # meta outputs (sibling to the zip file)
  meta_dir <- file.path(parent_dir, paste0("meta_", hash_stub, "_", ts))
  dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)
  excluded_list_path <- file.path(meta_dir, "EXCLUDED_LIST.txt")
  large_list_path    <- file.path(meta_dir, "LARGE_FILES.txt")
  
  # Defaults: exclude typical dev folders/files; we will add back required ones explicitly
  default_exclude_dirs <- c(".git", ".Rproj.user", "dist", "release")
  default_exclude_glob <- c("*.zip")  # keep as before
  
  # Always include these top-level paths (relative to project_root)
  include_dirs_exact <- c(
    "scripts",          # whole scripts/ (legacy behavior)
    "data",             # special handling for >1MB
    "scripts/bench",
    "bench",
    "test",
    "output"
  )
  
  # Build candidate file list from project_root
  all_files <- list_files_recursive(project_root)
  
  # Filter out default excludes first (but not our explicit include dirs)
  rel_all <- gsub(paste0("^", gsub("\\\\", "/", project_root), "/?"), "", gsub("\\\\", "/", all_files))
  keep <- rep(TRUE, length(all_files))
  
  # exclude by directory unless inside our explicit include roots
  for (i in seq_along(all_files)) {
    rel <- rel_all[i]
    top <- strsplit(rel, "/")[[1]]
    top <- if (length(top)) top[1] else ""
    # If path is under an explicit include root, we keep it (handle later)
    under_include <- any(startsWith(rel, paste0(include_dirs_exact, "/")) | rel %in% include_dirs_exact)
    if (!under_include) {
      if (top %in% default_exclude_dirs) keep[i] <- FALSE
      # glob excludes
      if (keep[i] && length(default_exclude_glob)) {
        for (g in default_exclude_glob) {
          if (grepl(glob2rx(g), basename(rel), ignore.case = TRUE)) {
            keep[i] <- FALSE
            break
          }
        }
      }
    }
  }
  
  all_files <- all_files[keep]
  rel_all   <- rel_all[keep]
  
  # From the kept list, now compose the final include set:
  # 1) Everything under scripts/ (legacy)
  # 2) data/ (with >1MB stubbing)
  # 3) scripts/bench/, bench/, test/, output/ (unrestricted)
  is_under <- function(rel, root) {
    rel == root || startsWith(rel, paste0(root, "/"))
  }
  
  want <- logical(length(all_files))
  for (j in seq_along(all_files)) {
    rel <- rel_all[j]
    if (is_under(rel, "scripts") ||
        is_under(rel, "data") ||
        is_under(rel, "scripts/bench") ||
        is_under(rel, "bench") ||
        is_under(rel, "test") ||
        is_under(rel, "output")) {
      want[j] <- TRUE
    }
  }
  
  # Additionally: include root-level files like README.md, LICENSE, DESCRIPTION if present
  root_level_keep <- c("README.md", "LICENSE", "DESCRIPTION")
  root_paths <- file.path(project_root, root_level_keep)
  root_exists <- file.exists(root_paths)
  root_files <- root_paths[root_exists]
  root_rel   <- basename(root_files)
  
  final_src  <- c(all_files[want], root_files)
  final_rel  <- c(rel_all[want],   root_rel)
  
  # Deduplicate
  o <- order(final_rel)
  final_src <- final_src[o]
  final_rel <- final_rel[o]
  dupe <- duplicated(final_rel)
  if (any(dupe)) {
    final_src <- final_src[!dupe]
    final_rel <- final_rel[!dupe]
  }
  
  # Prepare meta logs
  writeLines(character(), con = excluded_list_path)
  writeLines(c(
    "# Files larger than threshold are stubbed as empty files in the zip.",
    paste0("# Threshold: ", data_large_threshold, " bytes"),
    "# Format: <size_bytes> <relative_path_from_project_root>",
    ""
  ), con = large_list_path)
  
  # Copy into staging, with special handling for data/ large files
  info <- file.info(final_src, extra_cols = FALSE)
  for (k in seq_along(final_src)) {
    src <- final_src[k]
    rel <- final_rel[k]
    
    # Determine if this is under data/
    if (is_under(rel, "data")) {
      sz <- info$size[k]
      if (!is.na(sz) && sz > data_large_threshold) {
        # create empty stub
        create_empty_stub(staging, rel)
        # record large file entry
        cat(sprintf("%d %s\n", as.integer(sz), rel), file = large_list_path, append = TRUE)
        next
      }
    }
    copy_with_dirs(src, staging, rel)
  }
  
  # Build list of excluded files (for transparency)
  # Anything under project_root not copied to staging and not a dir becomes "excluded"
  staged_all <- list.files(staging, recursive = TRUE, all.files = TRUE, full.names = FALSE)
  # filter out dirs in staged listing by mapping back to file.info on staging
  staged_full <- file.path(staging, staged_all)
  staged_all <- staged_all[file.info(staged_full)$isdir %in% FALSE]
  
  # all relative files from original
  all_project_files <- list_files_recursive(project_root)
  all_project_rel <- gsub(paste0("^", gsub("\\\\", "/", project_root), "/?"), "", gsub("\\\\", "/", all_project_files))
  all_project_rel <- all_project_rel[file.info(all_project_files)$isdir %in% FALSE]
  
  excluded_rel <- setdiff(all_project_rel, staged_all)
  if (length(excluded_rel)) {
    writeLines(sort(excluded_rel), con = excluded_list_path, useBytes = TRUE)
  }
  
  # Create zip in parent_dir
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(staging)
  zip_path <- file.path(parent_dir, paste0(zip_base, ".zip"))
  
  # Use utils::zip (platform independent)
  # Include everything under staging (relative paths)
  files_to_zip <- list.files(".", recursive = TRUE, all.files = TRUE, full.names = FALSE)
  utils::zip(zipfile = zip_path, files = files_to_zip)
  
  message("[OK] Wrote zip: ", safe_show(zip_path))
  message("[OK] Wrote meta: ", safe_show(meta_dir))
  invisible(list(zip = zip_path, meta = meta_dir, staging = staging))
}

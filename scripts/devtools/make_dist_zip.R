# scripts/devtools/make_dist_zip.R
# Build a distribution zip one level above the project root.
# Adds data/ (with >1MB files as zero-byte placeholders), scripts/bench/, test/ (incl test/bench), output/
# No multibyte characters in code/comments.

safe_show <- function(path) {
  p <- tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE),
                error = function(e) path)
  if (is.na(p) || p == "") path else p
}

make_dist_zip <- function(project_root = getwd(),
                          zip_prefix   = "DMP_dist",
                          parent_dir   = dirname(getwd()),  # one level above project root
                          head_lines   = 50,                # number of lines for head capture
                          small_data_threshold_bytes = 1024 * 1024  # 1MB
) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)
  
  timestamp <- format(Sys.time(), "%Y%m%d%H%M")
  # Use short git hash if available
  git_hash <- tryCatch({
    x <- system2("git", c("rev-parse", "--short=7", "HEAD"), stdout = TRUE, stderr = FALSE)
    if (length(x) == 0) "nohash" else gsub("\\s+", "", x[1])
  }, error = function(e) "nohash")
  
  zip_name <- sprintf("%s_%s_%s.zip", zip_prefix, git_hash, timestamp)
  zip_path <- file.path(parent_dir, zip_name)
  
  meta_dir     <- file.path(parent_dir, sprintf("meta_%s_%s", git_hash, timestamp))
  excluded_dir <- file.path(meta_dir, "excluded_heads")
  if (!dir.exists(meta_dir)) dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(excluded_dir)) dir.create(excluded_dir, recursive = TRUE, showWarnings = FALSE)
  
  excluded_list_path <- file.path(meta_dir, "EXCLUDED_LIST.txt")
  writeLines(character(0), excluded_list_path)  # reset
  
  # ---- base includes (existing behavior) ----
  # Keep your current include/exclude logic; then add the four paths requested.
  # Example of typical base excludes:
  base_exclude_dirs <- c(".git", ".Rproj.user", "dist", "release", "notes",
                         "bench",  # historically excluded; we will re-include scripts/bench explicitly
                         "output", # historically excluded; we will include explicitly below
                         "test",   # historically excluded; we will include explicitly below
                         "meta")   # meta is generated alongside zip
  base_exclude_globs <- c("*.zip", "*.rds")    # historical; we will re-include where requested
  base_include <- c("scripts", "data", "README.md", "LICENSE", "DESCRIPTION")
  
  # Gather all files under project_root
  all_files <- list.files(".", recursive = TRUE, all.files = FALSE, include.dirs = FALSE, no.. = TRUE)
  # Filter out meta output paths
  all_files <- all_files[!grepl(sprintf("^%s", gsub("\\.", "\\\\.", basename(meta_dir))), all_files)]
  
  # Start from base include, then subtract base excludes
  want <- all_files[
    startsWith(all_files, "scripts") |
      startsWith(all_files, "data")    |
      all_files %in% c("README.md", "LICENSE", "DESCRIPTION")
  ]
  
  # Remove base excludes
  drop_dir_pat <- paste0("^(", paste(gsub("\\.", "\\\\.", base_exclude_dirs), collapse = "|"), ")(/|$)")
  want <- want[!grepl(drop_dir_pat, want)]
  for (g in base_exclude_globs) {
    want <- want[!grepl(glob2rx(g), basename(want))]
  }
  
  # ---- explicit additions per request ----
  must_add_dirs <- c("scripts/bench", "test", "output")
  for (ad in must_add_dirs) {
    if (dir.exists(ad)) {
      add_files <- list.files(ad, recursive = TRUE, all.files = FALSE, include.dirs = FALSE, no.. = TRUE)
      add_files <- file.path(ad, add_files)
      # ensure no duplicates
      want <- unique(c(want, add_files))
    }
  }
  
  # ---- data/ handling: small vs large ----
  # Small (<= threshold) : include as-is
  # Large (> threshold)  : replace with zero-byte placeholder (same relative path) and log to EXCLUDED_LIST
  data_small <- character(0)
  data_large <- character(0)
  if (dir.exists("data")) {
    data_files <- list.files("data", recursive = TRUE, all.files = FALSE, include.dirs = FALSE, no.. = TRUE)
    data_files <- file.path("data", data_files)
    sizes <- file.info(data_files)$size
    data_small <- data_files[which(sizes <= small_data_threshold_bytes)]
    data_large <- data_files[which(sizes >  small_data_threshold_bytes)]
  }
  
  # Ensure all small data files are included
  want <- unique(c(want, data_small))
  
  # Stage placeholders for large data files into a temp staging dir
  staging_root <- file.path(tempdir(), sprintf("dmp_dist_staging_%s_%s", git_hash, timestamp))
  if (!dir.exists(staging_root)) dir.create(staging_root, recursive = TRUE, showWarnings = FALSE)
  
  stage_paths <- character(0)
  if (length(data_large) > 0) {
    for (rel in data_large) {
      stage_abs <- file.path(staging_root, rel)
      dir.create(dirname(stage_abs), recursive = TRUE, showWarnings = FALSE)
      file.create(stage_abs)  # zero-byte placeholder
      stage_paths <- c(stage_paths, stage_abs)
    }
    # Log large files to meta/EXCLUDED_LIST.txt
    log_lines <- sprintf("[data-large] %s\t(size=%s bytes)",
                         vapply(data_large, safe_show, character(1)),
                         format(file.info(data_large)$size, scientific = FALSE))
    write(log_lines, file = excluded_list_path, append = TRUE)
  }
  
  # Build final include list:
  #  - project-root-relative paths in 'want'
  #  - absolute staged paths in 'stage_paths' (we will zip them with the same relative layout)
  # Zip needs relative paths; so we will zip in two steps:
  #  (A) from project_root for 'want'
  #  (B) from staging_root for 'stage_paths'
  # Then merge into single zip (utils::zip cannot append; we recreate with combined file list via a temp tree)
  
  # Create a final, unified staging tree to zip once
  final_stage <- file.path(tempdir(), sprintf("dmp_final_stage_%s_%s", git_hash, timestamp))
  if (dir.exists(final_stage)) unlink(final_stage, recursive = TRUE, force = TRUE)
  dir.create(final_stage, recursive = TRUE, showWarnings = FALSE)
  
  # Copy regular includes
  for (rel in want) {
    if (!file.exists(rel)) next
    dst <- file.path(final_stage, rel)
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(rel, dst, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    if (!ok) {
      write(sprintf("[copy-fail] %s", safe_show(rel)), file = excluded_list_path, append = TRUE)
    }
  }
  
  # Copy staged placeholders for large data
  if (length(stage_paths) > 0) {
    for (abs_p in stage_paths) {
      rel <- substring(abs_p, nchar(staging_root) + 2L)  # relative to staging_root
      dst <- file.path(final_stage, rel)
      dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
      file.copy(abs_p, dst, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    }
  }
  
  # Write meta heads for selected excluded files (optional: disabled for simplicity here)
  # You can add head capture for excluded texts if needed using 'head_lines'.
  
  # Create the zip from the final staging tree
  old2 <- getwd()
  setwd(final_stage)
  on.exit(setwd(old2), add = TRUE)
  files_to_zip <- list.files(".", recursive = TRUE, all.files = FALSE, include.dirs = FALSE, no.. = TRUE)
  if (file.exists(zip_path)) file.remove(zip_path)
  utils::zip(zipfile = zip_path, files = files_to_zip, flags = "-r9X")
  
  # Write meta README
  meta_readme <- file.path(meta_dir, "README_meta.txt")
  writeLines(c(
    "This folder contains auxiliary information for the distribution zip.",
    sprintf("Zip: %s", safe_show(zip_path)),
    sprintf("Project root: %s", safe_show(project_root)),
    sprintf("Created at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Git hash: %s", git_hash),
    sprintf("Large data threshold: %d bytes", small_data_threshold_bytes),
    "Large data files were replaced with zero-byte placeholders in the zip.",
    "See EXCLUDED_LIST.txt for the list of large data files."
  ), meta_readme)
  
  message(sprintf("[OK] Wrote zip: %s", safe_show(zip_path)))
  message(sprintf("[OK] Meta dir: %s", safe_show(meta_dir)))
  invisible(list(zip = zip_path, meta = meta_dir))
}

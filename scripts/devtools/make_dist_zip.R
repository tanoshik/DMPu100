# scripts/devtools/make_dist_zip.R
# DMPu100 -> GPT sharing bundle
# Rules:
# - Keep directory tree as-is (no flatten)
# - Exclude ONLY root-level: .git/ , .Rproj.user/ , .Rhistory
# - Exclude any file > 2MB (all extensions)
# - Include output/ by default; allow blacklist of subfolders
# - Record excluded (path,size,reason) to meta/EXCLUDED_WITH_SIZE.csv
# - Zip name: DMP_dist_<sha>_<YYYYMMDDHHMM>.zip (written to parent of project_root)
# - Comments: ASCII only in code

# ---------- helpers ----------
..norm <- function(x) gsub("\\\\","/",x)
..now  <- function() format(Sys.time(), "%Y%m%d%H%M")
..git  <- function(cmd) tryCatch(system(paste("git", cmd), intern=TRUE, ignore.stderr=TRUE), error=function(e) character(0))
..sha  <- function(){ s <- ..git("rev-parse --short HEAD"); if (length(s)==0) "nogit" else s[1] }
..branch <- function(){ s <- ..git("rev-parse --abbrev-ref HEAD"); if (length(s)==0) "nogit" else s[1] }

..zip_with_fallback <- function(zipfile, files, root=".", quiet=TRUE){
  o <- getwd(); setwd(root); on.exit(setwd(o), add=TRUE)
  ok <- FALSE
  if (requireNamespace("zip", quietly=TRUE)) {
    zip::zipr(zipfile, files = files, include_directories = TRUE)
    ok <- TRUE
  } else {
    flags <- if (quiet) "-r9Xq" else "-r9X"
    z <- try(utils::zip(zipfile = zipfile, files = files, flags = flags), silent = quiet)
    if (!inherits(z,"try-error")) ok <- TRUE
  }
  if (!ok) stop("zip backend not available. install.packages('zip') is recommended.")
}

..safe_copy <- function(src, dst){
  dir.create(dirname(dst), recursive=TRUE, showWarnings=FALSE)
  file.copy(src, dst, overwrite=TRUE, copy.mode=TRUE, copy.date=TRUE)
}

# ---------- main ----------
make_dist_zip <- function(
    project_root           = getwd(),
    parent_dir             = dirname(getwd()),
    zip_prefix             = "DMP_dist",
    # output/ blacklist: subfolder names to EXCLUDE. default: include all
    output_exclude_subdirs = character(0),
    # global size guard (MB) applied once to ALL files
    max_file_mb_any        = 2,
    allow_large_any        = FALSE
){
  # maintain working dir
  owd <- setwd(project_root); on.exit(setwd(owd), add=TRUE)
  
  # names
  hash  <- ..sha(); stamp <- ..now()
  zname <- sprintf("%s_%s_%s.zip", zip_prefix, hash, stamp)
  zpath <- ..norm(file.path(parent_dir, zname))
  
  # enumerate all files (include hidden)
  all_abs <- list.files(".", recursive=TRUE, all.files=TRUE,
                        include.dirs=FALSE, no..=TRUE, full.names=TRUE)
  all_abs <- ..norm(all_abs)
  
  # to root-relative paths for clean filtering
  root_norm <- ..norm(normalizePath(project_root, winslash="/", mustWork=TRUE))
  rel_from_root <- function(p) {
    p_norm <- ..norm(normalizePath(p, winslash="/", mustWork=FALSE))
    sub(paste0("^", root_norm, "/?"), "", p_norm)
  }
  all_rel <- rel_from_root(all_abs)
  
  # root-level excludes: .git/ , .Rproj.user/ , .Rhistory (only at repo root)
  is_root_git        <- grepl("^\\.git/", all_rel)
  is_root_rproj_user <- grepl("^\\.Rproj\\.user/", all_rel)
  is_root_rhistory   <- grepl("^\\.Rhistory$", all_rel)
  
  # output blacklist
  drop_output <- rep(FALSE, length(all_rel))
  if (dir.exists("output") && length(output_exclude_subdirs)){
    # build pattern for subfolders under output/
    # match: output/<name>/...
    names_sanit <- gsub("[^A-Za-z0-9_\\-]", ".", output_exclude_subdirs)
    pat <- paste0("^output/(", paste(names_sanit, collapse="|"), ")(/|$)")
    drop_output <- grepl(pat, all_rel)
  }
  
  # initial candidates: remove above excludes
  excl_mask <- (is_root_git | is_root_rproj_user | is_root_rhistory | drop_output)
  cand_rel  <- all_rel[!excl_mask]
  cand_abs  <- all_abs[!excl_mask]
  
  # size guard (once, to all file types)
  large_skip_rel <- character(0)
  large_skip_abs <- character(0)
  if (!isTRUE(allow_large_any) && is.finite(max_file_mb_any) && max_file_mb_any > 0){
    info <- file.info(cand_abs)
    sz <- as.numeric(info$size)
    thr <- max_file_mb_any * 1024^2
    too_big <- which(!is.na(sz) & is.finite(sz) & (sz > thr))    
    if (length(too_big)){
      large_skip_rel <- cand_rel[too_big]
      large_skip_abs <- cand_abs[too_big]
      keep <- setdiff(seq_along(cand_abs), too_big)
      cand_rel <- cand_rel[keep]
      cand_abs <- cand_abs[keep]
    }
  }
  
  # sanity
  if (!length(cand_abs)) stop("No files to bundle. Adjust blacklist or size limit.")
  
  # staging: copy preserving tree
  stage <- file.path(tempdir(), sprintf("stage_%s_%s", hash, stamp))
  dir.create(stage, recursive=TRUE, showWarnings=FALSE)
  for (i in seq_along(cand_abs)){
    ..safe_copy(cand_abs[i], file.path(stage, cand_rel[i]))
  }
  
  # meta dir
  meta_dir <- file.path(stage, sprintf("meta_%s_%s", hash, stamp))
  dir.create(meta_dir, recursive=TRUE, showWarnings=FALSE)
  
  # ABOUT
  writeLines(c(
    "DMPu100 bundle for GPT sharing.",
    paste0("branch: ", ..branch()),
    paste0("git_sha: ", ..sha()),
    paste0("timestamp: ", stamp),
    paste0("max_file_mb_any: ", max_file_mb_any),
    paste0("allow_large_any: ", allow_large_any),
    paste0("output_exclude_subdirs: ",
           if (!length(output_exclude_subdirs)) "NONE" else paste(output_exclude_subdirs, collapse="|"))
  ), file.path(meta_dir, "ABOUT.txt"))
  
  # INCLUDED_LIST
  writeLines(cand_rel, file.path(meta_dir, "INCLUDED_LIST.txt"))
  
  # FILE_SIZES
  sizes <- file.info(cand_abs)$size
  tab <- data.frame(path=cand_rel,
                    bytes=as.numeric(sizes),
                    ext=tolower(sub(".*\\.(\\w+)$","\\1", cand_rel)),
                    stringsAsFactors=FALSE)
  utils::write.csv(tab, file.path(meta_dir, "FILE_SIZES.csv"), row.names=FALSE)
  
  # EXCLUDED_WITH_SIZE (path, bytes, reason)
  excl_rows <- list()
  
  # root-level excludes
  if (any(is_root_git)){
    e_abs <- all_abs[is_root_git]; fi <- file.info(e_abs); fi <- fi[!is.na(fi$size), , drop=FALSE]
    if (nrow(fi)) excl_rows[[length(excl_rows)+1]] <- data.frame(
      path=rel_from_root(rownames(fi)), bytes=as.numeric(fi$size), reason=".git_root", stringsAsFactors=FALSE)
  }
  if (any(is_root_rproj_user)){
    e_abs <- all_abs[is_root_rproj_user]; fi <- file.info(e_abs); fi <- fi[!is.na(fi$size), , drop=FALSE]
    if (nrow(fi)) excl_rows[[length(excl_rows)+1]] <- data.frame(
      path=rel_from_root(rownames(fi)), bytes=as.numeric(fi$size), reason=".Rproj.user_root", stringsAsFactors=FALSE)
  }
  if (any(is_root_rhistory)){
    e_abs <- all_abs[is_root_rhistory]; fi <- file.info(e_abs); fi <- fi[!is.na(fi$size), , drop=FALSE]
    if (nrow(fi)) excl_rows[[length(excl_rows)+1]] <- data.frame(
      path=rel_from_root(rownames(fi)), bytes=as.numeric(fi$size), reason=".Rhistory_root", stringsAsFactors=FALSE)
  }
  
  # output blacklist excludes
  if (any(drop_output)){
    e_abs <- all_abs[drop_output]; fi <- file.info(e_abs); fi <- fi[!is.na(fi$size), , drop=FALSE]
    if (nrow(fi)) excl_rows[[length(excl_rows)+1]] <- data.frame(
      path=rel_from_root(rownames(fi)), bytes=as.numeric(fi$size), reason="output_blacklist", stringsAsFactors=FALSE)
  }
  
  # large file excludes
  if (length(large_skip_abs)){
    fi <- file.info(large_skip_abs); fi <- fi[!is.na(fi$size), , drop=FALSE]
    if (nrow(fi)) excl_rows[[length(excl_rows)+1]] <- data.frame(
      path=rel_from_root(rownames(fi)), bytes=as.numeric(fi$size),
      reason=paste0("size>", max_file_mb_any, "MB"), stringsAsFactors=FALSE)
  }
  
  if (length(excl_rows)){
    excl <- do.call(rbind, excl_rows)
    excl <- excl[order(excl$reason, excl$path), ]
    utils::write.csv(excl, file.path(meta_dir, "EXCLUDED_WITH_SIZE.csv"), row.names=FALSE)
  } else {
    writeLines(character(0), file.path(meta_dir, "EXCLUDED_WITH_SIZE.csv"))
  }
  
  # lightweight git meta
  writeLines(..git("status --short"), file.path(meta_dir, "GIT_STATUS.txt"))
  writeLines(..git("log --oneline -n 50"), file.path(meta_dir, "GIT_LOG_1line.txt"))
  
  # zip: preserve tree (no flatten)
  ..zip_with_fallback(zipfile = zpath, files = ".", root = stage, quiet = TRUE)
  message("[OK] wrote: ", zpath)
  invisible(zpath)
}

# ---------- preset ----------
z2g <- function(){
  make_dist_zip(
    output_exclude_subdirs = character(0),  # include all under output/
    max_file_mb_any        = 2,
    allow_large_any        = FALSE
  )
}

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
..sha  <- function() {
  x <- try(system2("git", c("rev-parse", "--short", "HEAD"), stdout=TRUE, stderr=TRUE), silent=TRUE)
  if (inherits(x,"try-error") || length(x)==0) "nohash" else trimws(x[1])
}
..git <- function(args) {
  x <- try(system2("git", args, stdout=TRUE, stderr=TRUE), silent=TRUE)
  if (inherits(x,"try-error")) return(character(0))
  x
}

..zip_with_fallback <- function(zipfile, files, root=".", quiet=TRUE){
  owd <- getwd(); on.exit(setwd(owd), add=TRUE)
  setwd(root)
  # prefer base::zip (Windows Rtools zip also ok)
  res <- try(utils::zip(zipfile=zipfile, files=files, flags="-r9Xq"), silent=TRUE)
  if (inherits(res,"try-error")) {
    # fallback: shell via system2 (7z/zip if available)
    # do minimal; environment dependent
    cmd <- if (.Platform$OS.type=="windows") "powershell" else "sh"
    warning("zip fallback path used; please verify archive integrity.")
  }
}
..resolve_onedrive_root <- function(){
  # 1) Official env var if present
  od <- Sys.getenv("OneDrive", unset = "")
  if (nzchar(od) && dir.exists(od)) return(gsub("\\\\","/",od))
  # 2) Typical personal OneDrive under USERPROFILE
  up <- Sys.getenv("USERPROFILE", unset = "")
  if (nzchar(up)) {
    cand <- file.path(up, "OneDrive")
    if (dir.exists(cand)) return(gsub("\\\\","/",cand))
    # 3) Japanese Windows naming: "OneDrive - 個人"
    cand2 <- file.path(up, "OneDrive - 個人")
    if (dir.exists(cand2)) return(gsub("\\\\","/",cand2))
    # 4) Company tenant naming: "OneDrive - *"
    pats <- list.dirs(up, full.names=TRUE, recursive=FALSE)
    hits <- grep("OneDrive - ", basename(pats), fixed=TRUE, value=TRUE)
    if (length(hits) > 0) {
      cand3 <- file.path(up, hits[1])
      if (dir.exists(cand3)) return(gsub("\\\\","/",cand3))
    }
  }
  ""
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
    # allow single large file exceptions by extension (FALSE = hard cap)
    allow_large_any        = FALSE
){
  project_root <- ..norm(project_root)
  parent_dir   <- ..norm(parent_dir)
  
  # stage temp
  stage <- file.path(tempdir(), sprintf("dmp_stage_%s", as.integer(runif(1,1e6,9e6))))
  on.exit(unlink(stage, recursive=TRUE, force=TRUE), add=TRUE)
  dir.create(stage, recursive=TRUE, showWarnings=FALSE)
  
  # copy tree (exclude a few)
  root_files <- list.files(project_root, all.files=TRUE, no..=TRUE, full.names=TRUE, include.dirs=TRUE)
  root_files <- ..norm(root_files)
  exclude_root <- c(".git", ".Rproj.user", ".Rhistory")
  keep <- !basename(root_files) %in% exclude_root
  keep_files <- root_files[keep]
  
  # copy to stage
  for (src in keep_files) {
    dst <- file.path(stage, basename(src))
    if (dir.exists(src)) {
      dir.create(dst, recursive=TRUE, showWarnings=FALSE)
      fsrc <- list.files(src, recursive=TRUE, all.files=TRUE, full.names=TRUE, include.dirs=FALSE, no..=TRUE)
      fsrc <- ..norm(fsrc)
      for (p in fsrc) {
        rel <- substring(p, nchar(..norm(src))+2)
        .dst <- file.path(dst, rel)
        ..safe_copy(p, .dst)
      }
    } else {
      ..safe_copy(src, dst)
    }
  }
  
  # output/ subfolder blacklist
  if (length(output_exclude_subdirs) > 0) {
    for (nm in output_exclude_subdirs) {
      p <- file.path(stage, "output", nm)
      if (dir.exists(p)) unlink(p, recursive=TRUE, force=TRUE)
    }
  }
  
  # global size guard (MB)
  info <- file.info(list.files(stage, recursive=TRUE, all.files=TRUE, full.names=TRUE, include.dirs=FALSE, no..=TRUE))
  info$path <- rownames(info)
  too_big <- which(!is.na(info$size) & info$size > (max_file_mb_any * 1024^2))
  meta_dir <- file.path(dirname(stage), "meta")
  dir.create(meta_dir, recursive=TRUE, showWarnings=FALSE)
  write.csv(info[,c("path","size")], file.path(meta_dir, "ALL_WITH_SIZE.csv"), row.names=FALSE)
  if (length(too_big) > 0 && !isTRUE(allow_large_any)) {
    excl <- info[too_big, c("path","size")]
    write.csv(excl, file.path(meta_dir, "EXCLUDED_WITH_SIZE.csv"), row.names=FALSE)
    # remove large files from stage
    for (p in excl$path) unlink(p, force=TRUE)
  }
  
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
    sub(paste0("^", gsub("([.])","\\.", root_norm), "/?"), "", p_norm)
  }
  
  # zip build
  ..zip_with_fallback(zipfile = zpath, files = ".", root = stage, quiet = TRUE)
  message("[OK] wrote: ", zpath)
  invisible(zpath)
}

# ---------- preset ----------
z2g <- function(){
  # 1) build primary zip to parent_dir (default behavior)
  z <- make_dist_zip(
    output_exclude_subdirs = character(0),  # include all under output/
    max_file_mb_any        = 2,
    allow_large_any        = FALSE
  )
  # 2) mirror to OneDrive/projects if resolvable
  od_root <- ..resolve_onedrive_root()
  if (nzchar(od_root)) {
    onedrive_projects <- file.path(od_root, "projects")
    dir.create(onedrive_projects, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(z) && file.exists(z)) {
      dst <- file.path(onedrive_projects, basename(z))
      file.copy(z, dst, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
      message("[OK] mirrored: ", dst)
    } else {
      warning("[WARN] primary zip not found; skip OneDrive mirror.")
    }
  } else {
    warning("[WARN] OneDrive root not found; skip OneDrive mirror.")
  }
  invisible(z)
}

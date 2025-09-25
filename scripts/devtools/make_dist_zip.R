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
# - EOL normalize (LF) for text files; binary untouched
# - Embed meta (sizes, excludes, manifest, build log) inside ZIP

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

..text_ext <- function(p) grepl("\\.(r|R|cpp|hpp|h|c|sh|md|txt|csv|json|ya?ml|Rmd|js|ts|css|html)$", p, ignore.case=TRUE)
..code_ext <- function(p) grepl("\\.(r|R|cpp|hpp|h|c|sh|js|ts|css|html|Rmd)$", p, ignore.case=TRUE)
..is_under <- function(p, dirprefix) startsWith(..norm(p), paste0(..norm(dirprefix),"/"))

# FAST EOL normalizer for text files only (O(n))
..normalize_eol_text <- function(path){
  sz <- tryCatch(file.info(path)$size, error=function(e) NA_real_)
  if (is.na(sz) || sz <= 0) return(list(changed=FALSE, msg="empty"))
  txt <- readChar(path, nchars=sz, useBytes=TRUE)
  if (!length(txt)) return(list(changed=FALSE, msg="empty"))
  new <- gsub("\r\n", "\n", txt, fixed=TRUE)
  new2 <- gsub("\r", "\n", new, fixed=TRUE)
  if (!identical(txt, new2)) {
    con <- file(path, "wb"); on.exit(close(con), add=TRUE)
    writeChar(new2, con, eos=NULL, useBytes=TRUE)
    return(list(changed=TRUE, msg="normalized"))
  }
  list(changed=FALSE, msg="as-is")
}

..has_non_ascii <- function(path){
  sz <- tryCatch(file.info(path)$size, error=function(e) NA_real_)
  if (is.na(sz) || sz <= 0) return(FALSE)
  if (sz > 50e6) return(FALSE)
  con <- file(path, "rb"); on.exit(close(con), add=TRUE)
  rawv <- readBin(con, what="raw", n=sz)
  any(as.integer(rawv) > 127L)
}

..md5_manifest <- function(paths) {
  suppressWarnings(as.data.frame(
    cbind(path=paths, md5=as.character(tools::md5sum(paths))),
    stringsAsFactors=FALSE
  ))
}

..write_lines <- function(lines, path){
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  con <- file(path, "wb"); on.exit(close(con), add=TRUE)
  for (s in lines) writeBin(charToRaw(paste0(s, "\n")), con, useBytes=TRUE)
}

..zip_with_fallback <- function(zipfile, root=".", quiet=TRUE){
  owd <- getwd(); on.exit(setwd(owd), add=TRUE)
  setwd(root)
  ok <- FALSE
  res <- try(utils::zip(zipfile=zipfile, files=".", flags="-r9Xq"), silent=TRUE)
  if (!inherits(res, "try-error") && file.exists(zipfile)) ok <- TRUE
  if (!ok) {
    z7 <- Sys.which("7z")
    if (nzchar(z7)) {
      args <- c("a", "-tzip", "-mx=9", shQuote(zipfile), shQuote("."))
      suppressWarnings(system2(z7, args=args, stdout=TRUE, stderr=TRUE))
      if (file.exists(zipfile)) ok <- TRUE
    }
  }
  if (!ok) {
    zipbin <- Sys.which("zip")
    if (nzchar(zipbin)) {
      args <- c("-r9Xq", shQuote(zipfile), ".")
      suppressWarnings(system2(zipbin, args=args, stdout=TRUE, stderr=TRUE))
      if (file.exists(zipfile)) ok <- TRUE
    }
  }
  if (!ok) stop("Failed to create zip by all methods (utils::zip, 7z, zip).")
  invisible(TRUE)
}

..resolve_onedrive_root <- function(){
  od <- Sys.getenv("OneDrive", unset = "")
  if (nzchar(od) && dir.exists(od)) return(gsub("\\\\","/",od))
  up <- Sys.getenv("USERPROFILE", unset = "")
  if (nzchar(up)) {
    cand <- file.path(up, "OneDrive"); if (dir.exists(cand)) return(gsub("\\\\","/",cand))
    cand2 <- file.path(up, "OneDrive - 個人"); if (dir.exists(cand2)) return(gsub("\\\\","/",cand2))
    pats <- list.dirs(up, full.names=TRUE, recursive=FALSE)
    hits <- grep("OneDrive - ", basename(pats), fixed=TRUE, value=TRUE)
    if (length(hits) > 0) {
      cand3 <- file.path(up, hits[1]); if (dir.exists(cand3)) return(gsub("\\\\","/",cand3))
    }
  }
  ""
}

..safe_copy <- function(src, dst){
  dir.create(dirname(dst), recursive=TRUE, showWarnings=FALSE)
  file.copy(src, dst, overwrite=TRUE, copy.mode=TRUE, copy.date=TRUE)
}

# robust relative path (no regex)
..rel_from <- function(p, base){
  p <- ..norm(p); base <- ..norm(base)
  if (identical(p, base)) return(".")
  pref <- paste0(base, "/")
  # 先頭の '/' を確実に落とす（+1L）
  if (startsWith(p, pref)) substring(p, nchar(pref) + 1L) else p
}

# 追加ヘルパー：Windowsで違法になりうる名前（最後のパス要素のみ判定）
.is_illegal_win_basename <- function(bn){
  # 禁止文字:  < > : " | ? *   （/はディレクトリなので除外）
  if (grepl('[<>:"|?*]', bn)) return(TRUE)
  # 末尾スペース/ピリオドも不可
  if (grepl('[ \\.]$', bn)) return(TRUE)
  # 制御文字
  if (grepl('[[:cntrl:]]', bn)) return(TRUE)
  FALSE
}

.verify_zip_md5 <- function(zip_path){
  td <- file.path(tempdir(), sprintf("dmp_verify_%s", as.integer(runif(1,1e6,9e6))))
  on.exit(unlink(td, recursive=TRUE, force=TRUE), add=TRUE)
  dir.create(td, recursive=TRUE, showWarnings=FALSE)
  
  # 1) MANIFEST 抽出
  utils::unzip(zip_path, files = "meta/MANIFEST.md5", exdir = td)
  man <- try(read.csv(file.path(td, "meta", "MANIFEST.md5"), stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(man, "try-error")) stop("MANIFEST.md5 not found in zip")
  if (!all(c("path","md5") %in% names(man))) stop("MANIFEST.md5 schema invalid")
  
  # 2) zipのエントリ一覧（抽出前に存在照合：表記揺れを正規化）
  lst <- utils::unzip(zip_path, list = TRUE)
  if (!("Name" %in% names(lst))) stop("unzip list output missing 'Name'")
  
  .clean_names <- function(v){
    v <- ..norm(v)
    v <- sub("^\\./", "", v)  # 先頭 ./ を除去
    v <- sub("^/",    "", v)  # 先頭 /  を除去
    v
  }
  lst_names  <- .clean_names(lst$Name)
  man_paths  <- .clean_names(man$path)
  
  present <- man_paths %in% lst_names
  if (!all(present)) {
    miss <- man_paths[!present]
    stop(sprintf("Entries missing from zip listing: %d (e.g. %s)", length(miss), miss[1]))
  }

  
  bn <- basename(man_paths)
  illegal <- vapply(bn, .is_illegal_win_basename, logical(1))
  legal_paths <- man_paths[!illegal]
  
  ex <- file.path(td, "ex")
  dir.create(ex, recursive = TRUE, showWarnings = FALSE)
  if (length(legal_paths) > 0) {
    utils::unzip(zip_path, files = legal_paths, exdir = ex)
    paths <- file.path(ex, legal_paths)
    if (!all(file.exists(paths))) {
      stop("Some legal paths failed to extract; aborting md5 verification")
    }
    md5_now <- as.character(tools::md5sum(paths))
    # 参照は「正規化前の行」との対応を取るため、正規化版で一致を取ってから md5 を引く
    idx <- match(legal_paths, .clean_names(man$path))
    md5_ref <- man$md5[idx]
    bad <- which(md5_now != md5_ref)
    if (length(bad) > 0) {
      df <- data.frame(path = legal_paths[bad], manifest = md5_ref[bad], actual = md5_now[bad], stringsAsFactors = FALSE)
      write.csv(df, file.path(td, "MD5_MISMATCH.csv"), row.names = FALSE)
      stop(sprintf("MD5 mismatch for %d file(s); see %s", nrow(df), file.path(td, "MD5_MISMATCH.csv")))
    }
  }
  
  # 5) 違法名があった場合はスキップ報告（構造的存在は既に確認済み）
  if (any(illegal)) {
    msg <- sprintf("[VERIFY] skipped %d file(s) due to illegal Windows filename (e.g. %s)",
                   sum(illegal), man$path[which(illegal)[1]])
    message(msg)
  }
  TRUE
}

# ---------- main ----------
make_dist_zip <- function(
    project_root            = getwd(),
    parent_dir              = dirname(getwd()),
    zip_prefix              = "DMP_dist",
    output_exclude_subdirs  = character(0),
    max_file_mb_any         = 2,
    allow_large_any         = FALSE,
    normalize_text_eol      = TRUE,
    ascii_guard_code        = TRUE,
    ascii_guard_code_strict = FALSE,
    verbose                 = TRUE
){
  t_start <- proc.time()[["elapsed"]]
  project_root <- ..norm(project_root)
  parent_dir   <- ..norm(parent_dir)
  if (verbose) cat(sprintf("[make] root=%s\n", project_root))
  
  # stage temp
  stage <- file.path(tempdir(), sprintf("dmp_stage_%s", as.integer(runif(1,1e6,9e6))))
  on.exit(unlink(stage, recursive=TRUE, force=TRUE), add=TRUE)
  dir.create(stage, recursive=TRUE, showWarnings=FALSE)
  
  build_log <- character(0)
  
  # copy tree (exclude a few)
  root_files <- list.files(project_root, all.files=TRUE, no..=TRUE, full.names=TRUE, include.dirs=TRUE)
  root_files <- ..norm(root_files)
  exclude_root <- c(".git", ".Rproj.user", ".Rhistory")
  keep <- !basename(root_files) %in% exclude_root
  keep_files <- root_files[keep]
  if (verbose) cat(sprintf("[copy] items=%d\n", length(keep_files)))
  
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
  
  # Remove output/ subfolders in blacklist
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
  
  # meta dir (INSIDE stage)
  meta_dir <- file.path(stage, "meta")
  dir.create(meta_dir, recursive=TRUE, showWarnings=FALSE)
  write.csv(info[,c("path","size")], file.path(meta_dir, "ALL_WITH_SIZE.csv"), row.names=FALSE)
  
  if (length(too_big) > 0 && !isTRUE(allow_large_any)) {
    excl <- info[too_big, c("path","size")]
    excl$reason <- sprintf("size>%.1fMB", max_file_mb_any)
    write.csv(excl, file.path(meta_dir, "EXCLUDED_WITH_SIZE.csv"), row.names=FALSE)
    for (p in excl$path) unlink(p, force=TRUE)
    if (verbose) cat(sprintf("[size] excluded=%d (>%.1fMB)\n", nrow(excl), max_file_mb_any))
  } else {
    file.create(file.path(meta_dir, "EXCLUDED_WITH_SIZE.csv"))
  }
  
  # EOL normalize for text files & ASCII guard for code files
  files_all <- list.files(stage, recursive=TRUE, all.files=TRUE, full.names=TRUE, include.dirs=FALSE, no..=TRUE)
  files_all <- ..norm(files_all)
  non_ascii_rows <- data.frame(path=character(0), size=integer(0), stringsAsFactors=FALSE)
  eol_changed <- 0L
  
  text_files <- files_all[ vapply(files_all, function(fp) ..text_ext(substring(fp, nchar(stage)+2L)), logical(1)) ]
  n_txt <- length(text_files)
  if (verbose) cat(sprintf("[eol] text_files=%d\n", n_txt))
  
  if (isTRUE(normalize_text_eol) && n_txt > 0) {
    for (i in seq_along(text_files)) {
      fp <- text_files[i]
      res <- try(..normalize_eol_text(fp), silent=TRUE)
      if (!inherits(res, "try-error") && isTRUE(res$changed)) eol_changed <- eol_changed + 1L
      if (verbose && (i %% 200 == 0 || i == n_txt)) cat(sprintf("[eol] %d/%d (changed=%d)\n", i, n_txt, eol_changed))
    }
  }
  
  if (isTRUE(ascii_guard_code)) {
    code_files <- files_all[ vapply(files_all, function(fp) {
      rel <- substring(fp, nchar(stage)+2L)
      ..code_ext(rel) && !..is_under(rel, "test")
    }, logical(1)) ]
    for (fp in code_files) {
      rel <- substring(fp, nchar(stage)+2L)
      if (..has_non_ascii(fp)) {
        s <- tryCatch(file.info(fp)$size, error=function(e) NA_integer_)
        non_ascii_rows <- rbind(non_ascii_rows, data.frame(path=rel, size=as.integer(s), stringsAsFactors=FALSE))
      }
    }
  }
  
  if (nrow(non_ascii_rows) > 0) {
    write.csv(non_ascii_rows, file.path(meta_dir, "NON_ASCII_CODE.csv"), row.names=FALSE)
  } else {
    file.create(file.path(meta_dir, "NON_ASCII_CODE.csv"))
  }
  
  # names
  hash  <- ..sha(); stamp <- ..now()
  zname <- sprintf("%s_%s_%s.zip", zip_prefix, hash, stamp)
  zpath <- ..norm(file.path(parent_dir, zname))
  
  # manifest (md5) — NO REGEX, robust relative paths
  stage_files <- list.files(stage, recursive=TRUE, all.files=TRUE, full.names=TRUE, include.dirs=FALSE, no..=TRUE)
  stage_files_rel <- vapply(stage_files, ..rel_from, FUN.VALUE=character(1), base=stage)
  man <- ..md5_manifest(stage_files)
  names(man) <- c("path","md5")
  man$path <- stage_files_rel
  write.csv(man, file.path(meta_dir, "MANIFEST.md5"), row.names=FALSE)
  
  # build log
  t_copy <- proc.time()[["elapsed"]] - t_start
  build_log <- c(
    sprintf("[INFO] root=%s", project_root),
    sprintf("[INFO] files in stage=%d", length(stage_files)),
    sprintf("[INFO] text normalized=%d", eol_changed),
    sprintf("[INFO] non-ascii code files=%d", nrow(non_ascii_rows)),
    sprintf("[INFO] copy+prep sec=%.3f", t_copy),
    sprintf("[INFO] output zip=%s", zpath)
  )
  ..write_lines(build_log, file.path(meta_dir, "BUILD_LOG.txt"))
  
  # zip build (root="stage")
  if (verbose) cat("[zip] building...\n")
  ..zip_with_fallback(zipfile = zpath, root = stage, quiet = TRUE)
  if (!file.exists(zpath) || file.info(zpath)$size <= 0) stop("Zip not created or empty.")
  if (verbose) cat("[zip] done.\n")
  
  message("[OK] wrote: ", zpath)
  invisible(zpath)
}

# ---------- preset ----------
z2g <- function(verify_after_zip = TRUE){
  z <- make_dist_zip(
    output_exclude_subdirs  = character(0),
    max_file_mb_any         = 2,
    allow_large_any         = FALSE,
    normalize_text_eol      = TRUE,
    ascii_guard_code        = TRUE,
    ascii_guard_code_strict = FALSE,
    verbose                 = TRUE
  )
  if (isTRUE(verify_after_zip)) {
    ok <- .verify_zip_md5(z)
    if (isTRUE(ok)) message("[OK] md5 verified: ", z)
  }
  
  # 既存のOneDriveミラー（あなたの最新版）をそのまま踏襲
  od_root <- ..resolve_onedrive_root()
  if (nzchar(od_root)) {
    dst_dir <- file.path(od_root, "projects")
    dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(z) && file.exists(z)) {
      dst <- file.path(dst_dir, basename(z))
      okc <- file.copy(z, dst, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
      if (isTRUE(okc)) message("[OK] mirrored: ", dst) else warning("[WARN] failed to mirror: ", dst)
    } else {
      warning("[WARN] primary zip not found; skip OneDrive mirror.")
    }
  } else {
    warning("[WARN] OneDrive root not found; skip OneDrive mirror.")
  }
  invisible(z)
}

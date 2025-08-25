# scripts/devtools/make_dist_zip.R
# Build a distribution zip one level above the project root.
# Defaults: include scripts/, data/ (stub >1MB), scripts/bench, bench, test (incl. test/bench), output/
# Options to extend/override via arguments. Writes meta with EXCLUDED_LIST.txt / LARGE_FILES.txt.
# Windows/mac/Linux 対応（utils::zip 使用）。

# ---- helpers ---------------------------------------------------------------

bytes_1mb <- 1024L * 1024L

.safe_norm <- function(path) {
  tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE),
           error = function(e) path)
}

.is_under <- function(rel, root) {
  rel == root || startsWith(rel, paste0(root, "/"))
}

.list_files_recursive <- function(dir) {
  if (!dir.exists(dir)) return(character())
  f <- list.files(dir, recursive = TRUE, all.files = TRUE, full.names = TRUE)
  f[file.info(f)$isdir %in% FALSE]
}

.copy_with_dirs <- function(src, dst_root, rel) {
  dst <- file.path(dst_root, rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dst, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  if (!ok) stop("Failed to copy: ", src, " -> ", dst)
  dst
}

.create_empty_stub <- function(dst_root, rel) {
  dst <- file.path(dst_root, rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  con <- file(dst, open = "wb"); close(con); dst
}

# ---- main ------------------------------------------------------------------
# 既存の限定版と後方互換: 旧シグネチャだけで呼んでも動くように既定値を広めに設定
make_dist_zip <- function(
    project_root = getwd(),
    zip_prefix   = "DMP_dist",
    parent_dir   = dirname(getwd()),
    head_lines   = 50,                        # メタに使う場合の上部コメント行の長さ（互換用・任意）
    # 追加オプション（汎用化）
    extra_include       = character(),        # 例: c("bench","scripts/bench","output","test","data")
    extra_exclude_dirs  = character(),        # 例: c("notes","debug")
    extra_exclude_globs = character(),        # 例: c("*.rproj","*.tmp")
    include_root_files  = c("README.md","LICENSE","DESCRIPTION"),
    default_exclude_dirs = c(".git",".Rproj.user","dist","release"),
    default_exclude_glob = c("*.zip"),
    include_data        = TRUE,               # data/ を含めるか
    data_large_threshold = bytes_1mb,         # data/ のスタブ化閾値
    stub_large_data     = TRUE,               # TRUE: data/ で閾値超過を空ファイルで stub
    dry_run             = FALSE               # TRUE: コピー＆zipせず、メタ一覧だけ作る
) {
  project_root <- .safe_norm(project_root)
  parent_dir   <- .safe_norm(parent_dir)
  
  # 生成名
  hash_stub <- substr(as.character(as.integer(as.numeric(Sys.time()))), 1, 8)
  ts <- format(Sys.time(), "%Y%m%d%H%M")
  zip_base <- sprintf("%s_%s_%s", zip_prefix, hash_stub, ts)
  
  # ステージング先（zip内容の実体）
  staging <- file.path(tempdir(), paste0("staging_", zip_base))
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  
  # メタ出力（zipと同階層）
  meta_dir <- file.path(parent_dir, paste0("meta_", hash_stub, "_", ts))
  dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
  excluded_list_path <- file.path(meta_dir, "EXCLUDED_LIST.txt")
  large_list_path    <- file.path(meta_dir, "LARGE_FILES.txt")
  included_list_path <- file.path(meta_dir, "INCLUDED_LIST.txt")
  
  # 既定の「必ず含める候補」
  include_roots <- c(
    "scripts",         # 既定で scripts/ 全体
    if (include_data) "data",
    "scripts/bench",
    "bench",
    "test",
    "output"
  )
  # ユーザー追加分（重複は後でユニーク化）
  include_roots <- unique(c(include_roots, extra_include))
  
  # プロジェクト内の全ファイル候補
  all_files <- .list_files_recursive(project_root)
  rel_all   <- gsub(paste0("^", gsub("\\\\","/", project_root), "/?"), "", gsub("\\\\","/", all_files))
  
  # まず「既定の除外」を適用。ただし include_roots 配下は除外しない（必ず残す）
  keep <- rep(TRUE, length(all_files))
  for (i in seq_along(all_files)) {
    rel <- rel_all[i]
    top <- strsplit(rel, "/")[[1]]; top <- if (length(top)) top[1] else ""
    under_include <- any(startsWith(rel, paste0(include_roots, "/")) | rel %in% include_roots)
    if (!under_include) {
      if (top %in% default_exclude_dirs) keep[i] <- FALSE
      if (keep[i] && length(default_exclude_glob)) {
        for (g in default_exclude_glob) {
          if (grepl(glob2rx(g), basename(rel), ignore.case = TRUE)) { keep[i] <- FALSE; break }
        }
      }
      if (keep[i] && length(extra_exclude_dirs) && top %in% extra_exclude_dirs) keep[i] <- FALSE
      if (keep[i] && length(extra_exclude_globs)) {
        for (g in extra_exclude_globs) {
          if (grepl(glob2rx(g), basename(rel), ignore.case = TRUE)) { keep[i] <- FALSE; break }
        }
      }
    }
  }
  all_files <- all_files[keep]; rel_all <- rel_all[keep]
  
  # 最終的に含める集合: include_roots配下＋ルート直下の README.md 等
  want <- logical(length(all_files))
  for (j in seq_along(all_files)) {
    rel <- rel_all[j]
    if (any(vapply(include_roots, function(rt) .is_under(rel, rt), logical(1)))) want[j] <- TRUE
  }
  root_paths <- file.path(project_root, include_root_files)
  root_files <- root_paths[file.exists(root_paths)]
  root_rel   <- basename(root_files)
  
  final_src <- c(all_files[want], root_files)
  final_rel <- c(rel_all[want],   root_rel)
  
  # 重複除去
  o <- order(final_rel)
  final_src <- final_src[o]; final_rel <- final_rel[o]
  dupe <- duplicated(final_rel)
  if (any(dupe)) { final_src <- final_src[!dupe]; final_rel <- final_rel[!dupe] }
  
  # メタ初期化
  writeLines(character(), excluded_list_path)
  writeLines(c(
    "# Files larger than threshold are stubbed in the zip (data/ only when enabled).",
    paste0("# Threshold: ", data_large_threshold, " bytes"),
    "# Format: <size_bytes> <relative_path>",
    ""
  ), large_list_path)
  writeLines(c("# Included relative paths:", ""), included_list_path)
  
  # ステージングへコピー（data/ 大きいファイルは stub 可）
  info <- file.info(final_src, extra_cols = FALSE)
  for (k in seq_along(final_src)) {
    src <- final_src[k]; rel <- final_rel[k]
    if (stub_large_data && .is_under(rel, "data")) {
      sz <- info$size[k]
      if (!is.na(sz) && sz > data_large_threshold) {
        .create_empty_stub(staging, rel)
        cat(sprintf("%d %s\n", as.integer(sz), rel), file = large_list_path, append = TRUE)
        cat(rel, "\n", file = included_list_path, append = TRUE)
        next
      }
    }
    .copy_with_dirs(src, staging, rel)
    cat(rel, "\n", file = included_list_path, append = TRUE)
  }
  
  # 除外一覧（プロジェクトにあるが zip に入っていない実ファイル）
  staged_all <- list.files(staging, recursive = TRUE, all.files = TRUE, full.names = FALSE)
  staged_full <- file.path(staging, staged_all)
  staged_all <- staged_all[file.info(staged_full)$isdir %in% FALSE]
  
  all_proj_files <- .list_files_recursive(project_root)
  all_proj_rel <- gsub(paste0("^", gsub("\\\\","/", project_root), "/?"), "", gsub("\\\\","/", all_proj_files))
  all_proj_rel <- all_proj_rel[file.info(all_proj_files)$isdir %in% FALSE]
  
  excluded_rel <- setdiff(all_proj_rel, staged_all)
  if (length(excluded_rel)) writeLines(sort(excluded_rel), excluded_list_path, useBytes = TRUE)
  
  # dry-run はここまで（zip作成せずにメタだけ出す）
  if (isTRUE(dry_run)) {
    message("[DRY RUN] No zip created.")
    message("[META] Included: ", .safe_norm(included_list_path))
    message("[META] Large   : ", .safe_norm(large_list_path))
    message("[META] Excluded: ", .safe_norm(excluded_list_path))
    return(invisible(list(zip = NULL, meta = meta_dir, staging = staging)))
  }
  
  # zip 作成
  old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
  setwd(staging)
  zip_path <- file.path(parent_dir, paste0(zip_base, ".zip"))
  files_to_zip <- list.files(".", recursive = TRUE, all.files = TRUE, full.names = FALSE)
  utils::zip(zipfile = zip_path, files = files_to_zip)
  
  message("[OK] Wrote zip : ", .safe_norm(zip_path))
  message("[OK] Wrote meta: ", .safe_norm(meta_dir))
  invisible(list(zip = zip_path, meta = meta_dir, staging = staging))
}


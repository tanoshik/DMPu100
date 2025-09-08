# scripts/devtools/make_dist_zip.R
# 日本語コメントOK / 開発・評価共有向け。GPT参照に便利な meta を厚めに同梱。
# 依存: 'zip' パッケージがあればそれを優先（警告が減ります）。無ければ utils::zip を使用。

# ===== 内部ヘルパ =====
..norm <- function(x) gsub("\\\\","/",x)
..now  <- function() format(Sys.time(), "%Y%m%d%H%M")
..git  <- function(cmd) tryCatch(system(paste("git", cmd), intern=TRUE, ignore.stderr=TRUE), error=function(e) character(0))
..sha  <- function() { s <- ..git("rev-parse --short HEAD"); if (length(s)==0) "nogit" else s[1] }
..branch <- function(){ s <- ..git("rev-parse --abbrev-ref HEAD"); if (length(s)==0) "nogit" else s[1] }

..zip_with_fallback <- function(zipfile, files, root=".", quiet=TRUE){
  o <- getwd(); setwd(root); on.exit(setwd(o), add=TRUE)
  ok <- FALSE
  if (requireNamespace("zip", quietly=TRUE)) {
    # zip::zipr は静かで高速。Windowsでも安定。
    zip::zipr(zipfile, files = files, include_directories = TRUE)
    ok <- TRUE
  } else {
    # utils::zip は Windows では外部 zip.exe 依存で警告が出ることあり
    flags <- if (quiet) "-r9Xq" else "-r9X"
    z <- try(utils::zip(zipfile = zipfile, files = files, flags = flags), silent = quiet)
    if (!inherits(z, "try-error")) ok <- TRUE
  }
  if (!ok) stop("zip backend not available. install.packages('zip') を推奨。")
}

..sha256 <- function(path){
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("openssl", quietly=TRUE)) return(as.character(openssl::sha256(file(path, "rb"))))
  if (requireNamespace("digest",   quietly=TRUE)) return(digest::digest(file=path, algo="sha256"))
  "NA (install 'openssl' or 'digest')"
}

..safe_copy <- function(src, dst){
  dir.create(dirname(dst), recursive=TRUE, showWarnings=FALSE)
  file.copy(src, dst, overwrite=TRUE, copy.mode=TRUE, copy.date=TRUE)
}

# ===== メイン：配布ZIP作成 =====
make_dist_zip <- function(
    project_root   = getwd(),
    parent_dir     = dirname(getwd()),
    zip_prefix     = "DMP_dist",
    mode           = c("dev","review","release"),
    # RDSの取り扱い：none/s1000/s100000/s1000+s100000/all
    include_rds    = c("s1000+s100000","none","all","s1000","s100000"),
    # outputの取り扱い：sample/all/none
    include_output = c("sample","all","none"),
    include_test   = TRUE,      # test/ を含む（GPT検証・再現で有用）
    include_bench  = TRUE,      # bench/ を含む（ベンチログ共有）
    include_notes  = FALSE,     # notes/ を含む（内部メモも共有したければ TRUE）
    include_hidden = FALSE,     # .で始まる隠しを拾う（通常 FALSE）
    bundle_meta    = TRUE,      # meta_*/ をZIP内に同梱
    meta_level     = c("full","basic"),
    head_lines     = 120,       # head 抜粋行数
    dry_run        = FALSE,
    extra_include  = character(0),
    extra_exclude  = character(0),
    # 追加：事前に manifest / SESSIONINFO を自動更新
    auto_update_manifest = TRUE,
    auto_update_session  = TRUE
){
  mode           <- match.arg(mode)
  include_rds    <- match.arg(include_rds)
  include_output <- match.arg(include_output)
  meta_level     <- match.arg(meta_level)
  
  owd <- setwd(project_root); on.exit(setwd(owd), add=TRUE)
  
  # ---- 事前更新：manifest / SESSIONINFO ----
  if (isTRUE(auto_update_manifest) && file.exists("scripts/devtools/write_manifest.R")){
    try({
      message("[INFO] update manifest via scripts/devtools/write_manifest.R")
      system2("Rscript", c("scripts/devtools/write_manifest.R"), stdout=TRUE, stderr=TRUE)
    }, silent=TRUE)
  }
  if (isTRUE(auto_update_session)){
    try({
      message("[INFO] write SESSIONINFO to output/dev/SESSIONINFO.txt")
      dir.create("output/dev", recursive=TRUE, showWarnings=FALSE)
      writeLines(capture.output(sessionInfo()), "output/dev/SESSIONINFO.txt")
    }, silent=TRUE)
  }
  
  # ---- ZIP/Meta 名 ----
  hash  <- ..sha(); stamp <- ..now()
  zname <- sprintf("%s_%s_%s.zip", zip_prefix, hash, stamp)
  zpath <- ..norm(file.path(parent_dir, zname))
  metan <- sprintf("meta_%s_%s", hash, stamp)
  meta_tmp <- file.path(tempdir(), metan)
  dir.create(meta_tmp, showWarnings=FALSE, recursive=TRUE)
  
  # ---- 除外ポリシ（最小限） ----
  base_exclude <- c(
    "^\\.git($|/)", "^\\.Rproj\\.user($|/)",
    "^dist($|/)", "^release($|/)",
    "\\.zip$",
    "^debug_.*\\.R$", "^scratch_.*\\.R$"
  )
  if (!include_bench) base_exclude <- c(base_exclude, "^bench($|/)")
  if (!include_test)  base_exclude <- c(base_exclude, "^test($|/)")
  if (!include_notes) base_exclude <- c(base_exclude, "^notes($|/)")
  if (include_output=="none") base_exclude <- c(base_exclude, "^output($|/)")
  if (!include_hidden) base_exclude <- c(base_exclude, "^\\.")
  if (length(extra_exclude)) base_exclude <- c(base_exclude, extra_exclude)
  # 原則 .rds は除外（ホワイトリストで戻す）
  base_exclude <- c(base_exclude, "\\.rds$")
  
  list_all <- function(root="."){
    p <- list.files(root, recursive=TRUE, all.files=include_hidden, include.dirs=TRUE, no..=TRUE)
    ..norm(p)
  }
  all_paths <- list_all(".")
  is_excl <- function(path, pats) any(vapply(pats, function(rx) grepl(rx, path, perl=TRUE), logical(1)))
  
  # ---- 初期ホワイトリスト（網羅寄り） ----
  whitelist <- c("scripts","data","README.md","LICENSE","DESCRIPTION",".gitattributes",".Rprofile")
  if (include_output=="sample") whitelist <- c(whitelist, "output/sample", "output/manifest", "output/dev/SESSIONINFO.txt")
  if (include_output=="all")    whitelist <- c(whitelist, "output")
  whitelist <- unique(c(whitelist, extra_include))
  
  # RDSを戻す
  rds_keep <- character(0)
  if (include_rds %in% c("s1000","s1000+s100000","all")) rds_keep <- c(rds_keep, "data/virtual_db_u100_S1000_seed123.rds")
  if (include_rds %in% c("s100000","s1000+s100000","all")) rds_keep <- c(rds_keep, "data/virtual_db_u100_S100000_seed123.rds")
  if (include_rds=="all" && dir.exists("data")){
    rds_keep <- unique(c(rds_keep, file.path("data", grep("\\.rds$", list.files("data"), value=TRUE))))
  }
  rds_keep <- rds_keep[file.exists(rds_keep)]
  
  # ---- 収集 ----
  cand <- character(0)
  add_path <- function(p){
    if (dir.exists(p)) {
      x <- list.files(p, recursive=TRUE, all.files=include_hidden, include.dirs=FALSE, no..=TRUE, full.names=TRUE)
      cand <<- c(cand, ..norm(x))
    } else if (file.exists(p)) {
      cand <<- c(cand, ..norm(p))
    }
  }
  for (w in whitelist) add_path(w)
  cand <- unique(cand)
  
  keep <- vapply(cand, function(p) !is_excl(p, base_exclude), logical(1))
  cand <- cand[keep]
  cand <- unique(c(cand, ..norm(rds_keep)))
  
  if (length(cand)==0) stop("収集対象が空です。whitelist/フラグを調整してください。")
  
  # ---- ステージング ----
  stage <- file.path(tempdir(), paste0("stage_", hash, "_", stamp))
  dir.create(stage, recursive=TRUE, showWarnings=FALSE)
  for (src in cand) {
    ..safe_copy(src, file.path(stage, src))
  }
  
  # ---- meta 生成（GPT向けに充実） ----
  wmeta <- function(rel, lines){
    out <- file.path(meta_tmp, rel)
    dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
    writeLines(lines, out, useBytes=TRUE)
  }
  # ABOUT
  wmeta("ABOUT.txt", c(
    "Bundled meta for development/review/GPT reference.",
    paste0("branch: ", ..branch()),
    paste0("git_sha: ", ..sha()),
    paste0("timestamp: ", stamp),
    paste0("mode: ", mode),
    paste0("include_rds: ", include_rds),
    paste0("include_output: ", include_output)
  ))
  # INCLUDED / EXCLUDED
  included_rel <- ..norm(sub(paste0("^", ..norm(project_root), "/?"), "", ..norm(cand)))
  wmeta("INCLUDED_LIST.txt", included_rel)
  excluded_all <- all_paths[!dir.exists(all_paths)]
  excluded_rel <- setdiff(excluded_all, included_rel)
  wmeta("EXCLUDED_LIST.txt", excluded_rel)
  # TREE (ざっくり)
  wmeta("TREE_DIRS.txt", paste(sort(unique(dirname(included_rel))), collapse="\n"))
  # FILE_SIZES
  sizes <- file.info(cand)$size
  size_tab <- data.frame(path=included_rel, bytes=as.numeric(sizes),
                         ext=tolower(sub(".*\\.(\\w+)$","\\1", included_rel)), stringsAsFactors=FALSE)
  write.csv(size_tab, file.path(meta_tmp,"FILE_SIZES.csv"), row.names=FALSE)
  # CHECKSUMS（上限2000）
  ch <- c("SHA256  PATH")
  nmax <- min(2000L, length(cand))
  for (i in seq_len(nmax)) ch <- c(ch, sprintf("%s  %s", ..sha256(cand[i]), included_rel[i]))
  wmeta("CHECKSUMS_SHA256.txt", ch)
  # Git 情報
  wmeta("GIT_STATUS.txt", ..git("status --short"))
  if (meta_level=="full"){
    wmeta("GIT_DIFF_SUMMARY.txt", ..git("diff --stat"))
    wmeta("GIT_LOG_1line.txt", ..git("log --oneline -n 50"))
    wmeta("GIT_REMOTE_-v.txt", ..git("remote -v"))
  }
  # SESSIONINFO / manifest を meta にもコピー（存在すれば）
  if (file.exists("output/dev/SESSIONINFO.txt")) {
    ..safe_copy("output/dev/SESSIONINFO.txt", file.path(meta_tmp,"dev/SESSIONINFO.txt"))
  }
  if (file.exists("output/manifest/virtual_db_manifest.csv")) {
    ..safe_copy("output/manifest/virtual_db_manifest.csv", file.path(meta_tmp,"manifest/virtual_db_manifest.csv"))
  }
  # 代表テキストの先頭抜粋
  if (meta_level=="full"){
    head_dir <- file.path(meta_tmp,"heads"); dir.create(head_dir, showWarnings=FALSE)
    text_like <- included_rel[grepl("\\.(r|R|md|txt|csv|tsv|yaml|yml|json|Rmd)$", included_rel)]
    for (p in utils::head(text_like, 200L)){
      if (file.exists(p)) {
        ln <- try(readLines(p, n=head_lines, warn=FALSE), silent=TRUE)
        if (!inherits(ln,"try-error")) writeLines(ln, file.path(head_dir, paste0(gsub("[/:\\\\]", "__", p), ".head.txt")), useBytes=TRUE)
      }
    }
  }
  
  # ---- meta を stage へ同梱 ----
  if (isTRUE(bundle_meta)){
    mf <- list.files(meta_tmp, recursive=TRUE, all.files=TRUE, include.dirs=TRUE, no..=TRUE, full.names=FALSE)
    for (f in mf){
      src <- file.path(meta_tmp, f)
      dst <- file.path(stage, metan, f)
      if (dir.exists(src)) dir.create(dst, recursive=TRUE, showWarnings=FALSE) else ..safe_copy(src,dst)
    }
  }
  
  # ---- dry-run ----
  if (isTRUE(dry_run)){
    message("[DRY RUN] staged files: ", length(list.files(stage, recursive=TRUE)))
    message("[DRY RUN] stage dir: ", stage)
    return(invisible(list(zip=zpath, stage_dir=stage)))
  }
  
  # ---- ZIP化 ----
  files_for_zip <- list.files(stage, recursive=TRUE, include.dirs=TRUE, all.files=TRUE, no..=TRUE)
  ..zip_with_fallback(zipfile = zpath, files = files_for_zip, root = stage, quiet = TRUE)
  message("[OK] wrote: ", zpath)
  
  invisible(zpath)
}

# ===== プリセット（人間はこれだけ覚えればOK） =====

# 1) “いつもの”：GPT共有向け・最小操作・無引数
z2g <- function(){
  # S1000/S100000 RDS・sample出力・test/bench 同梱、metaは full
  make_dist_zip(mode="dev",
                include_rds    = "s1000+s100000",
                include_output = "all",
                include_test   = TRUE,
                include_bench  = TRUE,
                include_notes  = FALSE,
                bundle_meta    = TRUE,
                meta_level     = "full",
                auto_update_manifest = TRUE,
                auto_update_session  = TRUE)
}

# 2) これまでの make_dev_bundle（可変だけど安全な既定値）
make_dev_bundle <- function(include_rds=c("s1000+s100000","none","all"),
                            include_output=c("sample","all","none"),
                            ...){
  include_rds    <- match.arg(include_rds)
  include_output <- match.arg(include_output)
  make_dist_zip(mode="dev",
                include_rds=include_rds,
                include_output=include_output,
                bundle_meta=TRUE, meta_level="full",
                auto_update_manifest=TRUE, auto_update_session=TRUE,
                ...)
}

# 3) 何も省かない“全部入り”のバックアップ（除外ごく最小）
backup_all <- function(){
  make_dist_zip(mode="dev",
                include_rds    = "all",
                include_output = "all",
                include_test   = TRUE,
                include_bench  = TRUE,
                include_notes  = TRUE,
                include_hidden = TRUE,    # 隠しも含める
                bundle_meta    = TRUE,
                meta_level     = "full",
                auto_update_manifest = TRUE,
                auto_update_session  = TRUE)
}
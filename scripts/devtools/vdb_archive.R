# scripts/devtools/vdb_archive.R
suppressWarnings(suppressMessages({
  library(jsonlite); library(digest); library(tools)
}))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  store = "D:/DMP_vdb_store",
  legacy_tag = "legacy",
  also_keep_query_in_legacy = TRUE
)
for (a in args) {
  if (grepl("^--store=", a)) opt$store <- sub("^--store=", "", a)
  if (grepl("^--legacy_tag=", a)) opt$legacy_tag <- sub("^--legacy_tag=", "", a)
  if (grepl("^--also_keep_query_in_legacy=", a)) {
    opt$also_keep_query_in_legacy <- as.logical(sub("^--also_keep_query_in_legacy=", "", a))
  }
}

dir.create(opt$store, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$store, opt$legacy_tag), showWarnings = FALSE, recursive = TRUE)

hw <- list(cpu="Ryzen 7 5700G", ram_gb=32, gpu="RTX 4060 Ti")
commit <- tryCatch(system("git rev-parse --short HEAD", intern=TRUE), error=function(e) "unknown")
stamp  <- format(Sys.time(), "%Y%m%d%H%M")

vdbs <- Sys.glob("data/virtual_db_u100_S*_seed*.rds")
if (length(vdbs) == 0) {
  message("[INFO] no VDB under data/. nothing to do.")
  quit(save="no")
}

for (rds in vdbs) {
  bn <- basename(rds)
  m <- regexec("S([0-9]+)_seed([0-9]+)\\.rds$", bn); mm <- regmatches(bn, m)[[1]]
  size <- as.numeric(mm[2]); seed <- as.integer(mm[3])
  
  # legacy 判定（名前に *_legacy_* を含む or legacy配下）
  is_legacy <- grepl("_legacy_", bn, fixed=TRUE) || grepl("legacy", rds)
  target_dir <- if (is_legacy) file.path(opt$store, opt$legacy_tag) else opt$store
  dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)
  
  # メタ生成（data/ に出す→ ZIP後に掃除）
  meta <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    host = Sys.info()[["nodename"]],
    os   = sprintf("%s %s", Sys.info()[["sysname"]], Sys.info()[["release"]]),
    cpu  = hw$cpu, ram_gb = hw$ram_gb, gpu = hw$gpu,
    commit = commit, seed = seed, size = size,
    db_path = normalizePath(rds, winslash="/", mustWork = FALSE),
    db_sha256 = digest(rds, "sha256", file=TRUE)
  )
  meta_path <- sub("\\.rds$", ".meta.json", rds)
  jsonlite::write_json(meta, meta_path, pretty = TRUE, auto_unbox = TRUE)
  message(sprintf("[META] %s", meta_path))
  
  # 参照クエリ候補（今回の“誤当て”も含める）
  cand <- unique(c(
    sprintf("data/query_profile_seed%d_S1000.csv", seed),
    Sys.glob(sprintf("data/queries/query_seed%d_*", seed)),
    "data/query_profile_seed123_S1000.csv" # ←“誤当て”証跡
  ))
  cand <- cand[file.exists(cand)]
  
  # ZIP 作成（data/上の実体を束ねて store 側へ）
  zip_name <- sprintf("vdb_seed%d_S%d_%s.zip", seed, size, stamp)
  zip_path <- file.path(target_dir, zip_name)
  if (!requireNamespace("zip", quietly = TRUE)) {
    utils::zip(zipfile = zip_path, files = c(rds, meta_path, cand), flags = "-9")
  } else {
    zip::zipr(zipfile = zip_path, files = c(rds, meta_path, cand))
  }
  message(sprintf("[ZIP] %s", zip_path))
  
  # “誤当てクエリ”は legacy/ にも複製保持
  if (opt$also_keep_query_in_legacy && any(basename(cand) == "query_profile_seed123_S1000.csv")) {
    file.copy("data/query_profile_seed123_S1000.csv",
              file.path(opt$store, opt$legacy_tag, "query_profile_seed123_S1000.csv"),
              overwrite = TRUE)
  }
  
  # store 側へ生コピーはしない方針に変更（ZIPのみ残す）
  # data/ 側の原本は削除
  unlink(c(rds, meta_path), force = TRUE)
  if (length(cand)) {
    # 誤当てを legacy に残す以外は data/ から削除
    keep <- if (opt$also_keep_query_in_legacy) "data/query_profile_seed123_S1000.csv" else character(0)
    unlink(setdiff(cand, keep), force = TRUE)
  }
}
message("[OK] archived to store (ZIP only): ", opt$store)

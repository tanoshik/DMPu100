# --- test/test_prepare_profile.R ---

# プロジェクト直下へ寄せる（source() 実行時でも、RStudio実行時でも動く）
locate_this <- function() {
  # source() で呼ばれた場合
  pf <- parent.frame()
  of <- tryCatch(pf$ofile, error = function(e) NULL)
  if (!is.null(of)) return(normalizePath(of, winslash = "/"))
  
  # RStudio から直接実行された場合
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(p)) return(normalizePath(p, winslash = "/"))
  }
  
  stop("テストファイルの場所が特定できません。source() で呼ぶか、Rprojを開いてください。")
}

this_file <- locate_this()
proj_dir  <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
setwd(proj_dir)

source("scripts/utils_profile.R")

# test/test_prepare_profile.R がやっている内容と同等
test_profile <- list(D8S1179 = c("13", "15"))
print(prepare_profile(test_profile))
# 例）     locus allele1 allele2
#       D8S1179    1300    1500


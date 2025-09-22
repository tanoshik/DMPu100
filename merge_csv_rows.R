# 必要パッケージ
library(readr)
library(dplyr)

# 入力フォルダのパス（プロジェクトルートからの相対パス）
input_dir <- "output/smoke_20250922192042/1k_normal_off/_raw20"

# CSVファイル一覧を取得
file_list <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# 各ファイルを読み込んで結合（ヘッダなし）
merged_data <- lapply(file_list, function(file) {
  read_csv(file, col_names = FALSE)
}) |> bind_rows()

# 出力ファイル名（同じフォルダに保存）
output_file <- file.path(input_dir, "merged_output.csv")

# 書き出し（ヘッダなし）
write_csv(merged_data, output_file, col_names = FALSE)
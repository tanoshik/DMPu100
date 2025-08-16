## scripts/utils_profile.R
## u100/ANY でプロファイルを正規化する最小実装
## - 入力: named list 例 list(D8S1179 = c("13","15")) / data.frame でも可
## - 出力: data.frame(locus, allele1, allele2) すべて integer（100..9900 or 9999）
## - 仕様ポイント:
##    * "ANY","*" , NA,"" は ANY_CODE=9999
##    * 数値/文字は u100: floor(x*100 + 0.5)
##    * a1<=a2 に並べ替え（同値はそのまま）
##    * 片側しか無い場合は、strict=FALSEならホモ補完（a2<-a1）
##      strict=TRUEならエラー（必要に応じ変更）

ANY_CODE <- 9999L

# 1点を u100/ANY に変換
encode_u100_scalar <- function(x, strict = FALSE) {
  # ANY系
  if (is.null(x) || length(x) == 0) return(ANY_CODE)
  if (is.na(x) || isTRUE(x == "") || isTRUE(toupper(as.character(x)) %in% c("ANY","*"))) {
    return(ANY_CODE)
  }
  # 数値/文字 → 数値
  xv <- suppressWarnings(as.numeric(x))
  if (is.na(xv)) {
    if (strict) stop(sprintf("非数値を検出: [%s]", as.character(x)))
    return(ANY_CODE)
  }
  # 範囲ざっくりチェック（<0 または >=100 は異常とみなす）
  if (xv < 0 || xv >= 100) {
    if (strict) stop(sprintf("生値が範囲外: %s (期待: [0,100))", xv))
    return(ANY_CODE)
  }
  # 0.5以上は切り上げ（ties away from zero）
  as.integer(floor(xv * 100 + 0.5))
}

# ベクトル版（高速）
encode_u100_vec <- function(v, strict = FALSE) {
  sapply(v, encode_u100_scalar, strict = strict, USE.NAMES = FALSE)
}

# メイン: プロファイル正規化
# mode: "lenient"=寛容（不足は補完/ANY化） / "strict"=厳格（不足や不正でエラー）
prepare_profile <- function(x,
                            mode = c("lenient","strict"),
                            homo_to_any = FALSE) {
  mode <- match.arg(mode)
  strict <- identical(mode, "strict")
  
  # 入力を data.frame に寄せる
  df <- NULL
  if (is.data.frame(x)) {
    df <- x
  } else if (is.list(x) && !is.null(names(x))) {
    # named list: locus -> c(a1,a2) or length 1
    loci  <- names(x)
    a1 <- integer(length(loci))
    a2 <- integer(length(loci))
    for (i in seq_along(loci)) {
      alle <- x[[i]]
      if (length(alle) == 0) {
        if (strict) stop(sprintf("ローカス[%s]のアレルが空", loci[i]))
        alle <- c(NA, NA)
      } else if (length(alle) == 1) {
        # 片側のみ
        if (homo_to_any) {
          alle <- c(alle[1], NA)  # 片側は ANY 化（後で9999）
        } else {
          alle <- c(alle[1], alle[1])  # ホモ補完
        }
      } else {
        alle <- alle[1:2]
      }
      a1[i] <- encode_u100_scalar(alle[1], strict = strict)
      a2[i] <- encode_u100_scalar(alle[2], strict = strict)
    }
    df <- data.frame(locus = loci, allele1 = a1, allele2 = a2, stringsAsFactors = FALSE)
  } else {
    stop("prepare_profile: サポート外の入力。named list か data.frame を渡してください。")
  }
  
  # 列名を合わせる
  names(df) <- tolower(names(df))
  # locus 列 同定
  if (!("locus" %in% names(df))) {
    # 代替: marker を受けたら locus に
    if ("marker" %in% names(df)) {
      names(df)[names(df) == "marker"] <- "locus"
    } else {
      stop("prepare_profile: 'locus' 列が見つかりません（marker でも可）")
    }
  }
  # アレル列 同定
  # 許容: allele1|allele1_id, allele2|allele2_id
  if (!("allele1" %in% names(df))) {
    if ("allele1_id" %in% names(df)) names(df)[names(df) == "allele1_id"] <- "allele1"
  }
  if (!("allele2" %in% names(df))) {
    if ("allele2_id" %in% names(df)) names(df)[names(df) == "allele2_id"] <- "allele2"
  }
  if (!all(c("locus","allele1","allele2") %in% names(df))) {
    stop("prepare_profile: 必須列(locus, allele1, allele2)が不足しています")
  }
  
  # u100/ANY へエンコード（既に整数でも一応通す）
  df$allele1 <- encode_u100_vec(df$allele1, strict = strict)
  df$allele2 <- encode_u100_vec(df$allele2, strict = strict)
  
  # a1<=a2 に並べ替え（ANY=9999 なので通常は後ろへ流れる）
  swap <- df$allele1 > df$allele2
  if (any(swap)) {
    tmp <- df$allele1[swap]
    df$allele1[swap] <- df$allele2[swap]
    df$allele2[swap] <- tmp
  }
  
  # 値域チェック（全件）— 高速ベクトル判定
  ok <- (df$allele1 %in% c(ANY_CODE) | (df$allele1 >= 100L & df$allele1 <= 9900L)) &
    (df$allele2 %in% c(ANY_CODE) | (df$allele2 >= 100L & df$allele2 <= 9900L))
  if (!all(ok)) {
    bad_n <- sum(!ok)
    if (strict) stop(sprintf("u100/ANY 域外データを検出: %d 行", bad_n))
    # lenient: 域外は ANY に落とす
    df$allele1[!ok] <- ANY_CODE
    df$allele2[!ok] <- ANY_CODE
  }
  
  # 出力は data.frame(locus, allele1, allele2) / integer
  df$allele1 <- as.integer(df$allele1)
  df$allele2 <- as.integer(df$allele2)
  df[, c("locus","allele1","allele2")]
}

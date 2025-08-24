# scripts/gui/server_match_logic.R
# Result tab server logic (u100): run real matcher first; mock only if impossible.
# No multibyte chars in code/comments.

source("scripts/scoring_fast.R", local = TRUE)
source("scripts/matcher_fast.R", local = TRUE)

ANY_CODE <- 9999L

# ---- simple file logger ----
.logf <- function(...) {
  dir.create("output", showWarnings = FALSE, recursive = TRUE)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse=" "), "\n",
      file = "output/shiny_debug.log", append = TRUE)
}

# ---- tiny local helpers ----
.map_any <- function(x, any_code = ANY_CODE) {
  x <- as.character(x); x[is.na(x) | x == ""] <- "any"
  x[tolower(x) == "any"] <- as.character(any_code)
  suppressWarnings(as.integer(x))
}

.prepare_query_min <- function(q_raw, locus_order = NULL, any_code = ANY_CODE) {
  n <- names(q_raw)
  names(q_raw) <- sub("(?i)^locus$","Locus", n, perl=TRUE)
  names(q_raw) <- sub("(?i)^allele1$","allele1", names(q_raw), perl=TRUE)
  names(q_raw) <- sub("(?i)^allele2$","allele2", names(q_raw), perl=TRUE)
  stopifnot(all(c("Locus","allele1","allele2") %in% names(q_raw)))
  q <- q_raw[, c("Locus","allele1","allele2"), drop = FALSE]
  q$Locus   <- as.character(q$Locus)
  q$allele1 <- .map_any(q$allele1, any_code)
  q$allele2 <- .map_any(q$allele2, any_code)
  if (!is.null(locus_order) && length(locus_order) > 0L) {
    ord <- match(q$Locus, locus_order)
    q <- q[order(ord), , drop = FALSE]
  }
  rownames(q) <- NULL
  q
}

.build_db_index_v1 <- function(db_prep, locus_order, any_code = ANY_CODE) {
  req <- c("SampleID","Locus","Allele1","Allele2")
  stopifnot(is.data.frame(db_prep), all(req %in% names(db_prep)))
  SIDs <- unique(db_prep$SampleID)
  LIDs <- as.character(locus_order %||% unique(as.character(db_prep$Locus)))
  S <- length(SIDs); L <- length(LIDs)
  A1 <- matrix(any_code, nrow = S, ncol = L, dimnames = list(NULL, LIDs))
  A2 <- matrix(any_code, nrow = S, ncol = L, dimnames = list(NULL, LIDs))
  sid2i <- setNames(seq_len(S), SIDs)
  loc2j <- setNames(seq_len(L), LIDs)
  for (k in seq_len(nrow(db_prep))) {
    sid <- db_prep$SampleID[k]
    loc <- as.character(db_prep$Locus[k])
    i <- sid2i[[sid]]; j <- loc2j[[loc]]
    if (is.na(i) || is.na(j)) next
    A1[i, j] <- .map_any(db_prep$Allele1[k], any_code)
    A2[i, j] <- .map_any(db_prep$Allele2[k], any_code)
  }
  list(sample_ids = SIDs, locus_ids = LIDs, A1 = A1, A2 = A2)
}

## 追加：安全な取り出し・正規化
.as_int1 <- function(x) {
  if (is.null(x) || length(x) < 1L) return(NA_integer_)
  v <- suppressWarnings(as.integer(x[1]))
  if (is.na(v)) NA_integer_ else v
}
.norm_filter_type <- function(x) {
  if (is.null(x) || length(x) < 1L) return("all")
  v <- tolower(as.character(x[1]))
  # よくある値を正規化
  if (v %in% c("all","none","no_filter","all_items")) return("all")
  if (v %in% c("top_n","topn","top","top-n"))        return("top_n")
  if (v %in% c("score","score_min","score_ge","ge")) return("score")
  "all"
}

## 置き換え：単一選択仕様のフィルタ適用
.apply_summary_filters <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  
  # Confirm → rv に入っている想定（最新ZIP準拠）
  ft <- .norm_filter_type(rv$filter_type)
  tn <- .as_int1(rv$top_n)
  ms <- .as_int1(rv$score_min)
  
  # Score列の同定＆数値化
  scn <- "Score"
  if (!(scn %in% names(df))) {
    ix <- which(tolower(names(df)) == "score")
    if (length(ix) == 1L) scn <- names(df)[ix]
  }
  if (scn %in% names(df)) {
    suppressWarnings(df[[scn]] <- as.integer(df[[scn]]))
  }
  
  # 並びは一貫して Score 降順 → SampleID 昇順（安定化）
  if (scn %in% names(df)) {
    sid <- if ("SampleID" %in% names(df)) "SampleID" else names(df)[1L]
    o <- order(-df[[scn]], df[[sid]], na.last = TRUE)
    df <- df[o, , drop = FALSE]
  }
  
  # --- 単一選択で分岐（ANDにはしない） ---
  if (ft == "all") {
    # 何もしない（TopN/Score入力は無視）
    return(df)
  }
  
  if (ft == "top_n") {
    if (is.na(tn) || tn <= 0L) return(df)            # 入力無効なら素通し
    if (nrow(df) > tn) df <- df[seq_len(tn), , drop = FALSE]
    rownames(df) <- NULL
    return(df)
  }
  
  if (ft == "score") {
    if (is.na(ms)) return(df)                        # 入力無効なら素通し
    if (scn %in% names(df)) {
      df <- df[df[[scn]] >= ms, , drop = FALSE]
      # Score降順→SampleID昇順は維持（上で並べ直しているのでそのまま）
    }
    rownames(df) <- NULL
    return(df)
  }
  
  # 不明値はALL扱い
  df
}

server_match_logic <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    
    # get locus order
    .get_locus_order <- function() {
      if (!is.null(rv$locus_order)) return(rv$locus_order)
      if (exists("locus_order", envir = .GlobalEnv, inherits = FALSE)) return(get("locus_order", envir = .GlobalEnv))
      if (file.exists("data/locus_order.rds")) return(tryCatch(readRDS("data/locus_order.rds"), error=function(e) character(0)))
      character(0)
    }
    
    # build db_index on-demand
    .ensure_db_index <- function() {
      if (!is.null(rv$db_index)) return(rv$db_index)
      if (is.null(rv$db_profile)) return(NULL)
      
      dbp <- rv$db_profile
      n <- names(dbp)
      names(dbp) <- sub("(?i)^sample[_ ]?id$","SampleID", n, perl=TRUE)
      names(dbp) <- sub("(?i)^locus$","Locus", names(dbp), perl=TRUE)
      names(dbp) <- sub("(?i)^allele1$","Allele1", names(dbp), perl=TRUE)
      names(dbp) <- sub("(?i)^allele2$","Allele2", names(dbp), perl=TRUE)
      dbp <- dbp[, c("SampleID","Locus","Allele1","Allele2"), drop=FALSE]
      dbp$Allele1 <- .map_any(dbp$Allele1)
      dbp$Allele2 <- .map_any(dbp$Allele2)
      
      locs <- .get_locus_order()
      if (length(locs) == 0L) locs <- unique(dbp$Locus)
      
      idx <- .build_db_index_v1(dbp, locs)
      rv$db_index <- idx
      idx
    }
    
    # prepare query from rv (robust: show/std/df/list に対応)
    .prepare_query_from_rv <- function() {
      norm_cols <- function(df) {
        n <- names(df)
        names(df) <- sub("(?i)^locus$","Locus", n, perl=TRUE)
        names(df) <- sub("(?i)^allele1$","allele1", names(df), perl=TRUE)
        names(df) <- sub("(?i)^allele2$","allele2", names(df), perl=TRUE)
        df
      }
      # 1) show（Confirm表示用）を最優先
      if (!is.null(rv$query_profile_show) && is.data.frame(rv$query_profile_show)) {
        qdf <- norm_cols(rv$query_profile_show)
        if (all(c("Locus","allele1","allele2") %in% names(qdf))) {
          return(.prepare_query_min(qdf, .get_locus_order()))
        }
        # Allele1/Allele2 大文字のまま来る場合
        if (all(c("Locus","Allele1","Allele2") %in% names(qdf))) {
          qdf <- data.frame(
            Locus   = as.character(qdf$Locus),
            allele1 = qdf$Allele1,
            allele2 = qdf$Allele2,
            stringsAsFactors = FALSE
          )
          return(.prepare_query_min(qdf, .get_locus_order()))
        }
      }
      # 2) std（内部用）も見る
      if (!is.null(rv$query_profile_std) && is.data.frame(rv$query_profile_std)) {
        qdf <- norm_cols(rv$query_profile_std)
        if (all(c("Locus","allele1","allele2") %in% names(qdf))) {
          return(.prepare_query_min(qdf, .get_locus_order()))
        }
        if (all(c("Locus","Allele1","Allele2") %in% names(qdf))) {
          qdf <- data.frame(
            Locus   = as.character(qdf$Locus),
            allele1 = qdf$Allele1,
            allele2 = qdf$Allele2,
            stringsAsFactors = FALSE
          )
          return(.prepare_query_min(qdf, .get_locus_order()))
        }
      }
      # 3) 既存の候補も踏襲
      if (!is.null(rv$query_profile_df) && is.data.frame(rv$query_profile_df)) {
        qdf <- norm_cols(rv$query_profile_df)
        if (all(c("Locus","allele1","allele2") %in% names(qdf))) {
          return(.prepare_query_min(qdf, .get_locus_order()))
        }
      }
      if (!is.null(rv$query_df) && is.data.frame(rv$query_df)) {
        qdf <- norm_cols(rv$query_df)
        if (all(c("Locus","allele1","allele2") %in% names(qdf))) {
          return(.prepare_query_min(qdf, .get_locus_order()))
        }
      }
      # 4) list形式（将来拡張）: list(Locus=..., alleles=list(c(a1,a2),...))
      if (!is.null(rv$query_list) && is.list(rv$query_list)) {
        loci <- as.character(rv$query_list$Locus)
        alle <- do.call(rbind, lapply(rv$query_list$alleles, function(p) {
          p <- as.integer(p); if (length(p) < 2L) p <- c(p, ANY_CODE); p[1:2]
        }))
        q_raw <- data.frame(Locus = loci, allele1 = alle[,1], allele2 = alle[,2], stringsAsFactors = FALSE)
        return(.prepare_query_min(q_raw, .get_locus_order()))
      }
      NULL
    }
    
    # ---- UI-side debug status (always visible) ----
    output$result_status <- renderText({
      paste0("[status] ", rv$status_msg %||% "(none)",
             " | has_index=", !is.null(rv$db_index),
             " | n_summary=", if (is.null(rv$match_scores)) NA_integer_ else nrow(rv$match_scores),
             " | n_detail=",  if (is.null(rv$match_detail))  NA_integer_ else nrow(rv$match_detail))
    })
    
    # ---- central runner ----
    .run_match_and_set <- function() {
      if (!exists("run_match_fast", mode = "function")) {
        rv$status_msg <- "run_match_fast() not available (source scripts/matcher_fast.R)."
        .logf("[run] abort: run_match_fast missing")
        return(invisible(FALSE))
      }
      
      .logf("[run] invoked")
      
      # build prerequisites
      idx <- .ensure_db_index()
      if (is.null(idx)) {
        rv$status_msg <- "No database loaded (rv$db_profile / rv$db_index not found)."
        rv$match_detail <- data.frame(); rv$match_scores <- data.frame()
        .logf("[run] abort: no db_index")
        return(invisible(FALSE))
      }
      
      q <- .prepare_query_from_rv()
      if (is.null(q) || nrow(q) == 0L) {
        rv$status_msg <- "No prepared query available (rv$query_* empty)."
        rv$match_detail <- data.frame(); rv$match_scores <- data.frame()
        .logf("[run] abort: no query")
        return(invisible(FALSE))
      }
      
      # run matcher
      t0 <- Sys.time()
      rs <- tryCatch(
        run_match_fast(q, idx, pre_add_for_any_any = 2L, include_bits_in_detail = TRUE),
        error = function(e) {
          rv$status_msg <- paste0("run_match_fast failed: ", conditionMessage(e))
          rv$match_detail <- data.frame(); rv$match_scores <- data.frame()
          .logf("[run] error:", conditionMessage(e))
          return(NULL)
        }
      )
      if (is.null(rs)) return(invisible(FALSE))
      
      # set results (apply Confirm filters here)
      rv$match_detail <- rs$detail
      df_sum <- .apply_summary_filters(rs$summary)
      rv$match_scores <- df_sum
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      rv$status_msg <- paste0("Run OK (", rs$used, "). rows=", nrow(df_sum), " (from ", nrow(rs$summary), "). elapsed=", sprintf("%.2fs", el))
      .logf("[run] OK rows=", nrow(df_sum), " from=", nrow(rs$summary), " elapsed=", sprintf("%.2fs", el))
      rv$result_ready <- Sys.time()  # optional trigger
      invisible(TRUE)
    }
    
    # trigger wiring
    observeEvent(rv$trigger_run_match, {
      ok <- .run_match_and_set()
      if (!isTRUE(ok)) { }  # status_msg already set
    }, ignoreInit = TRUE)
    
    observeEvent(input$btn_run_match_local, {
      ok <- .run_match_and_set()
      if (!isTRUE(ok)) { }
    }, ignoreInit = TRUE)
    
    # ---- Result table rendering (Summary / Detail) ----
    .empty_df <- function() data.frame()
    
    output$tbl_result <- DT::renderDT({
      view <- input$result_view %||% "Summary"
      if (identical(view, "Summary")) {
        df <- rv$match_scores
      } else {
        df <- rv$match_detail
      }
      if (is.null(df)) df <- .empty_df()
      DT::datatable(
        df,
        rownames = FALSE,
        options = list(pageLength = 20, scrollX = TRUE)
      )
    })
    
    # safe DF flattener for CSV
    .write_safe_df <- function(df) {
      if (is.null(df)) return(data.frame())
      if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
      for (j in seq_along(df)) {
        if (is.list(df[[j]])) {
          df[[j]] <- vapply(df[[j]], function(x) paste(as.character(x), collapse=","), character(1))
        } else if (is.factor(df[[j]])) {
          df[[j]] <- as.character(df[[j]])
        }
      }
      row.names(df) <- NULL
      df
    }
    
    # ---- Downloads with timestamp (Asia/Tokyo) ----
    .ts <- function() {
      # タイムスタンプはJSTで
      old <- Sys.getenv("TZ", unset = NA)
      on.exit(if (!is.na(old)) Sys.setenv(TZ = old) else Sys.unsetenv("TZ"), add = TRUE)
      Sys.setenv(TZ = "Asia/Tokyo")
      format(Sys.time(), "%Y%m%d_%H%M%S")
    }
    
    output$dl_scores <- shiny::downloadHandler(
      filename = function() paste0("match_scores_", .ts(), ".csv"),
      content = function(file) {
        x <- .write_safe_df(rv$match_scores)
        if (is.null(x) || !nrow(x)) {
          x <- data.frame(SampleID=character(0), Score=integer(0), stringsAsFactors = FALSE)
        }
        utils::write.table(x, file = file, sep = ",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")
      }
    )
    
    output$dl_detail <- shiny::downloadHandler(
      filename = function() paste0("match_log_", .ts(), ".csv"),
      content = function(file) {
        x <- rv$match_detail
        if (is.null(x) || !nrow(x)) {
          x <- data.frame(SampleID=character(0), Locus=character(0), DB_Allele1=integer(0),
                          DB_Allele2=integer(0), Score=integer(0),
                          Bits=character(0), Code=integer(0),
                          stringsAsFactors = FALSE)
        }
        # list/unknown列で落ちないよう型を落とし切ってから保存
        to_chr <- c("SampleID","Locus","Bits")
        to_int <- c("DB_Allele1","DB_Allele2","Score","Code")
        for (nm in intersect(names(x), to_chr)) x[[nm]] <- as.character(x[[nm]])
        for (nm in intersect(names(x), to_int)) suppressWarnings(x[[nm]] <- as.integer(x[[nm]]))
        utils::write.table(x, file = file, sep = ",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")
      }
    )
  })
}

# scripts/gui/server_input_logic.R
# Input server:
# - DMP-aligned layout with CSV multi-sample pager.
# - Robust CSV reader (single- or multi-sample).
# - Validation against freq_table per locus; invalid selectors highlighted in red.
# - Navigation gating: show "Go to Confirm" only if all valid; otherwise show "Prepare Query".
# - Prepare Query auto-fixes invalid loci: if one allele invalid but the other valid -> result "valid, any";
#   if both invalid -> "any, any". Then re-validate and toggle buttons.
# No multibyte characters in code/comments.

# --- helpers: safe enable/disable even without shinyjs ---
.safe_enable  <- function(id) { if (requireNamespace("shinyjs", quietly = TRUE)) shinyjs::enable(id) }
.safe_disable <- function(id) { if (requireNamespace("shinyjs", quietly = TRUE)) shinyjs::disable(id) }

server_input_logic <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    
    `%||%` <- function(x, y) if (is.null(x) || (is.character(x) && !nzchar(x))) y else x
    sanitize_id <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)
    
    # ---- load dictionaries ----
    lo_path <- file.path("data","locus_order.rds")
    if (!file.exists(lo_path)) stop("Missing data/locus_order.rds")
    locus_order <- readRDS(lo_path)
    
    ft_path <- file.path("data","freq_table.rds")
    freq_table <- if (file.exists(ft_path)) tryCatch(readRDS(ft_path), error = function(e) NULL) else NULL
    
    # map of allowed alleles per locus (character vector, without "any")
    build_allowed_map <- function(loci, ft) {
      out <- setNames(vector("list", length(loci)), loci)
      if (is.null(ft)) {
        for (lc in loci) out[[lc]] <- character(0)
        return(out)
      }
      nm <- names(ft)
      locus_col  <- if ("Locus" %in% nm) "Locus" else nm[1]
      allele_col <- if ("Allele" %in% nm) "Allele" else setdiff(nm, locus_col)[1]
      for (lc in loci) {
        sub <- tryCatch(ft[ft[[locus_col]] == lc, allele_col, drop = TRUE], error = function(e) character(0))
        alleles <- unique(na.omit(as.character(sub)))
        alleles <- trimws(alleles)   # <- add: trim for safety
        out[[lc]] <- alleles
      }
      out
    }
    
    # choices used for UI ("any" + allowed sorted)
    build_choices_map <- function(loci, allowed) {
      out <- setNames(vector("list", length(loci)), loci)
      for (lc in loci) {
        alleles <- allowed[[lc]]
        if (length(alleles)) {
          num <- suppressWarnings(as.numeric(alleles))
          ord <- order(is.na(num), num, alleles)
          alleles <- alleles[ord]
        }
        out[[lc]] <- c("any", alleles)
      }
      out
    }
    
    allowed_map <- build_allowed_map(locus_order, freq_table)
    choices_map <- build_choices_map(locus_order, allowed_map)
    
    # ---- layout (left/right)
    load_layout <- function(path, locus_order) {
      if (!file.exists(path)) return(NULL)
      obj <- tryCatch(readRDS(path), error = function(e) NULL)
      if (is.null(obj)) return(NULL)
      norm_side <- function(x) {
        x <- as.character(x); x <- trimws(tolower(x))
        ifelse(x %in% c("l","left"), "L", ifelse(x %in% c("r","right"), "R", NA_character_))
      }
      if (is.list(obj) && all(c("left","right") %in% names(obj))) {
        left_raw  <- as.character(obj$left)
        right_raw <- as.character(obj$right)
      } else if (is.data.frame(obj) && all(c("Locus","Side") %in% names(obj))) {
        side <- norm_side(obj$Side)
        left_raw  <- as.character(obj$Locus[side == "L"])
        right_raw <- as.character(obj$Locus[side == "R"])
      } else {
        return(NULL)
      }
      left  <- locus_order[locus_order %in% left_raw]
      right <- locus_order[locus_order %in% right_raw]
      if (length(right)) right <- right[!right %in% left]
      assigned <- c(left, right)
      missing  <- locus_order[!locus_order %in% assigned]
      if (length(missing)) {
        for (i in seq_along(missing)) {
          if (length(left) <= length(right)) left <- c(left, missing[i]) else right <- c(right, missing[i])
        }
      }
      list(left = left[!is.na(left)], right = right[!is.na(right)])
    }
    
    layout_path <- file.path("data","locus_layout.rds")
    layout <- load_layout(layout_path, locus_order)
    if (!is.null(layout)) {
      left_loci  <- layout$left
      right_loci <- layout$right
    } else {
      left_loci  <- locus_order[seq(1, length(locus_order), by = 2)]
      right_loci <- locus_order[seq(2, length(locus_order), by = 2)]
    }
    
    # ---- UI builders (no labels; narrow width)
    locus_row_ui <- function(locus) {
      lid <- sanitize_id(locus)
      ns  <- session$ns
      shiny::div(id = ns(paste0("row_", lid)), class = "dmp-locus-row",
                 shiny::div(class = "dmp-locus-label", locus),
                 shiny::div(id = ns(paste0("box_a1_", lid)), class = "dmp-a1",
                            shiny::selectizeInput(
                              ns(paste0("a1_", lid)), label = NULL,
                              choices = choices_map[[locus]], selected = "any",
                              options = list(placeholder = "A1"), width = "100px"
                            )
                 ),
                 shiny::div(id = ns(paste0("box_a2_", lid)), class = "dmp-a2",
                            shiny::selectizeInput(
                              ns(paste0("a2_", lid)), label = NULL,
                              choices = choices_map[[locus]], selected = "any",
                              options = list(placeholder = "A2"), width = "100px"
                            )
                 )
      )
    }
    
    output$allele_grid_left  <- shiny::renderUI({
      do.call(shiny::tagList, lapply(left_loci,  locus_row_ui))
    })
    output$allele_grid_right <- shiny::renderUI({
      do.call(shiny::tagList, lapply(right_loci, locus_row_ui))
    })
    
    # =========================================================================
    # CSV multi-sample pager (robust)
    # =========================================================================
    rv$query_samples       <- NULL   # named list of data.frames (Locus, Allele1, Allele2)
    rv$query_sample_names  <- NULL   # names vector
    rv$query_sample_index  <- 0L
    
    find_col <- function(nm, candidates) {
      nm_low <- tolower(nm); cand_low <- tolower(candidates)
      hit <- which(nm_low %in% cand_low)[1]
      if (is.na(hit)) return(NULL); nm[hit]
    }
    
    read_query_csv_multi <- function(path, display_name = NULL, fallback_name = "Query") {
      df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                     error = function(e) NULL)
      if (is.null(df) || !nrow(df)) return(NULL)
      
      for (j in seq_along(df)) if (is.character(df[[j]])) df[[j]] <- trimws(df[[j]])
      nm <- names(df)
      
      sample_col <- find_col(nm, c("Sample","SampleID","sample","sampleid"))
      locus_col  <- find_col(nm, c("Locus","locus"))
      a1_col     <- find_col(nm, c("Allele1","Allele_1","allele1","allele_1"))
      a2_col     <- find_col(nm, c("Allele2","Allele_2","allele2","allele_2"))
      if (is.null(locus_col) || is.null(a1_col) || is.null(a2_col)) return(NULL)
      
      df[[locus_col]] <- as.character(df[[locus_col]])
      df[[a1_col]]    <- as.character(df[[a1_col]])
      df[[a2_col]]    <- as.character(df[[a2_col]])
      
      # single-sample
      if (is.null(sample_col)) {
        sub <- df[, c(locus_col, a1_col, a2_col)]
        names(sub) <- c("Locus","Allele1","Allele2")
        sub <- sub[sub$Locus %in% locus_order, , drop = FALSE]
        sub$Locus <- factor(sub$Locus, levels = locus_order)
        sub <- sub[order(sub$Locus), ]
        sub$Locus <- as.character(sub$Locus)
        
        # single-sample 分岐内の base_name 箇所を置換
        base_name <- if (!is.null(display_name) && nzchar(display_name)) {
          tools::file_path_sans_ext(display_name)
        } else {
          tryCatch(tools::file_path_sans_ext(basename(path)), error = function(e) fallback_name)
        }
        lst <- list(sub)
        names(lst) <- base_name
        return(list(samples = lst))
      }
      
      # multi-sample (long)
      df[[sample_col]] <- as.character(df[[sample_col]])
      df <- df[, c(sample_col, locus_col, a1_col, a2_col)]
      names(df) <- c("Sample","Locus","Allele1","Allele2")
      df <- df[df$Locus %in% locus_order, , drop = FALSE]
      if (!nrow(df)) return(NULL)
      
      ids <- unique(df$Sample)
      out_list <- vector("list", length(ids)); names(out_list) <- ids
      for (i in seq_along(ids)) {
        sid <- ids[i]
        sub <- df[df$Sample == sid, c("Locus","Allele1","Allele2"), drop = FALSE]
        sub$Locus <- factor(sub$Locus, levels = locus_order)
        sub <- sub[order(sub$Locus), ]
        sub$Locus <- as.character(sub$Locus)
        out_list[[i]] <- sub
      }
      list(samples = out_list, names = ids)
    }
    
    apply_sample_to_selectors <- function(idx) {
      if (is.null(rv$query_samples) || length(rv$query_samples) == 0L) return(invisible(NULL))
      n <- length(rv$query_samples)
      idx <- max(1L, min(as.integer(idx), n))
      
      pf <- rv$query_samples[[idx]]
      nm <- names(rv$query_samples)[idx] %||% "Query"   # <- ここを厳密に
      
      shiny::updateTextInput(session, "txt_sample_name", value = nm)
      
      for (i in seq_len(nrow(pf))) {
        locus <- as.character(pf$Locus[i]); if (!locus %in% locus_order) next
        lid <- sanitize_id(locus)
        a1 <- pf$Allele1[i] %||% "any"; a2 <- pf$Allele2[i] %||% "any"
        a1 <- if (!nzchar(a1)) "any" else as.character(a1)
        a2 <- if (!nzchar(a2)) "any" else as.character(a2)
        
        ch <- choices_map[[locus]]
        if (!a1 %in% ch) ch <- c(ch, a1)
        if (!a2 %in% ch) ch <- c(ch, a2)
        
        shiny::updateSelectizeInput(session, paste0("a1_", lid), choices = ch, selected = a1, server = TRUE)
        shiny::updateSelectizeInput(session, paste0("a2_", lid), choices = ch, selected = a2, server = TRUE)
      }
      
      rv$query_sample_index <- idx
      output$txt_q_index <- shiny::renderText(sprintf("%d / %d", idx, n))
    }
    
    shiny::observeEvent(input$file_query_csv, {
      fi <- input$file_query_csv
      if (is.null(fi) || !nzchar(fi$datapath)) return(invisible(NULL))
      
      parsed <- read_query_csv_multi(
        path = fi$datapath,
        display_name = fi$name  # ← ここがポイント：元のファイル名
      )
      if (is.null(parsed) || length(parsed$samples) == 0L) return(invisible(NULL))
      
      rv$query_samples      <- parsed$samples
      rv$query_sample_names <- names(parsed$samples)  # ← list の名前をそのまま採用
      apply_sample_to_selectors(1L)
    }, ignoreInit = TRUE)
    
    shiny::observeEvent(input$btn_q_prev, {
      if (is.null(rv$query_samples) || length(rv$query_samples) == 0L) return(invisible(NULL))
      n  <- length(rv$query_samples)
      ix <- rv$query_sample_index; if (!is.integer(ix) || ix <= 0L) ix <- 1L
      new_ix <- if (ix <= 1L) n else (ix - 1L)
      apply_sample_to_selectors(new_ix)
    })
    shiny::observeEvent(input$btn_q_next, {
      if (is.null(rv$query_samples) || length(rv$query_samples) == 0L) return(invisible(NULL))
      n  <- length(rv$query_samples)
      ix <- rv$query_sample_index; if (!is.integer(ix) || ix <= 0L) ix <- 1L
      new_ix <- if (ix >= n) 1L else (ix + 1L)
      apply_sample_to_selectors(new_ix)
    })
    
    # =========================================================================
    # Validation + gating
    # =========================================================================
    # get current selectors as data.frame (Locus, Allele1, Allele2)
    current_selectors <- shiny::reactive({
      loc <- character(0); a1v <- character(0); a2v <- character(0)
      for (locus in locus_order) {
        lid <- sanitize_id(locus)
        v1 <- input[[paste0("a1_", lid)]] %||% "any"
        v2 <- input[[paste0("a2_", lid)]] %||% "any"
        loc <- c(loc, locus); a1v <- c(a1v, as.character(v1)); a2v <- c(a2v, as.character(v2))
      }
      data.frame(Locus = loc, Allele1 = a1v, Allele2 = a2v, stringsAsFactors = FALSE)
    })
    
    # validity map per locus/allele
    validate_selectors <- function(df) {
      invalid <- rep(FALSE, nrow(df) * 2L)
      dim(invalid) <- c(nrow(df), 2L) # [,1]=A1, [,2]=A2
      for (i in seq_len(nrow(df))) {
        lc <- df$Locus[i]
        a1 <- df$Allele1[i]; a2 <- df$Allele2[i]
        allowed <- allowed_map[[lc]]
        invalid[i, 1] <- !(identical(a1, "any") || a1 %in% allowed)
        invalid[i, 2] <- !(identical(a2, "any") || a2 %in% allowed)
      }
      invalid
    }
    
    # apply red highlight and toggle buttons
    # add/remove red border on actual selectize controls and show messages
    apply_validation_ui <- function(df, invalid) {
      any_invalid <- any(invalid)
      
      # buttons: emphasize Prepare when invalid
      if (any_invalid) {
        shinyjs::show(id = "btn_prepare")
        shinyjs::hide(id = "btn_goto_confirm")
        shinyjs::addClass(id = "btn_prepare", class = "btn-danger")
      } else {
        shinyjs::hide(id = "btn_prepare")
        shinyjs::show(id = "btn_goto_confirm")
        shinyjs::removeClass(id = "btn_prepare", class = "btn-danger")
      }
      
      # messages only (no selector highlighting)
      msgs <- character(0)
      for (i in seq_len(nrow(df))) {
        lc <- df$Locus[i]
        if (invalid[i,1]) msgs <- c(msgs, sprintf("%s - Allele1 is not in frequency table", lc))
        if (invalid[i,2]) msgs <- c(msgs, sprintf("%s - Allele2 is not in frequency table", lc))
      }
      if (length(msgs)) {
        output$val_msgs <- shiny::renderUI(HTML(paste(msgs, collapse = "<br/>")))
      } else {
        output$val_msgs <- shiny::renderUI(HTML(""))
      }
    }
    
    # react to any selector change
    shiny::observe({
      df <- current_selectors()
      inv <- validate_selectors(df)
      apply_validation_ui(df, inv)
    })
    
    # =========================================================================
    # Prepare (auto-fix invalid), Go to Confirm
    # =========================================================================
    # auto-fix rule:
    # - if A1 invalid and A2 valid -> (A1 <- A2, A2 <- "any")
    # - if A1 valid and A2 invalid -> (A2 <- "any")
    # - if both invalid -> (A1 <- "any", A2 <- "any")
    auto_fix <- function(df) {
      out <- df
      for (i in seq_len(nrow(df))) {
        lc <- df$Locus[i]
        a1 <- df$Allele1[i]; a2 <- df$Allele2[i]
        allowed <- allowed_map[[lc]]
        v1 <- identical(a1, "any") || a1 %in% allowed
        v2 <- identical(a2, "any") || a2 %in% allowed
        if (v1 && v2) next
        if (!v1 && v2) {
          out$Allele1[i] <- a2
          out$Allele2[i] <- "any"
        } else if (v1 && !v2) {
          out$Allele2[i] <- "any"
        } else {
          out$Allele1[i] <- "any"
          out$Allele2[i] <- "any"
        }
      }
      out
    }
    
    # update selectors from df (no change to choices_map)
    set_selectors_from_df <- function(df) {
      for (i in seq_len(nrow(df))) {
        lc <- df$Locus[i]; if (!lc %in% locus_order) next
        lid <- sanitize_id(lc)
        a1 <- df$Allele1[i] %||% "any"; a2 <- df$Allele2[i] %||% "any"
        ch <- choices_map[[lc]]
        if (!a1 %in% ch) ch <- c(ch, a1)
        if (!a2 %in% ch) ch <- c(ch, a2)
        shiny::updateSelectizeInput(session, paste0("a1_", lid), choices = ch, selected = a1, server = TRUE)
        shiny::updateSelectizeInput(session, paste0("a2_", lid), choices = ch, selected = a2, server = TRUE)
      }
    }
    
    shiny::observeEvent(input$btn_prepare, {
      df  <- current_selectors()
      df2 <- auto_fix(df)
      set_selectors_from_df(df2)
      # re-validate UI (observer will run automatically too)
      inv <- validate_selectors(df2)
      apply_validation_ui(df2, inv)
    })
    
    # auto enable/disable Go to Confirm based on prepared query
    observe({
      ns <- session$ns
      ok <- !is.null(rv$query_profile_show) &&
        is.data.frame(rv$query_profile_show) &&
        nrow(rv$query_profile_show) > 0
      if (ok) .safe_enable(ns("btn_goto_confirm")) else .safe_disable(ns("btn_goto_confirm"))
    })
    
    # build standardized profile and push to Confirm
    # server_input_logic.R （モジュール内）
    shiny::observeEvent(input$btn_goto_confirm, {
      ns <- session$ns
      .safe_disable(ns("btn_goto_confirm"))
      on.exit({ .safe_enable(ns("btn_goto_confirm")) }, add = TRUE)
      
      # --- ここから既存の処理そのまま ---
      # only proceed when valid (button hidden otherwise)
      df <- current_selectors()
      rv$query_profile_std <- prepare_profile_df(
        profile_df = df,
        locus_order = locus_order,
        homo_to_any = isTRUE(input$chk_homo_to_any_query),
        lenient = FALSE,
        fill_missing_any = TRUE
      )
      rv$query_profile_show <- format_profile_df(rv$query_profile_std)
      rv$sample_name        <- input$txt_sample_name %||% "Query"
      
      output$tbl_confirm_query <- DT::renderDT({
        rv$query_profile_show
      }, options = list(pageLength = 25, scrollX = TRUE))
      
      output$txt_locus_count <- shiny::renderText({
        sprintf("Loci in query: %d", if (!is.null(rv$query_profile_std)) nrow(rv$query_profile_std) else 0L)
      })
      output$txt_db_count <- shiny::renderText({
        n <- if (!is.null(rv$db_std)) nrow(rv$db_std) else 0L
        sprintf("Database rows (standardized): %d", n)
      })
      
      rv$nav_request <- "Confirm"
      # --- ここまで既存の処理そのまま ---
    })
  })
}

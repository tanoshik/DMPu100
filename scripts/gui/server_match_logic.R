# scripts/gui/server_match_logic.R
# Result tab server module (run_match_fast integration with safe fallback)
# No multibyte characters in code/comments.

# ---- Helpers: common schema ----
# Build summary from detail (SampleID/Score total)
.build_scores_from_detail <- function(detail_df) {
  if (is.null(detail_df) || !nrow(detail_df)) {
    return(data.frame(SampleID = character(), Score = integer(), stringsAsFactors = FALSE))
  }
  agg <- stats::aggregate(Score ~ SampleID, data = detail_df, FUN = sum)
  names(agg) <- c("SampleID", "Score")
  agg[order(-agg$Score, agg$SampleID), , drop = FALSE]
}

# Keep only rows with required columns and types
.normalize_detail_schema <- function(df) {
  keep <- c("SampleID", "Locus", "DB_Allele1", "DB_Allele2", "Score")
  if (is.null(df) || !all(keep %in% names(df))) {
    # Return empty frame with proper columns
    out <- data.frame(matrix(ncol = length(keep), nrow = 0))
    names(out) <- keep
    return(out)
  }
  df$SampleID   <- as.character(df$SampleID)
  df$Locus      <- as.character(df$Locus)
  df$DB_Allele1 <- as.character(df$DB_Allele1)
  df$DB_Allele2 <- as.character(df$DB_Allele2)
  df$Score      <- as.integer(df$Score)
  df[, keep, drop = FALSE]
}

# ---- Fallback mock (used when DB missing or run_match_fast not available) ----
.build_mock_detail <- function(sample_ids, loci) {
  if (length(sample_ids) == 0L || length(loci) == 0L) {
    return(.normalize_detail_schema(NULL))
  }
  rows <- vector("list", length(sample_ids) * length(loci))
  k <- 0L
  for (sid in sample_ids) {
    for (lc in loci) {
      k <- k + 1L
      rows[[k]] <- list(
        SampleID   = sid,
        Locus      = lc,
        DB_Allele1 = "12",
        DB_Allele2 = "14",
        Score      = sample(c(0L, 1L, 2L), size = 1L)
      )
    }
  }
  df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  .normalize_detail_schema(df)
}

# ---- Adapter: call run_match_fast and coerce to detail schema ----
# Expected outputs supported:
#  A) detail-like data.frame with required columns already present
#  B) per-locus log with columns that can be renamed to required ones
#  C) per-sample scores only (then we synthesize a minimal detail)
.run_match_and_adapt <- function(query_df, db_df, locus_order) {
  # Ensure function exists
  if (!exists("run_match_fast", mode = "function")) {
    return(list(
      detail  = .build_mock_detail(sample_ids = sprintf("S%03d", 1:12), loci = unique(query_df$Locus)),
      summary = NULL, used = "mock:no_run_match_fast"
    ))
  }
  # Defensive checks
  if (is.null(query_df) || !nrow(query_df)) {
    stop("Query profile (query_df) is empty.")
  }
  if (is.null(db_df) || !nrow(db_df)) {
    stop("Database (db_df) is empty.")
  }
  
  # Try to call run_match_fast with common signatures
  # NOTE: Adjust these calls if your scoring_fast.R signature differs.
  res <- tryCatch({
    # Pattern 1: run_match_fast(query_df, db_df, locus_order)
    run_match_fast(query_df, db_df, locus_order)
  }, error = function(e1) {
    tryCatch({
      # Pattern 2: run_match_fast(query_df = query_df, db_df = db_df, locus_order = locus_order)
      run_match_fast(query_df = query_df, db_df = db_df, locus_order = locus_order)
    }, error = function(e2) {
      structure(list(`.err` = paste("run_match_fast failed:", conditionMessage(e2))), class = "runmatch_error")
    })
  })
  
  # If failed, fallback to mock
  if (inherits(res, "runmatch_error") || (is.null(res) && !is.data.frame(res))) {
    return(list(
      detail  = .build_mock_detail(sample_ids = sprintf("S%03d", 1:12), loci = unique(query_df$Locus)),
      summary = NULL, used = "mock:run_match_fast_error"
    ))
  }
  
  # Cases:
  # 1) res is a data.frame with locus-level details
  if (is.data.frame(res)) {
    # Try to coerce column names
    nm <- names(res)
    # Common variations:
    #  - "Allele1"/"Allele2" vs "DB_Allele1"/"DB_Allele2"
    #  - "sampleid" vs "SampleID"
    map <- list(
      SampleID   = c("SampleID", "sampleid", "ID", "id"),
      Locus      = c("Locus", "locus", "marker", "Marker"),
      DB_Allele1 = c("DB_Allele1", "Allele1", "allele1", "A1"),
      DB_Allele2 = c("DB_Allele2", "Allele2", "allele2", "A2"),
      Score      = c("Score", "score", "points")
    )
    # Rename in a copy
    df <- res
    for (tgt in names(map)) {
      if (!(tgt %in% names(df))) {
        cand <- intersect(map[[tgt]], names(df))
        if (length(cand)) names(df)[names(df) == cand[1]] <- tgt
      }
    }
    detail <- .normalize_detail_schema(df)
    summary <- .build_scores_from_detail(detail)
    return(list(detail = detail, summary = summary, used = "run_match_fast:detail_df"))
  }
  
  # 2) res is a list-like
  if (is.list(res)) {
    # Try common keys
    detail <- NULL
    summary <- NULL
    # A: detail-like under "detail" / "log" / "match_log"
    for (k in c("detail", "log", "match_log", "locus_log")) {
      if (is.data.frame(res[[k]])) {
        detail <- .normalize_detail_schema(res[[k]])
        break
      }
    }
    # B: summary-like under "scores" / "summary" / "score_df"
    for (k in c("scores", "summary", "score_df", "scores_df")) {
      if (is.data.frame(res[[k]])) {
        summary <- res[[k]]
        if (!all(c("SampleID", "Score") %in% names(summary))) {
          # Try to coerce
          nm <- names(summary)
          if ("score" %in% nm) names(summary)[nm == "score"] <- "Score"
          if ("sampleid" %in% nm) names(summary)[nm == "sampleid"] <- "SampleID"
          summary <- summary[, intersect(c("SampleID", "Score"), names(summary)), drop = FALSE]
        }
      }
    }
    # If only summary exists, synthesize a minimal detail (one row per sample)
    if (is.null(detail) && !is.null(summary) && nrow(summary)) {
      loci_one <- unique(query_df$Locus)
      if (!length(loci_one)) loci_one <- "NA"
      detail <- do.call(rbind, lapply(seq_len(nrow(summary)), function(i) {
        data.frame(
          SampleID   = summary$SampleID[i],
          Locus      = loci_one[1],
          DB_Allele1 = "NA",
          DB_Allele2 = "NA",
          Score      = 0L,
          stringsAsFactors = FALSE
        )
      }))
    }
    # If still no detail, fallback mock
    if (is.null(detail)) {
      detail <- .build_mock_detail(sample_ids = sprintf("S%03d", 1:12), loci = unique(query_df$Locus))
    }
    # Ensure summary
    if (is.null(summary)) summary <- .build_scores_from_detail(detail)
    return(list(detail = detail, summary = summary, used = "run_match_fast:list_adapt"))
  }
  
  # Unknown shape -> fallback
  list(
    detail  = .build_mock_detail(sample_ids = sprintf("S%03d", 1:12), loci = unique(query_df$Locus)),
    summary = NULL, used = "mock:unknown_shape"
  )
}

# ---- Main module ----
server_match_logic <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Render result table whenever rv$trigger_run_match changes
    shiny::observeEvent(rv$trigger_run_match, {
      # Guards
      qdf <- rv$query_profile_std
      if (is.null(qdf) || !nrow(qdf)) {
        rv$status_msg <- "Query profile is empty."
        shiny::showNotification("Query profile is empty.", type = "error")
        return(invisible(NULL))
      }
      # Database prepared df must be provided in rv$db_profile (Settings will populate this).
      # For now, if not provided, try to read a default CSV (safe fallback), else use mock.
      dbdf <- rv$db_profile
      if (is.null(dbdf)) {
        if (file.exists("data/database_profile.csv")) {
          dbdf <- tryCatch(utils::read.csv("data/database_profile.csv", stringsAsFactors = FALSE), error = function(e) NULL)
        }
      }
      
      # Decide path: real matcher vs mock
      use_mock <- is.null(dbdf) || !nrow(dbdf)
      if (use_mock) {
        # Mock path (no DB)
        sample_ids <- sprintf("S%03d", 1:12)
        loci <- if ("Locus" %in% names(qdf)) unique(qdf$Locus) else character(0)
        detail <- .build_mock_detail(sample_ids, loci)
        summary <- .build_scores_from_detail(detail)
        used <- "mock:no_db"
      } else {
        # Try run_match_fast and adapt output
        # Expect locus_order is available in the global env (as in query_gui_app.R)
        locs <- if (exists("locus_order", inherits = TRUE)) get("locus_order", inherits = TRUE) else unique(qdf$Locus)
        adap <- tryCatch(.run_match_and_adapt(query_df = qdf, db_df = dbdf, locus_order = locs),
                         error = function(e) list(detail = NULL, summary = NULL, used = paste("mock:error:", conditionMessage(e))))
        detail  <- .normalize_detail_schema(adap$detail)
        summary <- if (is.null(adap$summary)) .build_scores_from_detail(detail) else adap$summary
        used    <- adap$used %||% "run_match_fast"
      }
      
      # Apply display limit (filters stored in rv by Confirm)
      ftype <- rv$filter_type %||% "top_n"
      if (!nrow(summary)) {
        filtered_summary <- summary
        filtered_detail  <- detail[0, , drop = FALSE]
      } else if (ftype == "top_n") {
        n <- max(1L, as.integer(rv$top_n %||% 10L))
        filtered_summary <- head(summary[order(-summary$Score, summary$SampleID), , drop = FALSE], n = n)
        keep <- unique(filtered_summary$SampleID)
        filtered_detail  <- detail[detail$SampleID %in% keep, , drop = FALSE]
      } else if (ftype == "score_min") {
        thr <- as.integer(rv$score_min %||% 0L)
        filtered_summary <- summary[summary$Score >= thr, , drop = FALSE]
        keep <- unique(filtered_summary$SampleID)
        filtered_detail  <- detail[detail$SampleID %in% keep, , drop = FALSE]
      } else {
        filtered_summary <- summary
        filtered_detail  <- detail
      }
      
      # Store into rv
      rv$match_detail <- filtered_detail
      rv$match_scores <- filtered_summary
      
      # Render table (detail)
      output$tbl_result <- DT::renderDT({
        if (!nrow(filtered_detail)) return(filtered_detail)
        DT::datatable(
          filtered_detail,
          rownames = FALSE,
          options = list(pageLength = 25, scrollX = TRUE)
        )
      })
      
      # Status
      rv$status_msg <- paste("Match executed via", used)
      shiny::showNotification(rv$status_msg, type = "message")
    })
    
    # Downloads
    output$dl_scores <- shiny::downloadHandler(
      filename = function() paste0("match_scores_", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv"),
      content  = function(file) {
        df <- rv$match_scores
        if (is.null(df)) df <- data.frame(SampleID = character(), Score = integer())
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$dl_detail <- shiny::downloadHandler(
      filename = function() paste0("match_log_", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv"),
      content  = function(file) {
        df <- rv$match_detail
        keep <- c("SampleID", "Locus", "DB_Allele1", "DB_Allele2", "Score")
        if (!is.null(df) && all(keep %in% names(df))) {
          df <- df[, keep, drop = FALSE]
        } else {
          df <- data.frame(matrix(ncol = length(keep), nrow = 0))
          names(df) <- keep
        }
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
  })
}

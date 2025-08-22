# scripts/gui/server_settings_logic.R
# Settings tab server: read DB (RDS/CSV), normalize, prepare, store to rv.
# No multibyte characters in code/comments.

# ---- helper ----
.assert_db_index <- function(x) {
  is_list <- is.list(x)
  has_names <- all(c("sample_ids","locus_ids","A1","A2") %in% names(x))
  ok_A1 <- is.matrix(x$A1) && is.integer(x$A1)
  ok_A2 <- is.matrix(x$A2) && is.integer(x$A2)
  ok_sid <- is.character(x$sample_ids)
  ok_lid <- is.character(x$locus_ids)
  ok_dim <- if (ok_A1 && ok_A2) {
    identical(dim(x$A1), c(length(x$sample_ids), length(x$locus_ids))) &&
      identical(dim(x$A2), c(length(x$sample_ids), length(x$locus_ids)))
  } else FALSE
  is_list && has_names && ok_A1 && ok_A2 && ok_sid && ok_lid && ok_dim
}

`%||%` <- function(x, y) if (is.null(x)) y else x

server_settings_logic <- function(id, rv, locus_order,
                                  prepare_database_df_fn = NULL,
                                  csv_normalizer_fn = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Spinner placeholder
    output$load_spinner <- shiny::renderUI({ NULL })
    
    # ---- Current database status (db_index 優先) ----
    output$db_file_name <- shiny::renderText({
      rv$db_file_name %||% "NA"
    })
    output$db_sample_count <- shiny::renderText({
      n <- NA_integer_
      if (!is.null(rv$db_index) && .assert_db_index(rv$db_index)) {
        n <- length(rv$db_index$sample_ids)
      } else {
        n <- tryCatch(length(unique(rv$db_profile$SampleID)), error = function(e) 0L)
      }
      as.integer(n)
    })
    output$db_loaded_at <- shiny::renderText({
      rv$db_loaded_at %||% "NA"
    })
    
    # ---- Load button ----
    shiny::observeEvent(input$btn_load_db, {
      file <- input$db_file
      if (is.null(file) || is.null(file$datapath) || !nzchar(file$datapath)) {
        shiny::showNotification("No file selected.", type = "error")
        return(invisible(NULL))
      }
      
      fpath <- file$datapath
      ftyp  <- input$db_file_type %||% "csv"   # ui側: choices = c("RDS"="rds","CSV"="csv")
      
      # RDS ルート: db_index を期待
      if (identical(ftyp, "rds")) {
        obj <- NULL
        ok  <- FALSE
        msg <- NULL
        try({
          obj <- readRDS(fpath)
          if (.assert_db_index(obj)) {
            ok <- TRUE
          } else {
            msg <- "RDS is not a valid db_index (list with sample_ids/locus_ids/A1/A2 integer matrices)."
          }
        }, silent = TRUE)
        
        if (!ok) {
          shiny::showNotification(msg %||% "Failed to read RDS.", type = "error")
          return(invisible(NULL))
        }
        
        # store: index 優先。profile はクリア。
        rv$db_index    <- obj
        rv$db_profile  <- NULL
        rv$db_file_name <- file$name %||% basename(fpath)
        rv$db_loaded_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        
        shiny::showNotification(sprintf("Loaded RDS index: S=%d, L=%d",
                                        length(obj$sample_ids), length(obj$locus_ids)),
                                type = "message")
        return(invisible(NULL))
      }
      
      # ---- CSV ルート（現状維持） ----
      raw_df <- NULL
      ok <- FALSE
      msg <- NULL
      try({
        raw_df <- utils::read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE)
        ok <- is.data.frame(raw_df) && nrow(raw_df) > 0L
      }, silent = TRUE)
      if (!ok) {
        shiny::showNotification("Failed to read CSV.", type = "error")
        return(invisible(NULL))
      }
      
      # optional normalizer for CSV（列名の標準化など）
      if (is.function(csv_normalizer_fn)) {
        raw_df <- tryCatch(csv_normalizer_fn(raw_df), error = function(e) raw_df)
      }
      
      # optional homo->any flag
      apply_h2a <- isTRUE(input$apply_homo_to_any)
      
      # default prepare if not provided
      if (!is.function(prepare_database_df_fn)) {
        prepare_database_df_fn <- function(d, locus_order, homo_to_any = FALSE) {
          n <- names(d)
          names(d) <- sub("(?i)^sample[_ ]?id$", "SampleID", n, perl = TRUE)
          names(d) <- sub("(?i)^locus$", "Locus", names(d), perl = TRUE)
          names(d) <- sub("(?i)^allele1$", "allele1", names(d), perl = TRUE)
          names(d) <- sub("(?i)^allele2$", "allele2", names(d), perl = TRUE)
          
          # keep long-form as-is; matcher_fast 側で map/int 化する
          d$Locus   <- as.character(d$Locus)
          d$allele1 <- as.character(d$allele1)
          d$allele2 <- as.character(d$allele2)
          
          if (isTRUE(homo_to_any)) {
            idx <- !is.na(d$allele1) & !is.na(d$allele2) & nzchar(d$allele1) & nzchar(d$allele2) & d$allele1 == d$allele2
            d$allele2[idx] <- "any"
          }
          d
        }
      }
      
      prep <- NULL
      ok <- FALSE
      msg <- NULL
      try({
        prep <- prepare_database_df_fn(raw_df, locus_order, homo_to_any = apply_h2a)
        ok <- is.data.frame(prep) && all(c("SampleID","Locus","allele1","allele2") %in% names(prep))
      }, silent = TRUE)
      
      if (!ok) {
        shiny::showNotification("Failed to prepare database frame.", type = "error")
        return(invisible(NULL))
      }
      
      # store: profile 優先（従来どおり）。index はクリア。
      rv$db_profile   <- prep
      rv$db_index     <- NULL
      rv$db_file_name <- file$name %||% basename(fpath)
      rv$db_loaded_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      
      shiny::showNotification("Database loaded.", type = "message")
    })
  })
}

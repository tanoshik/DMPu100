# scripts/gui/server_settings_logic.R
# Settings tab server: read DB (RDS/CSV), normalize, prepare, store to rv$db_profile.
# No multibyte characters in code/comments.

server_settings_logic <- function(id, rv, locus_order,
                                  prepare_database_df_fn = NULL,
                                  csv_normalizer_fn = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Spinner placeholder
    output$load_spinner <- shiny::renderUI({ NULL })
    
    # Status outputs
    output$db_file_name <- shiny::renderText({
      rv$db_file_name %||% "NA"
    })
    output$db_sample_count <- shiny::renderText({
      n <- tryCatch(length(unique(rv$db_profile$SampleID)), error = function(e) 0L)
      as.integer(n)
    })
    output$db_loaded_at <- shiny::renderText({
      rv$db_loaded_at %||% "NA"
    })
    
    # Load button
    shiny::observeEvent(input$btn_load_db, {
      file <- input$db_file
      if (is.null(file) || is.null(file$datapath) || !nzchar(file$datapath)) {
        shiny::showNotification("No file selected.", type = "error")
        return(invisible(NULL))
      }
      
      fpath <- file$datapath
      ftyp  <- input$db_file_type %||% "rds"
      apply_h2a <- isTRUE(input$apply_homo_to_any)
      
      # Read raw
      raw_df <- NULL
      if (ftyp == "rds") {
        raw_df <- tryCatch(readRDS(fpath), error = function(e) NULL)
      } else {
        raw_df <- tryCatch(utils::read.csv(fpath, stringsAsFactors = FALSE), error = function(e) NULL)
      }
      if (is.null(raw_df)) {
        shiny::showNotification("Failed to read file.", type = "error")
        return(invisible(NULL))
      }
      
      # Normalize CSV columns if needed
      if (ftyp == "csv") {
        # default normalizer: rename common variants to canonical
        if (!is.function(csv_normalizer_fn)) {
          csv_normalizer_fn <- function(d) {
            nm <- names(d)
            nm <- sub("^sampleid$", "SampleID", nm, ignore.case = TRUE)
            nm <- sub("^locus$",    "Locus",    nm, ignore.case = TRUE)
            nm <- sub("^allele1$",  "allele1",  nm, ignore.case = TRUE)
            nm <- sub("^allele2$",  "allele2",  nm, ignore.case = TRUE)
            names(d) <- nm
            d
          }
        }
        raw_df <- csv_normalizer_fn(raw_df)
      }
      
      # Prepare DB df (standardize order, fill missing loci if your function does so)
      # If prepare_database_df_fn not provided, fall back to a minimal pass-through.
      if (!is.function(prepare_database_df_fn)) {
        prepare_database_df_fn <- function(d, locus_order, homo_to_any = FALSE) {
          # Expect columns: SampleID, Locus, allele1, allele2 (case-insensitive)
          nm <- names(d)
          nm <- sub("^Allele1$", "allele1", nm, ignore.case = FALSE)
          nm <- sub("^Allele2$", "allele2", nm, ignore.case = FALSE)
          names(d) <- nm
          # Coerce types
          d$SampleID <- as.character(d$SampleID)
          d$Locus    <- as.character(d$Locus)
          d$allele1  <- as.character(d$allele1)
          d$allele2  <- as.character(d$allele2)
          # Optional homo->any
          if (isTRUE(homo_to_any)) {
            idx <- !is.na(d$allele1) & !is.na(d$allele2) & nzchar(d$allele1) & nzchar(d$allele2) & d$allele1 == d$allele2
            d$allele2[idx] <- "any"
          }
          # Return long df; run_match_fast adapter will rename as needed
          d
        }
      }
      
      prep <- NULL
      ok <- FALSE
      msg <- NULL
      try({
        prep <- prepare_database_df_fn(raw_df, locus_order, homo_to_any = apply_h2a)
        ok <- is.data.frame(prep) && all(c("SampleID", "Locus") %in% names(prep))
      }, silent = TRUE)
      
      if (!ok) {
        shiny::showNotification("Failed to prepare database frame.", type = "error")
        return(invisible(NULL))
      }
      
      # Store into rv
      rv$db_profile   <- prep
      rv$db_file_name <- file$name %||% basename(fpath)
      rv$db_loaded_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      
      shiny::showNotification("Database loaded.", type = "message")
    })
  })
}

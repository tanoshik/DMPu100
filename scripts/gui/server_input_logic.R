# scripts/gui/server_input_logic.R
# Server logic for Input and Confirm tabs.
# Assumes utils_profile.R is sourced/loaded.
# No multibyte characters in code/comments.

server_input_logic <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # helper: read CSV -> data.frame(Locus, Allele1, Allele2)
    read_query_csv <- function(path) {
      df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE),
                     error = function(e) NULL)
      if (is.null(df)) return(NULL)
      req <- c("Locus","Allele1","Allele2")
      if (!all(req %in% names(df))) return(NULL)
      df[, req]
    }
    
    # raw preview
    output$tbl_query_raw <- shiny::renderDataTable({
      if (!is.null(rv$query_profile_raw)) rv$query_profile_raw else NULL
    }, options = list(pageLength = 25))
    
    # status box
    output$txt_status <- shiny::renderText({
      rv$status_msg %||% ""
    })
    
    # Prepare button
    shiny::observeEvent(input$btn_prepare, {
      rv$status_msg <- "Preparing query..."
      lo_path <- file.path("data","locus_order.rds")
      if (!file.exists(lo_path)) {
        rv$status_msg <- "Missing data/locus_order.rds"
        return(invisible(NULL))
      }
      lo <- readRDS(lo_path)
      
      # raw profile: from CSV or empty template
      pf_raw <- NULL
      if (!is.null(input$file_query_csv) && nzchar(input$file_query_csv$datapath)) {
        pf_raw <- read_query_csv(input$file_query_csv$datapath)
      }
      if (is.null(pf_raw)) {
        # create empty any-any template for all loci
        pf_raw <- data.frame(
          Locus = lo,
          Allele1 = rep("any", length(lo)),
          Allele2 = rep("any", length(lo)),
          stringsAsFactors = FALSE
        )
      }
      
      rv$query_profile_raw <- pf_raw
      
      # standardize
      pf_std <- prepare_profile_df(
        profile_df = pf_raw,
        locus_order = lo,
        homo_to_any = isTRUE(input$chk_homo_to_any_query),
        lenient = isTRUE(input$chk_lenient_query),
        fill_missing_any = TRUE
      )
      rv$query_profile_std  <- pf_std
      rv$query_profile_show <- format_profile_df(pf_std)
      
      # update confirm tab outputs
      output$tbl_confirm_query <- shiny::renderDataTable({
        rv$query_profile_show
      }, options = list(pageLength = 25))
      
      output$txt_locus_count <- shiny::renderText({
        sprintf("Loci in query: %d", if (!is.null(rv$query_profile_std)) nrow(rv$query_profile_std) else 0L)
      })
      
      output$txt_db_count <- shiny::renderText({
        n <- if (!is.null(rv$db_std)) nrow(rv$db_std) else 0L
        sprintf("Database rows (standardized): %d", n)
      })
      
      rv$status_msg <- "Prepared."
    })
    
    # Run Match button (delegated to match logic if available)
    shiny::observeEvent(input$btn_run_match, {
      if (is.null(rv$query_profile_std)) {
        rv$status_msg <- "Prepare query first."
        return(invisible(NULL))
      }
      # fire a flag for match logic
      rv$trigger_run_match <- Sys.time()
    })
  })
}

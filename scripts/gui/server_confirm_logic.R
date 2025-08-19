# scripts/gui/server_confirm_logic.R
# Confirm tab server module
# No multibyte characters in code/comments.

server_confirm_logic <- function(id, rv, locus_order, db_count_reactive, freq_calc_fn = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Left: prepared query table
    output$tbl_confirm_query <- shiny::renderTable({
      df <- rv$query_profile_show
      if (is.null(df) || !nrow(df)) return(NULL)
      df
    }, striped = TRUE, hover = TRUE, bordered = TRUE, digits = 0)
    
    # DB count (dummy reactive now; wire to Settings later)
    output$db_count <- shiny::renderText({
      n <- tryCatch(as.integer(db_count_reactive()), error = function(e) NA_integer_)
      paste("Database Samples : N =", ifelse(is.na(n), 0L, n))
    })
    
    # Total Frequency (no placeholder)
    output$total_frequency <- shiny::renderText({
      prof <- rv$query_profile_show
      if (is.null(prof) || !nrow(prof)) return("NA")
      if (!is.function(freq_calc_fn))  return("NA")
      val <- tryCatch(freq_calc_fn(prof), error = function(e) NA_real_)
      format_total_frequency(val, digits = 3)
    })
    
    # Spinner placeholder
    output$run_match_spinner <- shiny::renderUI({ NULL })
    
    # Initialize score_min default = full score
    shiny::observeEvent(locus_order, {
      full <- as.integer(get_locus_count(locus_order) * 2L)
      shiny::updateNumericInput(session, "score_min", value = full)
    }, ignoreInit = FALSE, once = TRUE)
    
    # Run button -> set filters into rv and trigger run
    shiny::observeEvent(input$btn_run_match, {
      rv$filter_type <- input$filter_type
      rv$top_n       <- as.integer(input$top_n %||% 10L)
      rv$score_min   <- as.integer(input$score_min %||% 0L)
      
      # trigger for Result-side processing
      rv$trigger_run_match <- Sys.time()
      
      # navigate to Result
      rv$nav_request <- "Result"
      
      shiny::showNotification("Run requested.", type = "message")
    })
  })
}

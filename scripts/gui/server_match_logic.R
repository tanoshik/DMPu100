# scripts/gui/server_match_logic.R
# Server logic for Result tab.
# Placeholder until run_match_fast() is implemented/wired.
# No multibyte characters in code/comments.

server_match_logic <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # When rv$trigger_run_match changes, build a placeholder result.
    shiny::observeEvent(rv$trigger_run_match, {
      if (is.null(rv$query_profile_std)) {
        rv$status_msg <- "No query to match."
        return(invisible(NULL))
      }
      
      # Placeholder: show standardized query as pseudo "result".
      df <- rv$query_profile_std
      res <- data.frame(
        SampleID   = rep(if (!is.null(rv$sample_name)) rv$sample_name else "Query", nrow(df)),
        Locus      = df$Locus,
        DB_Allele1 = df$Allele1,
        DB_Allele2 = df$Allele2,
        Score      = 0L,
        stringsAsFactors = FALSE
      )
      output$tbl_result <- shiny::renderDataTable({
        res
      }, options = list(pageLength = 25))
      
      rv$status_msg <- "Match executed (placeholder). Wire matcher next."
    })
  })
}

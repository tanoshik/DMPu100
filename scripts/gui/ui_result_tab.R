# scripts/gui/ui_result_tab.R
# Result tab UI: shows match results.
# No multibyte characters in code/comments.

ui_result_tab <- function(id = "result_tab") {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h4("Match Results"),
    DT::DTOutput(ns("tbl_result"))
  )
}

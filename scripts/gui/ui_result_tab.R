# scripts/gui/ui_result_tab.R
# Result tab UI: shows match results (placeholder until matcher is wired).

ui_result_tab <- function(id = "result_tab") {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h4("Match Results"),
    shiny::dataTableOutput(ns("tbl_result"))
  )
}

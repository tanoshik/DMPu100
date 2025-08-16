# scripts/gui/ui_confirm_tab.R
# Confirm tab UI: shows standardized query profile and basic counts.
# No multibyte characters in code/comments.

ui_confirm_tab <- function(id = "confirm_tab") {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h4("Prepared Query Profile"),
    DT::DTOutput(ns("tbl_confirm_query")),
    shiny::hr(),
    shiny::fluidRow(
      shiny::column(6, shiny::verbatimTextOutput(ns("txt_locus_count"))),
      shiny::column(6, shiny::verbatimTextOutput(ns("txt_db_count")))
    )
  )
}

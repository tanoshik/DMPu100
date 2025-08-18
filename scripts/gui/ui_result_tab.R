# scripts/gui/ui_result_tab.R
# Result tab UI: shows match results.
# No multibyte characters in code/comments.

ui_result_tab <- function(id = "result_tab") {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h4("Match Results"),
    # Download buttons (added)
    shiny::div(
      style = "display: flex; gap: 10px; justify-content: flex-end; margin-bottom: 8px;",
      shiny::downloadButton(ns("dl_scores"), label = "Summary (match_scores.csv)"),
      shiny::downloadButton(ns("dl_detail"), label = "Detail (match_log.csv)")
    ),
    DT::DTOutput(ns("tbl_result"))
  )
}

render_result_tab <- function() {
  tabPanel(
    title = "Result",
    fluidRow(
      column(
        width = 8,
        h3("Match Results"),
        div(
          style = "display: flex; justify-content: flex-end; gap: 10px;",
          downloadButton("download_result", label = "Summary"),
          downloadButton("download_detail", label = "Detail")
        )
      )
    ),
    tableOutput("result_table")
  )
}

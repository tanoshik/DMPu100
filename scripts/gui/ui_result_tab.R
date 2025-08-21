# scripts/gui/ui_result_tab.R
# Result tab UI: shows match results.
# No multibyte characters in code/comments.

ui_result_tab <- function(id = "result_tab") {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::tags$style(shiny::HTML("
      .result-bar {
        display:flex; gap:10px; align-items:center; justify-content:space-between;
        margin-bottom:8px; flex-wrap:wrap;
      }
      .result-left { display:flex; gap:10px; align-items:center; }
      .result-right { display:flex; gap:10px; align-items:center; }
      .result-status {
        font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
        font-size:12px; color:#6b7280;
      }
    ")),
    
    shiny::h4("Match Results"),
    shiny::div(class = "result-bar",
               shiny::div(class = "result-left",
                          shiny::radioButtons(
                            inputId = ns("result_view"),
                            label   = NULL,
                            choices = c("Summary", "Detail"),
                            selected = "Summary",
                            inline = TRUE
                          )
               ),
               shiny::div(class = "result-right",
                          shiny::downloadButton(ns("dl_scores"), label = "Summary (match_scores.csv)"),
                          shiny::downloadButton(ns("dl_detail"), label  = "Detail (match_log.csv)"),
                          # NOTE: verbatimTextOutput に class は渡せないので div で包む
                          shiny::div(class = "result-status",
                                     shiny::verbatimTextOutput(ns("result_status"), placeholder = TRUE)
                          )
               )
    ),
    
    DT::DTOutput(ns("tbl_result"))
  )
}

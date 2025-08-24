# scripts/gui/ui_result_tab.R
# Result tab UI: Summary-only, server-side DT.
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
        font-size: 12px; color:#444; background:#f8f8f8; padding:6px 8px; border-radius:6px;
        border:1px solid #eee; white-space:pre-wrap; max-width: 520px;
      }
    ")),
    shiny::h4("Match Results"),
    shiny::div(class = "result-bar",
               shiny::div(class = "result-left",
                          shiny::span("View:"),
                          shiny::span("Summary")  # ← Summary固定（ラジオボタン撤去）
               ),
               shiny::div(class = "result-right",
                          shiny::downloadButton(ns("dl_scores"), label = "Summary (match_scores.csv)"),
                          shiny::downloadButton(ns("dl_detail"), label = "Detail (disabled)"),
                          shiny::div(class = "result-status",
                                     shiny::verbatimTextOutput(ns("result_status"), placeholder = TRUE)
                          )
               )
    ),
    DT::dataTableOutput(ns("tbl_result"))
  )
}

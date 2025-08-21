# scripts/gui/ui_confirm_tab.R
# Confirm tab UI (module). IDs are aligned with server_confirm_logic.R.
# No multibyte characters in code/comments.

ui_confirm_tab <- function(id = "confirm") {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    title = "Confirm",
    shiny::fluidPage(
      # --- layout styles (329版テイスト) ---
      shiny::tags$style(shiny::HTML("
        .kpi-box { display:flex; gap:16px; flex-wrap:wrap; margin-bottom:12px; }
        .kpi { border:1px solid #e5e7eb; border-radius:10px; padding:12px 14px; background:#fafafa; }
        .kpi .label { font-size:12px; color:#6b7280; }
        .kpi .value { font-size:18px; font-weight:600; }
        .controls-line { display:flex; gap:12px; align-items:flex-end; flex-wrap:wrap; margin-top:8px; }
        .twocol { display:grid; grid-template-columns: 1fr 1fr; gap:16px; }
        @media (max-width: 960px) { .twocol { grid-template-columns:1fr; } }
      ")),
      
      shiny::h3("Confirm"),
      
      # KPI row: DB Samples / Total Frequency
      shiny::div(class = "kpi-box",
                 shiny::div(class = "kpi",
                            shiny::span(class="label","DB Samples"),
                            shiny::div(class = "value", shiny::textOutput(ns("db_count"), inline = TRUE))
                 ),
                 shiny::div(class = "kpi",
                            shiny::span(class="label","Total Frequency"),
                            shiny::div(class = "value", shiny::textOutput(ns("total_frequency"), inline = TRUE))
                 )
      ),
      
      # Ranking controls + Run button（server側のIDに合わせる）
      shiny::div(class="controls-line",
                 shiny::selectInput(
                   inputId = ns("filter_type"),
                   label   = "Ranking condition",
                   choices = c("Top N" = "top_n", "Score >= n" = "score_min", "All" = "all"),
                   selected = "top_n", width = "180px"
                 ),
                 shiny::numericInput(ns("top_n"),     "Top N",      value = 50, min = 1, step = 1, width = "140px"),
                 shiny::numericInput(ns("score_min"), "Score >= n", value = 0,  min = 0, step = 1, width = "140px"),
                 shiny::div(style="flex-grow:1"),
                 shiny::actionButton(ns("btn_run_match"), "Run Match")
      ),
      
      shiny::hr(),
      
      # Prepared Query（左） + Frequency table（右）
      shiny::div(class = "twocol",
                 shiny::div(
                   shiny::h4("Prepared Query Profile"),
                   shiny::tableOutput(ns("tbl_confirm_query"))
                 ),
                 shiny::div(
                   shiny::h4("Query \u00D7 Frequency table"),
                   shiny::tableOutput(ns("confirm_freq_table"))
                 )
      ),
      
      # spinner placeholder（server側が renderUI で差し込む想定）
      shiny::uiOutput(ns("run_match_spinner"))
    )
  )
}

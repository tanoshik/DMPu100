# scripts/gui/ui_confirm_tab.R
# Confirm tab UI (minimal, DMP-conformant)
# No multibyte characters in code/comments.

render_confirm_tab <- function(visible_loci, locus_order) {
  tabPanel(
    title = "Confirm",
    fluidRow(
      # Left block: prepared query profile table
      column(
        width = 8,
        h3("Prepared Query Profile"),
        tableOutput("confirm_table")
      ),
      # Right block: ranking controls and indicators
      column(
        width = 4,
        div(
          style = "display: flex; flex-direction: column; gap: 16px; margin-top: 8px; min-width: 260px;",
          
          # Database sample count (dummy until Settings is wired)
          strong(textOutput("db_count")),
          
          # Total Frequency (display only; numeric shown in scientific notation)
          div(
            style = "display: flex; flex-direction: column; gap: 6px;",
            tags$label("Total Frequency"),
            verbatimTextOutput("total_frequency", placeholder = TRUE)
          ),
          
          # Display Limit controls
          div(
            style = "display: flex; flex-direction: column; gap: 6px;",
            radioButtons(
              inputId = "filter_type", label = "Display Limit",
              choices = c("Top N" = "top_n", "Score \u2265 n" = "score_min", "All" = "all"),
              selected = "top_n", inline = FALSE
            ),
            conditionalPanel(
              "input.filter_type == 'top_n'",
              numericInput("top_n", "Top N", value = 10L, min = 1L, step = 1L, width = "160px")
            ),
            conditionalPanel(
              "input.filter_type == 'score_min'",
              numericInput(
                inputId = "score_min",
                label   = "Score threshold (n)",
                value   = get_locus_count(locus_order) * 2L,  # full score default
                min     = 0L, step = 1L, width = "180px"
              )
            )
          ),
          
          # Run Match button
          div(
            style = "display:flex; gap: 12px; align-items:center;",
            actionButton("run_match", "Run Match"),
            uiOutput("run_match_spinner") # optional spinner placeholder
          )
        )
      )
    )
  )
}

# scripts/gui/ui_confirm_tab.R
# Confirm tab UI module
# No multibyte characters in code/comments.

ui_confirm_tab <- function(id = "confirm") {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::fluidRow(
      shiny::column(
        width = 8,
        shiny::h4("Prepared Query Profile"),
        shiny::tableOutput(ns("tbl_confirm_query"))
      ),
      shiny::column(
        width = 4,
        shiny::div(
          style = "display: flex; flex-direction: column; gap: 14px; min-width: 260px;",
          
          shiny::strong(shiny::textOutput(ns("db_count"))),
          
          shiny::div(
            style = "display: flex; flex-direction: column; gap: 6px;",
            shiny::tags$label("Total Frequency"),
            shiny::verbatimTextOutput(ns("total_frequency"), placeholder = TRUE)
          ),
          
          shiny::div(
            style = "display: flex; flex-direction: column; gap: 6px;",
            shiny::radioButtons(
              inputId = ns("filter_type"),
              label   = "Display Limit",
              choices = c("Top N" = "top_n", "Score \u2265 n" = "score_min", "All" = "all"),
              selected = "top_n",
              inline = FALSE
            ),
            shiny::conditionalPanel(
              sprintf("input['%s'] == 'top_n'", ns("filter_type")),
              shiny::numericInput(ns("top_n"), "Top N", value = 10L, min = 1L, step = 1L, width = "160px")
            ),
            shiny::conditionalPanel(
              sprintf("input['%s'] == 'score_min'", ns("filter_type")),
              shiny::numericInput(ns("score_min"), "Score threshold (n)", value = 0L, min = 0L, step = 1L, width = "180px")
            )
          ),
          
          shiny::div(
            style = "display:flex; gap: 12px; align-items:center;",
            shiny::actionButton(ns("btn_run_match"), "Run Match"),
            shiny::uiOutput(ns("run_match_spinner"))
          )
        )
      )
    )
  )
}

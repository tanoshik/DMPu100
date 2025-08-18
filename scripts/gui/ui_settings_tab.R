# scripts/gui/ui_settings_tab.R
# Settings tab UI: load DB (CSV/RDS) and show status.
# No multibyte characters in code/comments.

ui_settings_tab <- function(id = "settings") {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h4("Settings"),
    shiny::fluidRow(
      shiny::column(
        width = 6,
        shiny::wellPanel(
          shiny::strong("Load database file"),
          shiny::br(),
          shiny::radioButtons(
            inputId = ns("db_file_type"),
            label   = "File type",
            choices = c("RDS" = "rds", "CSV" = "csv"),
            selected = "rds",
            inline = TRUE
          ),
          shiny::fileInput(
            inputId = ns("db_file"),
            label   = "Choose .rds or .csv file",
            accept  = c(".rds", ".csv"),
            buttonLabel = "Browse...",
            placeholder = "No file selected"
          ),
          shiny::checkboxInput(
            inputId = ns("apply_homo_to_any"),
            label   = "Apply homo-to-any conversion (optional)",
            value   = FALSE
          ),
          shiny::div(
            style = "display:flex; gap:12px; align-items:center; margin-top:6px;",
            shiny::actionButton(ns("btn_load_db"), "Load DB"),
            shiny::uiOutput(ns("load_spinner"))
          )
        )
      ),
      shiny::column(
        width = 6,
        shiny::wellPanel(
          shiny::strong("Current database status"),
          shiny::br(),
          shiny::div(
            shiny::tags$label("Current DB File:"),
            shiny::verbatimTextOutput(ns("db_file_name"), placeholder = TRUE)
          ),
          shiny::div(
            shiny::tags$label("Sample Count:"),
            shiny::verbatimTextOutput(ns("db_sample_count"), placeholder = TRUE)
          ),
          shiny::div(
            shiny::tags$label("Last Loaded At:"),
            shiny::verbatimTextOutput(ns("db_loaded_at"), placeholder = TRUE)
          ),
          shiny::div(
            shiny::tags$label("Notes:"),
            shiny::helpText("RDS is preferred. CSV requires column normalization.")
          )
        )
      )
    )
  )
}

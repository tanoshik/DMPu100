# scripts/gui/ui_input_tab.R
# Minimal DMP-aligned Input tab + CSV multi-sample pager + validation gating:
# - Top: title + allele grid (left/right columns)
# - Bottom row: left (homo->any), right (Sample Name)
# - Actions: right-aligned "Prepare Query" (only when invalid) / "Go to Confirm" (only when valid)
# - Pager row: Query CSV + Prev / indicator / Next
# - Uses shinyjs to toggle button visibility and red highlight on invalid selectors.
# No multibyte characters in code/comments.

ui_input_tab <- function(id = "input_tab") {
  ns <- shiny::NS(id)
  
  css <- shiny::tags$style(HTML("
    .dmp-wrap { padding: 8px; }
    .dmp-title { margin: 4px 0 8px 0; }
    .dmp-grid { margin-bottom: 10px; }
    .dmp-col { padding-left: 6px; padding-right: 6px; }
    .dmp-card { border: 1px solid #e5e7eb; border-radius: 8px; padding: 8px; }
    .dmp-locus-row { display: flex; align-items: center; gap: 8px; margin-bottom: 6px; }
    .dmp-locus-label { width: 120px; font-weight: 600; white-space: nowrap; }
    .dmp-a1, .dmp-a2 { width: 100px; }
    .selectize-control { margin: 0; }
    .dmp-bottom { margin-top: 8px; }
    .dmp-actions { margin-top: 6px; display: flex; gap: 8px; justify-content: flex-end; }
    .dmp-pager { margin-top: 6px; }
    .pager-row { display: flex; align-items: center; gap: 8px; }
    .pager-indicator { min-width: 60px; text-align: center; }
    .invalid-box .selectize-control { border: 2px solid #ef4444 !important; border-radius: 6px; }
    .invalid-row { background-color: #fff1f2; border-radius: 6px; }
    .invalid-select { border: 2px solid #ef4444 !important; border-radius: 6px; }
    .val-msgs { color: #b91c1c; margin-top: 6px; white-space: pre-line; font-weight: 600; }
    /* invalid highlight */
    .invalid .selectize-control { border: 2px solid #ef4444 !important; border-radius: 6px; }
  "))
  
  shiny::fluidPage(
    shinyjs::useShinyjs(),
    css,
    shiny::div(class = "dmp-wrap",
               
               # title + allele grid (top)
               shiny::h4(class = "dmp-title", "Allele selectors by locus"),
               shiny::div(class = "dmp-grid",
                          shiny::fluidRow(
                            shiny::column(
                              width = 6, class = "dmp-col",
                              shiny::div(class = "dmp-card", shiny::uiOutput(ns("allele_grid_left")))
                            ),
                            shiny::column(
                              width = 6, class = "dmp-col",
                              shiny::div(class = "dmp-card", shiny::uiOutput(ns("allele_grid_right")))
                            )
                          )
               ),
               
               # bottom controls (1st row)
               shiny::div(class = "dmp-bottom",
                          shiny::fluidRow(
                            shiny::column(
                              width = 6,
                              shiny::checkboxInput(ns("chk_homo_to_any_query"),
                                                   "Convert homozygous to any (Query)", value = FALSE)
                            ),
                            shiny::column(
                              width = 6,
                              shiny::textInput(ns("txt_sample_name"), "Sample Name", value = "")
                            )
                          ),
                          # bottom actions (2nd row, right-aligned)
                          shiny::div(class = "dmp-actions",
                                     shiny::actionButton(ns("btn_prepare"),      "Prepare Query"),
                                     shiny::actionButton(ns("btn_goto_confirm"), "Go to Confirm")
                          ),
                          shiny::div(class = "val-msgs", shiny::htmlOutput(ns("val_msgs"))),  # <-- add this
                          # pager row (3rd row): CSV + Prev / indicator / Next
                          shiny::div(class = "dmp-pager",
                                     shiny::div(class = "pager-row",
                                                shiny::fileInput(ns("file_query_csv"),
                                                                 "Query CSV (columns: [Sample|SampleID], Locus, Allele1, Allele2)",
                                                                 accept = c(".csv", "text/csv")
                                                ),
                                                shiny::actionButton(ns("btn_q_prev"), "Prev"),
                                                shiny::div(class = "pager-indicator", shiny::textOutput(ns("txt_q_index"), inline = TRUE)),
                                                shiny::actionButton(ns("btn_q_next"), "Next")
                                     )
                          )
               )
    )
  )
}

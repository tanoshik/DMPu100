# scripts/gui/ui_input_tab.R
# DMP-aligned Input tab: allele grid (left/right), Sample Name + CSV under right grid,
# bottom: center pager, left=Prepare(+messages), right=Go(+checkbox).
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
    .dmp-pager { margin-top: 6px; display: flex; justify-content: center; }
    .pager-row { display: flex; align-items: center; gap: 8px; }
    .pager-indicator { min-width: 60px; text-align: center; }

    .dmp-bottom-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; align-items: start; }
    .dmp-bottom-left  { display: flex; flex-direction: column; align-items: center; gap: 6px; }
    .dmp-bottom-right { display: flex; flex-direction: column; align-items: flex-end; gap: 6px; }

    .val-msgs { color: #b91c1c; white-space: pre-line; font-weight: 600; }
  "))
  
  shiny::fluidPage(
    shinyjs::useShinyjs(),
    css,
    shiny::div(class = "dmp-wrap",
               
               # title + allele grid (top)
               shiny::h4(class = "dmp-title", "Allele selectors by locus"),
               shiny::div(class = "dmp-grid",
                          shiny::fluidRow(
                            # left column: left loci
                            shiny::column(
                              width = 6, class = "dmp-col",
                              shiny::div(class = "dmp-card",
                                         shiny::uiOutput(ns("allele_grid_left"))
                              )
                            ),
                            # right column: right loci + Sample Name + CSV
                            shiny::column(
                              width = 6, class = "dmp-col",
                              shiny::div(class = "dmp-card",
                                         shiny::uiOutput(ns("allele_grid_right")),
                                         shiny::div(style = "margin-top:8px;",
                                                    shiny::textInput(ns("txt_sample_name"), "Sample Name", value = "")
                                         ),
                                         shiny::div(style = "margin-top:6px;",
                                                    shiny::fileInput(ns("file_query_csv"),
                                                                     "Query CSV (columns: [Sample|SampleID], Locus, Allele1, Allele2)",
                                                                     accept = c(".csv", "text/csv")
                                                    )
                                         )
                              )
                            )
                          )
               ),
               
               # bottom area (3 rows)
               shiny::div(class = "dmp-bottom",
                          
                          # row-1: pager in the center
                          shiny::div(class = "dmp-pager",
                                     shiny::div(class = "pager-row",
                                                shiny::actionButton(ns("btn_q_prev"), "Prev"),
                                                shiny::div(class = "pager-indicator", shiny::textOutput(ns("txt_q_index"), inline = TRUE)),
                                                shiny::actionButton(ns("btn_q_next"), "Next")
                                     )
                          ),
                          
                          # row-2+3: two-column grid (left=centered prepare+messages, right=go+checkbox)
                          shiny::div(class = "dmp-bottom-grid",
                                     
                                     # left column
                                     shiny::div(class = "dmp-bottom-left",
                                                shiny::actionButton(ns("btn_prepare"), "Prepare Query"),
                                                shiny::div(class = "val-msgs", shiny::htmlOutput(ns("val_msgs")))
                                     ),
                                     
                                     # right column
                                     shiny::div(class = "dmp-bottom-right",
                                                shiny::actionButton(ns("btn_goto_confirm"), "Go to Confirm"),
                                                shiny::checkboxInput(ns("chk_homo_to_any_query"),
                                                                     "Convert homozygous to any (Query)", value = FALSE)
                                     )
                          )
               )
    )
  )
}

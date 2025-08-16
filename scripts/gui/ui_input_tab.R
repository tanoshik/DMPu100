# scripts/gui/ui_input_tab.R
# Input tab UI: CSV upload or manual editing in future, flags, and actions.
# No multibyte characters in code/comments.

ui_input_tab <- function(id = "input_tab") {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::fluidRow(
      shiny::column(
        width = 6,
        shiny::fileInput(ns("file_query_csv"), "Query CSV (columns: Locus, Allele1, Allele2)",
                         accept = c(".csv", "text/csv")),
        shiny::checkboxInput(ns("chk_homo_to_any_query"),
                             "Convert homozygous to any (Query)", value = FALSE),
        shiny::checkboxInput(ns("chk_lenient_query"),
                             "Lenient: map unknown tokens to any", value = FALSE),
        shiny::textInput(ns("txt_sample_name"), "Sample Name", value = "")
      ),
      shiny::column(
        width = 6,
        shiny::actionButton(ns("btn_prepare"), "Prepare Query"),
        shiny::actionButton(ns("btn_run_match"), "Run Match"),
        shiny::verbatimTextOutput(ns("txt_status"))
      )
    ),
    shiny::hr(),
    shiny::h4("Query (raw preview)"),
    shiny::dataTableOutput(ns("tbl_query_raw"))
  )
}

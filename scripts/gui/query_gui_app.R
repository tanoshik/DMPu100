# scripts/gui/query_gui_app.R
# Minimal 3-tab Shiny app that uses utils_profile.R.
# No multibyte characters in code/comments.
options(shiny.legacy.datatable = FALSE)
library(shiny)

# load required packages
if (!requireNamespace("shiny", quietly = TRUE)) stop("Package 'shiny' is required")

# source utils
source(file.path("scripts","utils_profile.R"), local = TRUE)

# source modules
source(file.path("scripts","gui","ui_input_tab.R"),   local = TRUE)
source(file.path("scripts","gui","ui_confirm_tab.R"), local = TRUE)
source(file.path("scripts","gui","ui_result_tab.R"),  local = TRUE)
source(file.path("scripts","gui","server_input_logic.R"),  local = TRUE)
source(file.path("scripts","gui","server_match_logic.R"),  local = TRUE)

# UI
ui <- shiny::navbarPage(
  title = "DMPu100",
  id = "nav",
  tabPanel("Input",   ui_input_tab("input")),
  tabPanel("Confirm", ui_confirm_tab("confirm")),
  tabPanel("Result",  ui_result_tab("result"))
)

# Server
server <- function(input, output, session) {
  # Close R app when the browser window is closed (single-user dev).
  session$onSessionEnded(function() {
    stopApp()
  })
  # shared state
  rv <- shiny::reactiveValues(
    query_profile_raw  = NULL,
    query_profile_std  = NULL,
    query_profile_show = NULL,
    db_std             = NULL,   # to be set by settings flow in future
    sample_name        = NULL,
    status_msg         = NULL,
    trigger_run_match  = NULL,
    nav_request        = NULL 
  )
  
  # navigate when requested by modules
  shiny::observeEvent(rv$nav_request, {
    req(rv$nav_request)
    shiny::updateNavbarPage(session, inputId = "nav", selected = rv$nav_request)
    rv$nav_request <- NULL
  })
  
  # wire modules (unchanged)
  server_input_logic("input", rv = rv)
  server_match_logic("result", rv = rv)
}

# Build app
app <- shiny::shinyApp(ui = ui, server = server)

# Launch in "app mode" (no tabs, no address bar) on Windows using Chrome/Edge.
launch_app_mode <- function(url) {
  # candidate executables (Chrome / Edge)
  cands <- c(
    "C:/Program Files/Google/Chrome/Application/chrome.exe",
    "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe",
    "C:/Program Files/Microsoft/Edge/Application/msedge.exe",
    "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
  )
  exe <- cands[file.exists(cands)][1]
  
  if (!is.na(exe) && nzchar(exe)) {
    # --app=<url> opens a minimal window (no browser chrome)
    args <- sprintf("--app=%s", shQuote(url))
    # do not block R; ignore output
    try(suppressWarnings(system2(exe, args = args, wait = FALSE)), silent = TRUE)
  } else {
    # fallback: normal browser if Chrome/Edge not found
    utils::browseURL(url)
  }
}

# Optional: quiet down DT deprecation messages globally (recommended)
options(shiny.legacy.datatable = FALSE)

# Run the app using the custom launcher
shiny::runApp(app, launch.browser = launch_app_mode)



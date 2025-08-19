# scripts/gui/query_gui_app.R
# Minimal 3-tab Shiny app wired with Confirm/Result modules (mock).
# No multibyte characters in code/comments.

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
})

# --- sources (keep as in repo; add new modules) ---
source("scripts/utils_profile.R", local = TRUE)

# Input tab (existing)
source("scripts/gui/ui_input_tab.R", local = TRUE)
source("scripts/gui/server_input_logic.R", local = TRUE)

# Confirm tab (newly added in this iteration)
source("scripts/gui/ui_confirm_tab.R", local = TRUE)
source("scripts/gui/server_confirm_logic.R", local = TRUE)

# Result tab (extended with downloads in this iteration)
source("scripts/gui/ui_result_tab.R", local = TRUE)
source("scripts/gui/server_match_logic.R", local = TRUE)

# Setting tab
source("scripts/gui/ui_settings_tab.R", local = TRUE)
source("scripts/gui/server_settings_logic.R", local = TRUE)

# --- data bootstrap (do not change existing behavior) ---
# Expect locus_order and visible_loci prepared elsewhere in the pipeline.
# To avoid altering existing flow, only set them if missing.
if (!exists("locus_order", inherits = FALSE)) {
  if (file.exists("data/locus_order.rds")) {
    locus_order <- readRDS("data/locus_order.rds")
  } else {
    locus_order <- character(0)
  }
}
if (!exists("visible_loci", inherits = FALSE)) {
  # Fallback: use locus_order if no explicit visible_loci is defined.
  visible_loci <- locus_order
}

# --- reactives container ---
rv <- reactiveValues(
  # Input/Confirm
  query_profile_show = NULL,   # prepared query shown in Confirm table
  query_profile_std  = NULL,   # standardized query used by matcher
  
  # Result
  match_detail       = NULL,   # detail table: SampleID/Locus/DB_Allele1/DB_Allele2/Score
  match_scores       = NULL,   # summary table: SampleID/Score
  
  # Filters (Confirm)
  filter_type        = "top_n",
  top_n              = 10L,
  score_min          = 0L,
  
  # Trigger
  trigger_run_match  = NULL,
  
  # Status
  status_msg         = NULL,
  
  locus_order   = locus_order,
  visible_loci  = visible_loci
)

# --- window sizing and centering (keep existing spec) ---
# Width=1000, Height=1050, placed with top margin = screen height * (1/3)
.center_and_resize_js <- HTML(
  "
  (function() {
    function centerWindow() {
      try {
        var w = 1000, h = 1050;
        if (window.resizeTo) { window.resizeTo(w, h); }
        var left = (window.screen.width  - w) / 2;
        var top  = (window.screen.height - h) / 3; // 1/3 upper
        if (window.moveTo) { window.moveTo(Math.max(left, 0), Math.max(top, 0)); }
      } catch (e) { /* no-op */ }
    }
    if (document.readyState === 'complete') {
      centerWindow();
    } else {
      window.addEventListener('load', centerWindow, { once: true });
    }
  })();
  "
)

db_count <- reactive({
  df <- rv$db_profile
  if (is.null(df)) return(0L)
  tryCatch(length(unique(df$SampleID)), error = function(e) 0L)
})

# --- UI ---
ui <- navbarPage(
  title = "DMPu100",
  id    = "main_nav",
  
  # Input は既存が tabPanel を返している想定だが、念のため明示ラップでも安全
  tabPanel("Input",    ui_input_tab("input")),
  tabPanel("Confirm",  ui_confirm_tab("confirm")),
  tabPanel("Result",   ui_result_tab("result")),
  tabPanel("Settings", ui_settings_tab("settings")),
  
  header = tags$head(singleton(tags$script(.center_and_resize_js)))
)

# --- Server ---
server <- function(input, output, session) {
  server_input_logic(
    id           = "input",
    rv           = rv
  )
  
  # ADD: Settings wiring (provides rv$db_profile)
  server_settings_logic(
    id = "settings",
    rv = rv,
    locus_order = locus_order
    # , prepare_database_df_fn = your_prepare_database_df  # if you have a project-specific function
  )
  
  server_confirm_logic(
    id                = "confirm",
    rv                = rv,
    locus_order       = locus_order,
    db_count_reactive = db_count
  )
  
  server_match_logic(
    id = "result",
    rv = rv
  )
  
  # ADD: react to nav requests from modules (Input -> Confirm/Result/Settings)
  observeEvent(rv$nav_request, {
    dest <- rv$nav_request
    if (is.character(dest) && nzchar(dest)) {
      # "Input","Confirm","Result","Settings" は tabPanel のラベルと一致させる
      updateNavbarPage(session, "main_nav", selected = dest)
    }
  }, ignoreInit = TRUE)
}

# --- Standard shinyApp footer (keep) ---
app <- shinyApp(ui = ui, server = server)
if (interactive()) runApp(app)
app

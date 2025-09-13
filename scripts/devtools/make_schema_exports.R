# scripts/devtools/make_schema_exports.R
# ASCII-only comments. Export schema/codebook/sample from data assets.
# Families: virtual_db_u100, freq_table, locus_order, locus_layout
# Output dir: output/exports/<family>/

suppressWarnings(suppressMessages({
  library(jsonlite)
  library(utils)
}))

`%||%` <- function(a,b){ if (is.null(a) || !length(a)) b else a }

# -------- Registry loader (robust) ------------------------------------------
load_registry <- function(reg_path = NULL) {
  cands <- c(reg_path, "schema_registry.json",
             file.path("scripts","devtools","schema_registry.json"))
  cands <- cands[!is.na(cands)]
  reg <- NULL
  for (p in cands) {
    if (!is.null(p) && file.exists(p)) {
      reg <- tryCatch(jsonlite::read_json(p, simplifyVector = TRUE),
                      error=function(e) NULL)
      if (!is.null(reg)) break
    }
  }
  if (is.null(reg)) stop("schema_registry.json not found")
  
  # Allow array-style (top-level array of families)
  if (is.null(reg$families) && is.list(reg) && length(reg) > 0L) {
    reg <- list(families = reg)
  }
  reg$on_mismatch <- reg$on_mismatch %||% "stop"
  
  fams <- reg$families
  if (is.null(fams) || !length(fams)) stop("No families in registry")
  
  norm <- list()
  push <- function(x) { norm[[length(norm)+1L]] <<- x }
  
  for (i in seq_along(fams)) {
    fi <- fams[[i]]
    # Case A: character vector (one or more names)
    if (is.character(fi)) {
      for (nm in fi) {
        pats <- switch(nm,
                       "virtual_db_u100" = "^data/virtual_db_u100_S\\d+_seed\\d+\\.rds$",
                       "freq_table"      = "^data/freq_table\\.rds$",
                       "locus_order"     = "^data/locus_order\\.rds$",
                       "locus_layout"    = "^data/locus_layout\\.rds$",
                       ".*"
        )
        push(list(
          family   = nm,
          patterns = pats,
          handler  = switch(nm,
                            "virtual_db_u100" = "export_virtual_db_schema",
                            "freq_table"      = "export_freq_table_schema",
                            "locus_order"     = "export_locus_order_schema",
                            "locus_layout"    = "export_locus_layout_schema",
                            ""
          )
        ))
      }
      next
    }
    # Case B: proper list object
    if (!is.list(fi)) stop("Family entry must be object or string")
    fi$family   <- fi$family   %||% paste0("family_", i)
    fi$patterns <- fi$patterns %||% character(0)
    fi$handler  <- fi$handler  %||% switch(fi$family,
                                           "virtual_db_u100" = "export_virtual_db_schema",
                                           "freq_table"      = "export_freq_table_schema",
                                           "locus_order"     = "export_locus_order_schema",
                                           "locus_layout"    = "export_locus_layout_schema",
                                           ""
    )
    push(fi)
  }
  
  reg$families <- norm
  reg
}

match_family <- function(filepath, reg) {
  for (fam in reg$families) {
    pats <- fam$patterns %||% character(0)
    if (!length(pats)) next
    rx <- paste0("(", paste(pats, collapse = ")|("), ")")
    if (any(grepl(rx, filepath, perl = TRUE))) return(fam)
  }
  NULL
}

# -------- Helpers ------------------------------------------------------------
safe_read_rds <- function(path){
  tryCatch(readRDS(path), error=function(e){
    stop("Failed to read RDS: ", path, " :: ", conditionMessage(e))
  })
}
safe_mkdir <- function(d){ if (!dir.exists(d)) dir.create(d, recursive=TRUE, showWarnings=FALSE) }
write_schema_json <- function(path_json, schema){
  safe_mkdir(dirname(path_json))
  jsonlite::write_json(schema, path_json, pretty=TRUE, auto_unbox=TRUE)
}
write_codebook_md <- function(path_md, family, src_path, notes=NULL){
  safe_mkdir(dirname(path_md))
  cat("# Codebook\n\n", file=path_md)
  cat(sprintf("- family: %s\n", family), file=path_md, append=TRUE)
  cat(sprintf("- source: %s\n\n", src_path), file=path_md, append=TRUE)
  if (!is.null(notes) && length(notes)){
    cat("## Notes\n", file=path_md, append=TRUE)
    for (ln in notes) cat(paste0("- ", ln, "\n"), file=path_md, append=TRUE)
  }
}
write_sample_csv <- function(path_csv, df){
  safe_mkdir(dirname(path_csv)); utils::write.csv(df, path_csv, row.names=FALSE)
}

infer_schema_df <- function(x){
  props <- lapply(names(x), function(nm) list(type = class(x[[nm]])[1]))
  list(title="inferred schema (data.frame)", version="0.1.0",
       type="object", properties=setNames(props, names(x)), nrow=nrow(x))
}
infer_schema_vector <- function(x, item_type="string"){
  list(title="inferred schema (vector)", version="0.1.0",
       type="array", items=list(type=item_type), length=length(x))
}
infer_schema_list <- function(x){
  list(title="inferred schema (list)", version="0.1.0",
       type="object", names=names(x), length=length(x))
}
is_indexed_db <- function(obj){
  is.list(obj) && all(c("sample_ids","locus_ids","A1","A2") %in% names(obj))
}

.emit <- function(out_dir, stem, family, schema, sample_df=NULL, notes=NULL, src_path=NA_character_){
  write_schema_json(file.path(out_dir, paste0(stem, ".schema.json")), schema)
  write_codebook_md(file.path(out_dir, paste0(stem, ".codebook.md")), family, src_path, notes)
  if (!is.null(sample_df)) write_sample_csv(file.path(out_dir, paste0(stem, ".sample.csv")), sample_df)
}

# -------- Exporters ----------------------------------------------------------
export_virtual_db_schema <- function(path, family="virtual_db_u100"){
  x <- safe_read_rds(path)
  if (!is_indexed_db(x)) stop("virtual_db_u100 expects indexed_db list")
  sample_ids <- x$sample_ids; locus_ids <- x$locus_ids
  A1 <- x$A1; A2 <- x$A2
  if (!is.integer(A1)) A1 <- as.integer(A1)
  if (!is.integer(A2)) A2 <- as.integer(A2)
  s_take <- min(length(sample_ids), 4L); l_take <- min(length(locus_ids), 5L)
  rows <- vector("list", s_take * l_take); k <- 1L
  for (si in seq_len(s_take)) for (li in seq_len(l_take)) {
    rows[[k]] <- list(SampleID=sample_ids[si], Locus=as.character(locus_ids[li]),
                      A1=A1[si,li], A2=A2[si,li]); k <- k + 1L
  }
  sample_df <- do.call(rbind.data.frame, rows)
  schema <- list(
    title="Virtual DB (indexed grid)", version="0.1.0", type="object",
    properties=list(
      sample_ids=list(type="array", items=list(type="string"), length=length(sample_ids)),
      locus_ids =list(type="array", items=list(type="string"), length=length(locus_ids)),
      A1=list(type="integer", shape=c(length(sample_ids), length(locus_ids))),
      A2=list(type="integer", shape=c(length(sample_ids), length(locus_ids)))
    )
  )
  out_dir <- file.path("output","exports", family)
  .emit(out_dir, tools::file_path_sans_ext(basename(path)), family, schema, sample_df,
        notes=c("ANY_CODE=9999 treated as wildcard.",
                "Traversal order should align with locus_order."),
        src_path=normalizePath(path, winslash="/", mustWork=FALSE))
}

export_freq_table_schema <- function(path, family="freq_table"){
  x <- safe_read_rds(path)
  if (!is.data.frame(x)) stop("freq_table.rds must be a data.frame")
  schema <- infer_schema_df(x)
  sample_df <- utils::head(x, 100L)
  out_dir <- file.path("output","exports", family)
  .emit(out_dir, "freq_table", family, schema, sample_df,
        notes=c("Frequencies are used for display for now.",
                "Decimal alleles allowed as character.",
                "ANY_CODE(9999) handling: see project docs."),
        src_path=normalizePath(path, winslash="/", mustWork=FALSE))
}

export_locus_order_schema <- function(path, family="locus_order"){
  x <- safe_read_rds(path)
  if (!is.vector(x) || !is.character(x)) stop("locus_order.rds must be a character vector")
  schema <- infer_schema_vector(x, "string")
  sample_df <- data.frame(Locus = x, stringsAsFactors = FALSE)
  out_dir <- file.path("output","exports", family)
  .emit(out_dir, "locus_order", family, schema, sample_df,
        notes=c("Order used by UI and scoring traversal.",
                "Keep in sync with GlobalFiler 21 loci (fixed for now)."),
        src_path=normalizePath(path, winslash="/", mustWork=FALSE))
}

export_locus_layout_schema <- function(path, family="locus_layout"){
  x <- safe_read_rds(path)
  out_dir <- file.path("output","exports", family)
  if (is.data.frame(x)) {
    schema <- infer_schema_df(x); sample_df <- utils::head(x, 50L)
    notes <- c("UI grouping (left/right panes, dye assignment).",
               "Recommended columns: Locus, Pane, Group, Order, Dye.")
  } else if (is.list(x)) {
    schema <- infer_schema_list(x)
    sample_rows <- list(); take <- 0L
    if (length(names(x))) for (nm in names(x)) {
      v <- x[[nm]]
      if (is.vector(v) && is.character(v)) {
        for (val in utils::head(v, 50L)) {
          sample_rows[[length(sample_rows)+1L]] <- list(Pane=nm, Locus=val)
          take <- take + 1L; if (take >= 100L) break
        }
      }
      if (take >= 100L) break
    }
    sample_df <- if (length(sample_rows)) do.call(rbind.data.frame, sample_rows) else NULL
    notes <- c("List form detected (e.g., left/right -> loci).",
               "Consider normalizing to a data.frame for downstream simplicity.")
  } else {
    stop("Unsupported locus_layout class: ", paste(class(x), collapse = ","))
  }
  .emit(out_dir, "locus_layout", family, schema, sample_df, notes,
        src_path=normalizePath(path, winslash="/", mustWork=FALSE))
}

dispatch_export <- function(fam, input_path){
  family <- fam$family %||% "unknown"
  handler_name <- fam$handler %||% ""
  if (nzchar(handler_name) && exists(handler_name, mode="function")) {
    return(do.call(handler_name, list(input_path, family)))
  }
  if (identical(family, "virtual_db_u100")) return(export_virtual_db_schema(input_path, family))
  if (identical(family, "freq_table"))      return(export_freq_table_schema(input_path, family))
  if (identical(family, "locus_order"))     return(export_locus_order_schema(input_path, family))
  if (identical(family, "locus_layout"))    return(export_locus_layout_schema(input_path, family))
  obj <- safe_read_rds(input_path)
  out_dir <- file.path("output","exports", family)
  if (is.data.frame(obj)) {
    .emit(out_dir, tools::file_path_sans_ext(basename(input_path)), family,
          infer_schema_df(obj), utils::head(obj, 50L),
          notes="fallback by class: data.frame",
          src_path=normalizePath(input_path, winslash="/", mustWork=FALSE)); return(invisible(TRUE))
  }
  if (is.vector(obj) && is.character(obj)) {
    .emit(out_dir, tools::file_path_sans_ext(basename(input_path)), family,
          infer_schema_vector(obj, "string"),
          data.frame(value=head(obj, 100L)),
          notes="fallback by class: character vector",
          src_path=normalizePath(input_path, winslash="/", mustWork=FALSE)); return(invisible(TRUE))
  }
  if (is.list(obj)) {
    .emit(out_dir, tools::file_path_sans_ext(basename(input_path)), family,
          infer_schema_list(obj), NULL,
          notes="fallback by class: list",
          src_path=normalizePath(input_path, winslash="/", mustWork=FALSE)); return(invisible(TRUE))
  }
  stop("Unsupported RDS content class for: ", input_path)
}

main <- function(reg_path = NULL){
  reg <- load_registry(reg_path)
  files <- character(0)
  if (dir.exists("data")) {
    files <- list.files("data", recursive = TRUE, full.names = TRUE)
    files <- files[grepl("\\.(rds)$", files, ignore.case = TRUE)]
  }
  if (!length(files)) { message("[EXPORT] No candidate RDS under data/"); return(invisible(TRUE)) }
  for (f in files) {
    fam <- match_family(f, reg)
    if (is.null(fam)) {
      msg <- paste0("[REGISTRY] No family matched: ", f)
      if (identical(reg$on_mismatch, "stop")) stop(msg) else { warning(msg); next }
    }
    dispatch_export(fam, f)
  }
}

if (identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  reg_path <- if (length(args)) args[1] else NULL
  main(reg_path)
}

# rds_schema_export.R
# Purpose: Export schema.json, codebook.md, and sample.csv from a single .rds input.
# Notes:
# - ASCII comments only.
# - Do not infer domain logic here. Only structural scan and lightweight typing.
# - Trust order: schema_json > codebook_md > sample_csv > chat_text.

suppressPackageStartupMessages({
  library(jsonlite)
})

infer_type <- function(x) {
  cls <- class(x)[1L]
  if (inherits(x, "integer")) return("integer")
  if (inherits(x, "numeric")) return("number")
  if (inherits(x, "double"))  return("number")
  if (inherits(x, "logical")) return("boolean")
  if (inherits(x, "character")) return("string")
  return("string")
}

export_rds_schema <- function(rds_path, out_dir, family_name = NULL) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  obj <- readRDS(rds_path)
  
  schema <- list()
  codebook <- c()
  sample_csv_path <- file.path(out_dir, "sample.csv")
  
  if (is.data.frame(obj)) {
    cols <- lapply(names(obj), function(nm) {
      list(name = nm, type = infer_type(obj[[nm]]))
    })
    schema <- list(
      title = family_name %||% basename(rds_path),
      type = "object",
      properties = list(
        columns = cols
      ),
      required = list("columns"),
      additionalProperties = FALSE
    )
    
    # write sample.csv (head 10 rows)
    utils::write.csv(utils::head(obj, 10), sample_csv_path, row.names = FALSE)
    
    # codebook
    codebook <- c("# codebook", "", "| name | type |", "|---|---|")
    for (c in cols) {
      codebook <- c(codebook, sprintf("| %s | %s |", c$name, c$type))
    }
  } else if (is.vector(obj)) {
    schema <- list(
      title = family_name %||% basename(rds_path),
      type = "object",
      properties = list(
        vector = list(
          type = "array",
          items = list(type = infer_type(obj))
        )
      ),
      required = list("vector"),
      additionalProperties = FALSE
    )
    # sample.csv with first 10 entries
    utils::write.csv(data.frame(value = utils::head(obj, 10)), sample_csv_path, row.names = FALSE)
    codebook <- c("# codebook", "", "Vector sample exported as one-column CSV named 'value'.")
  } else {
    stop("Unsupported RDS content class: ", paste(class(obj), collapse = ","))
  }
  
  # schema.json
  schema_path <- file.path(out_dir, "schema.json")
  writeLines(jsonlite::prettify(jsonlite::toJSON(schema, auto_unbox = TRUE)), schema_path)
  
  # codebook.md
  writeLines(codebook, file.path(out_dir, "codebook.md"))
  
  return(invisible(list(
    schema_path = schema_path,
    codebook_path = file.path(out_dir, "codebook.md"),
    sample_path = sample_csv_path
  )))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

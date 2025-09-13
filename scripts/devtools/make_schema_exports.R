# make_schema_exports.R
# Purpose: Batch export schema artifacts for all .rds under data/.
# Output under output/exports/<family>/<basename_without_ext>/
# Also writes output/exports/manifest.json

suppressPackageStartupMessages({
  library(jsonlite)
})

source(file.path("scripts","devtools","rds_schema_export.R"))

load_registry <- function(path = "schema_registry.json") {
  fromJSON(path, simplifyVector = TRUE)
}

match_family <- function(path, registry) {
  for (fam in registry$families) {
    pats <- fam$patterns
    for (p in pats) {
      # glob-like: convert '*' to '.*' for grepl
      rx <- paste0("^", gsub("\\*", ".*", p), "$")
      if (grepl(rx, path)) return(fam$name)
    }
  }
  return(NA_character_)
}

main <- function() {
  registry <- load_registry()
  rds_files <- list.files("data", pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  
  manifest <- list(version = "1.0.0", exports = list())
  
  for (rds in rds_files) {
    fam <- match_family(rds, registry)
    if (is.na(fam)) {
      stop("No family matched for: ", rds)
    }
    base <- tools::file_path_sans_ext(basename(rds))
    out_dir <- file.path("output","exports", fam, base)
    res <- export_rds_schema(rds, out_dir, fam)
    
    manifest$exports[[length(manifest$exports)+1L]] <- list(
      family = fam,
      rds = rds,
      out_dir = out_dir,
      schema = res$schema_path,
      codebook = res$codebook_path,
      sample = res$sample_path
    )
  }
  
  man_path <- file.path("output","exports","manifest.json")
  if (!dir.exists(dirname(man_path))) dir.create(dirname(man_path), recursive = TRUE)
  writeLines(prettify(toJSON(manifest, auto_unbox = TRUE)), man_path)
  message("Wrote manifest: ", man_path)
}

if (identical(environment(), globalenv())) {
  main()
}

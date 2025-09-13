# scripts/devtools/make_schema_exports.R
# ASCII-only comments. Dev utility to export schema/codebook/sample from .rds.
# Supports:
#  - data.frame directly
#  - list(table = data.frame)
#  - list of equal-length vectors
#  - indexed_db list(sample_ids,locus_ids,A1,A2)

# ---- robust null-coalescing
`%||%` <- function(a,b){ if (is.null(a) || !length(a)) b else a }

# ---- registry loader (full replacement)
load_registry <- function(reg_path = NULL) {
  # search order: explicit path -> repo root -> scripts/devtools
  cands <- c(reg_path,
             "schema_registry.json",
             "scripts/devtools/schema_registry.json")
  cands <- unique(cands[!is.na(cands) & nzchar(cands)])
  hit <- NULL
  for (p in cands) { if (file.exists(p)) { hit <- p; break } }
  if (is.null(hit)) stop("registry json not found: schema_registry.json")
  
  # NOTE: keep raw list structure to avoid unwanted vector/data.frame simplification
  reg_raw <- jsonlite::fromJSON(hit, simplifyVector = FALSE)
  
  # allow either {"families":[...], "on_mismatch":"stop"} or just [...]
  if (is.list(reg_raw) && !is.null(reg_raw$families)) {
    fams <- reg_raw$families
    default_on_mismatch <- reg_raw$on_mismatch %||% "stop"
  } else if (is.list(reg_raw) && is.null(names(reg_raw))) {
    # top-level array
    fams <- reg_raw
    default_on_mismatch <- "stop"
  } else {
    stop("registry JSON shape not supported (expect {families:[...]} or top-level array)")
  }
  if (!length(fams)) stop("no families in registry")
  
  # normalize each family item
  fams2 <- lapply(fams, function(f) {
    # if family item was simplified to atomic, coerce to empty list to fail fast
    if (!is.list(f)) stop("family entry must be a JSON object (not an atomic array/vector)")
    list(
      family       = as.character(f$family %||% NA_character_),
      patterns     = as.character(unlist(f$patterns %||% character(0), use.names = FALSE)),
      schema       = as.character(f$schema %||% NA_character_),
      outdir       = as.character(f$outdir %||% "output/exports"),
      on_mismatch  = as.character(f$on_mismatch %||% default_on_mismatch)
    )
  })
  
  list(path = hit, families = fams2)
}

# ---- glob/regex matcher per family -----------------------------------------
match_family <- function(fam_entry, all_paths) {
  pats <- fam_entry$patterns
  if (!length(pats)) return(character(0))
  keep <- rep(FALSE, length(all_paths))
  for (rx in pats) keep <- keep | grepl(rx, all_paths, perl = TRUE)
  all_paths[keep]
}

# ---- RDS reader and normalization ------------------------------------------
read_rds_df <- function(path) {
  obj <- readRDS(path)
  # case 1: data.frame
  if (inherits(obj, "data.frame")) return(obj)
  # case 2: list(table = data.frame)
  if (is.list(obj) && !is.null(obj$table) && inherits(obj$table, "data.frame")) return(obj$table)
  # case 3: list of equal-length vectors -> data.frame
  if (is.list(obj) && is.null(obj$table) && length(obj) > 0 && all(vapply(obj, function(x) is.atomic(x) || is.factor(x), logical(1)))) {
    lens <- vapply(obj, length, integer(1))
    if (length(unique(lens)) == 1) return(as.data.frame(obj, stringsAsFactors = FALSE))
  }
  # case 4: indexed_db (sample_ids, locus_ids, A1, A2)
  if (is.list(obj) && all(c("sample_ids","locus_ids","A1","A2") %in% names(obj))) {
    # We do NOT fully materialize S*L rows to avoid heavy mem; return a tiny sample df.
    sample_ids <- obj$sample_ids
    locus_ids  <- obj$locus_ids
    A1 <- obj$A1; A2 <- obj$A2
    if (!is.integer(A1)) A1 <- as.integer(A1)
    if (!is.integer(A2)) A2 <- as.integer(A2)
    S <- length(sample_ids); L <- length(locus_ids)
    s_take <- min(S, 4L); l_take <- min(L, 5L)
    rows <- vector("list", s_take * l_take)
    k <- 1L
    for (si in seq_len(s_take)) {
      for (li in seq_len(l_take)) {
        rows[[k]] <- list(
          SampleID = sample_ids[si],
          Locus    = as.character(locus_ids[li]),
          A1       = A1[si, li],
          A2       = A2[si, li]
        )
        k <- k + 1L
      }
    }
    tiny <- do.call(rbind.data.frame, rows)
    return(tiny)
  }
  stop(sprintf("Unsupported RDS content class. class=%s keys=%s",
               paste(class(obj), collapse=","), paste(names(obj), collapse=",")))
}

# ---- schema inference for indexed_db ---------------------------------------
infer_schema_for_indexed_db <- function(obj) {
  # Schema for long-form view: SampleID(int), Locus(enum with levels), A1(int), A2(int)
  levels <- as.character(obj$locus_ids)
  list(
    title   = "virtual_db_u100 (indexed_db long view)",
    version = "0.1.0",
    columns = list(
      list(name="SampleID", type="integer"),
      list(name="Locus",    type="string", levels = levels),
      list(name="A1",       type="integer"),
      list(name="A2",       type="integer")
    )
  )
}

# ---- schema inference generic (best-effort) --------------------------------
infer_schema_from_df <- function(df) {
  cols <- lapply(names(df), function(nm){
    x <- df[[nm]]
    tp <- if (is.integer(x)) "integer" else if (is.numeric(x)) "number" else if (is.logical(x)) "boolean" else "string"
    out <- list(name = nm, type = tp)
    if (is.factor(x)) out$levels <- levels(x)
    out
  })
  list(
    title   = "inferred schema",
    version = "0.1.0",
    columns = cols
  )
}

# ---- codebook writer --------------------------------------------------------
write_codebook <- function(path_md, schema, family, src_path) {
  dir.create(dirname(path_md), recursive = TRUE, showWarnings = FALSE)
  cat("# Codebook\n\n", file = path_md)
  cat(sprintf("- family: %s\n", family), file = path_md, append = TRUE)
  cat(sprintf("- source: %s\n", src_path), file = path_md, append = TRUE)
  cat(sprintf("- title:  %s\n", schema$title %||% ""), file = path_md, append = TRUE)
  cat(sprintf("\n\n## Columns\n\n"), file = path_md, append = TRUE)
  for (c in schema$columns) {
    cat(sprintf("- %s : %s", c$name, c$type), file = path_md, append = TRUE); 
    if (!is.null(c$levels)) cat(sprintf(" (levels: %s)", paste(head(c$levels, 20L), collapse=", ")), file = path_md, append = TRUE)
    cat("\n", file = path_md, append = TRUE)
  }
}

# ---- sample csv writer (tiny) ----------------------------------------------
write_sample_csv <- function(path_csv, df) {
  dir.create(dirname(path_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(head(df, 50L), path_csv, row.names = FALSE)
}

# ---- export one -------------------------------------------------------------
export_rds_schema <- function(input_rds, out_prefix, family, schema_path = NA_character_, on_mismatch = "stop") {
  obj <- readRDS(input_rds)
  is_indexed <- is.list(obj) && all(c("sample_ids","locus_ids","A1","A2") %in% names(obj))
  # choose schema
  if (!is.na(schema_path) && nzchar(schema_path) && file.exists(schema_path)) {
    schema <- jsonlite::fromJSON(schema_path, simplifyVector = TRUE)
  } else {
    schema <- if (is_indexed) infer_schema_for_indexed_db(obj) else infer_schema_from_df(read_rds_df(input_rds))
  }
  # very light mismatch check (column names only if df)
  ok <- TRUE
  if (!is_indexed) {
    df <- read_rds_df(input_rds)
    want <- vapply(schema$columns, function(c) c$name, character(1))
    have <- names(df)
    if (!all(want %in% have)) {
      msg <- sprintf("schema columns not all present. want=%s / have=%s", paste(want, collapse=","), paste(have, collapse=","))
      if (identical(on_mismatch, "stop")) stop(msg) else warning(msg)
      ok <- FALSE
    }
  }
  # write files
  schema_json <- paste0(out_prefix, ".schema.json")
  codebook_md <- paste0(out_prefix, ".codebook.md")
  sample_csv  <- paste0(out_prefix, ".sample.csv")
  dir.create(dirname(schema_json), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(schema, schema_json, pretty = TRUE, auto_unbox = TRUE)
  write_codebook(codebook_md, schema, family, input_rds)
  # sample
  smalldf <- if (is_indexed) read_rds_df(input_rds) else utils::head(read_rds_df(input_rds), 50L)
  write_sample_csv(sample_csv, smalldf)
  invisible(list(schema = schema_json, codebook = codebook_md, sample = sample_csv, ok = ok))
}

# ---- main -------------------------------------------------------------------
main <- function() {
  reg <- load_registry(NULL)
  message("[REGISTRY] ", normalizePath(reg$path, winslash = "/"))
  all_paths <- list.files("data", pattern = "\\.rds$", recursive = FALSE, full.names = TRUE)
  for (fam in reg$families) {
    hits <- match_family(fam, all_paths)
    if (!length(hits)) next
    message(sprintf("[FAMILY] %s  files=%d", fam$family, length(hits)))
    for (f in hits) {
      base <- sub("\\.rds$", "", basename(f))
      outdir <- file.path(fam$outdir, fam$family)
      out_prefix <- file.path(outdir, base)
      message(sprintf("  - exporting %s -> %s", f, normalizePath(out_prefix, winslash = "/", mustWork = FALSE)))
      export_rds_schema(
        input_rds   = f,
        out_prefix  = out_prefix,
        family      = fam$family,
        schema_path = fam$schema,
        on_mismatch = fam$on_mismatch %||% "stop"
      )
    }
  }
}

if (identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  main()
}

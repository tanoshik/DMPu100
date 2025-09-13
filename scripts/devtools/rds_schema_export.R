# ASCII only.
# Helpers for exporting schema/codebook/sample from RDS.

suppressWarnings(suppressMessages({
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Need jsonlite")
}))

# ---- tiny utils ----
.norm <- function(x) gsub("\\\\","/", x)

# read RDS into a normalized description
# Supported:
# - data.frame
# - list(table=<data.frame>)
# - list with keys: sample_ids (int vec), locus_ids (int/chr vec), A1 (int matrix), A2 (int matrix)
read_rds_df <- function(path) {
  x <- readRDS(path)
  # direct data.frame
  if (is.data.frame(x)) {
    return(list(kind="data.frame", table=x))
  }
  # list(table=df)
  if (is.list(x) && "table" %in% names(x) && is.data.frame(x$table)) {
    return(list(kind="data.frame", table=x$table))
  }
  # indexed_db: list(sample_ids, locus_ids, A1, A2)
  if (is.list(x) && all(c("sample_ids","locus_ids","A1","A2") %in% names(x))) {
    # very light checks
    if (!is.matrix(x$A1) || !is.matrix(x$A2)) stop("A1/A2 must be matrix")
    if (nrow(x$A1) != nrow(x$A2) || ncol(x$A1) != ncol(x$A2)) stop("A1/A2 dim mismatch")
    if (length(x$sample_ids) != nrow(x$A1)) stop("sample_ids length mismatch")
    if (length(x$locus_ids)  != ncol(x$A1)) stop("locus_ids length mismatch")
    return(list(
      kind="indexed_db",
      sample_ids = x$sample_ids,
      locus_ids  = x$locus_ids,
      A1_dim     = dim(x$A1),
      A2_dim     = dim(x$A2)
    ))
  }
  # list of equal-length vectors -> data.frame
  if (is.list(x) && length(x) > 0 && all(vapply(x, function(col) is.atomic(col) || is.factor(col), logical(1)))) {
    lens <- vapply(x, length, integer(1))
    if (length(unique(lens)) == 1L) {
      return(list(kind="data.frame", table=as.data.frame(x, stringsAsFactors = FALSE)))
    }
  }
  cls <- class(x)
  stop(sprintf("Unsupported RDS content class. class=%s", paste(cls, collapse=",")))
}

# build schema object (list) from the normalized description
build_schema_obj <- function(desc, family) {
  if (desc$kind == "data.frame") {
    df <- desc$table
    cols <- lapply(names(df), function(nm){
      v <- df[[nm]]
      list(
        name = nm,
        type = class(v)[1],
        nullable = any(is.na(v))
      )
    })
    return(list(
      family = family,
      kind   = "data.frame",
      nrow   = nrow(df),
      ncol   = ncol(df),
      columns = cols
    ))
  }
  if (desc$kind == "indexed_db") {
    return(list(
      family = family,
      kind   = "indexed_db",
      sample_ids = list(type=class(desc$sample_ids)[1], length=length(desc$sample_ids)),
      locus_ids  = list(type=class(desc$locus_ids)[1],  length=length(desc$locus_ids)),
      A1         = list(type="integer_matrix", dim=desc$A1_dim),
      A2         = list(type="integer_matrix", dim=desc$A2_dim),
      notes      = "ANY_CODE=9999; locus order fixed; integers encoded"
    ))
  }
  stop("unknown desc kind")
}

# write minimal codebook.md
write_codebook_md <- function(path_prefix, schema_obj) {
  p <- paste0(path_prefix, ".codebook.md")
  out <- c(
    "# Codebook",
    paste0("- family: ", schema_obj$family),
    paste0("- kind  : ", schema_obj$kind)
  )
  if (schema_obj$kind == "data.frame") {
    out <- c(out, "", "## columns")
    for (c in schema_obj$columns) {
      out <- c(out, sprintf("- %s: type=%s nullable=%s", c$name, c$type, c$nullable))
    }
  } else if (schema_obj$kind == "indexed_db") {
    out <- c(out, "",
             sprintf("- sample_ids: type=%s length=%d", schema_obj$sample_ids$type, schema_obj$sample_ids$length),
             sprintf("- locus_ids : type=%s length=%d",  schema_obj$locus_ids$type,  schema_obj$locus_ids$length),
             sprintf("- A1        : %s dim=%s", schema_obj$A1$type, paste(schema_obj$A1$dim, collapse="x")),
             sprintf("- A2        : %s dim=%s", schema_obj$A2$type, paste(schema_obj$A2$dim, collapse="x")))
  }
  dir.create(dirname(p), recursive=TRUE, showWarnings=FALSE)
  writeLines(out, p, useBytes = TRUE)
}

# write small sample.csv (data.frame only). For indexed_db, we emit a tiny index summary CSV.
write_sample_csv <- function(path_prefix, desc) {
  p <- paste0(path_prefix, ".sample.csv")
  dir.create(dirname(p), recursive=TRUE, showWarnings=FALSE)
  if (desc$kind == "data.frame") {
    n <- min(nrow(desc$table), 50L)
    utils::write.csv(desc$table[seq_len(n), , drop=FALSE], p, row.names = FALSE)
  } else if (desc$kind == "indexed_db") {
    # emit summary rows: just first 10 sample_ids and first 10 locus_ids
    sid <- head(desc$sample_ids, 10L)
    lid <- head(desc$locus_ids,  10L)
    df  <- data.frame(kind="indexed_db",
                      sample_ids=I(list(sid)),
                      locus_ids =I(list(lid)),
                      stringsAsFactors = FALSE)
    # Write as a simple CSV: we stringify the vectors
    df2 <- data.frame(kind=df$kind,
                      sample_ids=sprintf("[%s]", paste(sid, collapse=" ")),
                      locus_ids =sprintf("[%s]", paste(lid, collapse=" ")),
                      stringsAsFactors = FALSE)
    utils::write.csv(df2, p, row.names = FALSE)
  }
}

# single file export
export_rds_schema <- function(input_rds, out_prefix, family, schema_path = NULL, on_mismatch = c("stop","report")) {
  on_mismatch <- match.arg(on_mismatch)
  desc <- read_rds_df(input_rds)
  schema_obj <- build_schema_obj(desc, family)
  
  # if schema_path given, compare minimal keys; otherwise just write detected schema
  if (!is.null(schema_path) && nchar(schema_path) > 0 && file.exists(schema_path)) {
    given <- jsonlite::fromJSON(schema_path, simplifyVector = TRUE)
    # very light comparison on "kind"
    same <- is.list(given) && isTRUE(identical(given$kind, schema_obj$kind))
    if (!same) {
      msg <- sprintf("schema mismatch: detected.kind=%s vs registry.kind=%s",
                     schema_obj$kind, if (is.list(given)) as.character(given$kind) else "NA")
      if (on_mismatch == "stop") stop(msg) else message("[SCHEMA] ", msg)
    }
  }
  
  # write schema.json
  sp <- paste0(out_prefix, ".schema.json")
  dir.create(dirname(sp), recursive=TRUE, showWarnings=FALSE)
  jsonlite::write_json(schema_obj, sp, pretty = TRUE, auto_unbox = TRUE)
  
  # write codebook and sample
  write_codebook_md(out_prefix, schema_obj)
  write_sample_csv(out_prefix, desc)
  
  invisible(list(schema=sp))
}

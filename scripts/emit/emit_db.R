# scripts/emit/emit_db.R
# Emit DB slice by SampleID with order preserved.
suppressPackageStartupMessages({ library(jsonlite) })

detect_id_col <- function(df) {
  nms <- tolower(names(df))
  cand <- which(nms %in% c("sampleid","sample_id","id","sid"))
  if (length(cand) == 0) stop("ids_path CSV must contain SampleID-like column")
  cand[1]
}

decode_allele <- function(v, any_code = 9999L) {
  v <- as.integer(v)
  out <- ifelse(v == any_code, "any", sprintf("%g", v/100))
  tolower(out)
}

emit_db <- function(db_path, ids_path, out_dir, format = c("csv","rds"), any_code = 9999L) {
  format <- match.arg(format)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ids_df <- read.csv(ids_path, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- detect_id_col(ids_df)
  target_ids <- as.character(ids_df[[id_col]])
  
  if (grepl("\\.rds$", db_path, ignore.case = TRUE)) {
    db <- readRDS(db_path)  # schema: list(sample_ids, locus_ids, A1, A2)
    stopifnot(is.list(db), !is.null(db$sample_ids), !is.null(db$locus_ids))
    ord <- match(db$sample_ids, db$sample_ids)  # identity (kept for clarity)
    keep <- which(db$sample_ids %in% target_ids)
    keep <- keep[order(keep)]  # preserve DB order
    out <- list(
      sample_ids = db$sample_ids[keep],
      locus_ids  = db$locus_ids,
      A1 = db$A1[keep, , drop=FALSE],
      A2 = db$A2[keep, , drop=FALSE]
    )
    if (format == "rds") {
      saveRDS(out, file.path(out_dir, "db_slice.rds"))
      return(invisible(TRUE))
    } else {
      # readable CSV
      res <- data.frame(
        SampleID = rep(out$sample_ids, each = length(out$locus_ids)),
        Locus = rep(out$locus_ids, times = length(out$sample_ids)),
        A1 = as.vector(out$A1),
        A2 = as.vector(out$A2),
        stringsAsFactors = FALSE
      )
      res$A1 <- decode_allele(res$A1, any_code)
      res$A2 <- decode_allele(res$A2, any_code)
      write.csv(res, file.path(out_dir, "db_slice.csv"), row.names = FALSE)
      return(invisible(TRUE))
    }
  } else {
    # CSV DB: expect columns SampleID,Locus,A1,A2 in internal integer form
    db <- read.csv(db_path, stringsAsFactors = FALSE, check.names = FALSE)
    nms <- tolower(names(db))
    stopifnot(all(c("sampleid","locus","a1","a2") %in% nms))
    ord <- match(unique(db$SampleID), unique(db$SampleID))
    keep_ids <- unique(db$SampleID)[unique(db$SampleID) %in% target_ids]
    # preserve DB order
    keep_ids <- keep_ids[match(keep_ids, unique(db$SampleID))]
    out <- db[db$SampleID %in% keep_ids, c("SampleID","Locus","A1","A2")]
    if (format == "rds") stop("RDS output requires RDS input schema")
    out$A1 <- decode_allele(out$A1, any_code)
    out$A2 <- decode_allele(out$A2, any_code)
    write.csv(out, file.path(out_dir, "db_slice.csv"), row.names = FALSE)
    return(invisible(TRUE))
  }
}

# No multibyte characters in code/comments.

# Benchmark: split a large DB RDS into chunks and run run_match_fast() per chunk.
# Outputs a single CSV summary under test/bench/ with timestamped filename.

suppressWarnings(suppressMessages({
  library(data.table)
}))

safe_dir_create <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

timestamp_tag <- function() format(Sys.time(), "%Y%m%d%H%M%S")

# ---- User-configurable ----
db_rds_path      <- "data/virtual_db_u100_S1000_seed123.rds"   # path to the big DB RDS (data.frame long format preferred)
locus_order_rds  <- "data/locus_order.rds"
result_dir       <- "test/bench"
chunk_sizes      <- c(100000, 150000, 200000, 250000)  # number of samples per chunk
tmp_chunk_dir    <- "output/_db_chunks_rds"      # where intermediate chunk RDS will be written
query_mode       <- "from_db_first"              # "from_db_first" or "from_file"
query_file       <- NA                            # if query_mode == "from_file", set CSV or RDS path
query_sample_id  <- NA                            # if "from_db_first", will auto-pick the first sample in first chunk

# ---- Helpers (project-shape tolerant) ----

read_locus_order <- function(path) {
  lo <- readRDS(path)
  if (!is.character(lo)) stop("locus_order.rds must be a character vector.")
  lo
}

# Accept:
#  (A) long data.frame with columns: SampleID, Locus, allele1, allele2
#  (B) named list: names = SampleID; each element convertible to df with same columns
load_db_any <- function(path) {
  obj <- readRDS(path)
  if (is.data.frame(obj)) {
    # long df
    df <- as.data.frame(obj, stringsAsFactors = FALSE)
    names(df) <- tolower(names(df))
    # normalize columns
    if (!"sampleid" %in% names(df) && "SampleID" %in% names(obj)) {
      df$sampleid <- obj$SampleID
    }
    names(df) <- sub("^allele1$", "allele1", names(df), ignore.case = TRUE)
    names(df) <- sub("^allele2$", "allele2", names(df), ignore.case = TRUE)
    names(df) <- sub("^locus$",   "locus",   names(df), ignore.case = TRUE)
    required <- c("sampleid","locus","allele1","allele2")
    miss <- setdiff(required, names(df))
    if (length(miss) > 0) stop(sprintf("DB df missing columns: %s", paste(miss, collapse=", ")))
    # keep only required columns in fixed order
    df <- df[, required]
    return(list(type="df", data=df))
  } else if (is.list(obj)) {
    # list of per-sample dfs; keep as-is, normalize later
    return(list(type="list", data=obj))
  } else {
    stop("Unsupported DB RDS shape. Expect df (long) or named list.")
  }
}

unique_sample_ids <- function(db_any) {
  if (db_any$type == "df") {
    return(unique(db_any$data[["SampleID"]]))
  } else {
    return(names(db_any$data))
  }
}

subset_db_by_sample_ids <- function(db_any, ids) {
  required <- c("SampleID","Locus","allele1","allele2")
  if (db_any$type == "df") {
    df <- db_any$data
    # df came normalized to lower-case names in load_db_any(); convert to proper caps here
    names(df) <- sub("^sampleid$", "SampleID", names(df))
    names(df) <- sub("^locus$",    "Locus",    names(df))
    # subset and enforce order
    sub <- df[df$SampleID %in% ids, c("SampleID","Locus","allele1","allele2")]
    rownames(sub) <- NULL
    return(sub)
  } else {
    # list => normalize each element to the same 4 columns in the same order
    dd <- lapply(ids, function(sid){
      x <- db_any$data[[sid]]
      if (is.null(x)) return(NULL)
      x <- as.data.frame(x, stringsAsFactors = FALSE)
      
      # normalize names
      nms <- tolower(names(x))
      names(x) <- nms
      if (!"sampleid" %in% names(x)) x$sampleid <- sid
      names(x) <- sub("^allele1$", "allele1", names(x))
      names(x) <- sub("^allele2$", "allele2", names(x))
      names(x) <- sub("^locus$",   "locus",   names(x))
      names(x) <- sub("^sampleid$","SampleID", names(x))  # set final casing for id
      names(x) <- sub("^locus$",   "Locus",    names(x))
      
      # add any missing required columns as NA
      miss <- setdiff(required, names(x))
      for (m in miss) x[[m]] <- NA_character_
      
      # keep only required in fixed order
      x <- x[, required]
      # fill Locus/allele if they were missing entirely
      x$Locus   <- as.character(x$Locus)
      x$allele1 <- as.character(x$allele1)
      x$allele2 <- as.character(x$allele2)
      x$SampleID <- as.character(x$SampleID)
      
      x
    })
    dd <- Filter(Negate(is.null), dd)
    if (length(dd) == 0) {
      stop("No matching samples found in list DB for given ids.")
    }
    df <- do.call(rbind, dd)
    rownames(df) <- NULL
    return(df)
  }
}
make_chunks <- function(ids, chunk_size) {
  split(ids, ceiling(seq_along(ids) / chunk_size))
}

# Build a query df:
# - from_db_first: take the first SampleID present in the given df chunk
# - from_file: load CSV/RDS with columns Locus, allele1, allele2 (SampleID optional)
build_query_df <- function(db_chunk_df, mode="from_db_first", file=NA) {
  if (mode == "from_db_first") {
    sid <- db_chunk_df$SampleID[1]
    q <- subset(db_chunk_df, SampleID == sid, c("Locus","allele1","allele2"))
    # ensure column names
    names(q)[names(q)=="Allele1"] <- "allele1"
    names(q)[names(q)=="Allele2"] <- "allele2"
    return(q)
  } else {
    if (is.na(file)) stop("query_mode==from_file but query_file not set.")
    if (grepl("\\.rds$", file, ignore.case = TRUE)) {
      q <- readRDS(file)
      q <- as.data.frame(q, stringsAsFactors = FALSE)
    } else {
      q <- read.csv(file, stringsAsFactors = FALSE)
    }
    names(q) <- tolower(names(q))
    required <- c("locus","allele1","allele2")
    miss <- setdiff(required, names(q))
    if (length(miss) > 0) stop(sprintf("Query df missing columns: %s", paste(miss, collapse=", ")))
    return(q[, required])
  }
}

# Try to call run_match_fast and return elapsed + result rows
run_one_bench <- function(query_df, db_df, locus_order) {
  source("scripts/scoring_fast.R")  # ensure function is available in vanilla R session
  gc()
  t0 <- proc.time()[["elapsed"]]
  res <- run_match_fast(query_df = query_df, db_df = db_df, locus_order = locus_order)
  t1 <- proc.time()[["elapsed"]]
  elapsed <- as.numeric(t1 - t0)
  nrows <- if (is.data.frame(res)) nrow(res) else NA_integer_
  list(elapsed_sec = elapsed, result_nrow = nrows)
}

# ---- Main ----
safe_dir_create(result_dir)
safe_dir_create(tmp_chunk_dir)

cat("[BENCH] loading DB and locus_order...\n")
db_any <- load_db_any(db_rds_path)
ids <- unique_sample_ids(db_any)
total_ids <- length(ids)
lo <- read_locus_order(locus_order_rds)
tag <- timestamp_tag()

summary_rows <- list()

for (cs in chunk_sizes) {
  cat(sprintf("[BENCH] chunk_size=%d\n", cs))
  chunks <- make_chunks(ids, cs)
  cat(sprintf("[BENCH] total samples=%d, chunks=%d\n", total_ids, length(chunks)))
  idx <- 1L
  for (ids_chunk in chunks) {
    cat(sprintf("  - building chunk %d/%d ...\n", idx, length(chunks)))
    db_chunk_df <- subset_db_by_sample_ids(db_any, ids_chunk)
    chunk_file <- file.path(tmp_chunk_dir, sprintf("db_chunk_%s_cs%d_%03d.rds", tag, cs, idx))
    saveRDS(db_chunk_df, chunk_file)
    
    # pick/build query
    q_df <- build_query_df(db_chunk_df, mode = query_mode, file = query_file)
    
    # bench
    cat("    run_match_fast() ...\n")
    bench <- run_one_bench(q_df, db_chunk_df, lo)
    
    summary_rows[[length(summary_rows)+1L]] <- data.frame(
      ts = Sys.time(),
      dataset_id = basename(db_rds_path),
      chunk_size = cs,
      chunk_index = idx,
      samples_in_chunk = length(ids_chunk),
      rows_in_chunk = nrow(db_chunk_df),
      elapsed_sec = bench$elapsed_sec,
      result_nrow = bench$result_nrow,
      stringsAsFactors = FALSE
    )
    
    idx <- idx + 1L
    gc()
  }
}

summary_df <- do.call(rbind, summary_rows)
out_csv <- file.path(result_dir, sprintf("bench_run_match_fast_%s.csv", tag))
write.csv(summary_df, out_csv, row.names = FALSE)

cat("\n[BENCH DONE] summary => ", out_csv, "\n", sep="")

# scripts/devtools/slice_query_fromdb.R
# No multibyte chars. Slice 1 sample from DB to human-readable CSV.

# --- defaults (house style) ---
opt <- list(
  db             = "data/virtual_db_u100_S1000_seed123.rds",
  sample         = 1L,
  out_csv        = "data/query_profile_seed123.csv",
  allow_defaults = TRUE
)

# track whether user explicitly provided options (to avoid "same-as-default" false positives)
seen <- list(db=FALSE, sample=FALSE, out_csv=FALSE, allow_defaults=FALSE)

# --- tiny arg parser: supports --key=value and --key value ---
argv <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(argv)) {
  a <- argv[i]
  if (startsWith(a, "--")) {
    kv <- sub("^--", "", a)
    if (grepl("=", kv, fixed = TRUE)) {
      key <- sub("=.*$", "", kv)
      val <- sub("^[^=]*=", "", kv)
    } else {
      key <- kv
      if ((i + 1L) <= length(argv) && !startsWith(argv[i + 1L], "--")) {
        val <- argv[i + 1L]; i <- i + 1L
      } else {
        val <- TRUE
      }
    }
    if (key %in% names(opt)) {
      # coerce types
      if (key %in% c("db","out_csv"))       opt[[key]] <- as.character(val)
      else if (key %in% c("sample"))        opt[[key]] <- as.integer(val)
      else if (key %in% c("allow_defaults"))opt[[key]] <- as.logical(val)
      seen[[key]] <- TRUE
    } else {
      warning(sprintf("Unknown option ignored: --%s", key))
    }
  }
  i <- i + 1L
}

# --- validation: require explicit args if allow_defaults=FALSE ---
if (!isTRUE(opt$allow_defaults)) {
  need <- c("db","out_csv")
  lacking <- need[!unlist(seen[need])]
  if (length(lacking)) {
    stop(sprintf(
      "[FATAL] The following options were not explicitly passed: %s\nPass them (e.g. --db=... --out_csv=...) or add --allow_defaults=TRUE.",
      paste(lacking, collapse=", ")
    ))
  }
}

# --- load DB ---
if (!file.exists(opt$db)) stop(sprintf("DB not found: %s", opt$db))
db <- readRDS(opt$db)
if (is.null(db$A1) || is.null(db$A2) || is.null(db$locus_ids))
  stop("Invalid DB index RDS (missing A1/A2/locus_ids)")

S <- nrow(db$A1)
s <- as.integer(opt$sample)
if (is.na(s) || s < 1L || s > S) stop(sprintf("sample must be 1..%d (got %s)", S, as.character(opt$sample)))

# --- fetch one sample (x100 integer -> human-readable) ---
a1 <- if (is.matrix(db$A1)) db$A1[s,] else vapply(db$A1, `[`, integer(1), s)
a2 <- if (is.matrix(db$A2)) db$A2[s,] else vapply(db$A2, `[`, integer(1), s)

decode <- function(z, ANY=9999L){
  if (is.na(z)) return("")
  if (as.integer(z) == as.integer(ANY)) return("any")
  v <- as.numeric(z) / 100
  if (abs(v - round(v)) < 1e-9) as.character(as.integer(round(v))) else format(v, trim=TRUE, nsmall=1)
}

out_df <- data.frame(
  Locus   = db$locus_ids,
  allele1 = vapply(a1, decode, ""),
  allele2 = vapply(a2, decode, ""),
  stringsAsFactors = FALSE
)

# --- write ---
dir.create(dirname(opt$out_csv), recursive=TRUE, showWarnings=FALSE)
write.csv(out_df, opt$out_csv, row.names=FALSE)
cat("[OK] wrote query CSV (readable): ", opt$out_csv, " (sample=", s, ")\n", sep="")

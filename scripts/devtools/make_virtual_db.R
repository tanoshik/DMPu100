# No multibyte chars. Build indexed DB from long freq table.
# Accepts columns: (Locus, Allele, P) or (Locus, Allele, Freq) or (Locus, Allele, Count)
# Prefix-invariant: first N samples are identical for a fixed seed even if DB size S changes.

# ----------------- arg parse (both "--k=v" and "--k v") -----------------
defaults <- list(
  size     = 2000000L,
  seed     = 123L,
  out      = "data/virtual_db_u100_S2000000_seed123.rds",
  any_code = 9999L,
  freq_rds = "data/freq_table.rds"
)
opt <- defaults
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
      if (i + 1L <= length(argv) && !startsWith(argv[i + 1L], "--")) {
        val <- argv[i + 1L]; i <- i + 1L
      } else {
        val <- "TRUE"
      }
    }
    if (key %in% c("out","freq_rds")) {
      opt[[key]] <- as.character(val)
    } else if (key %in% c("size","seed","any_code")) {
      opt[[key]] <- as.integer(val)
    }
  }
  i <- i + 1L
}

# ----------------- sanity -----------------
if (!file.exists("data/locus_order.rds")) stop("missing data/locus_order.rds")
if (!file.exists(opt$freq_rds)) stop("missing freq table: ", opt$freq_rds)

# ----------------- load resources -----------------
locus_order <- readRDS("data/locus_order.rds")     # character vector (length L)
ft <- readRDS(opt$freq_rds)                        # data.frame long format

# normalize colnames (case-insensitive)
nn <- names(ft)
nn <- sub("(?i)^locus$","Locus", nn, perl=TRUE)
nn <- sub("(?i)^allele$","Allele", nn, perl=TRUE)
nn <- sub("(?i)^p$","P", nn, perl=TRUE)
nn <- sub("(?i)^freq$","Freq", nn, perl=TRUE)
nn <- sub("(?i)^count$","Count", nn, perl=TRUE)
names(ft) <- nn

need <- c("Locus","Allele")
if (!all(need %in% names(ft))) stop("freq table must have Locus and Allele")

ft$Locus  <- as.character(ft$Locus)
ft$Allele <- as.character(ft$Allele)

# ----------------- choose prob column, per locus normalize if necessary -----------------
make_prob <- function(tbl) {
  if ("P" %in% names(tbl)) {
    p <- as.numeric(tbl$P)
  } else if ("Freq" %in% names(tbl)) {
    p <- as.numeric(tbl$Freq)
  } else if ("Count" %in% names(tbl)) {
    p <- as.numeric(tbl$Count)
  } else stop("no prob source at locus: ", unique(tbl$Locus))
  p[!is.finite(p) | p < 0] <- 0
  s <- sum(p)
  if (s <= 0) stop("non-positive sum of probs at locus: ", unique(tbl$Locus))
  p / s
}
ft_by <- split(ft, ft$Locus)

# ----------------- tiny 31-bit hash & seed-for helpers -----------------
# 0..(2^31-1) mod (2^31-1) arithmetic（set.seed の有効域に収める）
.M31 <- 2147483647.0  # doubleで扱い、最後に整数化
hash32 <- function(x) {
  # simple polynomial rolling hash (deterministic, no deps)
  bytes <- as.integer(charToRaw(x))
  h <- 0.0
  for (b in bytes) {
    h <- (h * 131.0 + b) %% .M31
  }
  as.integer(h)
}
seed_for <- function(seed, loc, field) {
  # combine user seed (integer) with hashed (loc#field) into a 31-bit seed
  base <- as.numeric(abs(as.integer(seed)) %% .M31)
  h    <- as.numeric(abs(hash32(paste0(loc, "#", field))) %% .M31)
  z <- as.integer((base + h) %% .M31)
  if (z == 0L) z <- 1L
  z
}

# ----------------- encode allele "xx.y" → integer ×100 -----------------
encode <- function(a_char) as.integer(round(as.numeric(a_char) * 100))

# ----------------- build DB (prefix-invariant) -----------------
S <- as.integer(opt$size)
L <- length(locus_order)
sample_ids <- sprintf("S%07d", seq_len(S))
locus_ids  <- as.character(locus_order)
A1 <- matrix(NA_integer_, nrow=S, ncol=L)
A2 <- matrix(NA_integer_, nrow=S, ncol=L)

for (j in seq_len(L)) {
  loc <- locus_ids[j]
  tbl <- ft_by[[loc]]
  if (is.null(tbl) || nrow(tbl) == 0) stop("no freq rows for locus: ", loc)
  p <- make_prob(tbl)
  
  # --- A1: unique RNG stream per (loc, "A1") ---
  set.seed(seed_for(opt$seed, loc, "A1"))
  idx1 <- sample.int(nrow(tbl), S, replace=TRUE, prob=p)
  e1   <- encode(tbl$Allele[idx1])
  
  # --- A2: unique RNG stream per (loc, "A2") ---
  set.seed(seed_for(opt$seed, loc, "A2"))
  idx2 <- sample.int(nrow(tbl), S, replace=TRUE, prob=p)
  e2   <- encode(tbl$Allele[idx2])
  
  # order so that e1 <= e2
  swap <- e1 > e2
  if (any(swap)) { tmp <- e1[swap]; e1[swap] <- e2[swap]; e2[swap] <- tmp }
  
  A1[, j] <- e1
  A2[, j] <- e2
}

db_index <- list(sample_ids=sample_ids, locus_ids=locus_ids, A1=A1, A2=A2)
dir.create(dirname(opt$out), recursive=TRUE, showWarnings=FALSE)
saveRDS(db_index, file=opt$out)
cat("[OK] wrote DB: ", opt$out, " S=", S, " L=", L, "\n", sep="")

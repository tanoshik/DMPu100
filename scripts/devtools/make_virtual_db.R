# scripts/devtools/make_virtual_db.R
# No multibyte chars. Build indexed DB from long freq table.
# Accepts columns: (Locus, Allele, P) or (Locus, Allele, Freq) or (Locus, Allele, Count)

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  size     = 2000000L,
  seed     = 123L,
  out      = "data/virtual_db_u100_S2000000_seed123.rds",
  any_code = 9999L,
  freq_rds = "data/freq_table.rds"   # we now consume your existing file
)
for (a in args) {
  kv <- strsplit(a, "=", fixed=TRUE)[[1]]
  if (length(kv)==2) opt[[sub("^--?","",kv[1])]] <- kv[2]
}
opt$size     <- as.integer(opt$size)
opt$seed     <- as.integer(opt$seed)
opt$any_code <- as.integer(opt$any_code)

if (!file.exists("data/locus_order.rds")) stop("missing data/locus_order.rds")
if (!file.exists(opt$freq_rds)) stop("missing freq table: ", opt$freq_rds)

locus_order <- readRDS("data/locus_order.rds")
ft <- readRDS(opt$freq_rds)
if (!is.data.frame(ft)) stop("freq table must be a data.frame")

# normalize colnames
nn <- names(ft)
nn <- sub("(?i)^locus$","Locus", nn, perl=TRUE)
nn <- sub("(?i)^allele$","Allele", nn, perl=TRUE)
nn <- sub("(?i)^p$","P", nn, perl=TRUE)
nn <- sub("(?i)^freq$","Freq", nn, perl=TRUE)
nn <- sub("(?i)^count$","Count", nn, perl=TRUE)
names(ft) <- nn

need <- c("Locus","Allele")
if (!all(need %in% names(ft))) stop("freq table must have Locus and Allele")

# choose a probability source in priority: P > Freq > Count(normalized)
if ("P" %in% names(ft)) {
  prob <- as.numeric(ft$P)
} else if ("Freq" %in% names(ft)) {
  prob <- as.numeric(ft$Freq)
} else if ("Count" %in% names(ft)) {
  prob <- NA_real_  # mark to normalize per-locus later from counts
} else {
  stop("freq table needs one of: P, Freq, or Count")
}

ft$Locus  <- as.character(ft$Locus)
ft$Allele <- as.character(ft$Allele)
if (!is.null(prob)) ft$P <- prob

# split by locus and ensure we have probs
ft_by <- split(ft, ft$Locus)
make_prob <- function(tbl) {
  if ("P" %in% names(tbl)) {
    p <- tbl$P
  } else if ("Freq" %in% names(tbl)) {
    p <- tbl$Freq
  } else if ("Count" %in% names(tbl)) {
    p <- as.numeric(tbl$Count)
  } else stop("no prob source at locus: ", unique(tbl$Locus))
  p <- as.numeric(p); p[!is.finite(p) | p < 0] <- 0
  s <- sum(p)
  if (s <= 0) stop("non-positive sum of probs at locus: ", unique(tbl$Locus))
  p / s
}

encode <- function(a_char) as.integer(round(as.numeric(a_char) * 100))

set.seed(opt$seed)
S <- opt$size
L <- length(locus_order)
sample_ids <- sprintf("S%07d", seq_len(S))
locus_ids  <- as.character(locus_order)
A1 <- matrix(NA_integer_, nrow=S, ncol=L)
A2 <- matrix(NA_integer_, nrow=S, ncol=L)

for (j in seq_len(L)) {
  loc <- locus_ids[j]
  tbl <- ft_by[[loc]]
  if (is.null(tbl) || nrow(tbl)==0) stop("no freq rows for locus: ", loc)
  p <- make_prob(tbl)
  # sample 2 alleles independently, then order (no ANY in freq table)
  a1 <- tbl$Allele[sample.int(nrow(tbl), S, replace=TRUE, prob=p)]
  a2 <- tbl$Allele[sample.int(nrow(tbl), S, replace=TRUE, prob=p)]
  e1 <- encode(a1); e2 <- encode(a2)
  swap <- e1 > e2
  if (any(swap)) { tmp <- e1[swap]; e1[swap] <- e2[swap]; e2[swap] <- tmp }
  A1[,j] <- e1; A2[,j] <- e2
}

db_index <- list(sample_ids=sample_ids, locus_ids=locus_ids, A1=A1, A2=A2)
saveRDS(db_index, file=opt$out)
cat("[OK] wrote DB: ", opt$out, " S=", S, " L=", L, "\n", sep="")

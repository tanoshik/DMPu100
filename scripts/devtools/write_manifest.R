# scripts/devtools/write_manifest.R
# No multibyte chars. Write manifest for virtual DBs to output/manifest.

suppressWarnings({ options(stringsAsFactors = FALSE) })
library(tools)

dir.create("output/manifest", recursive = TRUE, showWarnings = FALSE)
out_csv <- "output/manifest/virtual_db_manifest.csv"

list_vdb <- function() {
  list.files("data", pattern="^virtual_db_u100_S[0-9]+_seed[0-9]+\\.rds$", full.names=TRUE)
}

parse_meta <- function(path) {
  m <- regexec("virtual_db_u100_S([0-9]+)_seed([0-9]+)\\.rds$", basename(path))
  g <- regmatches(basename(path), m)[[1]]
  if (length(g) < 3) return(list(S=NA_integer_, seed=NA_integer_))
  list(S=as.integer(g[2]), seed=as.integer(g[3]))
}

files <- list_vdb()
if (!length(files)) stop("no DB rds in output/data")

rows <- lapply(files, function(p){
  meta <- parse_meta(p)
  sz <- file.info(p)$size
  md <- if (!is.na(sz) && sz > 0) as.character(md5sum(p)) else NA_character_
  data.frame(path=p, S=meta$S, seed=meta$seed, size=sz,
             md5=md, placeholder=FALSE, created_at=format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
             stringsAsFactors=FALSE)
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$S, tab$seed),]
write.csv(tab, out_csv, row.names = FALSE)
cat("[OK] wrote: ", out_csv, "\n", sep = "")

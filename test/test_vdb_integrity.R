# test/test_vdb_integrity.R
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1] else "data/virtual_db_u100_S2000000_seed123.rds"
db <- readRDS(path)

stopifnot(is.character(db$sample_ids), length(db$sample_ids) >= 1L)
stopifnot(is.character(db$locus_ids),  length(db$locus_ids)  == 21L)

ok_A1 <- (is.matrix(db$A1) && all(dim(db$A1)[2] == 21L)) || (is.list(db$A1) && length(db$A1) == 21L)
ok_A2 <- (is.matrix(db$A2) && all(dim(db$A2)[2] == 21L)) || (is.list(db$A2) && length(db$A2) == 21L)
stopifnot(ok_A1, ok_A2)

# 先頭サンプルを抽出して型を確認（行=サンプル、列=ローカス想定）
s <- 1L
L <- if (!is.null(colnames(db$A1))) colnames(db$A1) else db$locus_ids
getS <- function(A, s, L) if (is.matrix(A)) A[s, ] else vapply(A, `[`, numeric(1), s)
a1 <- getS(db$A1, s, L); a2 <- getS(db$A2, s, L)
stopifnot(length(a1) == 21L, length(a2) == 21L, is.integer(a1) || is.numeric(a1), is.integer(a2) || is.numeric(a2))

cat("[OK] VDB integrity:", path, " samples=", nrow(db$A1), " loci=21\n")

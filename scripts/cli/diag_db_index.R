# scripts/cli/diag_db_index.R
# Minimal inspector for the DB RDS and matcher expectations.
# No multibyte chars in code/comments.

args <- commandArgs(trailingOnly = TRUE)
db_path <- if (length(args) >= 1) args[[1]] else "data/virtual_db_u100_S1000_seed123.rds"

cat("=== diag_db_index ===\n")
cat("[db_path]", db_path, "\n")

# Load matcher to get .assert_db_index and run_match_fast
source("scripts/matcher_fast.R")

# Load object
obj <- readRDS(db_path)

cat("\n--- OBJ CLASS ---\n")
print(class(obj))

cat("\n--- OBJ NAMES (head) ---\n")
if (is.list(obj)) {
  nm <- names(obj)
  print(head(nm, 50))
} else {
  cat("(not a list)\n")
}

cat("\n--- HAS .assert_db_index ---\n")
print(exists(".assert_db_index"))

if (exists(".assert_db_index")) {
  cat("\n--- TRY .assert_db_index ---\n")
  tryCatch({
    .assert_db_index(obj)
    cat("[OK] assert passed\n")
  }, error = function(e) {
    cat("[FAIL] assert: ", conditionMessage(e), "\n", sep = "")
  })
}

cat("\n--- run_match_fast signature ---\n")
print(names(formals(run_match_fast)))
cat("=== end diag ===\n")

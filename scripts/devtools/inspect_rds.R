# scripts/devtools/inspect_rds.R
# Simple RDS inspector for virtual DB objects.
# No multibyte chars in code/comments.

inspect_rds <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  obj <- readRDS(path)
  cat("== class(obj): ", paste(class(obj), collapse = ", "), "\n")
  cat("== typeof(obj): ", typeof(obj), "\n")
  
  if (is.data.frame(obj)) {
    cat("== data.frame dim: ", nrow(obj), " x ", ncol(obj), "\n", sep = "")
    cat("== names: ", paste(names(obj), collapse = ", "), "\n")
    print(utils::head(obj, 5))
  } else if (is.list(obj)) {
    cat("== list length: ", length(obj), "\n", sep = "")
    nm <- names(obj)
    if (!is.null(nm)) cat("== names[1:10]: ", paste(utils::head(nm, 10), collapse = ", "), "\n")
    cat("== element[[1]] class: ", paste(class(obj[[1]]), collapse = ", "), "\n")
    if (is.data.frame(obj[[1]])) {
      cat("-- element[[1]] is data.frame with names: ",
          paste(names(obj[[1]]), collapse = ", "), "\n")
      print(utils::head(obj[[1]], 5))
    } else if (is.list(obj[[1]])) {
      cat("-- element[[1]] is list; showing str(elem[[1]])...\n")
      print(utils::str(obj[[1]], max.level = 1))
    } else {
      cat("-- element[[1]] typeof: ", typeof(obj[[1]]), "\n")
      print(utils::str(obj[[1]], max.level = 1))
    }
  } else {
    cat("== Unsupported top-level type; showing str(obj) level 1...\n")
    print(utils::str(obj, max.level = 1))
  }
  invisible(obj)
}

# run
path <- "data/virtual_db_u100_S1000_seed123.rds"
obj <- inspect_rds(path)

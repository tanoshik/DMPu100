# scripts/matcher_onepass.R
# One-pass matcher CLI wrapper for dmp_match_cpp/dmp_hist_cpp (strict)
# Strict: 入口/出口は寛容にせず、補完はfix_query/fix_dbに委譲

suppressPackageStartupMessages({
  library(optparse)
  library(Rcpp)
  library(jsonlite)
  library(tools)
})

# SCORE_TABLE (index: 0..15)
SCORE_TABLE <- as.integer(c(
  0,1,1,1,
  1,1,2,2,
  1,2,1,2,
  1,2,2,2
))

# ---- CLI options ----
option_list <- list(
  make_option("--db",      type="character", help="path to db.rds"),
  make_option("--query",   type="character", help="path to query.rds"),
  make_option("--out",     type="character", default="output/scores.csv"),
  make_option("--report",  type="character", default="all"),       # all|top|hist
  make_option("--n_cap",   type="integer",   default=200L),
  make_option("--score_min", type="integer", default=0L),
  make_option("--any_code",  type="integer", default=9999L),
  make_option("--idf_csv",   type="character", default=NA),
  make_option("--force_all_loci", type="integer", default=0L),
  make_option("--h2a_on", type="integer", default=0L),
  make_option("--debug",   type="integer",   default=1L),
  make_option("--path",    type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$db) || is.null(opt$query)) stop("--db and --query are required")
if (!file.exists(opt$db)) stop(sprintf("db not found: %s", opt$db))
if (!file.exists(opt$query)) stop(sprintf("query not found: %s", opt$query))

# ---- output path helpers ----
ensure_csv <- function(p) { if (grepl("\\.csv$",p,ignore.case=TRUE)) p else paste0(p,".csv") }
ensure_json<- function(p) { if (grepl("\\.json$",p,ignore.case=TRUE)) p else paste0(p,".json") }

# ---- compile core ----
sourceCpp("src/dmp_match.cpp")

t0_total <- proc.time()[["elapsed"]]

# ---- load DB ----
db <- readRDS(opt$db)
S <- length(db$sample_ids); L <- length(db$locus_ids)
A1m <- as.matrix(db$A1); A2m <- as.matrix(db$A2)
storage.mode(A1m) <- "integer"; storage.mode(A2m) <- "integer"

# ---- h2a ----
if (as.integer(opt$h2a_on)==1L) {
  ANY <- as.integer(opt$any_code)
  mask <- (A1m==A2m)&(A1m!=ANY)
  A2m[mask]<-ANY
  if (opt$debug==1L) cat(sprintf("[DBG] h2a_on: %d masked\n",sum(mask)))
}

# ---- load query ----
q <- readRDS(opt$query)
if (!is.list(q) || !all(c("q1","q2") %in% names(q))) stop("query RDS must have list(q1,q2)")
if (!(length(q$q1)==L && length(q$q2)==L)) stop("query length mismatch")

# ---- idf mask ----
make_all_ones <- function(L){ bits<-0L; for(j in 0:(L-1)){ bits<-bitwOr(bits,bitwShiftL(1L,j)) }; bits }
idf_mask_bits <- make_all_ones(L)
mask_name <- "normal"
if (!is.na(opt$idf_csv) && nzchar(opt$idf_csv) && file.exists(opt$idf_csv)) {
  mask_name <- "gf_idf_mask"
  dfm <- read.csv(opt$idf_csv, stringsAsFactors=FALSE, check.names=FALSE)
  nms <- tolower(names(dfm)); names(dfm)<-nms
  lc <- which(nms %in% c("locus","marker","locus_id"))[1]
  bc <- which(nms %in% c("bit","mask","enabled"))[1]
  if (!is.na(lc)&&!is.na(bc)){
    idf_mask_bits<-0L; mm<-match(db$locus_ids,dfm[[lc]])
    for(j in 0:(L-1)){
      v<-1L; if(!is.na(mm[j+1])){
        vv<-suppressWarnings(as.integer(dfm[[bc]][mm[j+1]]))
        if(!is.na(vv)) v<-ifelse(vv!=0L,1L,0L)
      }
      if(v==1L) idf_mask_bits<-bitwOr(idf_mask_bits,bitwShiftL(1L,j))
    }
  }
}
if (as.integer(opt$force_all_loci)==1L) idf_mask_bits<-make_all_ones(L)

# ---- run core ----
rep_mode <- tolower(opt$report)
compute_scores <- rep_mode %in% c("all","top")
compute_hist   <- rep_mode %in% c("all","hist")

dir.create(dirname(opt$out),recursive=TRUE,showWarnings=FALSE)
scores_csv <- ensure_csv(opt$out)
opts_json  <- ensure_json(file.path(dirname(scores_csv),
                                    sprintf("opts_%s.json",tools::file_path_sans_ext(basename(scores_csv)))))

res <- data.frame()
if (compute_scores){
  res <- dmp_match_cpp(A1=A1m,A2=A2m,q1=q$q1,q2=q$q2,
                       score_table=SCORE_TABLE,
                       idf_mask_bits=as.integer(idf_mask_bits),
                       opts=list(any_code=as.integer(opt$any_code),
                                 score_min=as.integer(opt$score_min),
                                 n_cap=as.integer(opt$n_cap),
                                 force_all_loci=as.logical(opt$force_all_loci)),
                       sample_ids=db$sample_ids)
  write.csv(res,scores_csv,row.names=FALSE)
  if(opt$debug==1L) cat(sprintf("[OK] wrote %s (%d rows)\n",scores_csv,nrow(res)))
}

if (compute_hist){
  hist_csv<-ensure_csv(file.path(dirname(scores_csv),"hist"))
  hdf<-dmp_hist_cpp(A1=A1m,A2=A2m,q1=q$q1,q2=q$q2,
                    score_table=SCORE_TABLE,
                    idf_mask_bits=as.integer(idf_mask_bits),
                    opts=list(any_code=as.integer(opt$any_code),
                              force_all_loci=as.logical(opt$force_all_loci)))
  write.csv(hdf,hist_csv,row.names=FALSE)
  if(opt$debug==1L) cat(sprintf("[OK] wrote %s (%d rows)\n",hist_csv,nrow(hdf)))
}

opts_list <- list(db_path=normalizePath(opt$db,winslash="/"),
                  query_path=normalizePath(opt$query,winslash="/"),
                  mask_name=mask_name,
                  any_code=as.integer(opt$any_code),
                  force_all_loci=as.logical(opt$force_all_loci),
                  h2a_on=as.logical(opt$h2a_on))
writeLines(toJSON(opts_list,auto_unbox=TRUE,pretty=TRUE),opts_json)
if(opt$debug==1L) cat(sprintf("[OK] wrote %s\n",opts_json))

# ---- time.log ----
core_name<-sprintf("GF/%s/0.0/%s/%d",mask_name,if(opt$h2a_on==1L)"on" else "off",S)
total_sec<-proc.time()[["elapsed"]]-t0_total
header<-"LOAD_DB_SEC,LOAD_Q_SEC,COMP_SEC,TOTAL_SEC,PEAK_MiB,name"
line<-sprintf("NA,NA,NA,%.3f,NA,%s",total_sec,core_name)
logp<-file.path(dirname(scores_csv),"time.log")
if(!file.exists(logp)) writeLines(header,logp)
cat(paste0(line,"\n"),file=logp,append=TRUE)
if(opt$debug==1L) cat("[DBG] done\n")

# scripts/devtools/build_all_dbs_and_query.R
# Build 2M -> slice to 1M, 100K, 1K; make query; write manifest; all under output/.

system("Rscript scripts/devtools/make_virtual_db.R --size=2000000 --seed=123 --out=output/data/virtual_db_u100_S2000000_seed123.rds", intern=FALSE)
system("Rscript scripts/devtools/slice_virtual_db.R --src=output/data/virtual_db_u100_S2000000_seed123.rds --size=1000000 --out=output/data/virtual_db_u100_S1000000_seed123.rds", intern=FALSE)
system("Rscript scripts/devtools/slice_virtual_db.R --src=output/data/virtual_db_u100_S2000000_seed123.rds --size=100000  --out=output/data/virtual_db_u100_S100000_seed123.rds",  intern=FALSE)
system("Rscript scripts/devtools/slice_virtual_db.R --src=output/data/virtual_db_u100_S2000000_seed123.rds --size=1000    --out=output/data/virtual_db_u100_S1000_seed123.rds",   intern=FALSE)
system("Rscript scripts/devtools/make_query_seed123.R --seed=123 --out_csv=output/sample/query_profile_seed123.csv", intern=FALSE)
system("Rscript scripts/devtools/write_manifest.R", intern=FALSE)

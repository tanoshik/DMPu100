## Quick Start (CLI)

### 0) 前提
- R 4.5.0（Windows）/ RStudio 2024.12.1
- GlobalFiler 21 loci 固定
- `data/locus_order.rds`, `data/freq_table.rds` 在中
- Query CSV ヘッダは **Locus, allele1, allele2**（ANY可, 小数可）。CLI側で自動正規化し、ANY→9999で数値化します。

### 1) 仮想DBを作る（インデックス済みRDS）
```bash
Rscript scripts/cli/dmp_make_virtual_db.R --size=1000 --seed=123
# => data/virtual_db_u100_S1000_seed123.rds

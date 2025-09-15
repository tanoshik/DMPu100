# siyo_specs.md — DMP 指示書（仕様本文）

> 本文は「制約・規則・運用ルール」を **指示書から分離** し、  
> `schema_registry.json` と各 `*.schema.json` を **唯一の正** としたうえで、開発・配布・連携の基準を示す。

## 1. 序文

### 1.1 目的
- R＋Shinyで GlobalFiler 21ローカスSTRプロファイルの検索・マッチング・ランキングを行う DMPu100 の**仕様本文**。
- スキーマ（`schema_registry.json` / family別 `*.schema.json`）を唯一の正として、入出力の**契約**を固定する。

### 1.2 適用範囲
- DMPu100（実験ベース）／配布は dist（クリーン）。  
- Shiny UI、スコアリング、DB照合、I/Oに関する**契約事項**。

### 1.3 非目標（Out of Scope）
- 外部公開ドキュメント（後日整理）。  
- 大規模ベンチや内部テスト手順の詳細（test/ 配下で別管理）。

---

## 2. 開発環境

### 2.1 ツールチェーン
- R 4.5.0（Windows）  
- RStudio 2025.09.0+387 "Cucumberleaf Sunflower" Release (af5fc22a687c0f462ee27c6afeeee38ee46507b9, 2025-09-11) for windows  
- Quarto 1.5.57  
- 入出力：`.CSV` / `.RDS`。  
- **RDSの構造は `*.schema.json` を唯一の正**（`schema_registry.json` に従う）。

### 2.2 実行ポリシー
- Git Bash を標準。  
- RStudioのコンソールも併用（特にエラー確認などで積極的に使用）。  
- PowerShellは制約が多いためなるべく使用しない。  
- コード中は2バイト文字を使わない（コメントは可）。  
- 最終的な配布公開版はコメントもASCIIのみとするが、レビュワー配布や開発ツール内では日本語コメントを許容する。

---

## 3. データ契約（Registry / Schemas）

### 3.1 Source of Truth
- `schema_registry.json` と familyごとの `*.schema.json`  
- codebook.md（人向け補足） / sample.csv（最小再現）を参照出力として併置。

### 3.2 schema_registry.json
- families に `virtual_db_u100` / `freq_table` / `locus_order` / `locus_layout` / `scores_csv`（将来）を登録。  
- `patterns` と `handler` を明記（外部schemasが無くても推定で出力可能）。

### 3.3 family schemas（*.schema.json）
- `locus_order.schema.json`：`array<string>`（GlobalFiler 21固定／UI並び・走査順の基準）  
- `locus_layout.schema.json`：`list` or `data.frame`（UI左/右などのグルーピング。正規化推奨）  
- `freq_table.schema.json`：`data.frame`（表示用頻度。ランキング未使用）  
- `virtual_db_u100.*.schema.json`：`list(sample_ids, locus_ids, A1, A2)`（インデックス格子）

### 3.4 参照出力（codebook.md / sample.csv）
- codebook.md：人向け説明（ANY_CODEや小数アレル、UI依存など）。  
- sample.csv：最小再現。**大容量化を避ける**（代表抽出）。

### 3.5 変更管理（Versioning）
- スキーマ破壊的変更時は **schema の `version` を上げる**。  
- `on_mismatch` は整備時 `"stop"`、開発時 `"warn"` に切替可能。

---

## 4. 制約・規則（Core Rules）

### 4.1 ANY_CODE / アレル整形
- ANY_CODE = **9999**。  
- 2アレルは**数値昇順**＋ANYは**右詰め**（encode/decodeで表示変換）。

### 4.2 ローカス順
- **`locus_order` に準拠（GlobalFiler 21固定）**。

### 4.3 一致ビットと表示順
- b0=(q1 vs r1), b1=(q1 vs r2), b2=(q2 vs r1), b3=(q2 vs r2)  
- 表示は **[b3 b2 b1 b0]**（内部は b0..b3）

### 4.4 スコア定義と SCORE_TABLE
- `code = b0*1 + b1*2 + b2*4 + b3*8`  
- `score = SCORE_TABLE[code+1]`（表は DMP共通規約どおり）

### 4.5 any,any の扱い
- **既定**：常に一致扱い（1111）。  
- **切替（計測時）**：「比較スキップ＋事前加点」。切替時は **spec の version を上げる**。

### 4.6 未知アレル（lenient / logging）
- lenient では **ANY へ写像**しログ可（大規模時はログ抑制可）。

### 4.7 IDF% マスク
- 1=比較有効、0=無効。**DB側無効化=該当ローカス比較対象外**（Q側ANY化とは異なる）。

### 4.8 整合性と互換（I/F変更ルール）
- 引数・構造変更時は **呼び出し元も同時修正**（テスト必須）。  
- 並行関数は統合方針を明記し、一時併存を許容。

---

## 5. 出力ポリシー

### 5.1 常時出力（scores.csv）
- 基本は **`output/scores.csv`（SampleID, Score）** のみ。

### 5.2 詳細出力（オンデマンド）
- detail は**別モジュールの要求時のみ**生成（負荷・容量対策）。

### 5.3 一周走査（TopN / Threshold）

- すべてのモードは **単一パス走査＋SoA構造** で実行する。  
- 出力制御は `--report {all, top, hist}` と `--n_cap` / `--score_min` によって決定される。

#### 選抜アルゴリズム
- **TopN**  
  - サイズ = `n_cap` の最小ヒープで上位N件を維持。  
  - 走査終了後に **Score降順・SampleID昇順** で整列し出力。  

- **Threshold**  
  - `Score >= score_min` のすべてを候補に入れる。  
  - 出力は `n_cap` 件までに制限。  
  - `n_cap=0` の場合は scores 出力を抑制し、hist のみを出力（研究用途など）。  

- **Hist**  
  - 全件走査し、Score=0..MAX_SC の頻度をカウント。  
  - 出力は hist のみ。  

#### 引数仕様
- `--n_cap` : 出力件数の上限（TopN と Threshold 共通）。  
- `--score_min` : Threshold モードでの最小スコア閾値。  
- `--report` : 出力形式を制御（all = scores＋hist, top = scoresのみ, hist = histのみ）。  
- `--report` が指定されるため、従来の `--mode` は廃止。  

#### ファイル命名規則
- scores/hist/opts の接尾辞は **`<report>_<条件>`** を共通化する。  
  - 例: `scores_topn_all.csv`, `hist_t30_top.csv`, `opts_histonly.json`  
- opts は **`opts_*.json`** とし、scores/hist と同一の接尾辞を持つ。


---

## 6. 運用ルール

### 6.1 Git・Zip・ログ
- ルート直下の adhoc `.sh`/`.ps1` は ignore。  
- 長時間ジョブ前に **1k / 1M スモーク**。  
- **make_dist_zip.R** を用いてZipを生成。  
- 多重Zipは禁止。大きなログは**1段Zip**に集約。  
- 参照Zipの優先順位は以下のとおり：  
  ⓪ 明示的に参照指示されたもの  
  ① スレッド内でアップロードされたZip（最新のものから順）  
  ② project files（スレッド引継ぎ時に更新する）  

### 6.2 実行・監視
- `nohup + disown` 可。`ps -ef | grep -i 'Rscript'` 等で監視。

### 6.3 エラーハンドリング
- `gzfile` レース回避に `sleep`。ワンライナーは**スクリプト化**。

---

## 7. 更新フロー（Documentation）
- schema出力に変更があった場合、**siyo_specs.mdの該当章**と **Lineage/ChangeLog** を更新する。  
  - → この意味：schema.jsonの内容が変わったら、それに対応する本文（例えば locus_order の仕様説明）を書き換える。  
  - → さらに、変更をLineage/ChangeLogに必ず記録する。  
- これ以外の部分を触る必要はない。

---

## 8. Lineage / ChangeLog（テンプレート）

### 8.1 Lineage（最終版に付属）
```json
{
  "origin": {
    "instruction": "DMP 指示書（本スレ・直近）",
    "siyo_specs": "siyo_specs.md v0.1.0",
    "schema_registry": "schema_registry.json v1.0.0"
  },
  "time": "<fill-when-commit>",
  "files": [
    "指示書.md",
    "siyo_index.txt",
    "siyo_specs.md",
    "schema_registry.json",
    "virtual_db_u100.schema.json",
    "freq_table.schema.json",
    "locus_order.schema.json",
    "scores_csv.schema.json",
    "scripts/devtools/rds_schema_export.R",
    "scripts/devtools/make_schema_exports.R"
  ],
  "commit": "UNKNOWN"
}
```

### 8.2 ChangeLog（例）
**Changed**
- 指示書を「制約・規則・運用ルール」に特化し、仕様本文を `siyo_specs.md` へ分離

**Rationale**
- 再現性（schema を唯一の正）・軽量化（scores.csv 常時）・保守性（detail 後出し）

**Risk**
- 既存スクリプトとの I/F 差（scores の列名・型）。schema 不一致で停止する設定

**Test**
- 1k スモーク（RDS→exports）、既存 fast 経路の `SampleID,Score` と整合比較

**Rollback**
- `schema_registry.json` を旧版へ、devtools 追加分を削除

**差分アンカー（例）**
```
@@ scripts/devtools/make_schema_exports.R:L1-L40 @@
（新規追加：バッチエクスポータのヘッダ・main 入口）
```

---

## 9. ChangeLog テンプレート

変更が発生した場合、以下の形式で記録する。  
（※schema出力の変更やspec本文修正に必ず対応させること）

### Changed
- 具体的に何を変更したか

### Rationale
- なぜ変更したのか

### Risk
- 変更に伴うリスクや影響範囲

### Test
- 確認方法

### Rollback
- 元に戻す方法

---

**記述ルール**  
- 必ず **siyo_specs.md に追記**する。  
- schema出力が変わった場合は **対応する本文の章も同時に更新**する。  
- 1回の変更ごとに1ブロック追加し、過去ログは削除せず保持する。

---

---

## ChangeLog

### Changed
- 5.3 出力ポリシーを更新  
  - `--n_cap` に名称統一（旧: display_limit, top_n）  
  - `--score_min` を Threshold 用に導入（旧: threshold 引数と区別）  
  - `--report` に統合し、`--mode` を廃止  
  - ファイル命名を `opts_*.json` に統一

### Rationale
- 引数とファイル名の一貫性を確保  
- 出力の分岐（all/top/hist）の明確化  

### Risk
- 既存スクリプトで `--mode` を使用している場合は動作しない  

### Test
- 1k / 1M スモークで report=all/top/hist, n_cap=0 ケースを確認済み  

### Rollback
- `siyo_specs.md` を前版に戻し、スクリプト側で旧 `display_limit` / `--mode` を復活させる

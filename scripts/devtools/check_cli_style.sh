#!/usr/bin/env bash
set -euo pipefail

# scripts/devtools/check_cli_style.sh
# 目的: Rスクリプトの optparse 長旗が snake-case のみであることを検査。
# - kebab (--foo-bar) の長旗が make_option に含まれていないこと
# - "usage:" 行があれば表示（無くても WARN 止まり）

ROOT="${1:-.}"
echo "[INFO] scanning R scripts under ${ROOT}/scripts"

# 1) kebab を検出（NG）
KebabHits=$(grep -RnoE 'make_option\([^)]*--[a-z0-9]+-[a-z0-9]+' "${ROOT}/scripts" || true)
if [[ -n "$KebabHits" ]]; then
  echo "[ERR] kebab-case long flags found (optparse/getopt in R does not support):"
  echo "$KebabHits"
  exit 1
fi
echo "[OK] no kebab-case long flags in make_option"

# 2) usage の有無（参考表示）
UsageHits=$(grep -RnoE '^[[:space:]]*#?[[:space:]]*usage:' "${ROOT}/scripts" || true)
if [[ -n "$UsageHits" ]]; then
  echo "[INFO] usage headings:"
  echo "$UsageHits"
else
  echo "[WARN] no 'usage:' headings found; consider adding usage comment lines."
fi

echo "[DONE] CLI style check finished."

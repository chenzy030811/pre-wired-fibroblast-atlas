#!/bin/bash
# Run all 4 species subsets, each in a fresh R process
# so Mac's 24 GB per-process ceiling doesn't blow up on Mus/Acomys.
#
# Usage:  bash run_all_subsets.sh
#
# Set R_MAX_VSIZE = 64 GB to allow R to grow beyond default 24 GB
# (uses swap if physical RAM insufficient — will be slow but works).

export R_MAX_VSIZE=64000000000

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
R_SCRIPT="$SCRIPT_DIR/subset_one.R"

for sp in reindeer human mus acomys; do
  echo ""
  echo "########################################"
  echo "# Processing: $sp"
  echo "########################################"
  Rscript "$R_SCRIPT" "$sp"
  if [ $? -ne 0 ]; then
    echo "!!! Failed on: $sp — stopping"
    exit 1
  fi
done

echo ""
echo "########################################"
echo "# ALL DONE"
echo "########################################"
ls -lh /Users/apple/Downloads/研究生毕业论文/Cross_Species/subset/Fibroblast/
ls -lh /Users/apple/Downloads/研究生毕业论文/Cross_Species/subset/Keratinocyte/

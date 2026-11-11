#!/usr/bin/env bash
# Combined figure (histograms + violin plots + statistics) from the two species. Run after run_timococcum.sh and run_zhukovskyi.sh.
# Usage (from the repository root):  bash run_merged.sh
set -euo pipefail
AM=data/kmer_windows/Am_TA10622_windows_scored_with_At_kmers.tsv; AT=data/kmer_windows/At_TA2804_windows_scored_with_Am_kmers.tsv
mkdir -p results/merged
python code/hist_merged_species.py --timo-pairs results/timococcum/pairs/pairs_for_hist.tsv --timo-ranges results/timococcum/pairs/ranges_for_hist.tsv \
   --zhuk-pairs results/zhukovskyi/pairs/pairs_for_hist.tsv --zhuk-ranges results/zhukovskyi/pairs/ranges_for_hist.tsv \
   --am-windows $AM --at-windows $AT --out results/merged/Fig_HE_breakpoints_vs_kmer_variation
python code/check_reproduction.py merged

#!/usr/bin/env bash
# T. timococcum (88 synthetic plants): read counts -> dosage calls -> breakpoints -> HE pairs -> HE-active ranges -> k-mer test, histogram and violin plot -> per-plant plots.
# Usage (from the repository root):  bash run_timococcum.sh            (uses data/all_plants_timococcum_100Kb.tsv.gz; or give another count file as the first argument)
set -euo pipefail
CNT=${1:-data/all_plants_timococcum_100Kb.tsv.gz}
C=code; O=results/timococcum
AM=data/kmer_windows/Am_TA10622_windows_scored_with_At_kmers.tsv; AT=data/kmer_windows/At_TA2804_windows_scored_with_Am_kmers.tsv
echo "[1/6] dosage calls and breakpoints"; python $C/call_dosage.py --combined "$CNT" --outdir $O/calls --min-tract 3 --plots
mkdir -p $O/calls/per_plant; mv $O/calls/*.segments.tsv $O/calls/*.overview.png $O/calls/per_plant/ 2>/dev/null || true
echo "[2/6] HE classification and pairing"; python $C/he_pipeline.py --breakpoints $O/calls/breakpoints_for_pipeline.tsv --outdir $O/he --shared-min-plants 79 --name-regex '^timo_(?P<lineage>MVGB\d+-\d+T\d+)_(?P<generation>F\d+)-(?P<plant_no>\d+)$'
echo "[3/6] reliability filter"; python $C/filter_reliable_he.py --long $O/he/breakpoints_long.tsv --outdir $O/reliable --min-tract 3 --recurrent-min 20
echo "[4/6] unique exact reciprocal pairs and HE-active ranges"; python $C/make_unique_pairs.py --reliable $O/reliable/breakpoints_reliable.tsv --free-ranges $O/reliable/HE_free_ranges_used.tsv --outdir $O/pairs --am-windows $AM --at-windows $AT
echo "[5/6] k-mer test, histogram, violin plot"; mkdir -p $O/figure
python $C/hist_merged_species.py --timo-pairs $O/pairs/pairs_for_hist.tsv --timo-ranges $O/pairs/ranges_for_hist.tsv --am-windows $AM --at-windows $AT  --out $O/figure/Fig_timococcum_HE_vs_kmer_variation
echo "[6/6] per-plant plots (dosage axis + reliable HE breakpoints)"; python $C/plot_plant_overview.py --combined "$CNT" --calls $O/calls --reliable $O/reliable/breakpoints_reliable.tsv --outdir $O/overview_HE
python $C/check_reproduction.py timococcum

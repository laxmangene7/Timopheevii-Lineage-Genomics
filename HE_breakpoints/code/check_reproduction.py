#!/usr/bin/env python3
"""
check_reproduction.py -- compare freshly generated results with the numbers reported in the manuscript.
Usage (from the repository root):  python code/check_reproduction.py timococcum | zhukovskyi | merged
Counts must match exactly; means and effect sizes are allowed a small tolerance (permutation p-values are random and are not checked).
"""
import sys
import pandas as pd

ok = True
def check(name, got, want, tol=0):
    global ok
    good = abs(got - want) <= tol
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {name}: got {got:.4g}, expected {want} (tolerance {tol})")

EXPECT = {'timococcum': (492, {'HE': 1469, 'Non-HE structural': 2634, 'G-genome': 701}),
          'zhukovskyi': (129, {'HE': 923, 'Non-HE structural': 604, 'G-genome': 129})}
what = sys.argv[1] if len(sys.argv) > 1 else 'merged'
try:
    if what in EXPECT:
        n_unique, cats = EXPECT[what]
        u = pd.read_csv(f'results/{what}/pairs/HE_exact_reciprocal_pairs_unique.tsv', sep='\t')
        r = pd.read_csv(f'results/{what}/reliable/breakpoints_reliable.tsv', sep='\t')
        check(f'{what}: unique exact reciprocal Am/At pairs', len(u), n_unique)
        for c, v in cats.items(): check(f'{what}: reliable {c} breakpoints', int((r.category == c).sum()), v)
    elif what == 'merged':
        s = pd.read_csv('results/merged/Fig_HE_breakpoints_vs_kmer_variation_stats.tsv', sep='\t', index_col=0)
        for g, mean_bp, mean_bg, d in (('Am', 919.2, 1077.3, -0.70), ('At', 953.3, 1087.3, -0.66)):
            col = f'{g}: merged'
            check(f'merged {g}: pairs used', float(s.loc['n_pairs', col]), 616)
            check(f'merged {g}: mean k-mer variation at breakpoints', float(s.loc['mean_pairs', col]), mean_bp, 1.0)
            check(f'merged {g}: mean at non-breakpoint active windows', float(s.loc['mean_nonbreakpoint', col]), mean_bg, 1.0)
            check(f'merged {g}: Cohen d', float(s.loc['cohens_d', col]), d, 0.02)
    else:
        sys.exit('argument must be timococcum, zhukovskyi or merged')
except FileNotFoundError as e:
    print('FAIL  missing file:', e.filename); ok = False
print('\nALL CHECKS PASSED' if ok else '\nSOME CHECKS FAILED')
sys.exit(0 if ok else 1)

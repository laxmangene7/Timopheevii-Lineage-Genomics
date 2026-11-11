#!/usr/bin/env python3
"""
make_unique_pairs.py -- reliable HE calls -> unique exact reciprocal Am-At pairs + HE-active ranges

Input : breakpoints_reliable.tsv and HE_free_ranges_used.tsv written by filter_reliable_he.py
Output: HE_exact_reciprocal_pairs.xlsx / .tsv   Unique_pairs (one row per distinct Am/At site pair, with plant IDs) and Pair_per_plant
        pairs_for_hist.tsv                      input of hist_merged_species.py
        ranges_for_hist.tsv                     HE-free interior per chromosome (input of hist_merged_species.py)
        HE_active_ranges.tsv                    HE-active ranges (outside the HE-free interior), with chromosome lengths

Exact reciprocal pair = one Am breakpoint and one At breakpoint in the same plant, same chromosome group, whose dosage changes are
exactly equal and opposite (e.g. Am 2>3 with At 2>1). Pairs with compound (split) changes or opposite-sign-only matches are not used.
A site pair that occurs in several plants is counted once (the plants are listed).
"""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--reliable', required=True); ap.add_argument('--free-ranges', required=True); ap.add_argument('--outdir', required=True)
ap.add_argument('--am-windows', help='k-mer window table of Am (for chromosome lengths)'); ap.add_argument('--at-windows', help='k-mer window table of At')
a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)

d = pd.read_csv(a.reliable, sep='\t')
h = d[(d.category == 'HE') & d.subgenome.isin(['Am', 'At'])]
rows = []
for (pl, pid), g in h.groupby(['plant_id', 'pair_id']):
    am, at = g[g.subgenome == 'Am'], g[g.subgenome == 'At']
    if len(am) == 1 and len(at) == 1 and am.dosage_change.iloc[0] == -at.dosage_change.iloc[0] and str(am.match_type.iloc[0]).startswith('Exact'):
        rows.append(dict(plant=pl, chromosome=am.chromosome.iloc[0], Am_pos_Mb=am.position_Mb.iloc[0], At_pos_Mb=at.position_Mb.iloc[0],
                         Am_dosage=f'{am.dosage_left.iloc[0]}>{am.dosage_right.iloc[0]}', At_dosage=f'{at.dosage_left.iloc[0]}>{at.dosage_right.iloc[0]}',
                         dosage_change_Am=int(am.dosage_change.iloc[0])))
P = pd.DataFrame(rows)
U = (P.groupby(['chromosome', 'Am_pos_Mb', 'At_pos_Mb'])
       .agg(Am_dosage=('Am_dosage', lambda s: ', '.join(sorted(set(s)))), At_dosage=('At_dosage', lambda s: ', '.join(sorted(set(s)))),
            n_plants=('plant', 'nunique'), plants=('plant', lambda s: '; '.join(sorted(set(s))))).reset_index())
U.insert(0, 'unique_pair_id', ['UP%03d' % (i + 1) for i in range(len(U))])
U.to_csv(os.path.join(a.outdir, 'pairs_for_hist.tsv'), sep='\t', index=False)
U.to_csv(os.path.join(a.outdir, 'HE_exact_reciprocal_pairs_unique.tsv'), sep='\t', index=False)
P.to_csv(os.path.join(a.outdir, 'HE_exact_reciprocal_pairs_per_plant.tsv'), sep='\t', index=False)
with pd.ExcelWriter(os.path.join(a.outdir, 'HE_exact_reciprocal_pairs.xlsx')) as x:
    U.to_excel(x, sheet_name='Unique_pairs', index=False); P.to_excel(x, sheet_name='Pair_per_plant', index=False)

Z = pd.read_csv(a.free_ranges, sep='\t')
Zh = Z.rename(columns={'subgenome': 'genome', 'free_start_Mb': 'HE_free_start_Mb', 'free_end_Mb': 'HE_free_end_Mb'})
Zh.to_csv(os.path.join(a.outdir, 'ranges_for_hist.tsv'), sep='\t', index=False)
if a.am_windows and a.at_windows:
    L = {'Am': pd.read_csv(a.am_windows, sep='\t').groupby('seqname').end.max() / 1e6, 'At': pd.read_csv(a.at_windows, sep='\t').groupby('seqname').end.max() / 1e6}
    out = []
    for r in Z.itertuples():
        ln, s, e = L[r.subgenome][r.chromosome], r.free_start_Mb, r.free_end_Mb
        act = ([f'0-{s:.1f}'] if s > 0 else []) + ([f'{e:.1f}-{ln:.1f}'] if e < ln else [])
        out.append(dict(subgenome=r.subgenome, chromosome=r.chromosome, length_Mb=round(ln, 1), HE_free_start_Mb=s, HE_free_end_Mb=e,
                        active_ranges_Mb=' and '.join(act), pct_active=round(100 * (1 - (e - s) / ln), 1), n_HE_sites=r.n_HE_sites))
    pd.DataFrame(out).sort_values(['subgenome', 'chromosome']).to_csv(os.path.join(a.outdir, 'HE_active_ranges.tsv'), sep='\t', index=False)
print(f'unique exact reciprocal pairs: {len(U)} (plant-level pairs {len(P)}; found in one plant only: {(U.n_plants == 1).sum()})')

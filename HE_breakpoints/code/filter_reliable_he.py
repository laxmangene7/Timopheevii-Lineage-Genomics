#!/usr/bin/env python3
"""
filter_reliable_he.py -- keep only reliable breakpoints / homoeologous exchanges (HE)

Input : breakpoints_long table produced by he_pipeline.py (.tsv or .xlsx)
Output: <outdir>/breakpoints_long_filtered.tsv   every row + one flag column per filter
        <outdir>/breakpoints_reliable.tsv        rows that pass all technical filters
        <outdir>/filter_waterfall.tsv            how many calls each filter removes
        <outdir>/Table_S_reliable.xlsx           same, as sheets (HE / non-HE / G / waterfall / readme)
        <outdir>/segment_length_hist.{svg,png}   evidence for the minimum tract length

Filters (applied in this order; every one is a column, so nothing is silently deleted)
  F_centromeric : A-subgenome breakpoint lies inside the HE-free (centromeric/pericentromeric)
                  interior of its chromosome. Interior = largest breakpoint-free gap computed
                  from the *tract-filtered, paired* HE calls of all plants (or --ranges file).
  F_short_tract : breakpoint bounds a tract shorter than --min-tract Mb (default 1 Mb).
                  At ~0.2x skim-seq a 100-kb bin holds ~150 reads per copy, so a 1-3 bin tract
                  cannot be distinguished from Poisson/GC noise; 1 Mb = 10 bins = ~1,500 reads/copy.
  F_two_res     : (optional, --calls-1mb) no call of the same sign within the 1-Mb-bin run.
                  The 1-Mb run has 10x the reads per bin, so agreement = independent confirmation.
  F_orphan      : HE whose reciprocal partner was removed by any filter above.
  inherited     : (not an error, a different class) parent carries the same site, or the site is
                  shared by >= --recurrent-min plants (assembly / mapping artefact or lineage-inherited).
Final classes:  reliable (passes technical filters) and new_event (reliable and not inherited).
"""
import argparse, os, sys
import numpy as np, pandas as pd


def read_table(p):
    if p.lower().endswith(('.xlsx', '.xls')):
        x = pd.ExcelFile(p)
        sh = 'breakpoints_long' if 'breakpoints_long' in x.sheet_names else x.sheet_names[0]
        return x.parse(sh)
    return pd.read_csv(p, sep='\t')


def tract_flags(df, min_tract_bp):
    """flag both breakpoints bounding a tract shorter than min_tract_bp (per plant/chrom)."""
    df = df.sort_values(['plant_id', 'chrom', 'position']).copy()
    g = df.groupby(['plant_id', 'chrom'])
    df['gap_next'] = g.position.shift(-1) - df.position
    df['gap_prev'] = df.position - g.position.shift(1)
    df['F_short_tract'] = (df.gap_next < min_tract_bp) | (df.gap_prev < min_tract_bp)
    # blip = short tract that returns to the starting dosage (net change 0) -> pure noise
    nl = g.dosage_right.shift(-1)
    df['blip_next'] = (df.gap_next < min_tract_bp) & (nl == df.dosage_left)
    return df


def derive_active(df, ok_mask):
    """largest breakpoint-free interval per (subgenome, chromosome) from surviving paired HE calls."""
    rows = []
    he = df[ok_mask & (df.category == 'HE') & df.subgenome.isin(['Am', 'At'])]
    for (sg, ch), d in df[df.subgenome.isin(['Am', 'At'])].groupby(['subgenome', 'chromosome']):
        L = d.interval_end.max()
        pos = np.sort(he[(he.subgenome == sg) & (he.chromosome == ch)].position.unique())
        pts = np.concatenate([[0], pos, [max(L, pos.max() if len(pos) else 0)]])
        if len(pos) == 0:
            rows.append((sg, ch, np.nan, np.nan, 0)); continue
        gaps = np.diff(pts); i = gaps.argmax()
        rows.append((sg, ch, pts[i] / 1e6, pts[i + 1] / 1e6, len(pos)))
    return pd.DataFrame(rows, columns=['subgenome', 'chromosome', 'free_start_Mb', 'free_end_Mb', 'n_HE_sites'])


def load_ranges(p):
    r = pd.read_csv(p, sep='\t')
    return r.rename(columns={'genome': 'subgenome', 'HE_free_start_Mb': 'free_start_Mb',
                             'HE_free_end_Mb': 'free_end_Mb'})[['subgenome', 'chromosome', 'free_start_Mb', 'free_end_Mb']]


def two_resolution(df, p1mb, pad_bp=100000):
    m = pd.read_csv(p1mb, sep='\t')
    m['sg'] = m.chrom.str.extract(r'_(TA\d+)$')[0].map({'TA2804': 'At', 'TA10622': 'Am'})
    m['chr'] = m.chrom.str.extract(r'^(chr\d[A-Z])')[0]
    m['sign'] = np.sign(m.change)
    idx = {k: d for k, d in m.groupby(['sample', 'sg', 'chr'])}
    ok = np.zeros(len(df), bool)
    for n, (pid, sg, ch, a, s) in enumerate(zip(df.plant_id, df.subgenome, df.chromosome, df.position,
                                                 np.sign(df.dosage_change))):
        d = idx.get((pid, sg, ch))
        if d is not None:
            ok[n] = ((a >= d.interval_start - pad_bp) & (a <= d.interval_end + pad_bp) & (d.sign == s)).any()
    return ok, m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--long', required=True, help='breakpoints_long.tsv/.xlsx from he_pipeline.py')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--min-tract', type=float, default=1.0, help='minimum tract length in Mb (default 1.0)')
    ap.add_argument('--ranges', help='optional HE-free ranges TSV (genome, chromosome, HE_free_start_Mb, HE_free_end_Mb)')
    ap.add_argument('--calls-1mb', help='optional 1-Mb-bin breakpoint calls (same columns as breakpoints.tsv)')
    ap.add_argument('--recurrent-min', type=int, default=10, help='site shared by >= this many plants = recurrent')
    ap.add_argument('--no-active', action='store_true', help='skip the centromeric filter')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    df = read_table(a.long)
    df = df.loc[:, ~df.columns.astype(str).str.startswith('Unnamed')]
    for c in ['own_parents_with_site', 'controls_with_site', 'plants_with_same_site']:
        if c not in df: df[c] = 0
    n0 = len(df)
    T = int(a.min_tract * 1e6)

    # --- F_short_tract
    df = tract_flags(df, T)
    # --- F_two_res (optional)
    if a.calls_1mb:
        isA = df.subgenome.isin(['Am', 'At'])
        ok, _ = two_resolution(df, a.calls_1mb)
        df['F_two_res'] = isA & ~ok          # G rows have no 1-Mb call -> not tested
    else:
        df['F_two_res'] = False
    # --- F_centromeric
    if a.no_active:
        rng = None; df['F_centromeric'] = False
    else:
        if a.ranges:
            rng = load_ranges(a.ranges)
        else:   # iterate once: ranges from HE calls that survive the tract filter and are paired
            rng = derive_active(df, ~df.F_short_tract)
        df = df.merge(rng, on=['subgenome', 'chromosome'], how='left')
        pm = df.position / 1e6
        df['F_centromeric'] = df.subgenome.isin(['Am', 'At']) & (pm > df.free_start_Mb) & (pm < df.free_end_Mb)
        df = df.drop(columns=['free_start_Mb', 'free_end_Mb'])
        rng.to_csv(os.path.join(a.outdir, 'HE_free_ranges_used.tsv'), sep='\t', index=False)

    tech = ['F_centromeric', 'F_short_tract', 'F_two_res']
    df['fail_any'] = df[tech].any(axis=1)
    # --- F_orphan: HE whose partner failed
    df['F_orphan'] = False
    he = df.category.eq('HE') & df.pair_id.notna()
    bad_pairs = set(df.loc[he & df.fail_any, 'pair_id'])
    df.loc[he & ~df.fail_any & df.pair_id.isin(bad_pairs), 'F_orphan'] = True
    df['reliable'] = ~(df.fail_any | df.F_orphan)
    # the surviving partner of a failed HE is not a reliable HE; keep it visible via F_orphan
    df['inherited'] = (df.own_parents_with_site > 0) | (df.plants_with_same_site >= a.recurrent_min)
    df['new_event'] = df.reliable & ~df.inherited

    # reason string
    def reason(r):
        out = [f for f in ['F_centromeric', 'F_short_tract', 'F_two_res', 'F_orphan'] if r[f]]
        return ';'.join(out)
    df['removed_by'] = df.apply(reason, axis=1)

    # --- waterfall (sequential)
    cats = ['HE', 'Non-HE structural', 'G-genome']
    wf = []
    keep = pd.Series(True, index=df.index)
    def snap(name):
        wf.append([name] + [int((keep & (df.category == c)).sum()) for c in cats] + [int(keep.sum())])
    snap('All calls (input)')
    for f, lab in [('F_centromeric', '- inside HE-free (centromeric/pericentromeric) interior'),
                   ('F_short_tract', f'- bounds a tract < {a.min_tract:g} Mb'),
                   ('F_two_res', '- not confirmed in 1-Mb-bin run'),
                   ('F_orphan', '- HE whose reciprocal partner was removed')]:
        keep &= ~df[f]; snap(lab)
    wf.append(['= RELIABLE'] + [int(((df.reliable) & (df.category == c)).sum()) for c in cats] + [int(df.reliable.sum())])
    wf.append(['  of which NEW events (not in parents, < %d plants)' % a.recurrent_min] +
              [int(((df.new_event) & (df.category == c)).sum()) for c in cats] + [int(df.new_event.sum())])
    wf = pd.DataFrame(wf, columns=['step'] + cats + ['total'])
    wf.to_csv(os.path.join(a.outdir, 'filter_waterfall.tsv'), sep='\t', index=False)

    per = df.groupby('plant_id').agg(all_calls=('position', 'size'), reliable=('reliable', 'sum'),
                                     new_events=('new_event', 'sum')).reset_index()
    perHE = df[df.category == 'HE'].groupby('plant_id').agg(HE_all=('position', 'size'), HE_reliable=('reliable', 'sum'),
                                                          HE_new=('new_event', 'sum')).reset_index()
    per = per.merge(perHE, on='plant_id', how='left').fillna(0)
    per.to_csv(os.path.join(a.outdir, 'per_plant_counts.tsv'), sep='\t', index=False)

    df = df.drop(columns=['gap_next', 'gap_prev', 'blip_next'], errors='ignore')
    df.to_csv(os.path.join(a.outdir, 'breakpoints_long_filtered.tsv'), sep='\t', index=False)
    rel = df[df.reliable]
    rel.to_csv(os.path.join(a.outdir, 'breakpoints_reliable.tsv'), sep='\t', index=False)

    readme = pd.DataFrame({'Item': ['min tract', 'recurrent', 'centromeric', 'two-resolution', 'reliable', 'new_event'],
                           'Meaning': [f'{a.min_tract:g} Mb', f'>= {a.recurrent_min} plants share the site',
                                       'largest breakpoint-free interval per chromosome (HE-free interior)',
                                       'confirmed by 1-Mb-bin run' if a.calls_1mb else 'not used',
                                       'passes all technical filters', 'reliable and not inherited / recurrent']})
    with pd.ExcelWriter(os.path.join(a.outdir, 'Table_S_reliable.xlsx')) as xw:
        rel[rel.category == 'HE'].sort_values(['pair_id', 'subgenome']).to_excel(xw, sheet_name='HE_reliable', index=False)
        rel[rel.category == 'Non-HE structural'].to_excel(xw, sheet_name='NonHE_reliable', index=False)
        rel[rel.category == 'G-genome'].to_excel(xw, sheet_name='G_reliable', index=False)
        wf.to_excel(xw, sheet_name='Waterfall', index=False)
        per.to_excel(xw, sheet_name='Per_plant', index=False)
        readme.to_excel(xw, sheet_name='Readme', index=False)

    # --- evidence figure
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        d = df.dropna(subset=['gap_next']) if 'gap_next' in df else None
    except Exception:
        d = None
    try:
        t = tract_flags(read_table(a.long).loc[:, lambda x: ~x.columns.astype(str).str.startswith('Unnamed')], T)
        g = t.gap_next.dropna() / 1e6
        fig, ax = plt.subplots(figsize=(6.2, 3.6))
        ax.hist(g, bins=np.logspace(-1, 3, 60), color='#4C78A8')
        ax.set_xscale('log'); ax.axvline(a.min_tract, color='#D62728', ls='--')
        ax.set_xlabel('Distance to next breakpoint in the same plant / chromosome (Mb)')
        ax.set_ylabel('Number of tracts')
        ax.set_title(f'Short tracts (< {a.min_tract:g} Mb, dashed) are the noise tail', fontsize=10)
        fig.tight_layout()
        for e in ('svg', 'png'): fig.savefig(os.path.join(a.outdir, f'segment_length_hist.{e}'), dpi=200)
    except Exception as e:
        print('figure skipped:', e, file=sys.stderr)

    print(wf.to_string(index=False))


if __name__ == '__main__':
    main()

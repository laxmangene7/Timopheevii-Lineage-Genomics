#!/usr/bin/env python3
"""
plot_plant_overview.py -- per-plant read-depth overview with dosage axis and reliable HE breakpoints

Grey dots = reads per 100-kb bin; red = called dosage; right-hand axis labels the copy number (0x ... 5x, 1x = per-plant unit);
dotted vertical lines = reliable HE breakpoints (exact reciprocal Am/At pairs, from breakpoints_reliable.tsv), blue = dosage gain, orange = dosage loss.
Usage: plot_plant_overview.py --combined counts.tsv.gz --calls calls_dir --reliable breakpoints_reliable.tsv --outdir figs [--plants P1 P2 ...] [--all-he]
"""
import argparse, os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
BIN = 100000

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--combined', required=True); ap.add_argument('--calls', required=True); ap.add_argument('--reliable', required=True)
ap.add_argument('--outdir', required=True); ap.add_argument('--plants', nargs='*', help='sample names (default: all)')
ap.add_argument('--strip-prefix', default=r'^[0-9a-f]{8}-'); ap.add_argument('--dpi', type=int, default=150)
ap.add_argument('--all-he', action='store_true', help='mark every reliable HE breakpoint, not only exact reciprocal pairs')
a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)

D = pd.read_csv(a.combined, sep='\t', header=None, names=['sample', 'n', 'chrom', 'pos'])
D['sample'] = D['sample'].str.replace(a.strip_prefix, '', regex=True) if a.strip_prefix else D['sample']
QC = pd.read_csv(os.path.join(a.calls, 'qc_per_plant.tsv'), sep='\t').set_index('sample')
R = pd.read_csv(a.reliable, sep='\t'); R = R[R.category == 'HE']
if not a.all_he: R = R[R.match_type.astype(str).str.startswith('Exact')]
plants = a.plants or sorted(QC.index)
rank = lambda c: (int(re.match(r'chr(\d)', c).group(1)), 0 if re.match(r'chr\dA_TA2804', c) else 1 if 'TA10622' in c else 2)

for s in plants:
    d = D[D['sample'] == s]; u = QC.loc[s, 'unit']
    seg = pd.read_csv(os.path.join(a.calls, 'per_plant', f'{s}.segments.tsv'), sep='\t'); hb = R[R.plant_id == s]
    order = sorted(d.chrom.unique(), key=rank)
    fig, axs = plt.subplots(7, 3, figsize=(16, 16), squeeze=False)
    for i, c in enumerate(order):
        ax = axs[i // 3][i % 3]; g = d[d.chrom == c]
        ax.scatter(g.pos / 1e6, g.n, s=.3, c='#888', alpha=.5, rasterized=True)
        sg = seg[seg.chrom == c]
        for r in sg.itertuples(): ax.plot([r.start / 1e6, r.end / 1e6], [r.dosage * u] * 2, c='#d62728', lw=1.2, solid_capstyle='butt')
        for k in range(1, 6): ax.axhline(k * u, c='#ddd', lw=.4, zorder=0)
        for r in hb[hb.chrom == c].itertuples():
            ax.axvline(r.position / 1e6, c='#1f77b4' if r.dosage_change > 0 else '#ff7f0e', ls=':', lw=1.1, zorder=3)
        ax.set_title(c, fontsize=8); ax.set_ylim(0, 6 * u); ax.tick_params(labelsize=6)
        ax.set_xlim(0, g.pos.max() / 1e6 + .1)
        ax2 = ax.twinx(); ax2.set_ylim(0, 6 * u); ax2.set_yticks([k * u for k in range(0, 6)]); ax2.set_yticklabels([f'{k}x' for k in range(0, 6)], fontsize=6.5, color='#d62728')
        ax2.tick_params(length=2, colors='#d62728')
    for j in range(len(order), 21): axs[j // 3][j % 3].axis('off')
    fig.suptitle(f'{s}   1 copy = {u:.0f} reads / 100 kb   (red = called dosage; right axis = copy number; dotted = reliable HE breakpoint)', fontsize=10)
    fig.legend(handles=[Line2D([], [], c='#1f77b4', ls=':', label='HE breakpoint, dosage gain'), Line2D([], [], c='#ff7f0e', ls=':', label='HE breakpoint, dosage loss')], loc='upper right', fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, .98)); fig.savefig(os.path.join(a.outdir, f'{s}.overview_HE.png'), dpi=a.dpi); plt.close(fig)
    print(s, len(hb), 'HE breakpoints marked')

#!/usr/bin/env python3
"""
call_dosage.py -- calibrated, segmentation-based dosage calls from skim-seq read counts in 100-kb bins.

Why this replaces per-bin thresholding
  * The copy unit is estimated per plant from the genome-wide 2x mode (G + A pooled, so a few changed
    chromosomes do not move it), then refined from the bins the model itself calls 2x. Nothing assumes
    that a particular chromosome is 2x, because G chromosomes can carry non-homoeologous changes too.
  * Bins that are absent from `uniq -c` output are true zeros and are filled in before anything is computed.
  * A hidden Markov model (states 0..MAXCN copies, negative-binomial emissions with an outlier component,
    sticky transitions) calls whole tracts. Isolated bins cannot flip the state.
  * State 0 has its own background (cross-mapping between At and Am), estimated from the data.
  * Tracts shorter than --min-tract Mb are merged into the neighbour with the closest copy number.
  * Optional per-bin cohort correction (--cohort-norm) removes shared mappability / GC bias; it is applied
    only to bins whose cohort level is stable and is limited to +-30 %.

Input : one file per plant, lines `count<TAB>chrom<TAB>start_bp`  (name: <sample>_100Kb.txt)
Output: <outdir>/<sample>.segments.tsv     one row per called tract
        <outdir>/<sample>.breakpoints.tsv  sample, chrom, position, interval_start, interval_end, dosage_left, dosage_right, change, flags
        <outdir>/qc_per_plant.tsv          unit, background, dispersion, fraction at 2x, breakpoint counts
        <outdir>/<sample>.overview.png     one panel per chromosome (--plots)
"""
import argparse, os, re, sys, glob
import numpy as np, pandas as pd
from scipy.special import gammaln

BIN = 100_000


def read_counts(path):
    d = pd.read_csv(path, sep=r'\s+', header=None, names=['n', 'chrom', 'pos'], dtype={'n': float, 'pos': np.int64})
    return d


def chrom_lengths(files, lengths_tsv=None):
    if lengths_tsv:
        t = pd.read_csv(lengths_tsv, sep='\t', header=None, names=['chrom', 'length'])
        return {c: int(np.ceil(l / BIN)) for c, l in zip(t.chrom, t.length)}
    L = {}
    for f in files:
        d = read_counts(f)
        for c, m in d.groupby('chrom').pos.max().items():
            L[c] = max(L.get(c, 0), int(m // BIN) + 1)
    return L


def to_matrix(d, L):
    """dict chrom -> float array of reads per bin (absent bins = 0)."""
    out = {c: np.zeros(n) for c, n in L.items()}
    for c, g in d.groupby('chrom'):
        if c in out:
            idx = (g.pos.values // BIN).astype(int)
            ok = idx < len(out[c])
            out[c][idx[ok]] = g.n.values[ok]
    return out


# ---------------------------------------------------------------- emissions
def nb_logpmf(x, mu, r):
    mu = np.maximum(mu, 1e-3)
    return (gammaln(x + r) - gammaln(r) - gammaln(x + 1) + r * np.log(r / (r + mu)) + x * np.log(mu / (r + mu)))


def emission(x, mus, r, eps=0.01, xmax=None, mask=None):
    """log P(x | state) with a broad outlier component; shape (len(x), n_states)."""
    xmax = xmax or max(x.max(), 1.0)
    lo = np.log(eps / (xmax + 1.0))
    E = np.stack([nb_logpmf(x, m, r) for m in mus], axis=1)
    E = np.logaddexp(np.log1p(-eps) + E, lo)
    if mask is not None:
        E[mask] = 0.0                      # masked bins carry no information about the state
    return E


def viterbi(E, p_switch):
    n, k = E.shape
    ls, lw = np.log1p(-p_switch), np.log(p_switch / (k - 1))
    T = np.full((k, k), lw); np.fill_diagonal(T, ls)
    dp = np.empty((n, k)); bp = np.zeros((n, k), dtype=np.int8)
    dp[0] = E[0] - np.log(k)
    for i in range(1, n):
        s = dp[i - 1][:, None] + T
        bp[i] = s.argmax(0); dp[i] = s.max(0) + E[i]
    path = np.empty(n, dtype=np.int8); path[-1] = dp[-1].argmax()
    for i in range(n - 1, 0, -1):
        path[i - 1] = bp[i, path[i]]
    return path


# ---------------------------------------------------------------- calibration
def mode_2x(vals):
    v = vals[vals > 0]
    med = np.median(v)
    h, e = np.histogram(v, bins=np.linspace(0.3 * med, 2.2 * med, 160))
    h = np.convolve(h, np.ones(5) / 5, mode='same')
    return 0.5 * (e[h.argmax()] + e[h.argmax() + 1])


def estimate_dispersion(x2, mu2):
    """method of moments for var = mu + mu^2 / r on bins called 2x."""
    v = np.var(x2) - mu2
    return float(np.clip(mu2 ** 2 / v, 3, 500)) if v > 0 else 500.0


def calibrate_standard(counts, std=('chr3G', 'chr5G', 'chr7G'), tol=0.15):
    """copy unit from chromosomes assumed 2x (default 3G, 5G, 7G). If one disagrees with the others by > tol,
    it is dropped. Returns unit, dispersion r, flag text."""
    keys = [c for c in counts if c.split('_')[0] in std]
    med = {c: np.median(counts[c]) for c in keys}
    m = np.median(list(med.values()))
    good = [c for c in keys if abs(med[c] - m) / m <= tol]
    flag = ''
    if len(good) < len(keys):
        flag = 'dropped:' + ','.join(c.split('_')[0] for c in keys if c not in good)
    if len(good) < 2:
        good, flag = keys, 'standard chromosomes disagree'
    x = np.concatenate([counts[c] for c in good])
    mu = np.median(x)
    lo, hi = np.percentile(x, [5, 95]); xt = x[(x >= lo) & (x <= hi)]
    # trimmed dispersion: scale the trimmed variance back up for a normal-like core
    v = np.var(xt) * (1 / 0.55)          # variance of the central 90% of a normal is ~0.55 x full variance
    r = float(np.clip(mu ** 2 / max(v - mu, 1e-6), 5, 500))
    return mu / 2.0, r, flag


def cohort_mask(counts, cal, low=0.75, high=1.35, frac=0.5, close_bins=3, max_run_bins=80):
    """bins whose depth is far from 2x in >= frac of ALL plants (mappability, assembly gaps, differences fixed
    in the whole cohort) are uninformative for a segregating exchange and are masked."""
    out, rows = {}, []
    for c in next(iter(counts.values())):
        M = np.stack([counts[s][c] / (2 * cal[s][0]) for s in counts])
        m = ((M < low).mean(0) >= frac) | ((M > high).mean(0) >= frac)
        from scipy.ndimage import binary_closing, label
        m = binary_closing(m, structure=np.ones(2 * close_bins + 1), border_value=0) | m
        lab, n = label(m)                       # a masked run longer than max_run_bins is a real shared state (e.g. a whole lost arm), not an artefact dip
        for i in range(1, n + 1):
            if (lab == i).sum() > max_run_bins: m[lab == i] = False
        out[c] = m
        ch = np.flatnonzero(np.diff(np.r_[0, m.astype(int), 0]))
        for a_, b_ in zip(ch[::2], ch[1::2]): rows.append((c, a_ * BIN, b_ * BIN, (b_ - a_) * BIN / 1e6, float(np.median(M[:, a_:b_]))))
    return out, pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'length_Mb', 'median_level_vs_2x'])


def call_plant(counts, maxcn=6, p_switch=1e-3, min_tract_bins=10, iters=2, u_fixed=None, r_fixed=None, mask=None):
    mask = mask or {}
    chroms = list(counts)
    allv = np.concatenate([counts[c] for c in chroms])
    u = u_fixed if u_fixed else mode_2x(allv) / 2.0
    bg = 0.08 * u
    r = r_fixed if r_fixed else 30.0
    for _ in range(iters):
        mus = np.array([bg] + [k * u for k in range(1, maxcn + 1)])
        paths = {c: viterbi(emission(counts[c], mus, r, xmax=allv.max(), mask=mask.get(c)), p_switch) for c in chroms}
        P = np.concatenate([paths[c] for c in chroms])
        if not u_fixed and (P == 2).sum() > 1000:
            x2 = allv[P == 2]; u = np.median(x2) / 2.0; r = estimate_dispersion(x2, np.mean(x2))
        if (P == 0).sum() > 200:
            bg = float(np.clip(np.mean(allv[P == 0]), 0.01 * u, 0.3 * u))
    mus = np.array([bg] + [k * u for k in range(1, maxcn + 1)])
    paths = {c: merge_short(viterbi(emission(counts[c], mus, r, xmax=allv.max(), mask=mask.get(c)), p_switch), counts[c], mus, min_tract_bins) for c in chroms}
    qc = dict(unit=u, background=bg, dispersion_r=r, frac_bins_2x=float((np.concatenate(list(paths.values())) == 2).mean()))
    return paths, mus, qc


def merge_short(path, x, mus, min_bins):
    """merge tracts shorter than min_bins into the neighbour with the closest copy number (ends kept if >=3 bins)."""
    path = path.copy()
    while True:
        ch = np.flatnonzero(np.diff(path)) + 1
        st = np.r_[0, ch]; en = np.r_[ch, len(path)]; ln = en - st
        cand = [i for i in range(len(st)) if ln[i] < (3 if i in (0, len(st) - 1) else min_bins)]
        if not cand or len(st) == 1:
            return path
        i = min(cand, key=lambda j: ln[j])
        nb = [j for j in (i - 1, i + 1) if 0 <= j < len(st)]
        # choose neighbour whose level fits the tract mean best (fall back to the longer one)
        m = x[st[i]:en[i]].mean()
        j = min(nb, key=lambda j: (abs(mus[path[st[j]]] - m), -ln[j]))
        path[st[i]:en[i]] = path[st[j]]


def segments_from_path(sample, chrom, path, x, u, mus):
    ch = np.flatnonzero(np.diff(path)) + 1
    st = np.r_[0, ch]; en = np.r_[ch, len(path)]
    return [dict(sample=sample, chrom=chrom, start=int(s * BIN), end=int(e * BIN), n_bins=int(e - s), dosage=int(path[s]),
                 mean_reads=float(x[s:e].mean()), mean_copies=float(x[s:e].mean() / u)) for s, e in zip(st, en)]


def refine(counts_c, path, mus, r, ci_drop=1.92, maxshift=8, mask=None):
    """For every state change, scan the boundary +-maxshift bins and return (best_bin, lo_bin, hi_bin):
    the boundary position maximising the likelihood and the range within ci_drop log-lik units (~95 % CI)."""
    ch = np.flatnonzero(np.diff(path)) + 1
    out = {}
    for k, b in enumerate(ch):
        sa, sb = path[b - 1], path[b]
        left = ch[k - 1] if k > 0 else 0
        right = ch[k + 1] if k + 1 < len(ch) else len(path)
        lo_s, hi_s = max(left + 1, b - maxshift), min(right - 1, b + maxshift)
        cand = np.arange(lo_s, hi_s + 1)
        if len(cand) == 0:
            out[b] = (b, b, b); continue
        ea = nb_logpmf(counts_c[lo_s - 1:hi_s + 1], mus[sa], r); eb = nb_logpmf(counts_c[lo_s - 1:hi_s + 1], mus[sb], r)
        if mask is not None:
            mk = mask[lo_s - 1:hi_s + 1]; ea = np.where(mk, 0, ea); eb = np.where(mk, 0, eb)
        ll = np.array([ea[:c - (lo_s - 1)].sum() + eb[c - (lo_s - 1):].sum() for c in cand])
        ok = cand[ll >= ll.max() - ci_drop]
        out[b] = (int(cand[ll.argmax()]), int(ok.min()), int(ok.max()))
    return out


def breakpoints(seg, chrom_len_bp, end_mask_bp, ci=None):
    rows = []
    for (s, c), g in seg.groupby(['sample', 'chrom'], sort=False):
        g = g.sort_values('start')
        for a, b in zip(g.iloc[:-1].itertuples(), g.iloc[1:].itertuples()):
            L = chrom_len_bp[c]
            lo, hi = (ci or {}).get((c, b.start), (b.start - BIN, b.start))
            rows.append(dict(sample=s, chrom=c, position=b.start, interval_start=lo, interval_end=hi + BIN,
                             dosage_left=a.dosage, dosage_right=b.dosage, change=b.dosage - a.dosage,
                             left_tract_Mb=a.n_bins * BIN / 1e6, right_tract_Mb=b.n_bins * BIN / 1e6,
                             near_chrom_end=bool(b.start < end_mask_bp or L - b.start < end_mask_bp)))
    return pd.DataFrame(rows)


def cohort_factors(all_counts, units, maxrel=0.3):
    """per-bin multiplicative correction = cohort median of (reads / (2*unit)), clipped; only used for >=5 plants."""
    out = {}
    for c in next(iter(all_counts.values())):
        M = np.stack([all_counts[s][c] / (2 * units[s]) for s in all_counts])
        f = np.median(M, axis=0)
        out[c] = np.clip(f, 1 - maxrel, 1 + maxrel)
    return out


def overview(sample, counts, path, u, outpng, mus):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    rank = lambda c: (int(re.match(r'chr(\d)', c).group(1)), 0 if re.match(r'chr\dA_TA2804', c) else 1 if 'TA10622' in c else 2)
    order = sorted(counts, key=rank)
    fig, axs = plt.subplots(7, 3, figsize=(16, 16), squeeze=False)
    cols = ['#9e9e9e', '#9ecae1', '#4292c6', '#2171b5', '#08306b', '#6a51a3', '#54278f']
    for i, c in enumerate(order):
        ax = axs[i // 3][i % 3]; x = counts[c]; pos = np.arange(len(x)) * BIN / 1e6
        ax.scatter(pos, x, s=.3, c='#888', alpha=.5)
        ax.step(pos, mus[path[c]], where='post', c='#d62728', lw=1.1)
        for k in range(1, 6): ax.axhline(k * u, c='#ddd', lw=.4)
        ax.set_title(c, fontsize=8); ax.set_ylim(0, 6 * u); ax.tick_params(labelsize=6)
    fig.suptitle(f'{sample}   1 copy = {u:.0f} reads / 100 kb   (red = called dosage, grey lines = 1..5 copies)', fontsize=10)
    fig.tight_layout(); fig.savefig(outpng, dpi=90); plt.close(fig)


def read_combined(path):
    d = pd.read_csv(path, sep='\t', header=None, names=['sample', 'n', 'chrom', 'pos'], dtype={'n': float, 'pos': np.int64})
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('--counts', nargs='+', help='per-plant *_100Kb.txt files (or globs)')
    src.add_argument('--combined', help='one TSV(.gz): sample, count, chrom, start_bp  (all plants)')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--strip-suffix', default='_100Kb.txt')
    ap.add_argument('--strip-prefix', default=r'^[0-9a-f]{8}-', help='regex removed from the start of file names; "" to keep')
    ap.add_argument('--lengths', help='optional TSV chrom<TAB>length_bp (otherwise max position over all plants)')
    ap.add_argument('--min-tract', type=float, default=3.0, help='minimum tract length in Mb (default 3.0)')
    ap.add_argument('--p-switch', type=float, default=1e-8)
    ap.add_argument('--max-cn', type=int, default=6)
    ap.add_argument('--standard', default='chr3G,chr5G,chr7G', help='chromosomes assumed 2x for calibration ("" = genome-wide 2x mode)')
    ap.add_argument('--end-mask', type=float, default=3.0, help='flag breakpoints within this many Mb of a chromosome end')
    ap.add_argument('--cohort-mask', type=float, default=0.5, help='mask bins that are off 2x in >= this fraction of plants (0 = off)')
    ap.add_argument('--mask-max-run', type=float, default=8.0, help='masked runs longer than this (Mb) are real shared states and are NOT masked')
    ap.add_argument('--plots', action='store_true')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.combined:
        D = read_combined(a.combined)
        L = {}
        for c, m in D.groupby('chrom').pos.max().items(): L[c] = int(m // BIN) + 1
        if a.lengths: L = chrom_lengths([], a.lengths)
        counts = {s: to_matrix(g[['n', 'chrom', 'pos']], L) for s, g in D.groupby('sample', sort=True)}
        del D
    else:
        files = sorted(sum([glob.glob(p) for p in a.counts], []))
        L = chrom_lengths(files, a.lengths)
        def sname(f):
            b = os.path.basename(f); b = b[:-len(a.strip_suffix)] if b.endswith(a.strip_suffix) else b
            return re.sub(a.strip_prefix, '', b) if a.strip_prefix else b
        counts = {sname(f): to_matrix(read_counts(f), L) for f in files}
    std = tuple(x for x in a.standard.split(',') if x)
    min_bins = int(round(a.min_tract * 1e6 / BIN))
    bps, qcs = [], []
    cal = {s: calibrate_standard(cc, std) if std else (mode_2x(np.concatenate(list(cc.values()))) / 2.0, 30.0, '') for s, cc in counts.items()}
    mask = {}
    if a.cohort_mask > 0:
        mask, mt = cohort_mask(counts, cal, frac=a.cohort_mask, max_run_bins=int(a.mask_max_run * 1e6 / BIN))
        mt.to_csv(os.path.join(a.outdir, 'cohort_masked_regions.tsv'), sep='\t', index=False)
        print(f"cohort mask: {sum(int(m.sum()) for m in mask.values())} of {sum(len(m) for m in mask.values())} bins masked in {len(mt)} regions", flush=True)
    for s, cc in counts.items():
        u0, r0, flag = cal[s]
        paths, mus, qc = call_plant(cc, a.max_cn, a.p_switch, min_bins, u_fixed=u0, r_fixed=r0, mask=mask)
        u, r = qc['unit'], qc['dispersion_r']
        ci = {}
        for k in cc:                                     # refine every boundary to its maximum-likelihood bin + 95 % interval
            p = paths[k]; ref = refine(cc[k], p, mus, r, mask=mask.get(k))
            for b, (best, lo, hi) in sorted(ref.items()):
                if best < b: p[best:b] = p[b]
                elif best > b: p[b:best] = p[b - 1]
                ci[(k, best * BIN)] = (lo * BIN, hi * BIN)
        sg = pd.DataFrame([row for k in cc for row in segments_from_path(s, k, paths[k], cc[k], u, mus)])
        bp = breakpoints(sg, {k: len(v) * BIN for k, v in cc.items()}, a.end_mask * 1e6, ci)
        sg.to_csv(os.path.join(a.outdir, f'{s}.segments.tsv'), sep='\t', index=False)
        isG = bp.chrom.str.contains(r'^chr\dG_')
        # level check: A chromosomes in state 2 should have the same reads/copy as the standard chromosomes
        A2 = np.concatenate([cc[k][paths[k] == 2] for k in cc if 'G_' not in k]); G2 = np.concatenate([cc[k][paths[k] == 2] for k in cc if 'G_' in k])
        qc.update(sample=s, calibration_flag=flag, A_over_G_level=float(np.median(A2) / np.median(G2)) if len(A2) and len(G2) else np.nan,
                  n_breakpoints=len(bp), n_A=int((~isG).sum()), n_G=int(isG.sum()), n_A_not_near_end=int(((~isG) & ~bp.near_chrom_end).sum()))
        qcs.append(qc); bps.append(bp)
        if a.plots: overview(s, cc, paths, u, os.path.join(a.outdir, f'{s}.overview.png'), mus)
        print(f"{s}: unit {u:.1f}, bg {qc['background']:.1f}, r {r:.0f}, A/G level {qc['A_over_G_level']:.2f} {flag} | breakpoints {len(bp)} (A {qc['n_A']}, G {qc['n_G']})", flush=True)
    pd.DataFrame(qcs)[['sample', 'unit', 'background', 'dispersion_r', 'A_over_G_level', 'calibration_flag', 'frac_bins_2x', 'n_breakpoints', 'n_A', 'n_G', 'n_A_not_near_end']].to_csv(os.path.join(a.outdir, 'qc_per_plant.tsv'), sep='\t', index=False)
    pd.concat(bps).to_csv(os.path.join(a.outdir, 'breakpoints_all_plants.tsv'), sep='\t', index=False)
    pd.concat(bps)[['sample', 'chrom', 'position', 'interval_start', 'interval_end', 'dosage_left', 'dosage_right', 'change']].to_csv(os.path.join(a.outdir, 'breakpoints_for_pipeline.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()

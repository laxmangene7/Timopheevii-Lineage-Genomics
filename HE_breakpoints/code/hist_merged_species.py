#!/usr/bin/env python3
"""Merge the unique exact reciprocal Am-At HE pairs of T. timococcum and T. zhukovskyi and compare their breakpoint windows with
non-breakpoint windows of the HE-active range, on the same k-mer variation tracks. Histogram + violin plots (SVG/PDF/PNG) + statistics.
Active window = outside the HE-free interior of EITHER species (= windows where HE occurs in at least one species).
Unit = one unique pair (Am site, At site) per species; a pair is used only if both its Am and At windows have k-mer data and are active."""
import argparse, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"font.family": ["TeX Gyre Heros", "Liberation Sans", "DejaVu Sans"],   # TeX Gyre Heros is metric-identical to Helvetica
                      "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 8, "axes.linewidth": .6})
GREY, TIMO, ZHUK, COMB = "#808080", "#d95f02", "#1b9e77", "#219ebc"; rng = np.random.default_rng(11)

def load_windows(path, rangesT, rangesZ, g):
    w = pd.read_csv(path, sep="\t"); w = w[(w.total_kmers >= 99900) & (w.observed_kmers > 0)].reset_index(drop=True); w["mid"] = (w.start + w.end) / 2e6
    for tag, Z in (("T", rangesT), ("Z", rangesZ)):
        z = Z[Z.genome == g].set_index("chromosome")
        w["act" + tag] = (w.start / 1e6 < w.seqname.map(z.HE_free_start_Mb)) | (w.end / 1e6 > w.seqname.map(z.HE_free_end_Mb))
    w["active"] = w.actT | w.actZ; return w

def locate(w, ch, mb):
    d = w[w.seqname == ch]; i = np.searchsorted(d.start.values, mb * 1e6, side="right") - 1
    return d.index[i] if 0 <= i < len(d) and mb * 1e6 < d.end.values[i] and d.active.values[i] else None

def tests(w, idx, local_mb=20, n_perm=10000):
    s = set(idx); N = w[w.active & ~w.index.isin(s)]; y = N.variations.values; x = w.loc[idx, "variations"].values
    chs = w.loc[idx, "seqname"].values; grp = {c: N[N.seqname == c].variations.values for c in np.unique(chs)}
    obs = x.mean() - y.mean(); sims = np.array([np.concatenate([rng.choice(grp[c], (chs == c).sum(), replace=False) for c in grp]).mean() - y.mean() for _ in range(n_perm)])
    pooled = np.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / (len(x) + len(y) - 2)); diffs = []
    for i in idx:
        r = w.loc[i]; l = N[(N.seqname == r.seqname) & ((N.mid - r.mid).abs() <= local_mb)]
        if len(l) >= 20: diffs.append(r.variations - l.variations.mean())
    diffs = np.array(diffs); boots = [rng.choice(diffs, len(diffs)).mean() for _ in range(5000)]
    return dict(n_pairs=len(x), n_distinct_windows=len(s), n_nonbreakpoint_active_windows=len(y), mean_pairs=x.mean(), mean_nonbreakpoint=y.mean(), median_pairs=np.median(x), median_nonbreakpoint=np.median(y),
                mean_difference=obs, cohens_d=obs / pooled, mannwhitney_p=stats.mannwhitneyu(x, y).pvalue, permutation_p_by_chromosome=(1 + (np.abs(sims) >= abs(obs)).sum()) / (1 + n_perm),
                local_mean_difference=diffs.mean(), local_diff_95CI=f"{np.percentile(boots,2.5):.0f} to {np.percentile(boots,97.5):.0f}", local_wilcoxon_p=stats.wilcoxon(diffs).pvalue,
                frac_pairs_below_local_mean=(diffs < 0).mean(), sign_test_p=stats.binomtest(int((diffs < 0).sum()), len(diffs), .5).pvalue), N

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("timo-pairs", "timo-ranges", "am-windows", "at-windows"): ap.add_argument("--" + k, required=True)
    ap.add_argument("--first-name", default="T. timococcum", help="species name of the --timo-* files (use \"T. zhukovskyi\" to analyse that species alone)")
    for k in ("zhuk-pairs", "zhuk-ranges"): ap.add_argument("--" + k, help="optional: omit both to analyse T. timococcum alone")
    ap.add_argument("--out", default="hist_merged"); a = ap.parse_args()
    RT = pd.read_csv(a.timo_ranges, sep="\t"); RZ = pd.read_csv(a.zhuk_ranges, sep="\t") if a.zhuk_ranges else RT
    W = {"Am": load_windows(a.am_windows, RT, RZ, "Am"), "At": load_windows(a.at_windows, RT, RZ, "At")}
    P = pd.concat([pd.read_csv(a.timo_pairs, sep="\t").assign(species=a.first_name)] + ([pd.read_csv(a.zhuk_pairs, sep="\t").assign(species="T. zhukovskyi")] if a.zhuk_pairs else []), ignore_index=True)
    SPECIES = list(P.species.unique()); LAB = {"T. timococcum": "T. timococcum", "T. zhukovskyi": "T. zhukovskyi"}
    P["iAm"] = [locate(W["Am"], c, m) for c, m in zip(P.chromosome, P.Am_pos_Mb)]; P["iAt"] = [locate(W["At"], c, m) for c, m in zip(P.chromosome, P.At_pos_Mb)]
    n0 = P.groupby("species").size().to_dict(); P = P.dropna(subset=["iAm", "iAt"]).copy(); P[["iAm", "iAt"]] = P[["iAm", "iAt"]].astype(int)
    shared = P.duplicated(["chromosome", "Am_pos_Mb", "At_pos_Mb"], keep=False).sum() // 2
    print("pairs in:", n0, "| used:", P.groupby("species").size().to_dict(), "| identical site pairs present in both species:", shared)
    P.drop(columns=["iAm", "iAt", "dosage_Am", "dosage_At"], errors="ignore").to_csv(a.out + "_pairs_used.tsv", sep="\t", index=False)
    rows, allN, VT = {}, {}, []
    fig = plt.figure(figsize=(9.6, 8.4)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=.62, wspace=.7)
    for j, (g, tag) in enumerate((("Am", "a"), ("At", "b"))):
        w = W[g]; col = "i" + g; idx_all = P[col].values
        s, N = tests(w, idx_all); rows[(g, "merged")] = s; allN[g] = N; y = N.variations.values; x = w.loc[idx_all, "variations"].values
        for sp in P.species.unique():
            rows[(g, sp)], _ = tests(w, P.loc[P.species == sp, col].values)
        top = np.quantile(w[w.active].variations, .999); edges = np.linspace(0, top, 46)
        ax = fig.add_subplot(gs[0, j]); axr = ax.twinx()                                                                     # non-breakpoint windows on the RIGHT axis, drawn behind
        hg, _, _ = axr.hist(y, bins=edges, density=True, color=GREY, alpha=.65)
        hb, _, _ = ax.hist(x, bins=edges, density=True, color=COMB, alpha=.65)             # HE breakpoints on the LEFT axis, drawn in front
        ax.set_zorder(axr.get_zorder() + 1); ax.patch.set_visible(False)
        ax.set_ylim(0, hb.max() * 1.08); axr.set_ylim(0, hg.max() * 1.08)                  # each axis scaled to the same peak height
        ax.tick_params(axis="y", colors="#0b7285"); axr.tick_params(axis="y", colors="#555555"); ax.spines["left"].set_color("#0b7285")
        axr.spines["top"].set_visible(False); axr.spines["right"].set_color("#555555"); axr.set_ylabel("Density, non-breakpoint windows", color="#555555")
        ax.set_title(f"{tag}  {g}", loc="left", fontsize=9, fontweight="bold", pad=26)
        ax.text(0, 1.03, f"n = {s['n_pairs']} unique pairs{" (both species)" if len(SPECIES) > 1 else " (" + SPECIES[0] + ")"} vs {s['n_nonbreakpoint_active_windows']:,} non-breakpoint windows\nMann–Whitney p = {s['mannwhitney_p']:.1e};  permutation p {('< 1e-4' if s['permutation_p_by_chromosome'] <= 1.01e-4 else '= %.1e' % s['permutation_p_by_chromosome'])}", transform=ax.transAxes, fontsize=6.6, va="bottom")
        ax.set_xlabel("K-mer variations per 100-kb window"); ax.set_ylabel("Density, HE breakpoints", color="#0b7285")
        for sp_ in ("top", "right"): ax.spines[sp_].set_visible(False)
        # violin: non-breakpoint, timococcum pairs, zhukovskyi pairs
        ax2 = fig.add_subplot(gs[1, j]); groups = [("Non-breakpoint\nHE active\nregion", y, GREY)]
        for sp, c in dict([("T. zhukovskyi", ZHUK), (a.first_name, TIMO)]).items() if a.first_name == "T. zhukovskyi" else ((a.first_name, TIMO), ("T. zhukovskyi", ZHUK)):
            if sp in SPECIES: groups.append((f"{sp}\nHE breakpoint\nwindows", w.loc[P.loc[P.species == sp, col].values, "variations"].values, c))
        vp = ax2.violinplot([v for _, v, _ in groups], positions=range(len(groups)), widths=.8, showextrema=False)
        for body, (_, v, c) in zip(vp["bodies"], groups): body.set_facecolor(c); body.set_alpha(.7); body.set_edgecolor("none")
        for k, (_, v, c) in enumerate(groups):
            q1, q3 = np.percentile(v, [25, 75]); ax2.plot([k, k], [q1, q3], color="black", lw=2.2, solid_capstyle="butt", zorder=4)   # interquartile range
            m_, se_ = v.mean(), v.std(ddof=1) / np.sqrt(len(v))                                                                     # mean and standard error of the mean
            ax2.plot([k], [m_], "o", mfc="white", mec="black", ms=4.5, zorder=5)
            ax2.text(k, 0, f"n = {len(v):,}\nmean = {m_:,.0f}\nSE = ±{se_:,.1f}", ha="center", va="bottom", fontsize=6.3, transform=ax2.get_xaxis_transform(), linespacing=1.25)
        ax2.set_xticks(range(len(groups))); ax2.set_xticklabels([l for l, _, _ in groups], fontsize=6.8); ax2.set_ylabel("K-mer variations per 100-kb window")
        glabel = ["non-breakpoint active windows"] + [f"{sp} HE breakpoints" for sp in SPECIES]
        # pairwise tests between the three violins: two-sided Mann-Whitney U, Holm-adjusted over the 3 comparisons; Cliff's delta as effect size
        vals = [v for _, v, _ in groups]; ng = len(groups); cmp_ = [(0, 1), (0, 2), (1, 2)][:1 if ng == 2 else 3]; nc = len(cmp_); raw = []
        for p_, q_ in cmp_:
            u, pv = stats.mannwhitneyu(vals[p_], vals[q_], alternative="two-sided"); raw.append((pv, 2 * u / (len(vals[p_]) * len(vals[q_])) - 1))
        order = np.argsort([r[0] for r in raw]); adj = np.empty(nc); run = 0
        for rank, i_ in enumerate(order): run = max(run, min(1, (nc - rank) * raw[i_][0])); adj[i_] = run
        kw = stats.kruskal(*vals) if ng > 2 else None
        for (p_, q_), (pv, cd_), pa in zip(cmp_, raw, adj):
            VT.append(dict(subgenome=g, group_1=glabel[p_], group_2=glabel[q_], n_1=len(vals[p_]), n_2=len(vals[q_]), median_1=np.median(vals[p_]), median_2=np.median(vals[q_]),
                           mannwhitney_p=pv, holm_adjusted_p=pa, cliffs_delta_1_vs_2=cd_, kruskal_wallis_p=(kw.pvalue if kw else float("nan"))))
        ytop = np.quantile(w[w.active].variations, .9995); ax2.set_ylim(top=ytop * (1.22 if True else 1)); step = ytop * .075; base = ytop * 1.02
        for k_, ((p_, q_), pa) in enumerate(zip([(0, 1), (1, 2), (0, 2)][:nc] if nc == 3 else [(0, 1)], [adj[0], adj[2], adj[1]] if nc == 3 else [adj[0]])):
            yb = base + k_ * step; star = "***" if pa < 1e-3 else "**" if pa < 1e-2 else "*" if pa < .05 else "n.s."
            ax2.plot([p_, p_, q_, q_], [yb - step * .25, yb, yb, yb - step * .25], color="black", lw=.7)
            ax2.text((p_ + q_) / 2, yb, star if star != "n.s." else "n.s.", ha="center", va="bottom", fontsize=7)
        for sp_ in ("top", "right"): ax2.spines[sp_].set_visible(False)
        ax2.set_title(f"{'cd'[j]}  {g}", loc="left", fontsize=9, fontweight="bold")
    fig.legend(handles=[Patch(fc=COMB, alpha=.65, label="Unique HE breakpoints (" + " + ".join(SPECIES) + ")"), Patch(fc=GREY, alpha=.65, label="Non-breakpoint windows within HE active regions around telomere")], loc="upper center", bbox_to_anchor=(.5, .995), ncol=2, frameon=False, fontsize=7.8)
    fig.subplots_adjust(top=.88, bottom=.08, left=.09, right=.88); import re
    fig.savefig(a.out + ".png", dpi=600); fig.savefig(a.out + ".pdf")
    fig.savefig(a.out + ".svg")                                                 # editable text
    t = open(a.out + ".svg", encoding="utf-8").read()
    t = re.sub(r"font:([^;\"]*?)'?TeX Gyre Heros'?", lambda m: "font:" + m.group(1) + "Helvetica, Arial, 'Liberation Sans', sans-serif", t)
    t = re.sub(r"font-family:\s*'?TeX Gyre Heros'?", "font-family: Helvetica, Arial, 'Liberation Sans', sans-serif", t)
    open(a.out + ".svg", "w", encoding="utf-8").write(t)
    with plt.rc_context({"svg.fonttype": "path"}): fig.savefig(a.out + "_outlined.svg")   # text converted to outlines: looks identical everywhere
    S = pd.DataFrame(rows); S.columns = [f"{g}: {s}" for g, s in S.columns]; S.to_csv(a.out + "_stats.tsv", sep="\t")
    VTd = pd.DataFrame(VT); VTd.to_csv(a.out + "_violin_tests.tsv", sep="\t", index=False); print(VTd.to_string())
    pd.set_option("display.width", 250); print(S.to_string())

if __name__ == "__main__": main()

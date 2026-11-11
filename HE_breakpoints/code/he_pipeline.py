#!/usr/bin/env python3
"""Homoeologous exchange (HE) breakpoint pipeline for T. zhukovskyi.

Input : all_breakpoints.tsv  (accession, sample, chrom, position, interval_start, interval_end,
                              dosage_left, dosage_right, change)  one row per called breakpoint
        manual_corrections.tsv (optional; edits/additions made after checking the coverage plots)
Output: breakpoints_long.tsv   one row per breakpoint (the 'long' table)
        Table_S_breakpoints.xlsx  long table + per-plant/chromosome summary sheets

Conventions
  chrom   = chr<group><A|G>_<source accession>;  Am = chr?A_TA10622, At = chr?A_TA2804, G = chr?G_TA2804
  Breakpoint position = first bin (100 kb) of the new dosage state.
Classification
  A-chromosome breakpoints with a reciprocal At/Am partner -> 'HE'
  A-chromosome breakpoints without partner               -> 'Non-HE structural'
  G-chromosome breakpoints                               -> 'G-genome'
"""
import argparse, sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd

# ---------------- parameters (change here, then re-run) ----------------
SHARED_MIN_PLANTS = 17      # a site carried by >= this many plants is 'shared' (set with --shared-min-plants; ~90% of the plants)
N_PLANTS = 19               # overwritten from the data
FIXED_DOMINANT    = 0.90    # ...and 'Fixed' if one dosage transition covers >= this fraction of those plants
WIN_MIN_OFFSET_MB = -5.0    # Am position may be this much before the At position ...
WIN_BASE_MB, WIN_FRAC = 15.0, 0.12   # ... or up to base + frac*At position after it (Am chromosomes are longer)
SCORE_FRAC, SCORE_BASE = 0.07, 3.0   # expected offset = 0.07*At + 3 Mb (used only to rank candidate pairs)
CONSENSUS_TOL_MB  = 2.0     # tolerance when re-using a site pairing seen in >= 2 other plants
COMPOUND_SPAN_MB  = 15.0    # max distance between the two breakpoints of a split (compound) change
BIN_BP            = 100_000 # coverage bin size

SUB = {("A","TA10622"):"Am", ("A","TA2804"):"At", ("G","TA2804"):"G"}

# ---------------- load / correct / annotate ----------------
def load(bp_path, corr_path=None, accession_regex=None, strip_suffix="", name_regex=None):
    df = pd.read_csv(bp_path, sep="\t")
    if strip_suffix: df["sample"] = df["sample"].str.replace(strip_suffix+"$", "", regex=True)
    if name_regex:   # named groups (lineage, generation, plant_no) are parsed from the sample name
        ext = df["sample"].str.extract(name_regex)
        for c in ext.columns: df[c] = ext[c]
        if ext.isna().all(axis=1).any(): print("WARNING: --name-regex did not match some sample names, e.g.", df.loc[ext.isna().all(axis=1),"sample"].iloc[0])
        if "accession" not in df.columns and "lineage" in df.columns: df.insert(0, "accession", df["lineage"].fillna("NA"))
    if "accession" not in df.columns:   # e.g. the skim-seq file has no accession column
        df.insert(0, "accession", df["sample"].str.extract(accession_regex)[0].fillna("NA") if accession_regex else "NA")
    df = df.drop(columns=[c for c in ("change",) if c in df.columns])
    if corr_path:
        for r in pd.read_csv(corr_path, sep="\t").itertuples():
            m = (df["sample"]==r.sample)&(df.chrom==r.chrom)&(df.position==r.position)
            if r.action == "edit":
                assert m.sum()==1, f"edit target not found: {r}"
                for col in ("interval_start","interval_end","dosage_left","dosage_right"):
                    v = getattr(r,col)
                    if pd.notna(v): df.loc[m,col] = int(v)
            elif r.action == "add":
                assert m.sum()==0, f"row to add already exists: {r}"
                acc = df.loc[df["sample"]==r.sample,"accession"].iloc[0]
                df.loc[len(df)] = [acc,r.sample,r.chrom,int(r.position),int(r.interval_start),int(r.interval_end),int(r.dosage_left),int(r.dosage_right)]
            else: sys.exit(f"unknown action {r.action}")
    return df   # annotate() is called in main() after controls are split off

def add_controls(df, ctrl, tol=300_000):
    """For every breakpoint: number of parental/control samples with a breakpoint on the same chromosome within +-tol bp.
    controls_with_site counts all controls; own_parents_with_site counts only the plant's own parents, taken from its lineage
    (lineage 'MVGB530-1T1' -> parents '1T1' and MVGB530 / MVGB530-1 / MVGB530-2 ...). A site in an own parent is inherited, not new."""
    cp = {ch: g[["position","sample"]].values for ch, g in ctrl.groupby("chrom")}
    names = sorted(ctrl["sample"].unique())
    def own(lin):
        if not isinstance(lin, str): return set()
        p1 = lin.split("-")[0]; p2 = lin.split("-",1)[1] if "-" in lin else None
        return {c for c in names if c == p1 or c.startswith(p1+"-") or (p2 is not None and c == p2)}
    have_lin = "lineage" in df.columns
    own_map = {l: own(l) for l in df["lineage"].unique()} if have_lin else {}
    n, o = [], []
    lin = df["lineage"].values if have_lin else [None]*len(df)
    for ch, pos, l in zip(df.chrom, df.position, lin):
        g = cp.get(ch)
        hits = set() if g is None else set(g[np.abs(g[:,0].astype(np.int64)-pos) <= tol, 1])
        n.append(len(hits))
        if have_lin: o.append(len(hits & own_map[l]))
    df["controls_with_site"] = n; df["controls_total"] = len(names)
    if have_lin:
        df["own_parents_with_site"] = o; df["own_parents_total"] = [len(own_map[l]) for l in lin]
    return df

def annotate(df):
    df = df.copy(); global N_PLANTS; N_PLANTS = df["sample"].nunique()
    df[["chr","source"]] = df.chrom.str.split("_", expand=True)
    df["group"] = df.chr.str[3].astype(int); df["genome"] = df.chr.str[4]
    df["subgenome"] = [SUB[(g,s)] for g,s in zip(df.genome, df.source)]
    df["Mb"] = df.position/1e6; df["chg"] = df.dosage_right - df.dosage_left
    site = df.groupby(["chrom","position"])
    df["plants_with_site"] = site["sample"].transform("nunique")
    df["dominant_frac"] = site.apply(lambda d: d.groupby(["dosage_left","dosage_right"]).size().max()/len(d), include_groups=False)\
                              .reindex(pd.MultiIndex.from_frame(df[["chrom","position"]])).values
    df["site_class"] = np.where((df.plants_with_site>=SHARED_MIN_PLANTS)&(df.dominant_frac>=FIXED_DOMINANT), "Fixed", "Variable")
    return df.sort_values(["sample","chrom","position"]).reset_index(drop=True)

# ---------------- reciprocal pairing (A chromosomes only) ----------------
def in_window(a, m):
    d = m.Mb - a.Mb
    return WIN_MIN_OFFSET_MB <= d <= WIN_BASE_MB + WIN_FRAC*a.Mb
def score(a, m): return abs((m.Mb-a.Mb) - (SCORE_FRAC*a.Mb + SCORE_BASE))

def pair_breakpoints(df):
    """returns pairs (list of dict) and partners {breakpoint index: [partner index,...]}"""
    A = df[df.genome=="A"]
    pairs, partners = [], defaultdict(list)
    def add(ia, im, mtype):
        pairs.append(dict(at=ia, am=im, match_type=mtype))
        partners[ia].append(im); partners[im].append(ia)
    for (p,g), d in A.groupby(["sample","group"]):
        At = d[d.subgenome=="At"]; Am = d[d.subgenome=="Am"]; used=set()
        for mtype, cond in (("Exact (equal and opposite)", lambda a,m: a.chg==-m.chg),
                            ("Opposite sign only",        lambda a,m: a.chg*m.chg<0)):
            cand = sorted((score(a,m), ia, im) for ia,a in At.iterrows() for im,m in Am.iterrows() if cond(a,m) and in_window(a,m))
            for _, ia, im in cand:
                if ia in used or im in used: continue
                used |= {ia, im}; add(ia, im, mtype)
    # consensus pass: re-use the most common At<->Am site pairing (seen in >=2 plants) for breakpoints still unpaired
    cnt = defaultdict(Counter)
    for r in pairs: cnt[(df.at[r["at"],"group"], df.at[r["at"],"position"])][df.at[r["am"],"position"]] += 1
    cmap = {k: v.most_common(1)[0][0] for k,v in cnt.items() if sum(v.values())>=2}
    for (p,g), d in A.groupby(["sample","group"]):
        At = d[d.subgenome=="At"]; Am = d[d.subgenome=="Am"]
        for ia,a in At.iterrows():
            if ia in partners: continue
            tgt = cmap.get((g, a.position))
            if tgt is None: continue
            for im,m in Am.iterrows():
                if im in partners: continue
                if abs(m.position-tgt)/1e6 <= CONSENSUS_TOL_MB and a.chg*m.chg<0:
                    add(ia, im, "Consensus (site pairing seen in other plants)"); break
    # compound pass: a reciprocal change split over two nearby breakpoints on one subgenome
    for r in list(pairs):
        a, m = df.loc[r["at"]], df.loc[r["am"]]
        d = A[(A["sample"]==a["sample"])&(A.group==a.group)]
        for iu,u in d.iterrows():
            if iu in partners: continue
            if u.subgenome=="At" and abs(u.Mb-a.Mb)<=COMPOUND_SPAN_MB and a.chg+u.chg==-m.chg and a.chg*u.chg>0:
                add(iu, r["am"], "Compound (split breakpoint)")
            elif u.subgenome=="Am" and abs(u.Mb-m.Mb)<=COMPOUND_SPAN_MB and m.chg+u.chg==-a.chg and m.chg*u.chg>0:
                add(r["at"], iu, "Compound (split breakpoint)")
    return pairs, partners

# ---------------- long table ----------------
def long_table(df, pairs, partners):
    out = df.copy()
    out["category"] = np.where(out.genome=="G", "G-genome", np.where(out.index.isin(partners.keys()), "HE", "Non-HE structural"))
    pm = {}
    for i,r in enumerate(pairs): pm.setdefault(r["at"],[]).append(i); pm.setdefault(r["am"],[]).append(i)
    pid, ppos, pall, ptype = [], [], [], []
    for i in out.index:
        ps = partners.get(i, [])
        pid.append("; ".join(f"P{n+1:04d}" for n in pm.get(i,[])) or "")
        ppos.append(int(out.at[ps[0],"position"]) if ps else np.nan)
        pall.append("; ".join(f"{out.at[j,'position']/1e6:.2f}" for j in ps))
        ptype.append("; ".join(pairs[n]["match_type"] for n in pm.get(i,[])))
    out["pair_id"]=pid; out["partner_position_bp"]=ppos; out["partner_positions_Mb"]=pall; out["match_type"]=ptype
    out["partner_subgenome"] = np.where(out.partner_position_bp.notna(), np.where(out.subgenome=="Am","At","Am"), "")
    out["window_start_bp"] = (out.position//BIN_BP)*BIN_BP          # 100-kb bin containing the breakpoint (join key for window tables)
    out["window_end_bp"] = out.window_start_bp + BIN_BP
    out["breakpoint_id"] = [f"BP{n+1:05d}" for n in range(len(out))]
    extra = [c for c in ("lineage","generation","plant_no") if c in out.columns]
    cols = ["breakpoint_id","accession","sample"]+extra+["chrom","chr","group","genome","subgenome","source","position","interval_start","interval_end",
            "window_start_bp","window_end_bp","dosage_left","dosage_right","chg","category","site_class","plants_with_site",
            "pair_id","partner_subgenome","partner_position_bp","partner_positions_Mb","match_type"]+[c for c in ("controls_with_site","controls_total","own_parents_with_site","own_parents_total") if c in out.columns]
    out = out[cols].rename(columns={"accession":"parental_accession","sample":"plant_id","chr":"chromosome","genome":"genome_letter",
            "source":"source_accession","chg":"dosage_change","group":"homoeologous_group","plants_with_site":"plants_with_same_site"})
    out["position_Mb"] = (out.position/1e6).round(2)
    return out.sort_values(["plant_id","homoeologous_group","genome_letter","subgenome","position"]).reset_index(drop=True)

# ---------------- wide sheets (from the long table) ----------------
def wide(long_, category, order):
    rows=[]
    for (p,acc,g), d0 in long_[long_.category.isin(category)].groupby(["plant_id","parental_accession","homoeologous_group"]):
        for sub in order:
            d = d0[d0.subgenome==sub].sort_values("position")
            if d.empty: continue
            pos = ", ".join(f"{x:.2f}" for x in d.position_Mb)
            chg = ", ".join(f"{a}x>{b}x" for a,b in zip(d.dosage_left,d.dosage_right))
            rec = ", ".join((" + ".join(s.split("; "))+f" ({d.partner_subgenome.iloc[i]})") if s else "–" for i,s in enumerate(d.partner_positions_Mb))
            rows.append([p,acc,d.chromosome.iloc[0],d.source_accession.iloc[0],sub,len(d),pos,chg,rec,
                         ", ".join(str(int(x)) for x in d.plants_with_same_site)])
    return rows

def write_excel(long_, path):
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    hdr = PatternFill("solid", fgColor="1F3864")
    def style(ws, widths):
        for c in ws[1]: c.font=Font(name="Arial",bold=True,color="FFFFFF",size=10); c.fill=hdr; c.alignment=Alignment(horizontal="center",vertical="center",wrap_text=True)
        ws.row_dimensions[1].height=42
        for row in ws.iter_rows(min_row=2):
            for c in row: c.font=Font(name="Arial",size=10)
        for i,w in enumerate(widths): ws.column_dimensions[ws.cell(1,i+1).column_letter].width=w
        ws.freeze_panes="A2"; ws.auto_filter.ref=ws.dimensions
    wb = Workbook(); rd = wb.active; rd.title="Readme"
    for t in ["Homoeologous exchange breakpoints in T. zhukovskyi (generated by he_pipeline.py)",
              "Long_format: one row per breakpoint (all plants, all chromosomes). category = HE | Non-HE structural | G-genome. window_start_bp/window_end_bp give the 100-kb bin of the breakpoint, for joining to window-based tables.",
              "HE: A-chromosome breakpoints with a reciprocal At/Am partner (Am = TA10622 A chromosomes, At = TA2804 A chromosomes). Partner positions use each chromosome's own coordinates; At and Am are not collinear.",
              "Pairing: same plant and homoeologous group, opposite dosage change; exact (equal and opposite) first, then opposite sign only; one-to-one, nearest to the expected offset first; then consensus (site pairing seen in >=2 other plants); then compound (change split over two nearby breakpoints).",
              "Non-HE structural: A-chromosome breakpoints without a reciprocal partner. plants_with_same_site close to the number of plants suggests a discrepancy with the reference (possible assembly artifact); 1 = plant-specific.",
              "G-genome: all breakpoints on G chromosomes. Parameters and manual corrections are listed at the top of he_pipeline.py and in manual_corrections.tsv."]: rd.append([t])
    rd["A1"].font=Font(name="Arial",bold=True,size=12); rd.column_dimensions["A"].width=170
    ws = wb.create_sheet("Long_format"); ws.append(list(long_.columns))
    for r in long_.itertuples(index=False): ws.append([None if (isinstance(v,float) and np.isnan(v)) else v for v in r])
    style(ws, [11,11,13,14,11,9,9,10,10,10,14,14,14,14,12,12,10,16,10,10,14,10,16,18,30,10])
    H = ["Plant_ID","Parental_accession","Chromosome","Source_accession","Subgenome","N_breakpoints","Breakpoint_positions_Mb","Dosage_change_at_breakpoints","Reciprocal_breakpoint_Mb (partner subgenome)","Plants_with_same_site (of N)"]
    H[-1] = f"Plants_with_same_site (of {N_PLANTS})"
    ws = wb.create_sheet("HE_Am_At"); ws.append(H[:-1])
    for r in wide(long_, ["HE"], ["Am","At"]): ws.append(r[:-1])
    style(ws,[14,12,11,12,10,12,52,52,50])
    ws = wb.create_sheet("Non_HE_structural"); ws.append(H[:8]+[H[9]])
    for r in wide(long_, ["Non-HE structural"], ["Am","At"]): ws.append(r[:8]+[r[9]])
    style(ws,[14,12,11,12,10,12,48,48,24])
    ws = wb.create_sheet("G_genome_changes"); ws.append(H[:8]+[H[9]])
    for r in wide(long_, ["G-genome"], ["G"]): ws.append(r[:8]+[r[9]])
    style(ws,[14,12,11,12,10,12,48,48,24])
    wb.save(path)

def main():
    global SHARED_MIN_PLANTS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--breakpoints", required=True); ap.add_argument("--corrections", default=None); ap.add_argument("--outdir", default="results")
    ap.add_argument("--shared-min-plants", type=int, default=SHARED_MIN_PLANTS, help="a breakpoint site in at least this many plants is shared (use ~90%% of the plants; 17 for 19 plants, 85 for 95)")
    ap.add_argument("--accession-regex", default=None, help="regex with one group to take the parental accession from the sample name, e.g. '^(\\d+T\\d+)' (default: NA)")
    ap.add_argument("--strip-suffix", default="", help="remove this suffix from sample names, e.g. _100Kb")
    ap.add_argument("--name-regex", default=None, help="regex with named groups (?P<lineage>..)(?P<generation>..)(?P<plant_no>..) applied to the sample name (after --strip-suffix); adds lineage/generation/plant_no columns")
    ap.add_argument("--control-regex", default=None, help="regex matching sample names of parents/controls (e.g. '^(?!timo_)'); they are removed from the analysis, written to controls_breakpoints.tsv, and each breakpoint gets controls_with_site")
    a = ap.parse_args()
    import os; os.makedirs(a.outdir, exist_ok=True)
    SHARED_MIN_PLANTS = a.shared_min_plants
    df = load(a.breakpoints, a.corrections, a.accession_regex, a.strip_suffix, a.name_regex)
    ctrl = None
    if a.control_regex:
        isc = df["sample"].str.contains(a.control_regex, regex=True)
        ctrl = df[isc].copy(); df = df[~isc].copy()
        ctrl.to_csv(f"{a.outdir}/controls_breakpoints.tsv", sep="\t", index=False)
        print(f"{df['sample'].nunique()} individuals analysed; {ctrl['sample'].nunique()} controls set aside: {', '.join(sorted(ctrl['sample'].unique()))}")
    df = annotate(df)
    if ctrl is not None: df = add_controls(df, ctrl)
    print(f"shared threshold: {SHARED_MIN_PLANTS} of {N_PLANTS} individuals")
    pairs, partners = pair_breakpoints(df)
    long_ = long_table(df, pairs, partners)
    long_.to_csv(f"{a.outdir}/breakpoints_long.tsv", sep="\t", index=False)
    write_excel(long_, f"{a.outdir}/Table_S_breakpoints.xlsx")
    print(long_.category.value_counts().to_string()); print("pairs:", pd.Series([p['match_type'] for p in pairs]).value_counts().to_dict())
if __name__ == "__main__": main()

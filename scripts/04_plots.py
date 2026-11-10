#!/usr/bin/env python
"""Presentation figures.
Usage: python scripts/04_plots.py            # QC + prevalence (works on the pilot)
       python scripts/04_plots.py --gwas     # + Manhattan/QQ for every GEMMA result,
                                             #   gwas_hits_<trait>.tsv and gwas_summary.tsv
Figures -> $WORK/results/figures, tables -> $WORK/results/tables
"""
import glob, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORK, PROJ = os.environ["WORK"], os.environ["PROJ"]
T = f"{WORK}/results/tables"; F = f"{WORK}/results/figures"; os.makedirs(F, exist_ok=True)
qc = pd.read_csv(f"{T}/qc_summary.tsv", sep="\t", index_col=0)
tg = pd.read_csv(f"{T}/targets.tsv", sep="\t")

# 1. Where do the unmapped reads go? (stacked fractions per accession, sorted)
parts = pd.DataFrame({
    "mitochondrion": qc.mito, "chloroplast": qc.chloroplast, "PhiX": qc.phix,
    "Kraken-classified": qc.kraken_classified,
    "unclassified": (qc.clean_reads - qc.kraken_classified).clip(lower=0)})
frac = parts.div(parts.sum(axis=1), axis=0).sort_values("Kraken-classified")
ax = frac.plot.bar(stacked=True, width=1.0, figsize=(10, 4), linewidth=0)
ax.set_xticks([]); ax.set_xlabel(f"accessions (n={len(frac)})"); ax.set_ylabel("fraction of unmapped reads")
ax.legend(ncol=5, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), frameon=False)
plt.tight_layout(); plt.savefig(f"{F}/01_unmapped_composition.png", dpi=200); plt.close()

# 2. Top 25 species by mean RPM
rpm = pd.read_csv(f"{T}/species_rpm.tsv", sep="\t", index_col=0)
top = rpm.mean().sort_values(ascending=False).head(25)[::-1]
plt.figure(figsize=(7, 7)); plt.barh(top.index, top.values)
plt.xscale("log"); plt.xlabel("mean reads per million total reads")
plt.title("Most abundant non-rice species in discarded reads")
plt.tight_layout(); plt.savefig(f"{F}/02_top_species.png", dpi=200); plt.close()

# 3. Pathogen RPM by subpopulation (log scale, called-present in color)
taxa = tg.taxon.unique()
fig, axes = plt.subplots(1, len(taxa), figsize=(3 * len(taxa), 4), sharey=False)
for ax, t in zip(np.atleast_1d(axes), taxa):
    d = tg[tg.taxon == t]; subs = sorted(d.subpop_k9.dropna().unique())
    for i, s in enumerate(subs):
        x = d[d.subpop_k9 == s]; jit = i + np.random.uniform(-.3, .3, len(x))
        ax.scatter(jit, x.rpm + 1e-3, s=6, c=np.where(x.present == 1, "crimson", "grey"), alpha=.7)
    ax.set_yscale("log"); ax.set_xticks(range(len(subs))); ax.set_xticklabels(subs, rotation=90, fontsize=7)
    ax.set_title(t, fontsize=9, style="italic")
axes.flat[0].set_ylabel("reads per million (red = called present)")
plt.tight_layout(); plt.savefig(f"{F}/03_pathogens_by_subpop.png", dpi=200); plt.close()

# 4. Validation scatter: Kraken RPM vs direct-alignment breadth
v = tg.dropna(subset=["breadth_pct"])
plt.figure(figsize=(5, 4))
for t, d in v.groupby("taxon"):
    plt.scatter(d.rpm + 1e-3, d.breadth_pct + 1e-3, s=8, label=t, alpha=.7)
plt.xscale("log"); plt.yscale("log"); plt.xlabel("Kraken2 RPM"); plt.ylabel("genome breadth covered (%)")
plt.axhline(1.0, ls="--", c="k", lw=.7); plt.legend(fontsize=7)
plt.tight_layout(); plt.savefig(f"{F}/04_kraken_vs_alignment.png", dpi=200); plt.close()

# 5. Contamination check: does load track sequencing depth or read length (batch proxies)?
k = tg.groupby("accession").rpm.sum()
fig, ax = plt.subplots(1, 2, figsize=(8, 3.5))
ax[0].scatter(qc.total, k.reindex(qc.index), s=6); ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("total reads"); ax[0].set_ylabel("summed target RPM")
qc.boxplot(column="frac_classified_of_clean", by="read_len", ax=ax[1]); ax[1].set_title("")
ax[1].set_ylabel("fraction classified"); plt.suptitle("")
plt.tight_layout(); plt.savefig(f"{F}/05_batch_checks.png", dpi=200); plt.close()

# 6. GWAS
if "--gwas" in sys.argv:
    genes = pd.read_csv(f"{PROJ}/metadata/known_R_genes.tsv", sep="\t")
    files = sorted(glob.glob(f"{WORK}/gwas/output/lmm_*.assoc.txt"))
    print(f"{len(files)} GWAS result files")
    summ = []
    for f in files:
        name = os.path.basename(f)[4:-10]
        a = pd.read_csv(f, sep="\t", usecols=["chr", "ps", "p_wald"]).dropna()
        a["lp"] = -np.log10(a.p_wald); a = a.sort_values(["chr", "ps"])
        off = a.groupby("chr").ps.max().cumsum().shift(fill_value=0)
        a["x"] = a.ps + a.chr.map(off)
        fig, (m, q) = plt.subplots(1, 2, figsize=(14, 4), gridspec_kw={"width_ratios": [4, 1]})
        for c, d in a.groupby("chr"):
            m.scatter(d.x, d.lp, s=2, c="#1f4e79" if c % 2 else "#7fa7d0")
        bonf = -np.log10(0.05 / len(a)); m.axhline(bonf, c="crimson", lw=.8)
        top_y = max(a.lp.max(), bonf) * 1.35; m.set_ylim(0, top_y)   # headroom so gene labels stay inside
        for _, g in genes.iterrows():
            if g.chr in off.index:
                gx = g.start + off[g.chr]; m.axvline(gx, c="orange", lw=.6, ls=":")
                m.text(gx, top_y * 0.98, g.gene.split(" ")[0], rotation=90, fontsize=6, color="darkorange", va="top", ha="right")
        m.set_xticks([off[c] + a[a.chr == c].ps.max() / 2 for c in off.index]); m.set_xticklabels(off.index)
        m.set_xlabel("chromosome"); m.set_ylabel("-log10 p"); m.set_title(name)
        e = -np.log10((np.arange(1, len(a) + 1) - .5) / len(a)); o = np.sort(a.lp.values)[::-1]
        lam = np.median(np.square(__import__("scipy.stats").stats.norm.isf(a.p_wald / 2))) / 0.4549
        q.scatter(e, o, s=2); q.plot([0, e.max()], [0, e.max()], c="k", lw=.7)
        q.set_title(f"QQ  λ={lam:.2f}"); q.set_xlabel("expected"); q.set_ylabel("observed")
        plt.tight_layout(); plt.savefig(f"{F}/gwas_{name}.png", dpi=200); plt.close()
        hits = a[a.lp > bonf]
        hits.to_csv(f"{T}/gwas_hits_{name}.tsv", sep="\t", index=False)
        top = a.loc[a.lp.idxmax()]
        summ.append(dict(trait=name, n_snps=len(a), lambda_gc=round(lam, 3), n_bonf_hits=len(hits),
                         top_chr=int(top.chr), top_ps=int(top.ps), top_log10p=round(top.lp, 2), bonf_log10p=round(bonf, 2)))
    if summ:
        sm = pd.DataFrame(summ); sm.to_csv(f"{T}/gwas_summary.tsv", sep="\t", index=False)
        print(sm.to_string(index=False))
print(f"figures in {F}")

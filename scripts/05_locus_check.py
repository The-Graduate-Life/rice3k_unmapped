#!/usr/bin/env python
"""Is a GWAS hit for pathogen presence robust to sample batch AND subpopulation?

1. Reads tables/gwas_hits_<trait>.tsv, groups significant SNPs within 1 Mb into loci and keeps the
   best SNP per locus.
2. Extracts that SNP's genotype (0/1/2 copies) with plink.
3. Compares homozygotes (0 vs 2) only inside strata of the same ID batch x K9 subpopulation and pools
   them: Mantel-Haenszel odds ratio + Cochran-Mantel-Haenszel test (no continuity correction).

Usage: python scripts/05_locus_check.py [trait]        (default: present_Xanthomonas_oryzae)
Output: printed report + tables/locus_check_<trait>.tsv
"""
import os, subprocess, sys
import numpy as np, pandas as pd
from scipy.stats import chi2

W = os.environ["WORK"]
TRAIT = sys.argv[1] if len(sys.argv) > 1 else "present_Xanthomonas_oryzae"
TAXON = TRAIT.split("_", 1)[1].replace("_", " ")
T = f"{W}/results/tables"; G = f"{W}/gwas"; L = f"{G}/loci"; os.makedirs(L, exist_ok=True)

hits_f = f"{T}/gwas_hits_{TRAIT}.tsv"
if not os.path.exists(hits_f) or os.path.getsize(hits_f) == 0:
    sys.exit(f"{hits_f} missing: run scripts/04_plots.py --gwas first")
hits = pd.read_csv(hits_f, sep="\t")
if hits.empty:
    sys.exit(f"no genome-wide significant SNPs for {TRAIT}; nothing to check")

# group hits into loci (same chromosome, < 1 Mb apart) and keep the best SNP of each
hits = hits.sort_values(["chr", "ps"]).reset_index(drop=True)
locus = (hits.chr.ne(hits.chr.shift()) | hits.ps.diff().gt(1_000_000)).cumsum()
best = hits.loc[hits.groupby(locus).lp.idxmax()]

t = pd.read_csv(f"{T}/targets.tsv", sep="\t")
x = t[t.taxon == TAXON].set_index("accession")
acc = x.index.to_series()
x["batch"] = (acc.str.extract(r"IRIS_313-(\d+)", expand=False).astype(float).floordiv(1000)
              .map(lambda v: f"IRIS{int(v)}k", na_action="ignore").fillna(acc.str.replace(r"\d+", "", regex=True)))

print(f"== {TAXON}: prevalence by subpopulation ==")
print(x.groupby("subpop_k9").present.agg(["sum", "count", "mean"]).round(3).to_string())


def cmh(df):
    """df: geno (0/2), present (0/1), stratum -> (MH OR, CMH chi2, p, informative strata)."""
    num = den = a_sum = e_sum = v_sum = 0.0; used = 0
    for _, s in df.groupby("stratum"):
        a = ((s.geno == 2) & (s.present == 1)).sum(); b = ((s.geno == 2) & (s.present == 0)).sum()
        c = ((s.geno == 0) & (s.present == 1)).sum(); d = ((s.geno == 0) & (s.present == 0)).sum()
        n = a + b + c + d
        if n < 2 or 0 in (a + b, c + d, a + c, b + d):
            continue                                        # stratum carries no information
        used += 1
        num += a * d / n; den += b * c / n
        a_sum += a; e_sum += (a + b) * (a + c) / n
        v_sum += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
    if used == 0 or v_sum == 0:
        return np.nan, np.nan, np.nan, 0
    stat = (a_sum - e_sum) ** 2 / v_sum
    return (num / den if den > 0 else np.inf), stat, chi2.sf(stat, 1), used


rows = []
for _, h in best.iterrows():
    c, ps = int(h.chr), int(h.ps); out = f"{L}/chr{c}_{ps}"
    subprocess.run(["plink", "--bfile", f"{G}/pruned_v2.1", "--chr", str(c), "--from-bp", str(ps),
                    "--to-bp", str(ps), "--recode", "A", "--allow-no-sex", "--out", out],
                   check=True, stdout=subprocess.DEVNULL)
    g = pd.read_csv(f"{out}.raw", sep=r"\s+").set_index("IID").iloc[:, -1].rename("geno")
    y = x.join(g, how="inner"); y = y[y.geno.isin([0, 2])].copy()
    print(f"\n######## chr{c}:{ps:,}  (-log10 p = {h.lp:.2f}; {int((locus == locus[h.name]).sum())} SNPs in locus)")
    print("-- infection by genotype (all homozygous accessions) --")
    print(y.groupby("geno").present.agg(["sum", "count", "mean"]).round(3).to_string())
    print("-- within each subpopulation carrying both genotypes --")
    both = sorted(s for s in y.subpop_k9.unique() if y[y.subpop_k9 == s].geno.nunique() == 2)
    if both:
        print(y[y.subpop_k9.isin(both)].groupby(["subpop_k9", "geno"]).present.agg(["sum", "count", "mean"]).round(3).to_string())
    row = dict(trait=TRAIT, chr=c, ps=ps, gwas_log10p=round(h.lp, 2))
    for name, strat in [("batch", y.batch), ("subpop", y.subpop_k9), ("batch_x_subpop", y.batch + "|" + y.subpop_k9)]:
        orr, stat, p, k = cmh(y.assign(stratum=strat))
        print(f"   stratified by {name:15s}: MH OR = {orr:6.2f}   CMH chi2 = {stat:6.2f}   p = {p:.2e}   ({k} informative strata)")
        row.update({f"OR_{name}": round(orr, 2), f"p_{name}": p, f"strata_{name}": k})
    rows.append(row)

pd.DataFrame(rows).to_csv(f"{T}/locus_check_{TRAIT}.tsv", sep="\t", index=False)
print(f"\nwrote {T}/locus_check_{TRAIT}.tsv")

#!/usr/bin/env python
"""Merge per-accession outputs into analysis tables + GWAS phenotype files.

Usage: python scripts/03_aggregate.py  (reads $WORK/results, $PROJ/metadata)

Outputs ($WORK/results/tables):
  qc_summary.tsv            per-accession read accounting
  species_counts.tsv        accession x species (Bracken est. reads)
  species_rpm.tsv           same, reads per million TOTAL sequenced reads
  targets.tsv               target pathogens: reads, RPM, minimizers, breadth, identity, call
  target_hits.tsv           targets.tsv rows that are not "absent" (input for 02_identity_check.sh)
  prevalence_by_subpop.tsv / prevalence_by_country.tsv / prevalence_by_batch.tsv
  gwas_pheno.txt            GEMMA -p file, rows in pruned_v2.1.fam order, one column per trait
  gwas_pheno_names.txt      column index -> trait name (use with gemma -n)
  gwas_covar.txt            intercept + PC1-5 + log10(total reads)
  gwas_covar_batch.txt      gwas_covar.txt + sample-batch dummies
The GWAS files are skipped when $WORK/gwas/pruned_v2.1.fam is absent (e.g. the smoke test).
"""
import glob, os, sys
import numpy as np, pandas as pd
from scipy.stats import norm, rankdata

WORK = os.environ.get("WORK"); PROJ = os.environ.get("PROJ")
if not WORK or not PROJ:
    sys.exit("source config.sh first (needs WORK and PROJ)")
R = f"{WORK}/results"; OUT = f"{R}/tables"; os.makedirs(OUT, exist_ok=True)

meta = pd.read_csv(f"{PROJ}/metadata/accession_metadata.tsv", sep="\t").set_index("accession")
targets = pd.read_csv(f"{PROJ}/metadata/target_taxa.tsv", sep="\t", dtype={"validation_ref": str})
tot = pd.read_csv(f"{R}/qc/total_reads.tsv", sep="\t").set_index("accession")

# ---------- QC ----------
qcols = ["accession", "unmapped_reads", "chloroplast", "mito", "phix", "clean_reads",
         "kraken_classified", "read_len", "gc_clean"]
qc = pd.concat([pd.read_csv(f, sep="\t", header=None, names=qcols)
                for f in glob.glob(f"{R}/qc/*.qc.tsv")]).set_index("accession")
qc = qc.join(tot["total"]).join(meta[["country", "subpop_k9"]])
qc["frac_unmapped"] = qc.unmapped_reads / qc.total
qc["frac_classified_of_clean"] = qc.kraken_classified / qc.clean_reads.clip(lower=1)
# Flag (not drop) libraries far below the panel's depth: RPM is unstable with a tiny denominator.
qc["low_depth"] = (qc.total < 1e6).astype(int)
if qc.low_depth.sum():
    print("low-depth libraries (<1M total reads, flagged in qc_summary.tsv):", ", ".join(qc.index[qc.low_depth == 1]))
qc.to_csv(f"{OUT}/qc_summary.tsv", sep="\t")
accs = qc.index.tolist()
print(f"{len(accs)} accessions with results")

# ---------- Bracken species matrix ----------
rows = []
for a in accs:
    f = f"{R}/bracken/{a}.S.bracken"
    if os.path.exists(f) and os.path.getsize(f) > 0:
        b = pd.read_csv(f, sep="\t")
        rows.append(pd.Series(b.new_est_reads.values, index=b.name.values, name=a))
sp = pd.DataFrame(rows).fillna(0)
sp.to_csv(f"{OUT}/species_counts.tsv", sep="\t")
rpm = sp.div(qc.loc[sp.index, "total"], axis=0) * 1e6
rpm.to_csv(f"{OUT}/species_rpm.tsv", sep="\t")

# ---------- Kraken clade reads + distinct minimizers for targets & kingdoms ----------
def k2_clades(a, taxids):
    """Return {taxid: (clade_reads, distinct_minimizers)} from a --report-minimizer-data report."""
    out = {t: (0, 0) for t in taxids}
    f = f"{R}/kraken/{a}.k2report"
    if not os.path.exists(f):
        return out
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        if len(p) < 8:
            continue
        tid = int(p[6])
        if tid in out:
            out[tid] = (int(p[1]), int(p[4]))
    return out

KINGDOMS = {2: "Bacteria", 4751: "Fungi", 10239: "Viruses", 9606: "Homo sapiens"}
tids = list(targets.taxid) + list(KINGDOMS)
k2 = {a: k2_clades(a, tids) for a in accs}

# ---------- direct-alignment breadth ----------
cmap = pd.read_csv(f"{WORK}/refs/pathogen_contigs.tsv", sep="\t", header=None,
                   names=["contig", "genome"]).set_index("contig").genome
cov_cols = ["accession", "contig", "start", "end", "numreads", "covbases", "coverage",
            "meandepth", "meanbaseq", "meanmapq"]
cov_files = [f for f in glob.glob(f"{R}/patho/*.coverage.tsv") if os.path.getsize(f) > 0]
if cov_files:
    cov = pd.concat([pd.read_csv(f, sep="\t", header=None, names=cov_cols) for f in cov_files])
    cov["genome"] = cov.contig.map(cmap)
    g = cov.groupby(["accession", "genome"]).agg(reads=("numreads", "sum"), covbases=("covbases", "sum"),
                                                  length=("end", "sum"))
    g["breadth_pct"] = 100 * g.covbases / g.length
else:
    g = pd.DataFrame(columns=["reads", "covbases", "length", "breadth_pct"])

# ---------- per-read identity (scripts/02_identity_check.sh; optional) ----------
id_files = [f for f in glob.glob(f"{R}/identity/*.tsv") if os.path.getsize(f) > 0]
idt = (pd.concat([pd.read_csv(f, sep="\t", header=None,
                              names=["accession", "genome", "id_reads", "mean_identity", "frac_ge97", "frac_ge99"])
                  for f in id_files]).set_index(["accession", "genome"]) if id_files else None)
id_checked = set(idt.index.get_level_values(0)) if idt is not None else set()
print(f"identity results for {len(id_checked)} accessions")

# ---------- target table + calls ----------
# Revised after the 98-accession pilot (2026-11-08):
#  * Distinct-minimizer cut-off removed for taxa with a genome: an 8 kb virus (RTBV) cannot reach 200.
#  * Breadth is judged against the breadth EXPECTED if the aligned reads were spread uniformly:
#      expected = 1 - exp(-aligned_reads * READ_BP / genome_len)
#    A real genome gives ratio ~0.6-1; clumped hits from congeners on conserved genes give ~0.3-0.5.
#  call levels:
#    confirmed   : >=20 Kraken reads, >=50 aligned reads, breadth >=1%, breadth_ratio >= 0.6
#    congener    : Kraken + alignment signal but clumped (ratio < 0.6) -> related, non-target species
#    unvalidated : no validation genome; >=20 Kraken reads and >=200 distinct minimizers
#    absent      : otherwise
#  Revised after the full run (2026-11-09): breadth_ratio does not reject near relatives that align
#  genome-wide at ~85-90% identity (B. glumae: 6-15% breadth, >10k aligned reads but few Kraken reads).
#  When 02_identity_check.sh has been run for an accession, confirmed also needs frac_ge97 >= 0.5
#  (most aligned reads near-identical to the reference); failing that -> congener.
MIN_FRAC97 = 0.5
READ_BP = 83                      # 3K RG mean read length (B001: 11,277.88 Mb / 135.93 M reads)
MIN_READS, MIN_ALN, MIN_BREADTH, MIN_RATIO, MIN_MINIMIZERS = 20, 50, 1.0, 0.6, 200
trows = []
for a in accs:
    for _, t in targets.iterrows():
        reads, mins = k2[a][t.taxid]
        ref = t.validation_ref if isinstance(t.validation_ref, str) else None
        if ref:
            breadth = g.breadth_pct.get((a, ref), 0.0)
            aln = g.reads.get((a, ref), 0)
            glen = g.length.get((a, ref), np.nan)
            expct = 100 * (1 - np.exp(-aln * READ_BP / glen)) if aln and glen == glen else 0.0
            ratio = breadth / expct if expct > 0 else 0.0
            if a in id_checked:
                f97 = idt.frac_ge97.get((a, ref), 0.0); mid = idt.mean_identity.get((a, ref), np.nan)
            else:
                f97 = mid = np.nan
            id_ok = not (f97 == f97) or f97 >= MIN_FRAC97        # unchecked -> not penalised
            if reads >= MIN_READS and aln >= MIN_ALN and breadth >= MIN_BREADTH and ratio >= MIN_RATIO and id_ok:
                call = "confirmed"
            elif reads >= MIN_READS and aln >= MIN_ALN and breadth >= MIN_BREADTH:
                call = "congener"
            else:
                call = "absent"
        else:
            breadth = aln = expct = ratio = f97 = mid = np.nan
            call = "unvalidated" if reads >= MIN_READS and mins >= MIN_MINIMIZERS else "absent"
        trows.append(dict(accession=a, taxon=t["name"], disease=t.disease, k2_reads=reads,
                          rpm=1e6 * reads / qc.at[a, "total"], distinct_minimizers=mins,
                          aligned_reads=aln, breadth_pct=breadth, expected_breadth_pct=expct,
                          breadth_ratio=ratio, mean_identity=mid, frac_ge97=f97, call=call, present=int(call == "confirmed")))
tg = pd.DataFrame(trows).merge(meta[["country", "subpop_k9"]], left_on="accession", right_index=True)
tg.to_csv(f"{OUT}/targets.tsv", sep="\t", index=False)
for by in ["subpop_k9", "country"]:
    (tg.groupby([by, "taxon"]).present.agg(["sum", "count", "mean"]).rename(columns={"mean": "prevalence"})
       .reset_index().to_csv(f"{OUT}/prevalence_by_{by.split('_')[0]}.tsv", sep="\t", index=False))
print(pd.crosstab(tg.taxon, tg.call).to_string())
tg[tg.call != "absent"].sort_values(["taxon", "rpm"], ascending=[True, False]).to_csv(
    f"{OUT}/target_hits.tsv", sep="\t", index=False)

# Sample-batch proxy: ID series (B = 3K pilot, CX = CAAS/China, IRIS_313 by 1000-number block).
# Xoo prevalence differs strongly between IRIS blocks (seed grow-out / DNA-extraction batches),
# so a second covariate set adds batch dummies (largest batch = reference).
def id_batch(a):
    if a.startswith("IRIS_313-"):
        return f"IRIS{int(a.split('-')[1]) // 1000}k"
    return "".join(ch for ch in a if not ch.isdigit()) or "other"
tg["id_batch"] = tg.accession.map(id_batch)
bt = tg.pivot_table(index="id_batch", columns="taxon", values="present", aggfunc="mean").round(3)
bt.insert(0, "n", tg.groupby("id_batch").accession.nunique())
bt.to_csv(f"{OUT}/prevalence_by_batch.tsv", sep="\t")
print("\nprevalence by ID batch:\n" + bt.to_string())

# ---------- GWAS phenotypes (fam order) ----------
FAM = f"{WORK}/gwas/pruned_v2.1.fam"
if not os.path.exists(FAM):
    print(f"{FAM} not found: skipping GWAS phenotype/covariate files")
else:
    fam = pd.read_csv(f"{WORK}/gwas/pruned_v2.1.fam", sep=r"\s+", header=None)[0].tolist()

    def rint(x):  # rank-based inverse normal transform, NaN-safe
        x = pd.Series(x, dtype=float); m = x.notna()
        r = rankdata(x[m]); x[m] = norm.ppf((r - 0.5) / m.sum()); return x

    ph = pd.DataFrame(index=fam)
    for tid, name in KINGDOMS.items():
        ph[f"load_{name.replace(' ', '_')}"] = rint(np.log10(
            pd.Series({a: k2[a][tid][0] for a in accs}) / qc.total * 1e6 + 1)).reindex(fam)
    for _, t in targets.iterrows():
        sub = tg[tg.taxon == t["name"]].set_index("accession")
        key = t["name"].replace(" ", "_")
        ph[f"rpm_{key}"] = rint(np.log10(sub.rpm + 0.01)).reindex(fam)
        if sub.present.sum() >= 30:          # enough cases for a case/control scan
            ph[f"present_{key}"] = sub.present.reindex(fam)
    ph.to_csv(f"{OUT}/gwas_pheno.txt", sep="\t", header=False, index=False, na_rep="NA")
    with open(f"{OUT}/gwas_pheno_names.txt", "w") as fh:
        for i, c in enumerate(ph.columns, 1):
            fh.write(f"{i}\t{c}\n")

    # covariates: intercept, PC1-5 (SNP-Seek precomputed), log10 total reads (sequencing depth)
    pcs = pd.read_csv(f"{WORK}/gwas/pca_pruned_v2.1.eigenvec", sep=r"\s+")
    pcs = pcs.set_index("IID")[[f"PC{i}" for i in range(1, 6)]]
    cv = pd.DataFrame({"intercept": 1.0}, index=fam).join(pcs.reindex(fam))
    cv["log_depth"] = np.log10(tot.total.reindex(fam))
    cv = cv.fillna(cv.mean())
    cv.to_csv(f"{OUT}/gwas_covar.txt", sep="\t", header=False, index=False)

    batch = pd.Series({a: id_batch(a) for a in fam})
    dummies = pd.get_dummies(batch, dtype=float)
    dummies = dummies.drop(columns=batch.value_counts().idxmax())
    dummies = dummies.loc[:, dummies.sum() >= 20]          # tiny batches would make GEMMA unstable
    cv.join(dummies).to_csv(f"{OUT}/gwas_covar_batch.txt", sep="\t", header=False, index=False)
    print("phenotypes:", ", ".join(ph.columns))
print(f"wrote tables to {OUT}")

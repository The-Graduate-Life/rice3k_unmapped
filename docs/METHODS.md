# Methods

Detailed description of every step, its parameters, and why each choice was made. Script names refer
to `scripts/`; parameters live in `config.sh` and at the top of `scripts/03_aggregate.py`.

## 1. Data

| Item | Source | Notes |
|---|---|---|
| Read alignments | `s3://3kricegenome/Nipponbare/<ACC>.realigned.bam` + `.bai` | 3,024 accessions; reads aligned to Nipponbare IRGSP-1.0 **chr01–12 only** (no organelles), so organelle reads are among the unmapped reads; mean read length 83 bp |
| SNP genotypes | `s3://3kricegenome/snpseek-dl/3k-pruned-v2.1/` | LD-pruned set, ~1.0 M SNPs (496 k after MAF ≥ 0.05), PLINK bed/bim/fam |
| Principal components | same folder, `pca_pruned_v2.1.eigenvec` | PC1–5 used as covariates |
| Subpopulation | SNP-Seek K = 9 ADMIXTURE Q-matrix | assigned when the largest Q ≥ 0.6, otherwise `admix` (`metadata/accession_metadata.tsv`) |
| Country, variety name, SRA sample | `3K_list_sra_ids.txt` (same bucket) | |
| Known R/S genes | NCBI Gene, IRGSP-1.0 coordinates | `metadata/known_R_genes.tsv`; plotted on Manhattan plots only |

All IDs match 1:1 between BAMs, SNP `.fam`, Q-matrix and country list.

## 2. Retrieving only the discarded reads (`00_setup_refs.sh`, `01_process_accession.sh`)

A coordinate-sorted BAM stores unplaced unmapped reads in a final block that the index points to.
`samtools view -u <https URL> '*'` fetches **only that block** (≈ 7 s and 190 MB per accession;
575 GB for the panel instead of 27 TB). samtools caches each remote index in the working directory,
so every accession runs in its own folder that is deleted afterwards.

* **Total reads** per accession (the RPM denominator) come from `samtools idxstats` on the remote
  index alone (mapped + unmapped).
* **55 accessions** have zero unplaced-unmapped reads in the deposited BAM (unmapped pairs were
  stripped before upload; over-represented in Thailand 12/150 and China 17/473). They are logged to
  `logs/no_unmapped.txt` and given **missing** values, not zero, so they cannot pose as "healthy".
* **7 accessions** returned no data on four attempts and are treated the same way.
* Reads are converted to paired FASTQ with `samtools fastq -n`; singletons are dropped.

## 3. Removing host organelles and spike-in (`01_process_accession.sh`)

Reads are aligned with `minimap2 -ax sr` to a decoy set: rice chloroplast **NC_001320.1**, rice
mitochondrion **BA000029.3** and PhiX **NC_001422.1**. Pairs with both mates unmapped
(`samtools fastq -f 12 -F 0x900`) continue. Mitochondrial alignments are kept as
`results/mito/<ACC>.mito.bam` for follow-up work (cytoplasm / CMS typing). Chloroplast-like reads are
almost absent because they map to chloroplast insertions in the nuclear genome (NUPTs) and were
therefore never unmapped.

## 4. Taxonomic classification

* **Kraken2** with the prebuilt **PlusPF-16GB** index (`k2_pluspf_16gb_20250402`: RefSeq archaea,
  bacteria, viruses, plasmids, human, UniVec, protozoa, fungi; capped at 16 GB), `--paired
  --confidence 0.1 --report-minimizer-data`. The minimizer columns give *distinct minimizers*
  per taxon, the KrakenUniq-style guard against many reads hitting the same few k-mers.
* **Bracken** re-estimates species (`-l S`) and genus (`-l G`) abundance with the 75-mer
  distribution (read length is measured per accession and binned to 75/100/150).
* **Abundance** = reads per million *total* sequenced reads (RPM), so accessions sequenced to
  different depths are comparable.

## 5. Independent validation by alignment (`01_process_accession.sh`, `02_identity_check.sh`)

Cleaned reads are aligned with `minimap2 -ax sr` to a panel of reference genomes:

| Pathogen (disease) | Kraken2 taxid | Validation genome |
|---|---|---|
| *Xanthomonas oryzae* (bacterial leaf blight / streak) | 347 | PXO99A, CP000967.2 |
| *Burkholderia glumae* (bacterial panicle blight) | 337 | BGR1, GCF_000022645.2 |
| *Pyricularia oryzae* (rice blast) | 318829 | 70-15, GCF_000002495.2 |
| *Ustilaginoidea virens* (false smut) | 1159556 | UV-8b, GCF_000687475.1 |
| Rice tungro bacilliform virus | 3048298 | NC_001914.1 |
| *Rhizoctonia solani* (sheath blight) | 456999 | none, reported as *unvalidated* |
| *Fusarium fujikuroi* (bakanae) | 5127 | none, reported as *unvalidated* |

From the alignments (MAPQ ≥ 20, `samtools coverage`):

* **Breadth** = % of the genome covered by ≥ 1 read.
* **Expected breadth** if the same reads were spread uniformly: `1 − exp(−aligned_reads × 83 bp / genome_length)`.
* **Breadth ratio** = observed / expected. A real genome gives 0.6–1.0. Reads from a *related*
  species pile up on shared, conserved genes and give ~0.2–0.5.
* **Per-read identity** (`02_identity_check.sh`, run for every accession with a hit) =
  `1 − NM / (M + I + D)` for each primary alignment ≥ 50 bp. Reported per genome: mean identity and
  the fraction of reads at ≥ 97% and ≥ 99%. A real pathogen gives most reads ≥ 97%. A relative
  that resembles the reference along its whole genome passes the breadth ratio but fails here; for
  example, B001's *Xanthomonas*-like reads average 93.7%, and only 15% reach ≥ 97%.

## 6. Presence calls (`03_aggregate.py`)

| Call | Rule |
|---|---|
| **confirmed** | ≥ 20 Kraken2 clade reads **and** ≥ 50 aligned reads **and** breadth ≥ 1% **and** breadth ratio ≥ 0.6 **and** (if identity was measured) ≥ 50% of aligned reads at ≥ 97% identity |
| **congener** | Kraken2 and alignment signal (first three rules) but fails the ratio or identity rule, i.e. a related, non-target organism |
| **unvalidated** | taxon without a validation genome: ≥ 20 Kraken2 reads and ≥ 200 distinct minimizers |
| **absent** | everything else |

How the rules evolved (each rule was added because the data showed a failure):

1. **Initial rule** (Kraken2 reads + distinct minimizers + breadth ≥ 1%). The 98-accession pilot showed
   that an 8 kb virus can never reach 200 distinct minimizers (RTBV was missed), and that 1% breadth
   admitted *Xanthomonas* / *Burkholderia* relatives.
2. **Breadth ratio** replaced the minimizer rule for taxa with a genome. In the full panel this still
   let through 70 *B. glumae* "infections" with tens of thousands of aligned reads but only ~1% of
   them classified as *B. glumae* by Kraken2, i.e. a relative aligning along the whole genome.
3. **Identity** was added: 70 → 2 *B. glumae*, 275 → 218 *X. oryzae*, 33 → 3 tungro.

Low-depth libraries (< 1 M total reads) are flagged in `qc_summary.tsv` (`low_depth`) but not
removed. RPM from a tiny denominator is unstable; one such library exists in the panel.

## 7. Sample-batch structure

Accession IDs come in series (`IRIS_313-7xxx` … `-15xxx`, `CX…`, `B…`). Confirmed *X. oryzae*
prevalence differs by more than 20 percentage points between series (§ RESULTS), consistent with
seed grow-outs / DNA extractions done in different seasons or fields. The series (IRIS by
1000-number block; others by prefix) is used as a **batch proxy**:

* `prevalence_by_batch.tsv` reports prevalence per batch;
* `gwas_covar_batch.txt` adds batch dummy variables (largest batch = reference; batches with < 20
  members merged into the reference).

## 8. GWAS (`jobs/kinship.sh`, `jobs/gwas.sh`)

* **Genotypes:** LD-pruned v2.1 set; GEMMA filters MAF ≥ 0.05, missingness ≤ 0.1 → 495,922 SNPs.
* **Kinship:** centred relatedness matrix recomputed from the same bed (`gemma -gk 1 -maf 0.05`), so
  sample order is guaranteed to match.
* **Model:** linear mixed model, Wald test (`gemma -lmm 1`; the published base scans used `-lmm 4`,
  whose `p_wald` column is identical).
* **Covariates:** intercept, PC1–5, log10 total reads (sequencing depth); the *batch* set adds the
  batch dummies.
* **Phenotypes** (`gwas_pheno_names.txt`):
  * `load_<kingdom>`: Bacteria, Fungi, Viruses, *Homo sapiens* Kraken2 clade RPM;
    rank-inverse-normal transform of log10(RPM + 1).
  * `rpm_<pathogen>`: rank-inverse-normal of log10(RPM + 0.01), for all 7 target taxa.
  * `present_<pathogen>`: 0/1 confirmed call, only when ≥ 30 accessions are confirmed (here only
    *X. oryzae*, 218 cases).
  * Accessions without results are `NA`. Scans with < 200 non-missing values are skipped.
* **Threshold:** Bonferroni, 0.05 / 495,922 → −log10 p = 7.0. With 12 traits, a stricter bar is
  ≈ 8.1; both are reported. Genomic inflation λ = median(χ²) / 0.4549.

## 9. Locus robustness check (`05_locus_check.py`)

Significant SNPs within 1 Mb are grouped into loci; the best SNP per locus is extracted with
`plink --recode A` (0/1/2 copies of the A1 allele). Comparing homozygotes only (rice is inbred;
heterozygotes are rare), infection is cross-tabulated against genotype **within strata** of
batch, subpopulation, and batch × subpopulation. Strata are pooled with the **Mantel–Haenszel odds
ratio** and tested with the **Cochran–Mantel–Haenszel χ²** (1 df, no continuity correction).
Strata without both genotypes or both outcomes carry no information and are skipped.

This is a robustness check on the discovery data, not a replication: odds ratios are inflated by
the winner's curse.

## 10. Software

samtools ≥ 1.21, minimap2, seqkit, Kraken2 ≥ 2.1.3, Bracken ≥ 2.9, PLINK 1.90b6.21, GEMMA ≥ 0.98.5,
Python 3.11 (pandas, numpy, scipy, matplotlib). Exact versions of a run:
`results/reuse3k.lock.yml` and `results/RUN_INFO.txt`.

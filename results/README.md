# results/

Small, shareable outputs of the published run. Large intermediates stay in `$WORK` on the cluster.

| Path | Contents |
|---|---|
| `figures/01_unmapped_composition.png` | what the discarded reads consist of, per accession |
| `figures/02_top_species.png` | most abundant microbes (mean reads per million) |
| `figures/03_pathogens_by_subpop.png` | pathogen abundance by subpopulation; red = confirmed |
| `figures/04_kraken_vs_alignment.png` | Kraken2 abundance vs genome breadth (validation) |
| `figures/05_batch_checks.png` | abundance vs sequencing depth and read length |
| `figures/gwas_<trait>[_batch].png` | Manhattan + QQ plot per GWAS scan (`_batch` = with batch covariates) |
| `tables/`, `gwas/`, `RUN_INFO.txt`, `reuse3k.lock.yml` | written by `scripts/06_collect_results.sh` at the end of a run (see FOLLOWME §8) |

Column definitions are in `docs/METHODS.md`; the numbers are discussed in `docs/RESULTS.md`.

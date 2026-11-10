# FOLLOWME: reproduce the project step by step

This guide takes you from nothing to the published figures and tables. Every step tells
you **what to type**, **what you should see**, and **how long it takes**. If something goes wrong, see
[§9 Troubleshooting](#9-troubleshooting).

**Short version**

```bash
git clone git@github.com:The-Graduate-Life/rice3k_unmapped.git && cd rice3k_unmapped
bash install.sh
bash run_pipeline.sh
```

---

## 0. What you need

| Need | Details |
|---|---|
| A GitHub account with an SSH key | to clone the repository (step 1 shows how to set up the key) |
| A Linux machine with bash | an HPC cluster for the full run (PBS or SLURM); any Linux/macOS laptop for the smoke test |
| Internet from the compute nodes | to stream reads from `3kricegenome.s3.amazonaws.com` (no account or key needed). If your nodes are offline, see §9. |
| Disk | ~45 GB of scratch space (set by `WORK` in `config.sh`) |
| Memory | 24 GB per processing job, 16 GB for GWAS |
| Time | about one night on a cluster; 10 minutes for the smoke test |
| Nothing else | `install.sh` installs every tool, including conda if you have none |

All commands below are run from the repository folder.

---

## 1. Get the code (over SSH)

**1a. Check that your machine can talk to GitHub over SSH** (once per machine, e.g. once on the cluster
login node and once on your laptop):

```bash
ssh -T git@github.com
```

You want: `Hi <your-username>! You've successfully authenticated, but GitHub does not provide shell access.`
The first time, answer `yes` to the "authenticity of host" question.

**1b. If it says `Permission denied (publickey)`, create a key and register it with GitHub:**

```bash
ssh-keygen -t ed25519 -C "your_email@example.com"   # press Enter at every prompt (default file, optional passphrase)
cat ~/.ssh/id_ed25519.pub                            # copy the whole line it prints
```

On GitHub: **Settings → SSH and GPG keys → New SSH key**, paste the line, save. Then repeat 1a.

**1c. Clone the repository:**

```bash
git clone git@github.com:The-Graduate-Life/rice3k_unmapped.git
cd rice3k_unmapped
```

On a cluster, put the repository in your home or project space. The large data go to a separate
scratch folder (step 3).

> **Editing on Windows?** Keep LF line endings (the included `.gitattributes` does this for git).
> A script saved with Windows line endings fails with `$'\r': command not found`.

---

## 2. Install the software (≈ 5–15 min, once)

```bash
bash install.sh
```

What it does: finds conda (or installs Miniforge into `~/miniforge3` if there is none), then creates
the environment **`reuse3k`** from `env/reuse3k.yml`. It ends with a tool check:

```
>> tool check
   samtools  ok
   minimap2  ok
   ...
   gemma     ok
```

Every line must say `ok`. Re-running `install.sh` is safe; it updates the environment.

---

## 3. Configure (≈ 2 min; often nothing to change)

Open `config.sh`. There are only two things you may need to change:

1. **`WORK`**, the scratch folder for downloads and results. The default is `$SCRATCH/rice3k_work`
   (or `~/rice3k_work` if `$SCRATCH` is not set). It needs ~45 GB.
2. **The `RES_*` lines**, the queue (PBS) or partition (SLURM), CPUs, memory and wall time for each
   stage. The defaults fit the **Alabama Supercomputer**. Format: `"queue|cpus|mem_GB|walltime"`.

   Example for a SLURM cluster whose partition is called `compute`:
   ```bash
   export RES_PROCESS="compute|8|24|08:00:00"
   export RES_GWAS="compute|4|16|24:00:00"
   ```
   Leave the queue empty (`"|8|24|08:00:00"`) to use the cluster's default.

Any setting can also be overridden on the command line without editing the file:
`WORK=/big/disk/rice3k bash run_pipeline.sh`.

---

## 4. (Optional) Smoke test on 3 varieties (≈ 10 min, ~2 GB)

Proves the whole chain works before you spend cluster time. Runs anywhere, no scheduler needed.

```bash
bash tests/smoke_test.sh
```

It uses a small viral database and skips the GWAS. Expected ending:

```
== checks
  ok   3 accessions processed
  ok   targets.tsv has 3 x 7 rows
  ok   calls are valid labels
  ok   5 QC/prevalence figures
  ok   mito BAMs written
  ok   results collected
SMOKE TEST PASSED
```

It writes only to `./smoke_work/`, which you can delete afterwards.

---

## 5. (Optional) Check a compute node (≈ 1 min plus queue time)

```bash
bash run_pipeline.sh --check-node
cat logs/rice3k_check.log          # PBS writes the log when the job ends
```

You want `INTERNET_FROM_COMPUTE=yes`. If it says `no`, follow §9, "Compute nodes have no internet".

---

## 6. Run everything (one command)

```bash
bash run_pipeline.sh
```

It detects your scheduler and submits **all stages at once**, chained by job dependencies, so you
can log out:

```
scheduler: pbs | list: .../metadata/accessions_all.txt (3024 accessions) | WORK: ... | Kraken2 DB: pluspf16
setup     : 119700.asax-pbs1
process   : 121 jobs
aggregate : 119822.asax-pbs1
identity  : 6 jobs
aggregate : 119829.asax-pbs1
kinship   : 119830.asax-pbs1
gwas      : 12 jobs
report    : 119843.asax-pbs1
All jobs submitted. Check progress any time with:  bash run_pipeline.sh --status
```

What each stage does:

| # | Stage | Jobs | Typical time | What happens |
|---|---|---|---|---|
| 1 | setup | 1 | 30–90 min | downloads decoys, 5 pathogen genomes, Kraken2 DB (12 GB), SNP set (0.8 GB); counts reads per variety from the BAM indexes |
| 2 | process | 121 (25 varieties each) | ~1–2 h total | streams each variety's unmapped reads, removes organelle/PhiX reads, runs Kraken2/Bracken and the alignment check, deletes the reads |
| 3 | aggregate | 1 | 5 min | first presence calls, `target_hits.tsv` |
| 4 | identity | 6 | ~1 h | per-read identity for every hit (separates pathogens from relatives) |
| 5 | aggregate | 1 | 5 min | final calls, GWAS phenotype and covariate files |
| 6 | kinship + GWAS | 1 + 12 | 0.5 h + 1–6 h | 12 traits × 2 covariate sets = 24 scans |
| 7 | report | 1 | 10 min | figures, GWAS summary, locus check, copy into `results/` |

**Only want part of it?**

| Goal | Command |
|---|---|
| The 99-variety pilot | `bash run_pipeline.sh --list metadata/accessions_pilot100.txt --no-gwas` |
| Disease survey only, no GWAS | `bash run_pipeline.sh --no-gwas` |
| See the job graph without submitting | `bash run_pipeline.sh --dry-run` |
| Force a scheduler | `--scheduler pbs`, `--scheduler slurm` or `--scheduler local` |

---

## 7. Monitor

```bash
bash run_pipeline.sh --status
```

```
accessions in list        : 3024
processed (qc files)      : 2962
skipped, no unmapped block: 55
failed (logs/failed.txt)  : 7
identity checks done      : 565
GWAS scans done           : 24
figures                   : 31
```

**The numbers to expect at the end** (from the published run):

| Counter | Expected | Meaning |
|---|---|---|
| processed | **2,962** | varieties with results |
| skipped, no unmapped block | **55** | deposited BAMs whose unmapped reads were removed before upload; listed in `logs/no_unmapped.txt` |
| failed | **7** | IRIS_313-11643, -11657, -11669, -11671, -11723, -11801, -11812: S3 returns no data for them |
| identity checks | **565** | every variety with a confirmed or congener call |
| GWAS scans | **24** | 12 traits × (base, batch) |

2,962 + 55 + 7 = 3,024. Any other failures are usually temporary network errors: just run
`bash run_pipeline.sh --from process` again. Finished varieties are skipped.

Logs are in `logs/<job name>.log`. Under PBS they appear only when each job ends.

---

## 8. Check the results

When the `report` job has finished, the repository's `results/` folder contains:

```
results/
  RUN_INFO.txt              date, host, database, tool versions
  reuse3k.lock.yml          exact conda package versions used
  figures/                  01–05 QC/prevalence figures + one Manhattan/QQ per GWAS scan
  tables/
    qc_summary.tsv          per variety: total, unmapped, organelle, classified reads, low_depth flag
    targets.tsv             per variety x pathogen: reads, RPM, breadth, identity, call
    target_hits.tsv         the non-absent rows of targets.tsv
    prevalence_by_*.tsv     by subpopulation, country, sample batch
    gwas_summary.tsv        per scan: lambda, number of significant SNPs, top SNP
    species_rpm.tsv.gz      variety x species abundance matrix
  gwas/
    gwas_hits_<trait>.tsv   SNPs above the Bonferroni line
    locus_check_present_Xanthomonas_oryzae.tsv   robustness of each blight locus
```

Compare with the published numbers:

```bash
conda activate reuse3k
python -c "import pandas as pd; t=pd.read_csv('results/tables/targets.tsv', sep='\t'); print(pd.crosstab(t.taxon, t.call))"
column -t results/tables/gwas_summary.tsv | grep present_Xan
column -t results/gwas/locus_check_present_Xanthomonas_oryzae.tsv
```

| Check | Published value |
|---|---|
| *X. oryzae* confirmed / congener | 218 / 296 |
| *B. glumae* confirmed / congener | 2 / 170 |
| Tungro confirmed / congener | 3 / 30 |
| *P. oryzae* confirmed | 7 |
| `present_Xanthomonas_oryzae` λ | 0.99 |
| significant loci (base scan) | chr5 ~7.06 Mb (4 SNPs), chr6 4,641,440, chr7 29,223,762 |
| MH odds ratio within batch × subpop | ≈ 3.0 / 3.6 / 2.6 |

Small differences are possible if a public file changed or a download failed; the RUN_INFO and lock
files tell you what was used. More detail: [docs/RESULTS.md](docs/RESULTS.md).

---

## 9. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `git@github.com: Permission denied (publickey)` when cloning | no SSH key registered with GitHub on this machine | follow step 1b, then `ssh -T git@github.com` |
| `ssh: connect to host github.com port 22: Connection timed out` | the cluster/firewall blocks port 22 | use SSH over port 443: add `Host github.com` / `Hostname ssh.github.com` / `Port 443` to `~/.ssh/config` |
| `conda: command not found` | no conda on the cluster or not on PATH | `bash install.sh` (installs Miniforge); config.sh finds `~/miniforge3` automatically |
| `EnvironmentNameNotFound: reuse3k` inside jobs | a site `module load` switched to another conda | set `CONDA_ENV` in config.sh to the full path, e.g. `$HOME/miniforge3/envs/reuse3k` |
| `sbatch: command not found` / `qsub: command not found` | other scheduler | auto-detect normally handles it; or `--scheduler pbs` / `slurm` |
| `Job violates queue and/or server resource limits` (PBS) | queue limits: on ASC every queue has a 1 h minimum; medium ≤ 16 GB; small ≤ 4 GB | adjust the `RES_*` lines in config.sh (`qstat -Qf` lists the limits) |
| `PBS: job killed: walltime … exceeded` | GWAS scan too slow | raise the wall time in `RES_GWAS`, then `bash run_pipeline.sh --from gwas` (finished scans are skipped) |
| `[E::idx_find_and_load] Could not retrieve index` / `retrieval of region failed` for a few varieties | transient S3/network error, or a broken deposited file | re-run `--from process`; persistent failures are the 7 listed in §7 |
| a variety listed in `logs/no_unmapped.txt` | the deposited BAM has no unmapped block | expected for 55 varieties; they get NA, not zero |
| `$'\r': command not found` | Windows line endings | `sed -i 's/\r$//' <file>` (or `dos2unix`) |
| **Compute nodes have no internet** | site policy | run setup on a login node (`bash scripts/00_setup_refs.sh`), fetch reads there (`bash scripts/utils/fetch_on_login_node.sh metadata/accessions_all.txt 4`, ~575 GB in `$WORK/fastq`), then `bash run_pipeline.sh --from process --classify-only` |
| dependency list rejected (very long PBS `depend=` string) | site limit on dependency length | submit in two steps: `bash run_pipeline.sh --no-gwas` then, when finished, `bash run_pipeline.sh --from gwas` |
| `gwas_hits_*.tsv: No such file` | looked in `figures/` | hit lists are in `$WORK/results/tables/` and `results/gwas/` |

Re-running is always safe: every stage skips work that is already complete (finished varieties,
identity files, GWAS scans, kinship). To redo a stage from scratch, delete its outputs in `$WORK`
first (e.g. `rm $WORK/gwas/output/lmm_*` to redo all scans).

---

## 10. Where everything lives

| Location | Contents | Size |
|---|---|---|
| repository (`PROJ`) | code, metadata, small results | < 20 MB |
| `$WORK/db`, `$WORK/refs`, `$WORK/gwas` | Kraken2 DB, reference genomes, SNP set, kinship, GEMMA outputs | ~21 GB |
| `$WORK/results/kraken`, `bracken`, `patho`, `identity`, `qc` | per-variety outputs | ~2 GB |
| `$WORK/results/mito` | mitochondrial reads per variety (BAM), for follow-up work on cytoplasm types | ~20 GB |
| `$WORK/fastq` | reads in flight; each variety's folder is deleted when it finishes | small |

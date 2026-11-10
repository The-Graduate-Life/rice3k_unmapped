#!/bin/bash
# Stage 1 (one-time downloads; needs internet). Normally run by run_pipeline.sh as a job.
#   bash scripts/00_setup_refs.sh
# Env: LIST (accessions to count reads for; default all), SKIP_GWAS_DATA=1 (smoke test).
# Everything is skipped if already present, so re-running is cheap.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_env
mkdir -p "$WORK"/{db,refs,gwas,fastq,results/{kraken,bracken,qc,mito,patho,identity,tables,figures}}
cd "$WORK/refs"

efetch() { curl -fsS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=$1"; }
dsfetch() {  # NCBI Datasets genome FASTA
  curl -fsSL -o "$1.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$1/download?include_annotation_type=GENOME_FASTA"
  python -c "import zipfile,sys; z=zipfile.ZipFile(sys.argv[1]); sys.stdout.buffer.write(b''.join(z.read(n) for n in z.namelist() if n.endswith('.fna')))" "$1.zip" > "$1.fa" && rm "$1.zip"
}

# 1) Decoy set removed before Kraken2: rice organelles (absent from the 3K BAM reference) + PhiX spike-in
[ -s decoy.fa ] || { efetch NC_001320.1 > rice_cp.fa; efetch BA000029.3 > rice_mt.fa; efetch NC_001422.1 > phix.fa
                     cat rice_cp.fa rice_mt.fa phix.fa > decoy.fa; }

# 2) Pathogen genomes for hit validation (breadth of coverage)
[ -s CP000967.2.fa ]      || efetch CP000967.2 > CP000967.2.fa        # Xanthomonas oryzae pv. oryzae PXO99A
[ -s NC_001914.1.fa ]     || efetch NC_001914.1 > NC_001914.1.fa      # Rice tungro bacilliform virus
[ -s GCF_000002495.2.fa ] || dsfetch GCF_000002495.2                  # Pyricularia oryzae 70-15
[ -s GCF_000022645.2.fa ] || dsfetch GCF_000022645.2                  # Burkholderia glumae BGR1
[ -s GCF_000687475.1.fa ] || dsfetch GCF_000687475.1                  # Ustilaginoidea virens UV-8b
# Concatenated pathogen panel + contig->genome map (for samtools coverage aggregation)
if [ ! -s pathogens.fa ]; then
  : > pathogens.fa; : > pathogen_contigs.tsv
  for g in CP000967.2 NC_001914.1 GCF_000002495.2 GCF_000022645.2 GCF_000687475.1; do
    cat "$g.fa" >> pathogens.fa
    grep '>' "$g.fa" | sed 's/>//;s/ .*//' | awk -v g="$g" '{print $1"\t"g}' >> pathogen_contigs.tsv
  done
  minimap2 -d pathogens.mmi pathogens.fa
fi

# 3) Kraken2 database (prebuilt, Langmead lab index collection)
cd "$WORK/db"
if [ ! -s "$KDB/hash.k2d" ]; then
  mkdir -p "$KDB"
  case "$KDB_FLAVOR" in
    pluspf16) TGZ=k2_pluspf_16gb_20250402.tar.gz ;;   # ~12 GB download
    pluspf)   TGZ=k2_pluspf_20251015.tar.gz ;;        # ~83 GB download
    viral)    TGZ=k2_viral_20251015.tar.gz ;;         # 0.5 GB, plumbing tests only
  esac
  curl -fSL -C - -o "$TGZ" "https://genome-idx.s3.amazonaws.com/kraken/$TGZ"
  tar -xzf "$TGZ" -C "$KDB" && rm "$TGZ"
fi

# 4) GWAS inputs from 3K RG / SNP-Seek (LD-pruned SNP set v2.1, ~1.0M SNPs, ~0.8 GB + PCA)
#    Kinship is recomputed from this exact .bed in jobs/kinship.sh so sample order is guaranteed.
cd "$WORK/gwas"
[ "${SKIP_GWAS_DATA:-0}" = 1 ] || for f in pruned_v2.1.bed pruned_v2.1.bim pruned_v2.1.fam pruned_v2.1.pca/pca_pruned_v2.1.eigenvec; do
  [ -s "$(basename "$f")" ] || curl -fSL -o "$(basename "$f")" "$S3_BASE/snpseek-dl/3k-pruned-v2.1/$f"
done

# 5) Total read counts per accession (reads only the ~1 MB .bai index per accession) -> normalisation
#    denominator. Only accessions not yet in the table are queried, so a later, longer LIST extends it.
cd "$WORK/results/qc"
[ -s total_reads.tsv ] || echo -e "accession\tmapped\tunmapped_placed\tunmapped_unplaced\ttotal" > total_reads.tsv
cut -f1 total_reads.tsv | sort > .have.txt
sort "${LIST:-$PROJ/metadata/accessions_all.txt}" | comm -23 - .have.txt | while read -r acc; do
  [ -z "$acc" ] && continue
  samtools idxstats "$BAM_URL/$acc.realigned.bam" 2>/dev/null \
    | awk -v a="$acc" '{m+=$3; u+=$4; if($1=="*") s=$4} END{if(NR) print a"\t"m"\t"u-s"\t"s"\t"m+u}' >> total_reads.tsv
  rm -f "$acc.realigned.bam.bai"     # samtools caches the remote index in CWD
done
rm -f .have.txt
echo "read totals: $(( $(wc -l < total_reads.tsv) - 1 )) accessions"
echo "Setup done. Kraken2 DB: $KDB"

#!/usr/bin/env bash
set -eo pipefail
source /kaggle/working/project13/scripts/config.sh
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq_cli
D=$PROJ/data

# 1. ENA run metadata and sample table (all 18 runs)
curl -s --retry 3 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA263263&result=read_run&fields=run_accession,sample_accession,sample_title,experiment_title,library_layout,library_strategy,read_count,fastq_ftp&format=tsv" \
  -o $D/ena_runs_PRJNA263263.tsv
N_LINES=$(wc -l < $D/ena_runs_PRJNA263263.tsv)
if [ "$N_LINES" -ne 19 ]; then echo "expected 19 lines (header + 18 runs), got $N_LINES"; exit 1; fi
printf 'run\tsample_title\tstage\ttreatment\treplicate\tread_count\n' > $D/samples.tsv
awk -F'\t' 'NR>1 {split($3,a,"_"); t=a[3]; if (t=="un") t="untreated"; print $1"\t"$3"\t"a[2]"\t"t"\t"a[4]"\t"$7}' \
  $D/ena_runs_PRJNA263263.tsv | sort -k3,3 -k4,4 -k5,5 >> $D/samples.tsv
column -t -s $'\t' $D/samples.tsv

# 2. Reference: NCBI RefSeq GCF_000002765.6
cd $TMPD/ref
N=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/765
DIR=$(curl -s -m 60 $N/ | grep -o 'GCF_000002765\.6_[A-Za-z0-9_.-]*' | sort -u | head -n 1)
echo "NCBI directory: [$DIR]"
if [ -z "$DIR" ]; then echo "could not find the assembly directory"; exit 1; fi
curl -sL -m 300 -o genomic.fna.gz $N/$DIR/${DIR}_genomic.fna.gz
curl -sL -m 300 -o genomic.gff.gz $N/$DIR/${DIR}_genomic.gff.gz
gzip -t genomic.fna.gz
gzip -t genomic.gff.gz
zcat genomic.fna.gz > genome.fa
zcat genomic.gff.gz > annotation.gff
echo "$DIR" > $PROJ/results/reference_release.txt
cp annotation.gff $D/annotation.gff
gffread annotation.gff -T -o genes.gtf
sed -i 's/gene_id "gene-/gene_id "/' genes.gtf
seqkit stats genome.fa
echo "sequences: $(grep -c '^>' genome.fa) (expected 14)"
echo "genes with exons: $(awk -F'\t' '$3=="exon"' genes.gtf | grep -o 'gene_id "[^"]*"' | sort -u | wc -l) (Day 2: 5611)"
echo "kelch13 PF3D7_1343700 lines in GTF: $(grep -c 'PF3D7_1343700' genes.gtf)"

# 3. HISAT2 index with known splice sites
hisat2_extract_splice_sites.py genes.gtf > splice_sites.txt
hisat2_extract_exons.py genes.gtf > exons.txt
wc -l splice_sites.txt exons.txt
hisat2-build -p $NT --ss splice_sites.txt --exon exons.txt genome.fa $TMPD/index/pf3d7 2>&1 | tail -n 3
ls -lh $TMPD/index

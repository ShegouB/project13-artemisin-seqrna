#!/usr/bin/env bash
set -eo pipefail
source /kaggle/working/project13/scripts/config.sh
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq_cli
FC_OPTS="${1:--C}"
cd $TMPD
TSV=$PROJ/data/ena_runs_PRJNA263263.tsv
R=$PROJ/results
mkdir -p $R/hisat2 $R/fastp
tail -n +2 $PROJ/data/samples.tsv | cut -f1 | while read RUN; do
  if [ -f bam/${RUN}.done ]; then echo "skip $RUN"; continue; fi
  L=$(awk -F'\t' -v r=$RUN 'NR==1{for(i=1;i<=NF;i++) if($i=="fastq_ftp") c=i} $1==r{print $c}' $TSV)
  U1=$(echo "$L" | tr ';' '\n' | grep '_1.fastq.gz')
  U2=$(echo "$L" | tr ';' '\n' | grep '_2.fastq.gz')
  curl -sL --retry 3 -o fastq/${RUN}_1.fastq.gz https://$U1
  curl -sL --retry 3 -o fastq/${RUN}_2.fastq.gz https://$U2
  gzip -t fastq/${RUN}_1.fastq.gz; gzip -t fastq/${RUN}_2.fastq.gz
  fastp -i fastq/${RUN}_1.fastq.gz -I fastq/${RUN}_2.fastq.gz \
    -o fastq/${RUN}_t1.fastq.gz -O fastq/${RUN}_t2.fastq.gz \
    --detect_adapter_for_pe -w 4 --length_required 30 \
    -j $R/fastp/${RUN}.json -h $R/fastp/${RUN}.html 2> $R/fastp/${RUN}.log
  hisat2 -p $NT --rna-strandness RF --max-intronlen 5000 -x index/pf3d7 \
    -1 fastq/${RUN}_t1.fastq.gz -2 fastq/${RUN}_t2.fastq.gz \
    --summary-file $R/hisat2/${RUN}.summary 2> /dev/null \
    | samtools sort -@ $NT -o bam/${RUN}.sorted.bam -
  samtools index bam/${RUN}.sorted.bam
  touch bam/${RUN}.done
  rm -f fastq/${RUN}_*.fastq.gz
  echo "done $RUN"
done
featureCounts -T $NT -p --countReadPairs $FC_OPTS -s 2 -t exon -g gene_id \
  -a ref/genes.gtf -o $R/counts_all.txt \
  $(tail -n +2 $PROJ/data/samples.tsv | cut -f1 | sed 's#^#bam/#; s#$#.sorted.bam#')
multiqc $R -o $R/multiqc -f

#!/usr/bin/env bash
set -eo pipefail
mkdir -p /kaggle/temp/mm
cd /kaggle/temp/mm
if [ ! -x bin/micromamba ]; then
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
fi
export MAMBA_ROOT_PREFIX=/kaggle/temp/mm/root
export PATH=/kaggle/temp/mm/bin:$PATH
micromamba --version
if [ ! -d $MAMBA_ROOT_PREFIX/envs/rnaseq_cli ]; then
  micromamba create -y -n rnaseq_cli -c conda-forge -c bioconda --strict-channel-priority \
    hisat2 samtools subread fastp multiqc seqkit gffread
fi
if [ ! -d $MAMBA_ROOT_PREFIX/envs/rnaseq_r ]; then
  micromamba create -y -n rnaseq_r -c conda-forge -c bioconda --strict-channel-priority \
    r-base r-tidyverse r-rmarkdown r-knitr r-pheatmap r-ggrepel pandoc \
    bioconductor-deseq2 bioconductor-edger bioconductor-limma bioconductor-apeglm \
    bioconductor-clusterprofiler bioconductor-enrichplot bioconductor-fgsea \
    bioconductor-enhancedvolcano bioconductor-complexheatmap
fi
cat > /kaggle/working/project13/scripts/config.sh << 'CFG'
export PROJ=/kaggle/working/project13
export TMPD=/kaggle/temp/project13
export MAMBA_ROOT_PREFIX=/kaggle/temp/mm/root
export PATH=/kaggle/temp/mm/bin:$PATH
export NT=$(nproc)
mkdir -p $PROJ/data $PROJ/results $PROJ/report $TMPD/ref $TMPD/fastq $TMPD/bam $TMPD/index
CFG
source /kaggle/working/project13/scripts/config.sh
micromamba env export -n rnaseq_cli > $PROJ/results/env_rnaseq_cli.yml
micromamba env export -n rnaseq_r > $PROJ/results/env_rnaseq_r.yml
echo "=== versions ==="
micromamba run -n rnaseq_cli bash -c 'hisat2 --version | head -n 1; featureCounts -v 2>&1 | head -n 2; samtools --version | head -n 1; fastp --version 2>&1; multiqc --version; seqkit version; gffread --version'
micromamba run -n rnaseq_r Rscript -e 'for (p in c("DESeq2","edgeR","apeglm","clusterProfiler","enrichplot","EnhancedVolcano","ComplexHeatmap","rmarkdown","ggplot2")) cat(sprintf("%-18s %s\n", p, as.character(packageVersion(p)))); cat(R.version.string, "\n")'

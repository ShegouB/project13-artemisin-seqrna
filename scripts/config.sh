export PROJ=/kaggle/working/project13
export TMPD=/kaggle/temp/project13
export MAMBA_ROOT_PREFIX=/kaggle/temp/mm/root
export PATH=/kaggle/temp/mm/bin:$PATH
export NT=$(nproc)
mkdir -p $PROJ/data $PROJ/results $PROJ/report $TMPD/ref $TMPD/fastq $TMPD/bam $TMPD/index

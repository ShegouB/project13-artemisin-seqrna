# Project 13: transcriptomic response of *Plasmodium falciparum* to dihydroartemisinin

A reproducible re-analysis of a public RNA-seq dataset (GEO GSE62132, ENA project
PRJNA263263) in which *P. falciparum* K1 parasites were exposed to dihydroartemisinin
(DHA), to a DMSO vehicle, or to nothing, at three life-cycle stages.

Block 6 (Transcriptomics and Functional Genomics), Week 13, Saturday project.
Author: B. O. Djagou (GitHub: ShegouB).

## Question and context

Artemisinin-based drugs are the backbone of malaria treatment and resistance threatens
them. The project brief proposes that knowing which genes are dysregulated in
artemisinin-treated parasites can reveal resistance mechanisms before they appear as
clinical treatment failures, which matters for malaria programmes in West Africa. The
questions here are which genes and processes respond to DHA, how reliable those results
are with two replicates, and whether genes of the kelch13 pathway change in expression.

## Main findings

| Stage | DHA vs vehicle | Standing |
|---|---|---|
| Trophozoite | 1,703 of 4,497 tested genes at padj < 0.05 (923 up, 780 down); 973 consensus genes (593 up, 380 down) | Supported |
| Ring | 413 genes (4-sample model), 81 (3-group model), 7 consensus genes; depends on the control used | Tentative |
| Schizont | 0 genes | Negative control |

Consensus gene = significant against both the vehicle and the untreated control, with the
same sign.

- **Trophozoite, lower in DHA:** translation machinery (102 of 265 annotated down genes
  carry the term "translation", adjusted p 4.6e-31; ribosome terms strongest in the
  molecular-function and cellular-component ontologies) and proteasome core subunits.
- **Trophozoite, higher in DHA:** merozoite invasion genes (MTRAP, EBA181, EBA175, EBA140,
  AMA1, SUB2, RIPR, CyRPA), signalling kinases (PKG, CDPK2, CDPK5, CDPK7, STK2, TKL3, TKL4,
  PKB) and the exported-protein system (PRESAN domain proteins, FIKK kinases, PTEX150,
  SBP1, SEMP1, HSP101).
- **Stage profile:** DHA-up genes are enriched among genes that are low at the trophozoite
  stage in untreated parasites and higher in both ring and schizont (36% of up genes,
  enrichment 2.7); DHA-down genes are enriched among genes that peak at trophozoite
  (14%, enrichment 2.6). This is neither a simple delay nor a simple advance. It fits a
  flattening of the trophozoite-specific profile, a mixture of ring-like and schizont-like
  parasites, or a direct effect on stage-regulated genes, and bulk data with two
  replicates cannot tell these apart.
- **kelch13 pathway (16 genes: K13, KIC1 to KIC10, UBP1, EPS15, AP2-MU, MyoF, MCA2):** 12 of
  16 are significantly higher in trophozoites (K13: log2FC 0.94, padj 1e-6) and the set is
  shifted beyond expression-matched random sets (median Wald 3.97 vs 0.09). After removing
  the part of the DHA effect predicted by the stage axes, the set median log2FC falls from
  0.82 to 0.15 (random sets -0.04, p = 0.18). Ring and schizont show no signal. These data
  do not support a K13-specific transcriptional response, and they cannot address
  resistance, which is mostly mediated by K13 mutations acting at the protein level.

## Reproducibility

The whole analysis was rebuilt in a fresh notebook from public accessions with saved
scripts. Assigned fragment counts for all 18 samples were identical to an earlier,
independent run, although the environment resolved to newer tool versions.

## Data

- 18 paired-end, strand-specific runs: 3 stages (ring, trophozoite, schizont) x 3
  treatments (untreated, 0.1% DMSO vehicle, 500 nM DHA) x 2 replicates; strain K1; about
  1.1 to 2.0 million fragments counted per sample.
- Run table: `data/ena_runs_PRJNA263263.tsv`; derived sample table: `data/samples.tsv`.
- Reference: NCBI RefSeq assembly GCF_000002765.6 (14 nuclear chromosomes, 23,292,622 bp;
  5,611 genes with exons), fetched by script from the NCBI genomes FTP.
- GO annotation: UniProt, taxon 36329, fetched by REST (`data/uniprot_pf3d7_go2.tsv`).
  UniProt content changes over time, so this file is kept in the repository.

## Pipeline and run order

All scripts were run in a Kaggle notebook (Internet: On) and share `scripts/config.sh`.

| Order | Script | What it does |
|---|---|---|
| 1 | `00_setup.sh` | Installs micromamba and two environments (`rnaseq_cli`, `rnaseq_r`), writes `config.sh`, exports environment files |
| 2 | `01_prepare.sh` | ENA metadata and sample table, NCBI reference, GTF conversion (strips the `gene-` prefix), HISAT2 index with known splice sites |
| 3 | `02_align_count.sh` | Downloads FASTQ from ENA, trims (fastp), aligns (HISAT2), sorts, counts all runs together (featureCounts); a sample is skipped only if its completion marker exists |
| 4 | `03_deseq2.R` | Composition QC, PCA, per-stage DESeq2 (primary 2-group and 3-group models), consensus genes, ribosomal gene test, VST export |
| 5 | `04_kelch13.R` | Finds kelch13-pathway candidates by name in the annotation |
| 6 | `05_kelch13_sets.R` | Set-level tests against expression-matched random sets |
| 7 | `06_kelch13_stage_adjusted.R` | Stage-adjusted set test, ranked lists and top-20 table for Week 14 |
| 8 | `07_figures.R` | Volcano, MA and heatmap figures |
| 9 | `08_go_sources.R`, `09_go_sources2.R` | Exploratory comparison of GO annotation sources (coverage) |
| 10 | `10_go_enrichment.R` | GO over-representation analysis and dotplots |
| 11 | `11_stage_overlap.R` | Overlap of DHA gene lists with stage-defined gene sets; genes behind selected terms |
| 12 | `report/project13_report.Rmd` | Knits the HTML report from the result files |

Key parameters: HISAT2 `--rna-strandness RF --max-intronlen 5000`; featureCounts
`-p --countReadPairs -C -s 2 -t exon -g gene_id`; DESeq2 genes with at least 10 counts in at
least 2 samples, alpha 0.05, apeglm shrinkage; GO universe = the 4,497 genes tested in
trophozoite, terms of 10 to 500 genes, BH adjustment per ontology, redundancy reduction at
50% gene overlap.

### How to reproduce

1. Create a new Kaggle notebook with Internet on and run the scripts above in order
   (`bash scripts/00_setup.sh` and so on). The first run installs the environments, which
   takes several minutes.
2. Fast-changing external sources (UniProt, ENA, NCBI) are downloaded at run time. The NCBI
   assembly and ENA runs are fixed accessions; UniProt content may differ from the
   committed file.
3. Kaggle does not keep `/kaggle/temp` between sessions, so download the checkpoint
   archives (counts, tables, figures, logs) after each phase.

## Quality control

- Library strandedness is reverse (dUTP): on a 500,000-pair test, `-s 2` assigned 346,820
  pairs against 11,471 for `-s 1`; ambiguous assignments fell from 4,788 to 34 with
  strand-aware counting.
- Adapter read-through: trimming raised the overall alignment rate from 83.84% to 92.82%
  and concordant pairs from 65.6% to 86.9% on the test sample.
- Multi-mapped fragments are 3.5% to 12% per sample; unmapped fragments are 2% to 5% in ring
  and trophozoite and about 8% to 9.5% in schizont (cause not investigated).
- Ring vehicle libraries carry more rRNA (0.46% and 0.71%) than ring DHA (0.10% and 0.12%);
  GEO lists one extraction protocol for all samples, so the cause is unknown.
- DESeq2 size factors re-estimated without ribosomal and rRNA genes changed DEG counts by
  1% or less.

## Deviations from the brief

| Brief | This analysis | Reason |
|---|---|---|
| Treated vs untreated, 3 replicates | DHA vs vehicle (primary) with untreated as a second control, 2 replicates | The dataset has 2 replicates; DMSO is the matched control |
| Artemisinin | Dihydroartemisinin, 500 nM | The dataset's drug |
| PlasmoDB reference | NCBI RefSeq GCF_000002765.6 | PlasmoDB download pages are rendered by JavaScript and cannot be scripted; no organelle genomes in this reference |
| Heatmap of top 100 DEGs | Top 100 trophozoite consensus genes across six samples | Trophozoite is the supported result |
| GO enrichment with clusterProfiler | clusterProfiler with UniProt annotation propagated through GO.db; the deprecated `org.Pf.plasmo.db` as a sensitivity check | See below |
| Single comparison | Three stages analysed separately | Stage differences dominate any drug effect |

GO annotation coverage of the 4,497 tested genes: UniProt 66.9% (49.3% with a biological
process term); `org.Pf.plasmo.db` 84.9% (46.5%). The package states that its data are old
and that it is deprecated. UniProt was chosen as the primary source because the question
concerns biological processes and UniProt is current. Consensus-down BP terms agree
between sources (Jaccard 0.53); consensus-up terms agree poorly (0.17), so the host-cell
terms for the up genes exist only in the UniProt annotation.

## Problems encountered and how they were resolved

| Problem | Cause | Resolution |
|---|---|---|
| Scripted PlasmoDB download returned a web page saved as the genome file | JavaScript-rendered download pages | Switched to NCBI RefSeq; checks that genome and annotation are real files |
| An old test alignment was silently reused | A "skip if the file exists" rule matched a BAM with the same name | Found by comparing per-sample alignment rates; a sample is now skipped only if a `.done` marker exists |
| Annotation-based results empty (no ribosomal genes, composition shares zero, no kelch13 candidates) | Gene IDs carried a `gene-` prefix in the counts and not in the annotation | Prefix stripped when the GTF is built; the analysis stops with an error if IDs do not match |
| Low alignment of the R2 mate | Adapter read-through | Trimming with fastp |
| A stage regression gave misleading numbers | Stage axes were built from the same vehicle samples that form the denominator of the DHA fold change | Retracted and replaced by axes from the independent untreated samples |
| A partition-based calibration check looked clean | Only valid for schizont; in other stages the DHA effect entered the within-group variance | Documented as valid for schizont only |
| Early assumption that the Bioconductor annotation package was too poor to use | Coverage had not been measured | Measured: coverage is real, but the package is deprecated, so it is used for sensitivity only |

## Limitations

- Two replicates per condition (2 residual degrees of freedom per stage model); shallow
  depth; one dataset, one short exposure (about 1 h), one strain (K1) mapped to the 3D7
  reference; no organelle genomes in the reference.
- The vehicle differs from untreated in every stage (225 to 1,357 genes), and ring results
  depend on the control used.
- The stage-related component of the DHA response cannot be separated from drug-specific
  effects with this design. The stage-adjustment is noisy (two untreated samples per
  stage) and a p-value of 0.18 is a failure to detect a residual, not proof of none.
- GO term names describe annotation labels, not measured effects (for example, the
  adhesion term is driven by one exported protein family). Permutation p-values for gene
  sets ignore co-regulation and are optimistic.
- The exported-protein keyword enrichment used in exploratory work and the comparison set
  of other kelch-domain proteins were defined after part of the expression had been seen.
- Not investigated: the cause of the ring library composition differences, the excess of
  counted fragments over ENA read counts (about 2% to 7%), and the origin of unmapped reads
  (not BLASTed).

## Repository layout

```
project13/
  README.md
  report/        project13_report.Rmd, project13_report.html
  scripts/       config.sh, 00_setup.sh ... 11_stage_overlap.R
  data/          samples.tsv, ena_runs_PRJNA263263.tsv, uniprot_pf3d7_go2.tsv
  results/
    counts_all.txt, counts_all.txt.summary, alignment_qc.tsv, reference_release.txt
    env_rnaseq_cli.yml, env_rnaseq_r.yml, setup_log.txt
    deseq2/      per-stage tables, consensus lists, VST matrices, PCA
    kelch13/     candidate list, set tests, stage-adjusted tables
    go/          enrichment tables, annotation pairs, console logs
    figures/     volcano, MA, heatmap and GO figures (PNG and PDF)
    export/      ranked lists and top-20 table (see below)
    hisat2/, fastp/, multiqc/   per-sample logs and the MultiQC report
```

Not committed (regenerated by the scripts): BAM and FASTQ files, the HISAT2 index, and
`data/annotation.gff`.

## Outputs for Week 14

In `results/export/`: `ranked_by_wald_ring.tsv`, `ranked_by_wald_troph.tsv`,
`ranked_by_wald_schiz.tsv` (all tested genes ordered by the Wald statistic, with log2FC,
padj and consensus flag) and `troph_top20_consensus.csv` (the 20 trophozoite consensus genes
with the lowest adjusted p-values). These ranked lists suit pre-ranked GSEA, and the top-20
table suits manual annotation.

## Deliverables against the brief

- Rmarkdown HTML report: `report/project13_report.html`
- Volcano and MA plots, trophozoite top-100 heatmap, GO dotplot: `results/figures/`
- kelch13-pathway check: report section "kelch13-pathway genes" and `results/kelch13/`
- Reproducible pipeline: this repository

## Data sources

GEO GSE62132 and ENA project PRJNA263263 (reads and metadata); NCBI RefSeq GCF_000002765.6
(reference); UniProt (GO annotation); GO.db, Bioconductor 3.22 releases of DESeq2,
apeglm, clusterProfiler, EnhancedVolcano and ComplexHeatmap (R 4.5.3).

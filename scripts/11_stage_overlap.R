suppressPackageStartupMessages(library(DESeq2))
P <- "/kaggle/working/project13"; out <- file.path(P, "results/go")
set.seed(1)
fc <- read.delim(file.path(P, "results/counts_all.txt"), comment.char = "#", check.names = FALSE)
colnames(fc) <- sub("\\.sorted\\.bam$", "", sub(".*/", "", colnames(fc)))
cts <- as.matrix(fc[, 7:ncol(fc)]); rownames(cts) <- sub("^gene-", "", fc$Geneid); storage.mode(cts) <- "integer"
smp <- read.delim(file.path(P, "data/samples.tsv"), stringsAsFactors = FALSE)
cd <- smp[match(colnames(cts), smp$run), ]; rownames(cd) <- cd$run
ui <- cd$treatment == "untreated"
cdu <- cd[ui, ]; cdu$stage <- factor(cdu$stage, levels = c("troph", "ring", "schiz"))
du <- DESeqDataSetFromMatrix(cts[, ui], cdu, design = ~ stage)
du <- du[rowSums(counts(du) >= 10) >= 2, ]; du <- DESeq(du, quiet = TRUE)
rr <- as.data.frame(results(du, contrast = c("stage", "ring", "troph")))
rs <- as.data.frame(results(du, contrast = c("stage", "schiz", "troph")))
ok <- intersect(rownames(rr), rownames(rs)); rr <- rr[ok, ]; rs <- rs[ok, ]
sg <- function(r, dir) ((r$padj < 0.05) %in% TRUE) & (if (dir > 0) r$log2FoldChange > 1 else r$log2FoldChange < -1)
sets <- list(higher_in_ring_than_troph = ok[sg(rr, 1)], higher_in_schiz_than_troph = ok[sg(rs, 1)],
             trough_at_troph = ok[sg(rr, 1) & sg(rs, 1)], higher_in_troph_than_ring = ok[sg(rr, -1)],
             higher_in_troph_than_schiz = ok[sg(rs, -1)], peak_at_troph = ok[sg(rr, -1) & sg(rs, -1)])
t <- read.csv(file.path(P, "results/deseq2/deseq2_troph_all_genes.csv"), stringsAsFactors = FALSE)
univ <- intersect(t$locus[!is.na(t$padj)], ok)
isT <- function(x) x %in% TRUE
dl <- list(consensus_up = intersect(t$locus[isT(t$consensus) & isT(t$lfc_vs_vehicle > 0)], univ),
           consensus_down = intersect(t$locus[isT(t$consensus) & isT(t$lfc_vs_vehicle < 0)], univ))
cat(sprintf("untreated samples %d | genes with both stage contrasts %d | universe (also tested in trophozoite DE) %d\n", sum(ui), length(ok), length(univ)))
cat("stage-defined sets (|log2FC| > 1, padj < 0.05, untreated samples only), size in universe:\n")
for (n in names(sets)) cat(sprintf("  %-28s %d\n", n, length(intersect(sets[[n]], univ))))
cat(sprintf("DHA lists in universe: consensus_up %d | consensus_down %d\n", length(dl$consensus_up), length(dl$consensus_down)))
cat("\n== overlap of DHA lists with stage-defined sets (hypergeometric; enrichment = observed / expected) ==\n")
for (d in names(dl)) for (s in names(sets)) {
  a <- dl[[d]]; b <- intersect(sets[[s]], univ); N <- length(univ); k <- length(intersect(a, b))
  e <- length(a) * length(b) / N
  pu <- phyper(k - 1, length(b), N - length(b), length(a), lower.tail = FALSE)
  pl <- phyper(k, length(b), N - length(b), length(a), lower.tail = TRUE)
  cat(sprintf("%-15s x %-28s overlap %4d of %4d (%.0f%%) | expected %6.1f | enrichment %.2f | p(enriched) %.2g | p(depleted) %.2g\n",
      d, s, k, length(a), 100 * k / length(a), e, ifelse(e > 0, k / e, NA), pu, pl))
}
anyset <- unique(unlist(sets[c("higher_in_ring_than_troph", "higher_in_schiz_than_troph", "higher_in_troph_than_ring", "higher_in_troph_than_schiz")]))
cat(sprintf("\nshare of DHA genes that differ between trophozoite and ring or schizont in untreated samples: consensus_up %.0f%% | consensus_down %.0f%% | all universe genes %.0f%%\n",
    100 * mean(dl$consensus_up %in% anyset), 100 * mean(dl$consensus_down %in% anyset), 100 * mean(univ %in% anyset)))

show_genes <- function(file, id) {
  r <- read.csv(file.path(out, file), stringsAsFactors = FALSE); row <- r[r$ID == id, ]
  if (!nrow(row)) { cat(sprintf("\n%s not found in %s\n", id, file)); return(invisible()) }
  g <- strsplit(row$geneID, "/")[[1]]; tt <- t[match(g, t$locus), ]
  cat(sprintf("\n%s: %s (%d genes; up to 25 by lowest padj)\n", id, row$Description, length(g)))
  tt <- tt[order(tt$padj), ]; tt <- tt[seq_len(min(25, nrow(tt))), c("locus", "gene", "product", "log2FC", "padj")]
  tt$product <- substr(tt$product, 1, 42); tt$log2FC <- round(tt$log2FC, 2); tt$padj <- signif(tt$padj, 2)
  print(tt, row.names = FALSE)
}
cat("\n== genes behind selected terms (UniProt, consensus lists) ==\n")
show_genes("ora_uniprot_consensus_up_BP.csv", "GO:0016310")
show_genes("ora_uniprot_consensus_up_BP.csv", "GO:0044409")
show_genes("ora_uniprot_consensus_up_BP.csv", "GO:0020035")
show_genes("ora_uniprot_consensus_up_CC.csv", "GO:0020009")
show_genes("ora_uniprot_consensus_down_CC.csv", "GO:0005839")
cat("\nDONE\n")

suppressPackageStartupMessages({
  library(EnhancedVolcano); library(ComplexHeatmap); library(circlize); library(ggplot2); library(grid)
})
P  <- "/kaggle/working/project13"
dd <- file.path(P, "results/deseq2")
fg <- file.path(P, "results/figures"); dir.create(fg, showWarnings = FALSE, recursive = TRUE)
stages <- c("ring", "troph", "schiz")
stage_label <- c(ring = "Ring", troph = "Trophozoite", schiz = "Schizont")
setA <- c(K13 = "PF3D7_1343700", KIC1 = "PF3D7_0606000", KIC2 = "PF3D7_1227700", KIC3 = "PF3D7_0914400",
          KIC4 = "PF3D7_1246300", KIC5 = "PF3D7_1138700", KIC6 = "PF3D7_0609700", KIC7 = "PF3D7_0813000",
          KIC8 = "PF3D7_1014900", KIC9 = "PF3D7_1442400", KIC10 = "PF3D7_1122500", UBP1 = "PF3D7_0104300",
          EPS15 = "PF3D7_1025000", AP2MU = "PF3D7_1218300", MyoF = "PF3D7_1329100", MCA2 = "PF3D7_1438400")
tabs <- lapply(stages, function(st) read.csv(file.path(dd, paste0("deseq2_", st, "_all_genes.csv")), stringsAsFactors = FALSE))
names(tabs) <- stages
mk_label <- function(locus, gene) ifelse(is.na(gene) | gene == "" | gene == locus, locus, gene)
save_gg <- function(p, base, w, h) {
  for (ext in c("png", "pdf")) ggsave(file.path(fg, paste0(base, ".", ext)), p, width = w, height = h, dpi = 300)
}

volcano <- function(t, labels_ids, ttl, base, floor_p = 1e-80) {
  r <- t[!is.na(t$padj) & is.finite(t$log2FC), ]
  r$label <- mk_label(r$locus, r$gene)
  big <- r[r$padj < 0.05 & abs(r$log2FC) >= 1, ]
  sel <- if (is.null(labels_ids)) head(big$label[order(big$padj)], 15) else r$label[r$locus %in% labels_ids]
  n_floor <- sum(r$padj < floor_p)
  r$padj <- pmax(r$padj, floor_p)
  p <- EnhancedVolcano(r, lab = r$label, x = "log2FC", y = "padj", pCutoff = 0.05, FCcutoff = 1,
         selectLab = sel, title = ttl,
         subtitle = sprintf("n = 2/group | %d genes: padj < 0.05, |log2FC| >= 1", nrow(big)),
         caption = if (n_floor > 0) sprintf("%d genes with padj < 1e-80 are drawn at 1e-80", n_floor) else NULL,
         pointSize = 1.5, labSize = 3.5, drawConnectors = TRUE, max.overlaps = Inf,
         xlab = "log2 fold change (unshrunk)", ylab = "-log10 adjusted p-value",
         col = c("grey70", "#E69F00", "#56B4E9", "#D55E00"))
  save_gg(p, base, if (is.null(labels_ids)) 7 else 9, if (is.null(labels_ids)) 7 else 9)
  cat(sprintf("%-26s plotted %d | pass both cutoffs %d | labelled %d | floored %d\n",
              base, nrow(r), nrow(big), length(sel), n_floor))
}
ma <- function(t, st, base) {
  r <- t[is.finite(t$log2FC_shrunk) & !is.na(t$baseMean) & t$baseMean > 0, ]
  r$sig <- ifelse((r$padj < 0.05) %in% TRUE, "padj < 0.05", "padj >= 0.05")
  a <- r[r$locus %in% setA, ]
  p <- ggplot(r, aes(baseMean, log2FC_shrunk, colour = sig)) + geom_point(size = 0.8, alpha = 0.6) +
    scale_x_log10() + geom_hline(yintercept = 0) +
    geom_smooth(aes(group = 1), colour = "blue", method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.6) +
    geom_point(data = a, shape = 21, fill = NA, colour = "black", size = 2.2, stroke = 0.7) +
    scale_colour_manual(values = c("padj < 0.05" = "#D55E00", "padj >= 0.05" = "grey65"), name = NULL) +
    labs(title = paste0(stage_label[[st]], ": MA plot"), subtitle = "open circles: K13-compartment genes",
         x = "mean of normalised counts", y = "apeglm-shrunk log2 fold change") + theme_bw()
  save_gg(p, base, 7, 5)
  cat(sprintf("%-26s plotted %d | K13-compartment genes drawn %d\n", base, nrow(r), nrow(a)))
}
heat <- function(st, ids, ttl, base) {
  v <- read.csv(file.path(dd, paste0("vst_", st, "_6samples.csv")), row.names = 1, check.names = FALSE)
  keep <- ids[ids %in% rownames(v)]
  cat(sprintf("%-26s requested %d | in VST matrix %d\n", base, length(ids), length(keep)))
  if (length(keep) < 3) { cat("  too few genes, skipped\n"); return(invisible(NULL)) }
  ord <- c("untreated_rep1", "untreated_rep2", "vehicle_rep1", "vehicle_rep2", "DHA_rep1", "DHA_rep2")
  stopifnot(all(ord %in% colnames(v)))
  z <- t(scale(t(as.matrix(v[keep, ord, drop = FALSE]))))
  z <- z[complete.cases(z), , drop = FALSE]
  tb <- tabs[[st]]; lab <- mk_label(tb$locus, tb$gene)[match(rownames(z), tb$locus)]
  rownames(z) <- make.unique(ifelse(is.na(lab), rownames(z), lab))
  tr <- factor(sub("_rep[12]$", "", ord), levels = c("untreated", "vehicle", "DHA"))
  ann <- HeatmapAnnotation(treatment = tr, annotation_name_gp = gpar(fontsize = 8),
           col = list(treatment = c(untreated = "#CCCCCC", vehicle = "#56B4E9", DHA = "#D55E00")))
  ht <- Heatmap(z, name = "row z-score (VST)", col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
         top_annotation = ann, cluster_rows = TRUE, cluster_columns = FALSE,
         row_names_gp = gpar(fontsize = if (nrow(z) > 60) 6 else 8), column_names_gp = gpar(fontsize = 8),
         column_title = ttl, column_title_gp = gpar(fontsize = 8), heatmap_legend_param = list(at = c(-2, 0, 2)))
  hh <- 3 + 0.11 * nrow(z)
  png(file.path(fg, paste0(base, ".png")), width = 7 * 300, height = hh * 300, res = 300); draw(ht); dev.off()
  pdf(file.path(fg, paste0(base, ".pdf")), width = 7, height = hh); draw(ht); dev.off()
  cat(sprintf("  drawn: %d genes\n", nrow(z)))
}

cat("== volcano plots ==\n")
for (st in stages) volcano(tabs[[st]], NULL, paste0(stage_label[[st]], ": DHA vs vehicle"), paste0("volcano_", st))
for (st in c("troph", "ring"))
  volcano(tabs[[st]], unname(setA), paste0(stage_label[[st]], ": K13-compartment genes labelled"), paste0("volcano_", st, "_K13set"))
cat("\n== MA plots ==\n")
for (st in stages) ma(tabs[[st]], st, paste0("maplot_", st))
cat("\n== heatmaps (6 samples: untreated, vehicle, DHA) ==\n")
ct <- read.csv(file.path(dd, "consensus_troph.csv"), stringsAsFactors = FALSE)
cr <- read.csv(file.path(dd, "consensus_ring.csv"), stringsAsFactors = FALSE)
heat("troph", head(ct$locus, 100), "Trophozoite: top 100 consensus genes", "heatmap_troph_top100_consensus")
heat("troph", unname(setA), "Trophozoite: K13-compartment genes", "heatmap_troph_K13set")
heat("ring", cr$locus, "Ring: consensus genes (tentative)", "heatmap_ring_consensus")
cat("\nfiles:\n"); print(sort(list.files(fg)))
cat("\nDONE\n")

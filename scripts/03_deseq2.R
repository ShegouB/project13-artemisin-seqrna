suppressPackageStartupMessages({library(DESeq2); library(apeglm); library(ggplot2)})
P   <- "/kaggle/working/project13"
dd  <- file.path(P, "data")
out <- file.path(P, "results", "deseq2")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

## 1. data
fc <- read.delim(file.path(P, "results/counts_all.txt"), comment.char = "#", check.names = FALSE)
colnames(fc) <- sub("\\.sorted\\.bam$", "", sub(".*/", "", colnames(fc)))
cts <- as.matrix(fc[, 7:ncol(fc)]); rownames(cts) <- sub("^gene-", "", fc$Geneid)
storage.mode(cts) <- "integer"
smp <- read.delim(file.path(dd, "samples.tsv"), stringsAsFactors = FALSE)
coldata <- smp[match(colnames(cts), smp$run), ]
stopifnot(!anyNA(coldata$run), all(coldata$run == colnames(cts)))
rownames(coldata) <- coldata$run
coldata$stage <- factor(coldata$stage, levels = c("ring", "troph", "schiz"))
cat("count matrix:", dim(cts), "\n"); print(table(coldata$stage, coldata$treatment))

## 2. annotation from the NCBI GFF
gff <- read.delim(file.path(dd, "annotation.gff"), header = FALSE, comment.char = "#",
                  quote = "", stringsAsFactors = FALSE)
at <- gff$V9
getkey <- function(key) ifelse(grepl(paste0("(^|;)", key, "="), at),
  sub(paste0(".*(^|;)", key, "=([^;]*).*"), "\\2", at), NA_character_)
ann0 <- data.frame(locus = getkey("locus_tag"), gene = getkey("gene"),
                   product = getkey("product"), stringsAsFactors = FALSE)
ann0 <- ann0[!is.na(ann0$locus), ]
first <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA_character_ }
ann <- do.call(rbind, lapply(split(ann0, ann0$locus), function(d)
  data.frame(locus = d$locus[1], gene = first(d$gene), product = first(d$product), stringsAsFactors = FALSE)))
ann$product <- vapply(ann$product, function(x) if (is.na(x)) NA_character_ else URLdecode(x), "", USE.NAMES = FALSE)
ann$is_rp   <- !is.na(ann$product) & grepl("ribosomal protein", ann$product, ignore.case = TRUE)
ann$is_rrna <- !is.na(ann$product) & grepl("ribosomal RNA", ann$product, ignore.case = TRUE) &
               !grepl("protein|process|methyl", ann$product, ignore.case = TRUE)
am <- function(l, col) ann[[col]][match(l, ann$locus)]
cat("annotated loci:", nrow(ann), "| ribosomal protein genes:", sum(ann$is_rp), "| rRNA loci:", sum(ann$is_rrna), "\n")

## 3. QC: composition and PCA (all 18 samples)
if (mean(rownames(cts) %in% ann$locus) < 0.99) stop("gene IDs do not match the annotation: check the Geneid format in counts_all.txt")
rp <- am(rownames(cts), "is_rp") %in% TRUE; rr <- am(rownames(cts), "is_rrna") %in% TRUE
comp <- data.frame(coldata[, c("stage", "treatment", "replicate")], assigned = colSums(cts),
  rrna_share = round(colSums(cts[rr, , drop = FALSE]) / colSums(cts), 4),
  rp_share = round(colSums(cts[rp, , drop = FALSE]) / colSums(cts), 3))
comp <- comp[order(comp$stage, comp$treatment, comp$replicate), ]
cat("\n== library composition ==\n"); print(comp)
write.csv(comp, file.path(out, "composition_by_sample.csv"))

filt <- function(d) d[rowSums(counts(d) >= 10) >= 2, ]
dall <- filt(DESeqDataSetFromMatrix(cts, coldata, design = ~ 1))
vall <- assay(vst(dall, blind = TRUE))
pcs <- prcomp(t(vall[order(apply(vall, 1, var), decreasing = TRUE)[1:500], ]))
ve <- round(100 * pcs$sdev^2 / sum(pcs$sdev^2), 1)
sc <- data.frame(coldata[, c("stage", "treatment", "replicate")], round(pcs$x[, 1:3], 1))
cat("\n== PCA, 18 samples, top 500 variable genes ==\nvariance explained PC1-4:", ve[1:4], "\n")
print(sc); write.csv(sc, file.path(out, "pca_scores_18samples.csv"))
pp <- ggplot(sc, aes(PC1, PC2, colour = stage, shape = treatment)) + geom_point(size = 3.5) +
  labs(x = paste0("PC1: ", ve[1], "% variance"), y = paste0("PC2: ", ve[2], "% variance")) + theme_bw()
ggsave(file.path(out, "pca_18samples.png"), pp, width = 6.5, height = 5, dpi = 300)

## 4. helpers
n_sig <- function(padj, lfc) { s <- (padj < 0.05) %in% TRUE
  sprintf("%d (up %d, down %d)", sum(s), sum(s & lfc > 0, na.rm = TRUE), sum(s & lfc < 0, na.rm = TRUE)) }
ribo_test <- function(stat, bm, rpv, B = 5000) {
  g <- is.finite(stat) & !is.na(bm); stat <- stat[g]; bm <- bm[g]; rpv <- rpv[g]
  bins <- cut(bm, unique(quantile(bm, probs = seq(0, 1, 0.1))), include.lowest = TRUE)
  need <- table(bins[rpv]); need <- need[need > 0]; pools <- split(seq_along(stat), bins)
  null <- replicate(B, { s <- unlist(lapply(names(need), function(b) {
    pool <- pools[[b]]; pool <- pool[!rpv[pool]]; sample(pool, min(need[[b]], length(pool))) })); median(stat[s]) })
  c(n = sum(rpv), median_rp = median(stat[rpv]), null_mean = mean(null), share_neg = mean(stat[rpv] < 0))
}
ref <- list(ring  = c(prim = 413,  big = 101, veh_un = 1357, dha_un = 1256, dha_veh3 = 81,   cons = 7),
            troph = c(prim = 1703, big = 676, veh_un = 317,  dha_un = 1566, dha_veh3 = 1312, cons = 973),
            schiz = c(prim = 0,    big = 0,   veh_un = 225,  dha_un = 340,  dha_veh3 = 0,    cons = 0))

## 5. per-stage models
summ <- list()
for (st in levels(coldata$stage)) {
  cat("\n==========", st, "==========\n")
  i6 <- coldata$stage == st
  cd6 <- droplevels(coldata[i6, ]); c6 <- cts[, i6]
  cd6$trt <- factor(cd6$treatment, levels = c("untreated", "vehicle", "DHA"))
  i4 <- cd6$treatment %in% c("DHA", "vehicle")
  cd4 <- cd6[i4, ]; cd4$treatment <- factor(cd4$treatment, levels = c("vehicle", "DHA"))

  d2 <- DESeq(filt(DESeqDataSetFromMatrix(c6[, i4], cd4, design = ~ treatment)), quiet = TRUE)
  r2 <- results(d2, contrast = c("treatment", "DHA", "vehicle"), alpha = 0.05)
  sh <- lfcShrink(d2, coef = "treatment_DHA_vs_vehicle", type = "apeglm", quiet = TRUE)
  d3 <- DESeq(filt(DESeqDataSetFromMatrix(c6, cd6, design = ~ trt)), quiet = TRUE)
  rv <- results(d3, contrast = c("trt", "DHA", "vehicle"), alpha = 0.05)
  ru <- results(d3, contrast = c("trt", "DHA", "untreated"), alpha = 0.05)
  rx <- results(d3, contrast = c("trt", "vehicle", "untreated"), alpha = 0.05)

  tab <- data.frame(locus = rownames(cts), gene = am(rownames(cts), "gene"),
                    product = am(rownames(cts), "product"), stringsAsFactors = FALSE)
  m2 <- match(tab$locus, rownames(r2))
  tab$baseMean <- r2$baseMean[m2]; tab$log2FC <- r2$log2FoldChange[m2]; tab$lfcSE <- r2$lfcSE[m2]
  tab$log2FC_shrunk <- sh$log2FoldChange[match(tab$locus, rownames(sh))]
  tab$stat <- r2$stat[m2]; tab$pvalue <- r2$pvalue[m2]; tab$padj <- r2$padj[m2]
  mv <- match(tab$locus, rownames(rv)); mu <- match(tab$locus, rownames(ru)); mx <- match(tab$locus, rownames(rx))
  tab$baseMean3 <- rv$baseMean[mv]
  tab$lfc_vs_vehicle <- rv$log2FoldChange[mv]; tab$stat_vs_vehicle <- rv$stat[mv]; tab$padj_vs_vehicle <- rv$padj[mv]
  tab$lfc_vs_untreated <- ru$log2FoldChange[mu]; tab$stat_vs_untreated <- ru$stat[mu]; tab$padj_vs_untreated <- ru$padj[mu]
  tab$lfc_veh_vs_untr <- rx$log2FoldChange[mx]; tab$padj_veh_vs_untr <- rx$padj[mx]
  tab$consensus <- (tab$padj_vs_vehicle < 0.05) %in% TRUE & (tab$padj_vs_untreated < 0.05) %in% TRUE &
                   (sign(tab$lfc_vs_vehicle) == sign(tab$lfc_vs_untreated)) %in% TRUE

  cat("size factors, primary model:", paste(round(sizeFactors(d2), 3), collapse = " "), "\n")
  cat(sprintf("primary DHA vs vehicle (4 samples): tested %d | padj<0.05 %s | |log2FC|>=1: %d\n",
      sum(!is.na(tab$padj)), n_sig(tab$padj, tab$log2FC),
      sum(((tab$padj < 0.05) %in% TRUE) & ((abs(tab$log2FC) >= 1) %in% TRUE))))
  cat(sprintf("three-group (6 samples): vehicle vs untreated %s | DHA vs untreated %s | DHA vs vehicle %s\n",
      n_sig(tab$padj_veh_vs_untr, tab$lfc_veh_vs_untr), n_sig(tab$padj_vs_untreated, tab$lfc_vs_untreated),
      n_sig(tab$padj_vs_vehicle, tab$lfc_vs_vehicle)))
  cat(sprintf("consensus genes: %d (up %d, down %d)\n", sum(tab$consensus),
      sum(tab$consensus & tab$lfc_vs_vehicle > 0), sum(tab$consensus & tab$lfc_vs_vehicle < 0)))
  cur <- c(prim = sum((tab$padj < 0.05) %in% TRUE),
           big = sum(((tab$padj < 0.05) %in% TRUE) & ((abs(tab$log2FC) >= 1) %in% TRUE)),
           veh_un = sum((tab$padj_veh_vs_untr < 0.05) %in% TRUE), dha_un = sum((tab$padj_vs_untreated < 0.05) %in% TRUE),
           dha_veh3 = sum((tab$padj_vs_vehicle < 0.05) %in% TRUE), cons = sum(tab$consensus))
  cat("computed :", paste(names(cur), cur, collapse = " | "), "\n")
  cat("Day 3 ref:", paste(names(ref[[st]]), ref[[st]], collapse = " | "), "\n")
  cat("matches Day 3:", all(cur == ref[[st]]), "\n")

  rpv <- am(tab$locus, "is_rp") %in% TRUE
  for (nm in c("stat", "stat_vs_vehicle", "stat_vs_untreated")) {
    r <- ribo_test(tab[[nm]], tab[[if (nm == "stat") "baseMean" else "baseMean3"]], rpv)
    cat(sprintf("ribosomal protein genes, %-18s n=%d | median Wald %6.2f | matched-null mean %6.2f | share negative %.2f\n",
        nm, r["n"], r["median_rp"], r["null_mean"], r["share_neg"]))
  }
  write.csv(tab, file.path(out, paste0("deseq2_", st, "_all_genes.csv")), row.names = FALSE)
  cons <- tab[tab$consensus, ]
  cons <- cons[order(pmax(cons$padj_vs_vehicle, cons$padj_vs_untreated)), ]
  write.csv(cons, file.path(out, paste0("consensus_", st, ".csv")), row.names = FALSE)
  v6 <- assay(vst(d3, blind = TRUE)); colnames(v6) <- paste(cd6$treatment, cd6$replicate, sep = "_")
  write.csv(v6, file.path(out, paste0("vst_", st, "_6samples.csv")))
  summ[[st]] <- data.frame(stage = st, tested_primary = sum(!is.na(tab$padj)), t(cur), row.names = NULL)
}
sm <- do.call(rbind, summ); print(sm, row.names = FALSE)
write.csv(sm, file.path(out, "summary_by_stage.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_deseq2.txt"))
cat("\nDONE\n")

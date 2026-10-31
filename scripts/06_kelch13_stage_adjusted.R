suppressPackageStartupMessages(library(DESeq2))
P   <- "/kaggle/working/project13"
out <- file.path(P, "results", "kelch13"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(1)
setA <- c(K13 = "PF3D7_1343700", KIC1 = "PF3D7_0606000", KIC2 = "PF3D7_1227700", KIC3 = "PF3D7_0914400",
          KIC4 = "PF3D7_1246300", KIC5 = "PF3D7_1138700", KIC6 = "PF3D7_0609700", KIC7 = "PF3D7_0813000",
          KIC8 = "PF3D7_1014900", KIC9 = "PF3D7_1442400", KIC10 = "PF3D7_1122500", UBP1 = "PF3D7_0104300",
          EPS15 = "PF3D7_1025000", AP2MU = "PF3D7_1218300", MyoF = "PF3D7_1329100", MCA2 = "PF3D7_1438400")

## stage axes from the six untreated samples
fc <- read.delim(file.path(P, "results/counts_all.txt"), comment.char = "#", check.names = FALSE)
colnames(fc) <- sub("\\.sorted\\.bam$", "", sub(".*/", "", colnames(fc)))
cts <- as.matrix(fc[, 7:ncol(fc)]); rownames(cts) <- sub("^gene-", "", fc$Geneid); storage.mode(cts) <- "integer"
smp <- read.delim(file.path(P, "data/samples.tsv"), stringsAsFactors = FALSE)
cd <- smp[match(colnames(cts), smp$run), ]; rownames(cd) <- cd$run
ui <- cd$treatment == "untreated"
cdu <- cd[ui, ]; cdu$stage <- factor(cdu$stage, levels = c("ring", "troph", "schiz"))
du <- DESeqDataSetFromMatrix(cts[, ui], cdu, design = ~ stage)
du <- du[rowSums(counts(du) >= 10) >= 2, ]
du <- DESeq(du, quiet = TRUE)
a1 <- as.data.frame(results(du, contrast = c("stage", "troph", "ring")))
a2 <- as.data.frame(results(du, contrast = c("stage", "schiz", "troph")))
axes <- data.frame(locus = rownames(a1), a1 = a1$log2FoldChange, a2 = a2$log2FoldChange, stringsAsFactors = FALSE)
cat("untreated samples used for the stage axes:", sum(ui), "| genes with axes:", nrow(axes), "\n")

null_med <- function(v, inset, poolmask, bm, B = 5000) {
  bins <- cut(bm, unique(quantile(bm[poolmask | inset], probs = seq(0, 1, 0.1))), include.lowest = TRUE)
  need <- table(bins[inset]); need <- need[need > 0]
  pools <- split(which(poolmask), bins[poolmask])
  sims <- replicate(B, median(v[unlist(lapply(names(need), function(b) {
    p <- pools[[b]]; if (is.null(p)) return(integer(0)); p[sample.int(length(p), min(need[[b]], length(p)))] }))]))
  obs <- median(v[inset])
  c(obs = obs, null = mean(sims), p = (1 + sum(sims >= obs)) / (B + 1))
}
for (st in c("troph", "ring")) {
  cat("\n==========", st, "==========\n")
  t <- read.csv(file.path(P, "results/deseq2", paste0("deseq2_", st, "_all_genes.csv")), stringsAsFactors = FALSE)
  m <- merge(t[, c("locus", "gene", "product", "baseMean", "log2FC", "stat", "padj")], axes, by = "locus")
  m <- m[is.finite(m$log2FC) & is.finite(m$a1) & is.finite(m$a2) & m$baseMean >= 50, ]
  inA <- m$locus %in% unname(setA)
  rpf <- !is.na(m$product) & grepl("ribosomal protein|ribosomal RNA", m$product, ignore.case = TRUE)
  cat(sprintf("genes used: %d | set A genes among them: %d\n", nrow(m), sum(inA)))
  cat(sprintf("stage axes (log2FC, untreated): set A median a1 (troph vs ring) %.2f vs other genes %.2f | median a2 (schiz vs troph) %.2f vs %.2f\n",
      median(m$a1[inA]), median(m$a1[!inA]), median(m$a2[inA]), median(m$a2[!inA])))
  for (lab in c("all non-set genes", "excluding ribosomal genes")) {
    train <- !inA & (lab == "all non-set genes" | !rpf)
    fit <- lm(log2FC ~ a1 + a2, data = m[train, ])
    resid <- m$log2FC - predict(fit, newdata = m)
    cf <- summary(fit)$coefficients
    cat(sprintf("fit on %-26s genes=%d R2=%.3f slope a1=%.3f slope a2=%.3f\n", lab, sum(train), summary(fit)$r.squared, cf[2, 1], cf[3, 1]))
    r_raw <- null_med(m$log2FC, inA, train, m$baseMean)
    r_res <- null_med(resid, inA, train, m$baseMean)
    cat(sprintf("   set A median log2FC            %.3f | matched random %.3f | p=%.4f\n", r_raw["obs"], r_raw["null"], r_raw["p"]))
    cat(sprintf("   set A median stage-adjusted    %.3f | matched random %.3f | p=%.4f\n", r_res["obs"], r_res["null"], r_res["p"]))
    if (lab == "all non-set genes") m$resid <- resid
  }
  g <- merge(data.frame(name = names(setA), locus = unname(setA)), m[, c("locus", "baseMean", "log2FC", "a1", "a2", "resid")], by = "locus", all.x = TRUE)
  g[, -(1:2)] <- round(g[, -(1:2)], 2)
  cat("per-gene (a1 = troph vs ring, a2 = schiz vs troph in untreated samples; resid = DHA log2FC not explained by stage axes):\n")
  print(g[order(g$name), ], row.names = FALSE)
  write.csv(g, file.path(out, paste0("setA_stage_adjusted_", st, ".csv")), row.names = FALSE)
}

## exports for Week 14
ex <- file.path(P, "results/export")
for (st in c("ring", "troph", "schiz")) {
  t <- read.csv(file.path(P, "results/deseq2", paste0("deseq2_", st, "_all_genes.csv")), stringsAsFactors = FALSE)
  r <- t[!is.na(t$stat), c("locus", "gene", "product", "baseMean", "log2FC", "log2FC_shrunk", "stat", "pvalue", "padj", "consensus")]
  r <- r[order(-r$stat), ]
  write.table(r, file.path(ex, paste0("ranked_by_wald_", st, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("exported ranked list %s: %d genes\n", st, nrow(r)))
}
cons <- read.csv(file.path(P, "results/deseq2/consensus_troph.csv"), stringsAsFactors = FALSE)
top20 <- head(cons, 20)[, c("locus", "gene", "product", "baseMean", "log2FC", "log2FC_shrunk", "padj", "padj_vs_untreated")]
write.csv(top20, file.path(ex, "troph_top20_consensus.csv"), row.names = FALSE)
cat("\ntop 20 trophozoite consensus genes (for annotation):\n"); print(top20, row.names = FALSE)
cat("\nDONE\n")

P   <- "/kaggle/working/project13"
out <- file.path(P, "results", "kelch13"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
stages <- c("ring", "troph", "schiz")
tabs <- lapply(stages, function(st) read.csv(file.path(P, "results/deseq2", paste0("deseq2_", st, "_all_genes.csv")), stringsAsFactors = FALSE))
names(tabs) <- stages
base <- tabs$ring[, c("locus", "gene", "product")]

sym  <- "^(K13|KELCH13|PFK13|KIC[0-9]+|UBP1|EPS15|AP2MU|AP-2MU|AP2-MU|MYOF|MCA2|CORONIN)$"
prod <- "kelch|ubiquitin carboxyl-terminal hydrolase 1|eps15|AP-2 complex subunit mu|AP2 complex subunit mu|myosin F|metacaspase|coronin"
anchors <- c(K13 = "PF3D7_1343700", KIC1 = "PF3D7_0606000")
hit_sym  <- !is.na(base$gene) & grepl(sym, base$gene, ignore.case = TRUE)
hit_prod <- !is.na(base$product) & grepl(prod, base$product, ignore.case = TRUE)
base$match <- ifelse(base$locus %in% anchors, "anchor", ifelse(hit_sym, "symbol", ifelse(hit_prod, "product", NA)))
cand <- base[!is.na(base$match), ]
cat("candidate genes:", nrow(cand), "| anchor:", sum(cand$match == "anchor"),
    "| symbol:", sum(cand$match == "symbol"), "| product:", sum(cand$match == "product"), "\n")
cat("anchors present in the annotation:", paste(names(anchors), anchors %in% base$locus, collapse = ", "), "\n\n")

res <- cand
for (st in stages) {
  t <- tabs[[st]]; m <- match(cand$locus, t$locus)
  Fn <- ecdf(abs(t$stat[is.finite(t$stat)]))
  res[[paste0(st, "_bm")]]   <- round(t$baseMean[m], 0)
  res[[paste0(st, "_lfc")]]  <- round(t$log2FC[m], 2)
  res[[paste0(st, "_padj")]] <- signif(t$padj[m], 2)
  res[[paste0(st, "_pct")]]  <- round(100 * Fn(abs(t$stat[m])), 0)
  res[[paste0(st, "_cons")]] <- t$consensus[m]
}
write.csv(res, file.path(out, "kelch13_candidates_all_columns.csv"), row.names = FALSE)
show <- res; show$product <- substr(show$product, 1, 48)
cat("== candidates (curate by gene symbol and product; ignore the DE columns when deciding membership) ==\n")
print(show[, c("locus", "gene", "product", "match")], row.names = FALSE)
cat("\n== DE columns: baseMean, log2FC, padj (primary), |Wald| percentile within stage ==\n")
print(show[, c("locus", "gene", "ring_bm", "ring_lfc", "ring_padj", "ring_pct", "troph_bm", "troph_lfc",
               "troph_padj", "troph_pct", "schiz_bm", "schiz_lfc", "schiz_padj", "schiz_pct")], row.names = FALSE)
cat("\n== context: genome-wide fraction significant vs candidates ==\n")
for (st in stages) {
  t <- tabs[[st]]; sg <- (t$padj < 0.05) %in% TRUE; tested <- !is.na(t$padj)
  ps <- (res[[paste0(st, "_padj")]] < 0.05) %in% TRUE
  cat(sprintf("%-6s genome: %d of %d tested significant (%.1f%%) | candidates: %d of %d significant, %d of %d tested | consensus %d\n",
      st, sum(sg), sum(tested), 100 * sum(sg) / sum(tested), sum(ps), nrow(res),
      sum(!is.na(res[[paste0(st, "_padj")]])), nrow(res), sum(res[[paste0(st, "_cons")]] %in% TRUE)))
}
cat("\nDONE\n")

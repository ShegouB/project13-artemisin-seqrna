suppressPackageStartupMessages({library(GO.db); library(org.Pf.plasmo.db)})
P <- "/kaggle/working/project13"; out <- file.path(P, "results/go"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
t <- read.csv(file.path(P, "results/deseq2/deseq2_troph_all_genes.csv"), stringsAsFactors = FALSE)
univ <- t$locus[!is.na(t$padj)]; degs <- t$locus[(t$padj < 0.05) %in% TRUE]; cons <- t$locus[t$consensus %in% TRUE]
roots <- c("GO:0008150", "GO:0003674", "GO:0005575")
cat(sprintf("universe %d | primary DEGs %d | consensus %d\n", length(univ), length(degs), length(cons)))

ontology_of <- function(ids) { x <- AnnotationDbi::select(GO.db, keys = unique(ids), keytype = "GOID", columns = "ONTOLOGY"); setNames(x$ONTOLOGY, x$GOID) }
cov_line <- function(lbl, loci) cat(sprintf("  %-34s universe %4d (%.1f%%) | primary DEGs %4d (%.1f%%) | consensus %4d (%.1f%%)\n",
  lbl, sum(univ %in% loci), 100 * mean(univ %in% loci), sum(degs %in% loci), 100 * mean(degs %in% loci), sum(cons %in% loci), 100 * mean(cons %in% loci)))
summarise <- function(name, pr) {
  cat(sprintf("\n%s\n", name))
  pr <- unique(pr[!is.na(pr$go) & !(pr$go %in% roots), c("locus", "go")])
  om <- ontology_of(pr$go); pr$ont <- om[pr$go]
  cov_line("any non-root term", unique(pr$locus))
  for (o in c("BP", "MF", "CC")) cov_line(paste(o, "non-root term"), unique(pr$locus[pr$ont %in% o]))
  sz <- table(pr$go[pr$locus %in% univ])
  cat(sprintf("  terms with 10 to 500 universe genes: %d (of %d terms touching the universe)\n", sum(sz >= 10 & sz <= 500), length(sz)))
  invisible(pr)
}

## UniProt, second fetch
f <- file.path(P, "data/uniprot_pf3d7_go2.tsv")
up <- tryCatch(read.delim(f, quote = "", stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
if (is.null(up) || ncol(up) < 3) { cat("\nUniProt table not usable:\n"); print(head(readLines(f), 5)) } else {
  cat("\nUniProt columns:", paste(names(up), collapse = " | "), "\n")
  txt <- apply(up[, -c(1, ncol(up)), drop = FALSE], 1, paste, collapse = " ")
  loc <- regmatches(txt, gregexpr("PF3D7_[0-9]{7}", txt))
  cat("rows with a PF3D7 ID:", sum(lengths(loc) > 0), "of", nrow(up), "\n")
  go <- strsplit(ifelse(is.na(up[[ncol(up)]]), "", up[[ncol(up)]]), ";\\s*")
  pr <- do.call(rbind, lapply(seq_along(loc), function(i) { g <- go[[i]][nzchar(go[[i]])]
    if (length(loc[[i]]) && length(g)) expand.grid(locus = unique(loc[[i]]), go = g, stringsAsFactors = FALSE) }))
  if (is.null(pr)) cat("no PF3D7 IDs with GO terms\n") else {
    p1 <- summarise("UniProt (direct annotations, not propagated)", pr)
    write.table(p1, file.path(out, "go_pairs_uniprot_direct.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

## Bioconductor package: evidence codes, root terms, propagated annotations
db <- org.Pf.plasmo.db
ak <- keys(db, keytype = "ALIAS")
cat(sprintf("\nuniverse IDs found as ALIAS keys: %d of %d\n", sum(univ %in% ak), length(univ)))
cat("ORF key examples:", paste(head(keys(db, keytype = "ORF"), 3), collapse = ", "), "\n")
k <- intersect(univ, ak)
x <- AnnotationDbi::select(db, keys = k, keytype = "ALIAS", columns = c("GO", "EVIDENCE", "ONTOLOGY"))
cat("\nevidence codes (rows, all genes in universe):\n"); print(sort(table(x$EVIDENCE, useNA = "ifany"), decreasing = TRUE))
cat("root-term rows:", sum(x$GO %in% roots), "| root-term rows with evidence ND:", sum(x$GO %in% roots & x$EVIDENCE %in% "ND"), "\n")
pr_all <- data.frame(locus = x$ALIAS, go = x$GO, stringsAsFactors = FALSE)
summarise("org.Pf.plasmo.db direct, roots removed", pr_all)
xn <- x[!(x$EVIDENCE %in% "ND"), ]
p2 <- summarise("org.Pf.plasmo.db direct, roots and ND removed", data.frame(locus = xn$ALIAS, go = xn$GO, stringsAsFactors = FALSE))
y <- AnnotationDbi::select(db, keys = k, keytype = "ALIAS", columns = "GOALL")
p3 <- summarise("org.Pf.plasmo.db GOALL (propagated), roots removed", data.frame(locus = y$ALIAS, go = y$GOALL, stringsAsFactors = FALSE))
write.table(p3, file.path(out, "go_pairs_bioc_goall.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDONE\n")

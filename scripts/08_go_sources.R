P <- "/kaggle/working/project13"; out <- file.path(P, "results/go"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
t <- read.csv(file.path(P, "results/deseq2/deseq2_troph_all_genes.csv"), stringsAsFactors = FALSE)
univ <- t$locus[!is.na(t$padj)]; degs <- t$locus[(t$padj < 0.05) %in% TRUE]; cons <- t$locus[t$consensus %in% TRUE]
cat(sprintf("trophozoite universe %d | primary DEGs %d | consensus %d\n", length(univ), length(degs), length(cons)))
cat("GO.db version:", as.character(packageVersion("GO.db")), "\n")

cov_report <- function(name, pairs) {
  g <- unique(pairs$locus)
  cat(sprintf("\n%s: %d loci with >=1 GO term\n", name, length(g)))
  cat(sprintf("  universe: %d of %d (%.1f%%) | primary DEGs: %d of %d (%.1f%%) | consensus: %d of %d (%.1f%%)\n",
      sum(univ %in% g), length(univ), 100 * mean(univ %in% g), sum(degs %in% g), length(degs), 100 * mean(degs %in% g),
      sum(cons %in% g), length(cons), 100 * mean(cons %in% g)))
  per <- table(pairs$locus[pairs$locus %in% univ])
  cat(sprintf("  median terms per annotated gene in the universe: %.0f\n", median(per)))
  ont <- tryCatch({ suppressPackageStartupMessages(library(GO.db))
    AnnotationDbi::select(GO.db, keys = unique(pairs$go), keytype = "GOID", columns = "ONTOLOGY") }, error = function(e) NULL)
  if (!is.null(ont)) {
    m <- merge(pairs[pairs$locus %in% univ, ], ont, by.x = "go", by.y = "GOID")
    for (o in c("BP", "MF", "CC")) cat(sprintf("  %s: %d universe genes\n", o, length(unique(m$locus[m$ONTOLOGY == o]))))
  }
}
res_pairs <- list()

## 1. UniProt
f <- file.path(P, "data/uniprot_pf3d7_go.tsv")
up <- tryCatch(read.delim(f, quote = "", stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
if (is.null(up) || ncol(up) < 3) {
  cat("\nUniProt table not usable; first lines of the file:\n"); print(head(readLines(f), 5))
} else {
  cat("\nUniProt columns:", paste(names(up), collapse = " | "), "| rows:", nrow(up), "\n")
  loc <- regmatches(up[[2]], gregexpr("PF3D7_[0-9]{7}", up[[2]]))
  go <- strsplit(ifelse(is.na(up[[3]]), "", up[[3]]), ";\\s*")
  pr <- do.call(rbind, lapply(seq_along(loc), function(i) {
    g <- go[[i]][nzchar(go[[i]])]
    if (length(loc[[i]]) && length(g)) expand.grid(locus = loc[[i]], go = g, stringsAsFactors = FALSE) }))
  if (is.null(pr)) cat("no PF3D7 ordered locus names with GO terms found\n") else {
    pr <- unique(pr); res_pairs$uniprot <- pr
    write.table(pr, file.path(out, "go_pairs_uniprot.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cov_report("UniProt", pr)
  }
}

## 2. Bioconductor package
ok <- requireNamespace("org.Pf.plasmo.db", quietly = TRUE)
cat("\norg.Pf.plasmo.db installed:", ok, "\n")
if (ok) {
  suppressPackageStartupMessages(library(org.Pf.plasmo.db))
  d <- packageDescription("org.Pf.plasmo.db"); cat("version:", d$Version, "| date:", d$Date, "\n")
  cat("keytypes:", paste(keytypes(org.Pf.plasmo.db), collapse = ", "), "\n")
  cat("columns:", paste(columns(org.Pf.plasmo.db), collapse = ", "), "\n")
  x <- tryCatch({
    kt <- keytypes(org.Pf.plasmo.db); kty <- if ("GID" %in% kt) "GID" else kt[1]
    cat("key type used:", kty, "\n")
    s <- AnnotationDbi::select(org.Pf.plasmo.db, keys = keys(org.Pf.plasmo.db, keytype = kty), keytype = kty, columns = "GO")
    cat("example keys:", paste(head(unique(s[[kty]]), 3), collapse = ", "), "\n"); s
  }, error = function(e) { cat("select failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(x) && ncol(x) >= 2) {
    pr <- unique(data.frame(locus = x[[1]], go = x[["GO"]], stringsAsFactors = FALSE)); pr <- pr[!is.na(pr$go), ]
    res_pairs$bioc <- pr
    write.table(pr, file.path(out, "go_pairs_bioc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cov_report("org.Pf.plasmo.db", pr)
  }
}
if (length(res_pairs) == 2) {
  a <- unique(res_pairs$uniprot$locus); b <- unique(res_pairs$bioc$locus)
  cat(sprintf("\nloci with GO in both sources: %d | UniProt only: %d | Bioconductor only: %d | in universe, covered by either: %d of %d\n",
      length(intersect(a, b)), length(setdiff(a, b)), length(setdiff(b, a)), sum(univ %in% union(a, b)), length(univ)))
}
cat("\nDONE\n")

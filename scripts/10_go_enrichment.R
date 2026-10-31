suppressPackageStartupMessages({library(clusterProfiler); library(GO.db); library(org.Pf.plasmo.db); library(ggplot2); library(stringr)})
P   <- "/kaggle/working/project13"
out <- file.path(P, "results/go"); fg <- file.path(P, "results/figures")
dir.create(out, showWarnings = FALSE, recursive = TRUE); dir.create(fg, showWarnings = FALSE, recursive = TRUE)
set.seed(1)
t <- read.csv(file.path(P, "results/deseq2/deseq2_troph_all_genes.csv"), stringsAsFactors = FALSE)
univ <- t$locus[!is.na(t$padj)]
isT <- function(x) x %in% TRUE
lists <- list(
  consensus_up   = t$locus[isT(t$consensus) & isT(t$lfc_vs_vehicle > 0)],
  consensus_down = t$locus[isT(t$consensus) & isT(t$lfc_vs_vehicle < 0)],
  primary_up     = t$locus[isT(t$padj < 0.05) & isT(t$log2FC > 0)],
  primary_down   = t$locus[isT(t$padj < 0.05) & isT(t$log2FC < 0)])
cat(sprintf("universe %d\n", length(univ)))
for (n in names(lists)) cat(sprintf("  %-15s %d genes\n", n, length(lists[[n]])))

## annotation sources
roots <- c("GO:0008150", "GO:0003674", "GO:0005575", "all")
anc <- c(as.list(GOBPANCESTOR), as.list(GOMFANCESTOR), as.list(GOCCANCESTOR))
propagate <- function(pr) {
  n0 <- nrow(pr); pr <- pr[pr$go %in% names(anc) & pr$locus %in% univ, ]
  lst <- split(pr$go, pr$locus)
  res <- do.call(rbind, lapply(names(lst), function(l) {
    g <- setdiff(unique(c(lst[[l]], unlist(anc[lst[[l]]]))), roots)
    data.frame(go = g, locus = l, stringsAsFactors = FALSE) }))
  cat(sprintf("  direct pairs: %d | kept (known GO ID, in universe): %d | propagated pairs: %d | genes: %d\n",
              n0, nrow(pr), nrow(res), length(lst)))
  res
}
up <- read.delim(file.path(P, "data/uniprot_pf3d7_go2.tsv"), quote = "", stringsAsFactors = FALSE, check.names = FALSE)
txt <- apply(up[, -c(1, ncol(up)), drop = FALSE], 1, paste, collapse = " ")
loc <- regmatches(txt, gregexpr("PF3D7_[0-9]{7}", txt))
go  <- strsplit(ifelse(is.na(up[[ncol(up)]]), "", up[[ncol(up)]]), ";\\s*")
pr_up <- unique(do.call(rbind, lapply(seq_along(loc), function(i) { g <- go[[i]][nzchar(go[[i]])]
  if (length(loc[[i]]) && length(g)) expand.grid(locus = unique(loc[[i]]), go = g, stringsAsFactors = FALSE) })))
cat("UniProt:\n"); t2g_up <- propagate(pr_up)
y <- AnnotationDbi::select(org.Pf.plasmo.db, keys = intersect(univ, keys(org.Pf.plasmo.db, keytype = "ALIAS")),
                           keytype = "ALIAS", columns = "GOALL")
t2g_bc <- unique(data.frame(go = y$GOALL, locus = y$ALIAS, stringsAsFactors = FALSE))
t2g_bc <- t2g_bc[!is.na(t2g_bc$go) & !(t2g_bc$go %in% roots), ]
cat(sprintf("org.Pf.plasmo.db GOALL: pairs %d | genes %d\n", nrow(t2g_bc), length(unique(t2g_bc$locus))))
tn <- AnnotationDbi::select(GO.db, keys = unique(c(t2g_up$go, t2g_bc$go)), keytype = "GOID", columns = c("TERM", "ONTOLOGY"))
tn <- tn[!duplicated(tn$GOID), ]

## enrichment
ora <- function(genes, t2g, ont) {
  g <- t2g[t2g$go %in% tn$GOID[tn$ONTOLOGY == ont], ]
  r <- tryCatch(enricher(gene = genes, universe = univ, TERM2GENE = g,
        TERM2NAME = data.frame(go = tn$GOID, name = tn$TERM), pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500), error = function(e) NULL)
  if (is.null(r)) NULL else as.data.frame(r@result)
}
reduce_terms <- function(r, thr = 0.5, cut = 0.05) {
  r <- r[r$p.adjust < cut, ]; if (!nrow(r)) return(r)
  r <- r[order(r$p.adjust, -r$Count), ]; keep <- character(0)
  for (i in seq_len(nrow(r))) {
    gi <- strsplit(r$geneID[i], "/")[[1]]; dup <- FALSE
    for (k in keep) { gk <- strsplit(r$geneID[r$ID == k], "/")[[1]]
      if (length(intersect(gi, gk)) / length(union(gi, gk)) >= thr) { dup <- TRUE; break } }
    if (!dup) keep <- c(keep, r$ID[i])
  }
  r[r$ID %in% keep, ]
}
res <- list(); red <- list()
for (src in c("uniprot", "bioc")) {
  t2g <- if (src == "uniprot") t2g_up else t2g_bc
  for (ln in names(lists)) for (ont in c("BP", "MF", "CC")) {
    r <- ora(lists[[ln]], t2g, ont); key <- paste(src, ln, ont, sep = "|")
    if (is.null(r) || !nrow(r)) next
    r <- r[order(r$p.adjust), ]; res[[key]] <- r; red[[key]] <- reduce_terms(r)
    write.csv(r, file.path(out, sprintf("ora_%s_%s_%s.csv", src, ln, ont)), row.names = FALSE)
  }
}
summ <- do.call(rbind, lapply(names(res), function(k) { p <- strsplit(k, "\\|")[[1]]
  data.frame(source = p[1], list = p[2], ont = p[3], tested = nrow(res[[k]]),
             significant = sum(res[[k]]$p.adjust < 0.05), reduced = nrow(red[[k]]), stringsAsFactors = FALSE) }))
cat("\n=== terms tested, significant (BH < 0.05), and after redundancy reduction ===\n")
print(summ, row.names = FALSE); write.csv(summ, file.path(out, "ora_summary.csv"), row.names = FALSE)

show <- function(k, n) { r <- red[[k]]
  if (is.null(r) || !nrow(r)) { cat("  none\n"); return(invisible()) }
  r <- head(r, n); r$Description <- substr(r$Description, 1, 52); r$p.adjust <- signif(r$p.adjust, 2)
  print(r[, c("ID", "Description", "Count", "GeneRatio", "BgRatio", "p.adjust")], row.names = FALSE) }
for (ln in c("consensus_up", "consensus_down")) {
  cat(sprintf("\n=== UniProt, %s: top reduced terms ===\n", ln))
  cat("-- BP --\n"); show(paste("uniprot", ln, "BP", sep = "|"), 10)
  cat("-- MF --\n"); show(paste("uniprot", ln, "MF", sep = "|"), 5)
  cat("-- CC --\n"); show(paste("uniprot", ln, "CC", sep = "|"), 5)
}
cat("\n=== sensitivity: significant BP terms, UniProt vs org.Pf.plasmo.db ===\n")
sig_ids <- function(k) if (is.null(res[[k]])) character(0) else res[[k]]$ID[res[[k]]$p.adjust < 0.05]
jac <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
for (ln in names(lists)) { a <- sig_ids(paste("uniprot", ln, "BP", sep = "|")); b <- sig_ids(paste("bioc", ln, "BP", sep = "|"))
  cat(sprintf("%-15s UniProt %3d | Bioc %3d | shared %3d | Jaccard %s\n", ln, length(a), length(b), length(intersect(a, b)),
      ifelse(is.na(jac(a, b)), "NA", sprintf("%.2f", jac(a, b))))) }
cat("\n=== consensus vs primary lists, significant BP terms (UniProt) ===\n")
for (d in c("up", "down")) { a <- sig_ids(paste0("uniprot|consensus_", d, "|BP")); b <- sig_ids(paste0("uniprot|primary_", d, "|BP"))
  cat(sprintf("%-5s consensus %3d | primary %3d | shared %3d | Jaccard %s\n", d, length(a), length(b), length(intersect(a, b)),
      ifelse(is.na(jac(a, b)), "NA", sprintf("%.2f", jac(a, b))))) }
cat("\n=== org.Pf.plasmo.db, consensus lists: top 5 reduced BP terms ===\n")
for (ln in c("consensus_up", "consensus_down")) { cat(sprintf("-- %s --\n", ln)); show(paste("bioc", ln, "BP", sep = "|"), 5) }

mkplot <- function(src) {
  d <- do.call(rbind, lapply(c("consensus_up", "consensus_down"), function(ln) {
    r <- red[[paste(src, ln, "BP", sep = "|")]]; if (is.null(r) || !nrow(r)) return(NULL)
    data.frame(direction = ifelse(ln == "consensus_up", "higher in DHA", "lower in DHA"), head(r, 10), stringsAsFactors = FALSE) }))
  if (is.null(d)) { cat("no BP terms to plot for", src, "\n"); return(invisible(NULL)) }
  d$ratio <- sapply(strsplit(d$GeneRatio, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
  d$Description <- stringr::str_wrap(d$Description, 38)
  p <- ggplot(d, aes(ratio, reorder(Description, ratio), size = Count, colour = p.adjust)) + geom_point() +
    facet_wrap(~ direction, ncol = 1, scales = "free_y") + scale_colour_viridis_c(direction = -1) +
    labs(x = "gene ratio (DEGs in term / DEGs annotated)", y = NULL,
         title = "Trophozoite: GO biological process\n(redundancy-reduced, consensus genes)") + theme_bw()
  for (ext in c("png", "pdf")) ggsave(file.path(fg, sprintf("go_bp_dotplot_%s.%s", src, ext)), p, width = 8, height = 8, dpi = 300)
}
mkplot("uniprot"); mkplot("bioc")
cat("\nfiles in results/go:\n"); print(sort(list.files(out)))
cat("\nDONE\n")

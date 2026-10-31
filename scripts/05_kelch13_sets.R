P   <- "/kaggle/working/project13"
out <- file.path(P, "results", "kelch13"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(1)
stages <- c("ring", "troph", "schiz")
tabs <- lapply(stages, function(st) read.csv(file.path(P, "results/deseq2", paste0("deseq2_", st, "_all_genes.csv")), stringsAsFactors = FALSE))
names(tabs) <- stages

setA <- c(K13 = "PF3D7_1343700", KIC1 = "PF3D7_0606000", KIC2 = "PF3D7_1227700", KIC3 = "PF3D7_0914400",
          KIC4 = "PF3D7_1246300", KIC5 = "PF3D7_1138700", KIC6 = "PF3D7_0609700", KIC7 = "PF3D7_0813000",
          KIC8 = "PF3D7_1014900", KIC9 = "PF3D7_1442400", KIC10 = "PF3D7_1122500", UBP1 = "PF3D7_0104300",
          EPS15 = "PF3D7_1025000", AP2MU = "PF3D7_1218300", MyoF = "PF3D7_1329100", MCA2 = "PF3D7_1438400")
setB <- c("PF3D7_0724800", "PF3D7_1022600", "PF3D7_1125700", "PF3D7_1125800", "PF3D7_1131600",
          "PF3D7_1205400", "PF3D7_1210600", "PF3D7_1213400", "PF3D7_1466100")
sets <- list(A_K13_compartment = unname(setA), A_without_K13 = unname(setA[names(setA) != "K13"]),
             B_other_kelch_domain = setB)
stopifnot(all(unname(setA) %in% tabs$ring$locus), all(setB %in% tabs$ring$locus))

cat("== per-gene table, Set A (primary model DHA vs vehicle; cons = consensus gene) ==\n")
g <- data.frame(name = names(setA), locus = unname(setA), stringsAsFactors = FALSE)
for (st in stages) {
  t <- tabs[[st]]; m <- match(g$locus, t$locus)
  g[[paste0(st, "_bm")]] <- round(t$baseMean[m]); g[[paste0(st, "_lfc")]] <- round(t$log2FC[m], 2)
  g[[paste0(st, "_padj")]] <- signif(t$padj[m], 2); g[[paste0(st, "_cons")]] <- t$consensus[m]
}
print(g, row.names = FALSE); write.csv(g, file.path(out, "setA_per_gene.csv"), row.names = FALSE)

set_test <- function(tab, ids, B = 5000) {
  tab <- tab[!is.na(tab$stat) & !is.na(tab$baseMean), ]
  inset <- tab$locus %in% ids; n <- sum(inset)
  if (n == 0) return(NULL)
  bins <- cut(tab$baseMean, unique(quantile(tab$baseMean, probs = seq(0, 1, 0.1))), include.lowest = TRUE)
  need <- table(bins[inset]); need <- need[need > 0]; pools <- split(seq_len(nrow(tab)), bins)
  sig <- (tab$padj < 0.05) %in% TRUE; cons <- tab$consensus %in% TRUE
  f <- function(i) c(sum(sig[i]), sum(cons[i]), median(tab$stat[i]), mean(tab$stat[i] > 0))
  obs <- f(which(inset))
  sims <- replicate(B, f(unlist(lapply(names(need), function(b) {
    pool <- pools[[b]]; pool <- pool[!inset[pool]]; sample(pool, min(need[[b]], length(pool))) }))))
  list(n = n, obs = obs, null_mean = rowMeans(sims),
       p_sig = (1 + sum(sims[1, ] >= obs[1])) / (B + 1), p_cons = (1 + sum(sims[2, ] >= obs[2])) / (B + 1),
       p_med = (1 + sum(abs(sims[3, ] - mean(sims[3, ])) >= abs(obs[3] - mean(sims[3, ])))) / (B + 1))
}
cat("\n== set tests vs 5000 random sets matched on expression decile ==\n")
cat("(one-sided p for 'at least as many significant / consensus genes'; p floor 1/5001; sets are co-regulated, so p is optimistic)\n")
rows <- list()
for (st in stages) for (nm in names(sets)) {
  r <- set_test(tabs[[st]], sets[[nm]])
  if (is.null(r)) next
  cat(sprintf("%-6s %-22s n=%2d | significant %2d (random %.1f, p=%.4f) | consensus %2d (random %.1f, p=%.4f) | median Wald %5.2f (random %5.2f, p=%.4f) | share up %.2f (random %.2f)\n",
      st, nm, r$n, r$obs[1], r$null_mean[1], r$p_sig, r$obs[2], r$null_mean[2], r$p_cons,
      r$obs[3], r$null_mean[3], r$p_med, r$obs[4], r$null_mean[4]))
  rows[[length(rows) + 1]] <- data.frame(stage = st, set = nm, n = r$n, sig = r$obs[1], sig_random = r$null_mean[1], p_sig = r$p_sig,
      cons = r$obs[2], cons_random = r$null_mean[2], p_cons = r$p_cons, med_wald = r$obs[3], med_random = r$null_mean[3], p_med = r$p_med)
}
write.csv(do.call(rbind, rows), file.path(out, "set_tests.csv"), row.names = FALSE)
cat("\nDONE\n")

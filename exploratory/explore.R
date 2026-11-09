source("exploratory/paths.R")   # repo-relative inputs; run from the pancancer-deg repo root
if (file.exists(".Rprofile")) eval(parse(".Rprofile")[[2]])   # sandbox proxy fix (msigdbr download)
add_lib(BULK_LIB)
suppressMessages({ library(msigdbr) })
R <- file.path(overlap_results(), "full_signatures")   # two-sided overlap pipeline results
h <- msigdbr(species = "Homo sapiens", collection = "H")
hs <- split(h$gene_symbol, h$gs_name)
prolif <- unique(unlist(hs[c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
                             "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2")]))
dir.create("results/exploratory", recursive = TRUE, showWarnings = FALSE)
writeLines(prolif, "results/exploratory/prolif_genes_msigdbr.txt")   # frozen copy used by the analyses: inputs/other/prolif_genes.txt
sigs <- c("OIS", "RS", "SIPS", "CICQ", "SSCQ")
classes <- c("UP-UP", "UP-DOWN", "DOWN-UP", "DOWN-DOWN")
out <- list()
for (s in sigs) {
  cancers <- list.dirs(file.path(R, s), recursive = FALSE)
  for (cl in classes) {
    g <- lapply(cancers, function(d) read.csv(file.path(d, paste0("overlap_", cl, ".csv")))$gene_symbol)
    univ <- unique(unlist(lapply(cancers, function(d) read.csv(file.path(d, "common_tested_genes.csv"))[[1]])))
    tab <- sort(table(unlist(g)), decreasing = TRUE)
    core <- names(tab)[tab >= 12]
    allg <- unlist(g)
    frac_prolif <- mean(allg %in% prolif)
    # hallmark over-representation in core set vs universe
    hm <- sapply(hs, function(x) {
      x <- intersect(x, univ); a <- sum(core %in% x)
      if (length(core) == 0) return(c(a, 1))
      c(a, phyper(a - 1, length(x), length(univ) - length(x), length(core), lower.tail = FALSE))
    })
    hm <- data.frame(set = sub("HALLMARK_", "", colnames(hm)), n = hm[1, ], p = hm[2, ])
    hm$padj <- p.adjust(hm$p, "BH")
    top <- head(hm[order(hm$p), ][hm[order(hm$p), ]$padj < 0.05, ], 6)
    cat(sprintf("\n== %s %s: mean overlap/cancer %.0f, prolif genes %.0f%%, core(>=12/17) %d\n",
                s, cl, mean(lengths(g)), 100 * frac_prolif, length(core)))
    cat("  core examples:", head(core, 25), "\n")
    if (nrow(top)) cat("  hallmarks:", paste0(top$set, "(", top$n, ")"), "\n")
  }
}

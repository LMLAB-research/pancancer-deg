suppressMessages(library(GeneOverlap))
set.seed(1)
u <- paste0("g", 1:13000)
for (i in 1:5) {
  a <- sample(u, sample(500:3000, 1)); b <- sample(u, sample(500:3000, 1))
  if (i %% 2 == 0) b <- unique(c(b, sample(a, 300)))  # some enriched
  if (i == 5) b <- setdiff(b, a)                       # depleted (zero overlap)
  go <- testGeneOverlap(newGeneOverlap(a, b, genome.size = length(u)))
  o <- length(intersect(a, b))
  m <- matrix(c(o, length(a) - o, length(b) - o, length(u) - length(union(a, b))), 2)
  cat(sprintf("overlap %4d  GeneOverlap p=%.3g  fisher greater p=%.3g  two.sided p=%.3g  OR=%.2f\n",
              o, getPval(go), fisher.test(m, alternative = "greater")$p.value,
              fisher.test(m)$p.value, getOddsRatio(go)))
}

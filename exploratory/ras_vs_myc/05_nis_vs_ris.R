# Non-RAS oncogene senescence: NOTCH1-ICD-induced (NIS) vs RAS-induced (RIS) senescence in
# the same IMR90 ER:RAS cells (Hoare 2016, SRP062871; Control / NIS / RIS / RNIS, n=6).
# House find_degs_new() for each arm vs control. Run from the pancancer-deg repo root.
source(".claude/worktrees/two-sided-fisher/exploratory/ras_vs_ois/house_env.R")
# sandbox proxy fix from the branch .Rprofile (2nd expression; the 1st activates renv)
eval(parse(".claude/worktrees/two-sided-fisher/.Rprofile")[[2]])
suppressMessages(library(recount3))
options(recount3_url = "https://duffel.rail.bio/recount3")
in_dir <- "data/exploratory/ras_vs_ois"
out_dir <- "results/exploratory/ras_vs_ois/contrast_degs"
f <- file.path(in_dir, "SRP062871_gene_rse.rds")
if (!file.exists(f)) {
  p <- available_projects(organism = "human"); p <- p[p$project == "SRP062871", ]
  saveRDS(create_rse(p, type = "gene", verbose = FALSE), f)
}
rse <- readRDS(f)
cat("samples:\n"); print(table(sub("[ _-]*[0-9]+$", "", rse$sra.sample_title)))
assays(rse)$counts <- round(recount3::transform_counts(rse, by = "auc", targetSize = 40e6))
rd <- as.data.frame(rowData(rse))
rse <- rse[rd$gene_type == "protein_coding" & !is.na(rd$gene_name), ]
m <- rowsum(assay(rse, "counts"), rowData(rse)$gene_name)
# each biological replicate was sequenced as several runs: sum runs per sample title
m <- t(rowsum(t(m), rse$sra.sample_title))
se <- SummarizedExperiment(assays = list(counts = m),
                           colData = S4Vectors::DataFrame(sra.sample_title = colnames(m), row.names = colnames(m)))
se$external_id <- colnames(se)
grp <- sub("[ _-]*(rep)?[ _-]*[0-9]+$", "", se$sra.sample_title, ignore.case = TRUE)
cat("biological replicates per group:\n"); print(table(grp))
for (arm in setdiff(unique(grp), grep("control|vector|ctrl", unique(grp), ignore.case = TRUE, value = TRUE))) {
  ctrl <- grepl("control|vector|ctrl", grp, ignore.case = TRUE)
  s <- se[, c(which(grp == arm), which(ctrl))]
  s$grp <- factor(c(rep("G1", sum(grp == arm)), rep("G2", sum(ctrl))), levels = c("G2", "G1"))
  s <- s[rownames(filter_counts(assay(s, "counts"))), ]
  storage.mode(assay(s, "counts")) <- "integer"
  d <- find_degs_new(rse_obj = s, group_col = "grp", group_1 = "G1", group_2 = "G2")$degs
  name <- paste0("IMR90_N1ICD_", gsub("[^A-Za-z0-9]", "", arm))
  d$contrast <- name
  write.csv(d, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  cat(sprintf("%-28s n=%d vs %d up=%d down=%d\n", name, sum(grp == arm), sum(ctrl),
              sum(d$sig == "y" & d$log2FoldChange > 0), sum(d$sig == "y" & d$log2FoldChange < 0)))
}

# Download STAR counts for TCGA-PAAD and TCGA-SKCM (the most RAS/BRAF-driven cohorts;
# too few normals for the tumour-vs-normal pipeline) for tumour-only driver comparisons.
# Exploratory; run from the pancancer-deg repo root. Output in data/exploratory/tcga_extra/.
eval(parse(".claude/worktrees/two-sided-fisher/.Rprofile")[[2]])  # sandbox proxy fix only
suppressMessages({ library(TCGAbiolinks); library(SummarizedExperiment) })
out <- "data/exploratory/tcga_extra"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
old <- setwd(out); on.exit(setwd(old))
for (p in c("TCGA-PAAD", "TCGA-SKCM")) {
  f <- paste0(p, "_STAR_counts_RSE.rds")
  if (file.exists(f)) next
  q <- GDCquery(project = p, data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
  GDCdownload(q, method = "api", files.per.chunk = 20)
  r <- GDCprepare(q, summarizedExperiment = TRUE)
  saveRDS(r, f)
  cat(p, ncol(r), "samples; types:", paste(names(table(r$shortLetterCode)), table(r$shortLetterCode), collapse = " "), "\n")
}

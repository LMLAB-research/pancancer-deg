# RAS vs OIS (exploratory, local only): list every recount3 sample in the OIS studies
# used by the senescence meta-analysis, so bypass / escape / DDIS samples that the
# meta-analysis excluded can be identified. Run from the pancancer-deg repo root.
suppressMessages(library(recount3))
out_dir <- "data/exploratory/ras_vs_ois"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(recount3_url = "https://duffel.rail.bio/recount3")

studies <- c("SRP065206", "SRP123346", "SRP172671", "SRP062872", "SRP127037",
             "SRP066947", "SRP117883", "SRP066917", "SRP113329", "SRP113324",
             "SRP046254", "SRP060598", "SRP136727", "SRP070636", "SRP017142", "SRP040243")
proj <- available_projects(organism = "human")
proj <- proj[proj$project %in% studies, ]
cat("found in recount3:", nrow(proj), "of", length(studies), "\n")
meta <- do.call(rbind, lapply(seq_len(nrow(proj)), function(i) {
  rse <- create_rse(proj[i, ], type = "gene", verbose = FALSE)
  cd <- as.data.frame(colData(rse))
  keep <- c("external_id", "study", "sra.sample_title", "sra.sample_attributes", "sra.experiment_title")
  cd <- cd[, intersect(keep, colnames(cd))]
  saveRDS(rse, file.path(out_dir, paste0(proj$project[i], "_gene_rse.rds")))
  cd
}))
write.csv(meta, file.path(out_dir, "recount3_ois_studies_samples.csv"), row.names = FALSE)
print(table(meta$study))

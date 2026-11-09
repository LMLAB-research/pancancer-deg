# Pathway gene sets used alongside the RIS gene classes (P_*), metabolic score sets (P_GLYC/GAG/HEXOSAMINE) and the
# house Hallmark + lysosomal sets (HS_*, house_sets()). (TGF-beta, cGAS-STING, interferon, NF-kB, IL6-STAT3,
# p53, NOTCH, KRAS-up as a positive control). Sourced by 07, 09, 11, 13, 14, 15 and ras_vs_myc/02.
# msigdbr from the senescence_stress bulk library (read-only). Exploratory.
# cGAS-STING has no Hallmark set: KEGG "cytosolic DNA-sensing pathway" (sensor/adaptor components and
# their direct outputs) and Reactome "STING-mediated induction of host immune responses"; the type-I IFN
# output is read from Hallmark IFN-alpha response.
if (!exists("inp")) source("exploratory/paths.R")   # repo-relative inputs (run from the repo root)
pathway_sets <- function() {
  old <- .libPaths(); on.exit(.libPaths(old))
  add_lib(BULK_LIB)
  h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  k <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  r <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
  g <- function(d, n) { x <- unique(d$gene_symbol[d$gs_name == n]); stopifnot(length(x) > 0); x }
  list(
    P_TGFB = g(h, "HALLMARK_TGF_BETA_SIGNALING"),
    P_cGAS_STING_KEGG = g(k, "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY"),
    P_STING_Reactome = g(r, "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES"),
    P_IFNA = g(h, "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
    P_IFNG = g(h, "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    P_TNFA_NFKB = g(h, "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    P_IL6_STAT3 = g(h, "HALLMARK_IL6_JAK_STAT3_SIGNALING"),
    P_P53 = g(h, "HALLMARK_P53_PATHWAY"),
    P_NOTCH = g(h, "HALLMARK_NOTCH_SIGNALING"),
    P_KRAS_UP = g(h, "HALLMARK_KRAS_SIGNALING_UP")
  ) -> ps
  # metabolic score sets (glycolysis, GAG synthesis; metabolism/00_metabolic_sets.R), appended when present (2026-11-07)
  mf <- "results/exploratory/metabolism/metabolic_score_sets.csv"   # written by metabolism/00_metabolic_sets.R
  if (file.exists(mf)) { m <- read.csv(mf); ps <- c(ps, split(m$gene, m$set)[unique(m$set)]) }
  c(ps, house_sets())
}

# House pathway sets (2026-11-09): the 49 Hallmark sets + "Lysosomal Genes" used by the senescence_stress
# meta-analysis (meta_analysis/03 -> rerun_outputs/stress_response_pathways_RERUN.csv; read-only), prefixed HS_.
house_sets <- function() {
  h <- read.csv(SEN$hallmark_sets)
  nm <- paste0("HS_", gsub("[^A-Za-z0-9]+", "_", toupper(sub("^HALLMARK ", "", h$pathway))))
  stopifnot(length(unique(nm)) == 50)
  lapply(split(h$genes, nm), unique)
}

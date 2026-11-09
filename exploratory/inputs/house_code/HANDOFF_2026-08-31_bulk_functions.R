# functions.R
#
# All shared analysis functions for the bulk RNA-seq (recount3 meta-analysis +
# ERP021140 temporal) pipelines. Pure function definitions only — no
# hardcoded paths. Source config.R first (it sets SAVE_DIR_CSV etc.), then
# this file:
#
#   source("R/config.R")
#   source("R/functions.R")
#
# Restructured from the original Final/Scripts/for_github/general_functions.R
# (unmodified copy kept at archive_original/ for reference). Logic is
# byte-identical except: (1) the hardcoded save_dir/save_dir_csv/... path
# assignments that used to live at the top of this file were moved into
# config.R as PROJECT_DIR-relative paths, and (2) this header.

suppressPackageStartupMessages({
  library(ComplexUpset)
  library(data.table)
  library(recount3)
  library(biomaRt)
  library(ggplot2)
  library(ggpubr)
  library(edgeR)
  library(DESeq2)
  library(WGCNA)
  library(tibble)
  library(pheatmap)
  library(gplots)
  library(tidyr)
  library(dplyr)
  library(GeneOverlap)
  library(EnhancedVolcano)
  library(RColorBrewer)
  library(clusterProfiler)
  library(stringi)
  library(rrvgo)
  library(org.Hs.eg.db)
})


#General functions
#Save ggplot
# Exact, version-pinned protein-coding gene dictionary, matching what
# getBM(biotype='protein_coding') against the original apr2020 Ensembl-100
# archive would have returned -- but fetched from Ensembl's plain FTP file
# archive (still live) instead of the retired interactive
# archive+biomart web service. Verified: reproduces 13,681/13,681 (100%) of
# the gene symbols in the one surviving object from the original Ensembl-100
# pipeline run (cq_samples.rds). See config.R's ENSEMBL_PINNED_RELEASE note.
#
# Cached to disk after first download (the GTF is ~47MB) so reruns don't
# re-fetch it.
get_ensembl_release_pc <- function(release = ENSEMBL_PINNED_RELEASE, cache_dir = RERUN_DIR) {
  cache_rds <- file.path(cache_dir, sprintf("human_pc_ensembl%d_exact.rds", release))
  if (file.exists(cache_rds)) {
    return(readRDS(cache_rds))
  }
  gtf_path <- file.path(cache_dir, sprintf("Homo_sapiens.GRCh38.%d.gtf.gz", release))
  if (!file.exists(gtf_path)) {
    url <- sprintf("https://ftp.ensembl.org/pub/release-%d/gtf/homo_sapiens/Homo_sapiens.GRCh38.%d.gtf.gz",
                    release, release)
    message(sprintf("Downloading Ensembl release %d GTF (protein-coding dictionary source)...", release))
    utils::download.file(url, gtf_path, quiet = TRUE, mode = "wb")
  }
  gtf <- readLines(gzfile(gtf_path))
  gtf <- gtf[!grepl("^#", gtf)]
  fields <- strsplit(gtf, "\t")
  is_gene <- sapply(fields, function(x) length(x) >= 3 && x[3] == "gene")
  attr_col <- sapply(fields[is_gene], function(x) x[9])

  extract_attr <- function(attr_str, key) {
    pat <- paste0(key, ' "([^"]+)"')
    m <- regmatches(attr_str, regexpr(pat, attr_str))
    ifelse(nchar(m) > 0, sub(pat, "\\1", m), NA)
  }
  biotypes   <- extract_attr(attr_col, "gene_biotype")
  gene_names <- extract_attr(attr_col, "gene_name")
  gene_ids   <- extract_attr(attr_col, "gene_id")

  keep <- biotypes == "protein_coding" & !is.na(gene_names)
  df <- data.frame(external_gene_name = gene_names[keep], ensembl_gene_id = gene_ids[keep],
                    stringsAsFactors = FALSE)
  saveRDS(df, cache_rds)
  df
}

# Entrez ID cross-reference for the same pinned protein-coding gene list.
# The original pipeline pulled entrezgene_id from the same (now-dead)
# apr2020 biomart call. Ensembl's plain GTF doesn't carry Entrez IDs, so this
# maps gene symbol -> Entrez via org.Hs.eg.db instead -- NOT release-pinned
# the same exact way, but symbol<->Entrez mappings are far more stable across
# time than the protein-coding gene SET itself, so this is a much smaller
# approximation than falling back to a live Ensembl release would be.
get_ensembl_release_pc_entrez <- function(release = ENSEMBL_PINNED_RELEASE, cache_dir = RERUN_DIR) {
  pc <- get_ensembl_release_pc(release = release, cache_dir = cache_dir)
  entrez <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys = pc$external_gene_name,
    column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
  ))
  pc$entrezgene_id <- entrez[pc$external_gene_name]
  pc
}

save_p <- function(plot, file_name, save_dir, p_width = 7, p_height = 7) {
  if(!dir.exists(save_dir)){
    stop('Directory does not exist')
  }
  # check if .png is in the filename, if not add it
  if (!grepl("\\.png$", file_name)) {
    file_name <- paste0(file_name, ".png")
  }
  
  # construct full file path
  file_path <- file.path(save_dir, file_name)
  
  # save the plot
  ggsave(filename = file_path, plot = plot, width = p_width, height = p_height,
         device='png',dpi = 300)
  
  print(paste("Plot saved to:", file_path))
}

#save csv
save_csv <- function(data, file_name, path, row.names = FALSE,
                     merge_ensembl_symbol=FALSE,
                     gene_col='gene') {
  if(merge_ensembl_symbol){
    data=merge_gene_symbol_from_ensembl(df=data,
                                      ensembl_col=gene_col)
  }
  
  if(!dir.exists(path)){
    stop('Directory does not exist')
  }
  
  # check if .csv is in the filename, if not add it
  if (!grepl("\\.csv$", file_name)) {
    file_name <- paste0(file_name, ".csv")
  }
  
  # construct full file path
  file_path <- file.path(path, file_name)
  
  # write the data frame to a csv file
  write.csv(data, file = file_path, row.names = row.names)
  
  print(paste("Data saved to:", file_path))
}

#Function to download and collect studies into one object
download_studies=function(
    studies,
    sra_organism='human'){
  ##make sure organism makes sense
  if(sra_organism!='human'&sra_organism!='mouse'){
    message('Only human and mouse are supported')
    return(NULL)
  }
  projects <- available_projects(organism=sra_organism)
  proj_info <- subset(
    projects,
    project %in% studies & project_type == "data_sources"
  )
  
  rse=c()
  for(i in 1:nrow(proj_info)){
    per_row=proj_info[i,]
    message(which(per_row[['project']]==proj_info[['project']]),'/',nrow(proj_info))
    # temp_row=data.frame(matrix(per_row,ncol=length(per_row)))
    # colnames(temp_row)=colnames(proj_info)
    temp_return=tryCatch(create_rse(per_row, type = "gene"),
                         error=function(x){return(NULL)})
    if(is.null(temp_return)){
      message("Can't download ",per_row[['project']])
      next
    }else{
      if(is.null(rse)){
        rse=temp_return
      }else{
        if(ncol(colData(temp_return))==ncol(colData(rse))){
          rse=cbind(rse,temp_return)
        }else{
          next
        }
      }
    }
  }
  return(rse)
}

#function to calculate pi score to rank genes
rank_genelist <- function(genes,
                          gene_col='gene',
                          pval_col=NULL,
                          log2fc_col=NULL,
                          rank_type = c("pval", "fc", "both"),
                          background=NULL) {
  rank_type=match.arg(rank_type)
  if(!is.null(background)){
    genes=genes[genes[[gene_col]]%in%background,]
  }
  rownames(genes)=NULL
  rownames(genes)=genes[[gene_col]]
  
  if(!is.null(pval_col)){
    genes[[pval_col]]=as.numeric(as.character(genes[[pval_col]]))
    
    #add minimum point value for R
    genes[[pval_col]]=genes[[pval_col]]+.Machine$double.xmin
  }
  if(!is.null(log2fc_col)){
    genes[[log2fc_col]]=as.numeric(as.character(genes[[log2fc_col]]))
    # genes[[log2fc_col]]=genes[[log2fc_col]]+.Machine$double.xmin
  }
  
  # Initialize an empty vector for the ranking score
  ranking_score <- numeric(nrow(genes))
  
  if (rank_type == "pval") {
    # Rank based on p-value (in ascending order) and consider the sign of log2FC
    ranking_score <- -log10(genes[[pval_col]]) * sign(genes[[log2fc_col]])
    rank_label='pval'
  } else if (rank_type == "fc") {
    # Rank based on log2 fold change (considering both up and down regulation)
    ranking_score <- genes[[log2fc_col]]
    rank_label='fc'
  } else if (rank_type == "both") {
    # Rank based on a combination of p-value and log2 fold change
    ranking_score <- genes[[log2fc_col]] * -log10(genes[[pval_col]])
    rank_label='pi'
  }
  
  # Create a data frame with the ranking score
  ranked_genes <- data.frame(gene = rownames(genes), 
                             ranking_score = ranking_score)
  
  if(!is.null(log2fc_col)){
    ranked_genes[[log2fc_col]]=genes[[log2fc_col]]
  }
  log2FC = genes[[log2fc_col]]
  if(!is.null(pval_col)){
    ranked_genes[[pval_col]]=genes[[pval_col]]
  }
  ranked_genes[['rank_method']]=rank_label
  
  # Sort the data frame based on the ranking score
  ranked_genes <- ranked_genes[order(-ranking_score), ]
  return(ranked_genes)
}

#wrapper for rank_genelist
rank_genelist_wrapper=function(genes,
                               condition_col,
                               gene_col='gene',
                               rank_type=c("pval", "fc", "both"),
                               pval_col=NULL,
                               log2fc_col=NULL,
                               background=NULL){
  rank_type <- match.arg(rank_type)
  all_condition=unique(genes[[condition_col]])
  all_rank=do.call('rbind',lapply(all_condition,function(get_rank){
    per_condition_df=genes[genes[[condition_col]]==get_rank,]
    temp_rank=rank_genelist(genes = per_condition_df,
                            gene_col = gene_col,
                            rank_type=rank_type,
                            pval_col = pval_col,
                            log2fc_col = log2fc_col,
                            background = background)
    temp_rank[[condition_col]]=get_rank
    return(temp_rank)
  }))
  return(all_rank)
}

#Function to perform gsea on a single condition
gsea_func=function(rank,
                   rank_gene_col='gene',
                   rank_score_col='ranking_score',
                   minGSSize = 1,
                   maxGSSize = Inf,
                   gene_list,
                   gene_list_gene_col,
                   gene_list_col_subset=NULL){
  #prep gene list
  if(is.null(gene_list_col_subset)){
    gene_list$group='group'
  }else{
    gene_list[['group']]=gene_list[[gene_list_col_subset]]
  }
  gene_list$gene_col=gene_list[[gene_list_gene_col]]
  gene_list=gene_list%>%dplyr::select_('group','gene_col')
  
  #prep rankings
  temp_rankings=rank[[rank_score_col]]
  names(temp_rankings)=rank[[rank_gene_col]]
  temp_rankings=temp_rankings[order(-temp_rankings)]
  enrich_plot=GSEA(geneList = temp_rankings,
                   TERM2GENE=gene_list,
                   pvalueCutoff = 1,
                   pAdjustMethod = 'BH',
                   minGSSize = 1,
                   maxGSSize = Inf)
  enrich_plot=data.frame(enrich_plot)
  return(enrich_plot)
}

#wrapper for gsea_func
gsea_func_groups=function(rank,
                          rank_group_col,
                          rank_gene_col='gene',
                          rank_score_col='ranking_score',
                          minGSSize = 1,
                          maxGSSize = Inf,
                          gene_list,
                          gene_list_gene_col,
                          gene_list_col_subset=NULL){
  all_group=unique(rank[[rank_group_col]])
  do.call('rbind',lapply(all_group,function(per_condition){
    temp_df=rank[rank[[rank_group_col]]==per_condition,]
    
    temp_enrichment=gsea_func(rank = temp_df,
                              rank_gene_col = rank_gene_col,
                              rank_score_col=rank_score_col,
                              minGSSize=minGSSize,
                              maxGSSize=maxGSSize,
                              gene_list = gene_list,
                              gene_list_gene_col=gene_list_gene_col,
                              gene_list_col_subset=gene_list_col_subset)
    temp_enrichment[[rank_group_col]]=per_condition
    return(temp_enrichment)
  }))
}

#wrapper for gsea_func_groups and rank_genelist_wrapper
rank_gsea_wrapper=function(genes,
                           condition_col,
                           gene_col='gene',
                           pval_col=NULL,
                           log2fc_col=NULL,
                           background=NULL,
                           rank_type = c("pval", "fc", "both"),
                           gene_list,
                           gene_list_gene_col='gene',
                           gene_list_col_subset=NULL,
                           minGSSize=1,
                           maxGSSize=Inf){
  results=c()
  rank_type <- match.arg(rank_type)
  all_rankings=rank_genelist_wrapper(genes = genes,
                                     condition_col = condition_col,
                                     gene_col = gene_col,
                                     rank_type=rank_type,
                                     pval_col = pval_col,
                                     log2fc_col = log2fc_col,
                                     background=background)
  rank_groups=gsea_func_groups(rank=all_rankings,
                               rank_gene_col='gene',
                               rank_group_col=condition_col,
                               gene_list=gene_list,
                               gene_list_gene_col=gene_list_gene_col,
                               gene_list_col_subset=gene_list_col_subset,
                               minGSSize = minGSSize,
                               maxGSSize = maxGSSize)
  
  rank_groups$p.adjust=p.adjust(rank_groups$pvalue,'BH')
  results[['rankings']]=all_rankings
  results[['gsea_enrichment']]=rank_groups
  return(results)
}

#function to plot results from rank_gsea_wrapper$gsea_enrichment
plot_gsea_enrichment=function(enrichment_df,
                         condition_col,
                         pathway_col,
                         pval_col='p.adjust',
                         factor_condition_levels_x=NULL,
                         factor_condition_levels_facet=NULL,
                         tilt_x=90){
  enrichment_df[['pathway_col']]=enrichment_df[[pathway_col]]
  enrichment_df[['pval_col']]=log2(enrichment_df$p.adjust)*-1
  enrichment_df=enrichment_df%>%rstatix::add_significance(pval_col,
                                                          output.col = 'p.adjust.signif')
  enrichment_df$p.adjust.signif[enrichment_df$p.adjust.signif=='ns']=''
  if(!is.null(factor_condition_levels_x)){
    enrichment_df[[condition_col]]=factor(enrichment_df[[condition_col]],
                                          levels=factor_condition_levels_x)
  }
  if(!is.null(factor_condition_levels_facet)){
    enrichment_df[['pathway_col']]=factor(enrichment_df[['pathway_col']],
                                          levels=factor_condition_levels_facet)
  }
  
  plot_p=enrichment_df%>%ggplot(aes_string(y='enrichmentScore',
                                           x=condition_col,fill='pval_col'))+
    geom_bar(stat='identity',colour='black')+
    facet_wrap(~pathway_col)+
    scale_fill_gradient(expression('-log'[2]*'(p-val)'),low='white', high = "red2",
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks = TRUE, 
                                               ticks.colour = 'black',
                                               title.position = 'top',
                                               title.hjust=0.5))+
    geom_text(aes(label=p.adjust.signif),nudge_y = 0.1)+
    geom_hline(yintercept = 0)
  
  if(tilt_x==90){
    plot_p=plot_p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  }else if(tilt_x==45){
    plot_p=plot_p+theme(axis.text.x=element_text(angle=45,hjust=1))
  }
  
  return(plot_p)
}

#Function to make volcano plot
volcano_function=function(degs,
                          gene_col='external_gene_name',
                          custom_labs=NULL,
                          title=NULL,
                          not_sig_black=TRUE,
                          fc_cutoff=log2(1.5)){
  volcano_labs=degs[[gene_col]]
  if(is.null(custom_labs)){
    volcano_p=EnhancedVolcano(degs,x='log2FoldChange',y='padj',
                              lab=volcano_labs,subtitle=NULL,
                              title=title,
                              FCcutoff=fc_cutoff)
  }else{
    degs[[custom_labs]]=as.character(degs[[custom_labs]])
    temp_lab=degs[[custom_labs]]
    unique_group=unique(temp_lab)
    newCols <- colorRampPalette(grDevices::rainbow(length(unique_group)))
    mycolors <- newCols(length(unique_group))
    colour_key=data.frame(cbind(unique_group,mycolors))
    all_colours=do.call(c,lapply(temp_lab,function(relevant_colour){
      colour_key[colour_key[['unique_group']]==relevant_colour,][['mycolors']]
    }))
    if(not_sig_black){
      degs_temp=degs
      degs_temp$colour=all_colours
      degs_temp[degs_temp[['sig']]=='n',][['colour']]='lightgrey'
      colour_key=rbind(colour_key,c('NS','lightgrey'))
      all_colours=degs_temp$colour
    }
    for(i in 1:nrow(colour_key)){
      colours_temp=colour_key[i,]
      names(all_colours)[all_colours==colours_temp[['mycolors']]]=colours_temp[['unique_group']]
    }
    
    volcano_p=EnhancedVolcano(degs,x='log2FoldChange',y='padj',
                              lab=volcano_labs,subtitle=NULL,
                              title=title,colCustom=all_colours,
                              FCcutoff=fc_cutoff)
  }
}

# Function to remove transcript from ensembl ID
clean_id=function(dirty_id){
  do.call(c,lapply(dirty_id,function(per_id){
    unlist(strsplit(per_id,split='[.]'))[1]
  }))
}

replace_ensembl_with_gene_names <- function(rse, biomart_df, 
                                            ensembl_col_in_biomart="ensembl_gene_id", 
                                            gene_name_col_in_biomart="external_gene_names") {
  gene_duplicates=rownames(rse)[duplicated(rownames(rse))]
  rse=rse[!rownames(rse)%in%gene_duplicates,]
  # Extract the rowRanges from the RangedSummarizedExperiment object
  row_ranges <- rowRanges(rse)
  
  # Convert to a data.frame
  row_ranges_df <- as.data.frame(row_ranges)
  
  # Add a column for the row names (Ensembl IDs)
  row_ranges_df$ensembl_id <- row.names(row_ranges_df)
  
  # Left join with biomart dataframe
  merged_df <- left_join(row_ranges_df, biomart_df, 
                         by = setNames(ensembl_col_in_biomart, "ensembl_id"))
  
  # Remove the Ensembl column
  merged_df <- dplyr::select(merged_df, -one_of("ensembl_id"))
  
  # Rename gene name column to "ensembl_id"
  renamed_df <- rename(merged_df, ensembl_id = !!sym(gene_name_col_in_biomart))
  
  # Replace the rowRanges in the original RangedSummarizedExperiment object
  rowRanges(rse) <- makeGRangesFromDataFrame(renamed_df, keep.extra.columns = TRUE)
  
  # Set the rownames of the RangedSummarizedExperiment object to be the gene names
  row.names(rse) <- renamed_df$ensembl_id
  
  return(rse)
}

#function to filter RSE object for relevant samples and genes
process_rse=function(rse,
                     meta_add=NULL,
                     meta_id_col='sample_ID',
                     filter_pc=TRUE,
                     filter_meta=TRUE,
                     filter_meta_col=c('external_id','study','sra.study_title',
                                       'sra.library_layout','sra.platform_model',
                                       'recount_project.organism'),
                     ensembl_dictionary=NULL,
                     ensembl_gene_col='ensembl_gene_id',
                     ensembl_symbol_col='external_gene_name',
                     filter_duplicates=TRUE,
                     sample_filter=NULL,
                     low_expression_filter=c('all','per'),
                     scale_counts,
                     min_expression_filter=1,
                     min_perc_sample_filter=0.3,
                     method='cpm',
                     keep_filter=FALSE,
                     replace_ensembl=TRUE,
                     remove_cq=FALSE,
                     condition_col='cell_substate'){
  # filter for samples
  if(!is.null(sample_filter)){
    rse=rse[,colnames(rse)%in%sample_filter]
  }
  
  if((filter_pc|replace_ensembl)&is.null(ensembl_dictionary)){
    message('must give ensembl_dictionary')
    return(NULL)
  }
  #clean gene names
  rownames(rse)=clean_id(rownames(rse))
  
  #Filter for pc genes
  if(filter_pc){
    rse=rse[rownames(rse)%in%ensembl_dictionary[[ensembl_gene_col]],]
  }
  
  if(replace_ensembl){
    rse=replace_ensembl_with_gene_names(rse,
                                        biomart_df = ensembl_dictionary,
                                        ensembl_col_in_biomart = ensembl_gene_col,
                                        gene_name_col_in_biomart = ensembl_symbol_col)
  }
  
  if(!is.null(low_expression_filter)){
    if(low_expression_filter=='all'){
      counts_temp=assay(rse)
      counts_temp=filter_counts(dds=counts_temp,
                                min_expression_filter=min_expression_filter,
                                min_perc_sample_filter=min_perc_sample_filter,
                                method=method)
      genes_temp=rownames(counts_temp)
      samples_temp=colnames(counts_temp)
      rse=rse[rownames(rse)%in%genes_temp,]
      rse=rse[,colnames(rse)%in%samples_temp]
    }else if(low_expression_filter=='per'){
      meta_condition=meta_add[,colnames(meta_add)%in%c(meta_id_col,condition_col)]
      all_condition=unique(meta_condition[[condition_col]])
      counts_temp=assay(rse)
      genes_keep=do.call(c,lapply(all_condition,function(per_cond){
        temp_samples=meta_condition[meta_condition[['cell_substate']]==per_cond,][[meta_id_col]]
        
        counts_temp_assay=counts_temp[,colnames(counts_temp)%in%temp_samples]
        
        counts_temp_assay=filter_counts(dds=counts_temp_assay,
                                        min_expression_filter=min_expression_filter,
                                        min_perc_sample_filter=min_perc_sample_filter,
                                        method=method)
        genes_temp=rownames(counts_temp_assay)
      }))
      genes_keep=unique(genes_keep)
      
      rse=rse[rownames(rse)%in%genes_keep,]
    }
  }
  
  rse=rse[rowSums(assay(rse)) != 0, ]
  #remove duplicates
  if(filter_duplicates){
    gene_duplicates=unique(rownames(rse)[duplicated(rownames(rse))])
    rse=rse[!rownames(rse)%in%gene_duplicates,]
  }
  
  #scale counts
  if(scale_counts){
    temp_counts=transform_counts(rse)
    assay(rse)=temp_counts
  }
  
  if(filter_meta){
    colData(rse)=colData(rse)[,colnames(colData(rse))%in%filter_meta_col]
  }
  
  if(!is.null(meta_add)){
    rse=merge_meta(rse=rse,
                   meta_temp=meta_add,
                   to_merge_id_col=meta_id_col)
  }
  
  if(keep_filter){
    rse_keep=colData(rse)[colData(rse)[['keep']],][['external_id']]
    
    rse=rse[,colnames(rse)%in%rse_keep]
  }
  
  if(remove_cq){
    keep=colData(rse)$cell_state!='CQ'
    rse=rse[,keep]
  }
  
  return(rse)
}

#function to merge metadata into rse object
merge_meta=function(rse,
                    meta_temp,
                    to_merge_id_col='sample_ID'){
  temp_coldata=rownames(colData(rse))
  meta_temp=meta_temp[match(temp_coldata, meta_temp[[to_merge_id_col]]),]
  meta_temp=meta_temp[,colnames(meta_temp)!=to_merge_id_col]
  colData(rse)=cbind(colData(rse),meta_temp)
  return(rse)
}

#Function to filter counts by min counts and samples
filter_counts=function(dds,
                       min_expression_filter=1,
                       min_perc_sample_filter=0.3,
                       method='cpm'){
  if(method=='cpm'){
    keep <- rowSums(cpm(dds) >= min_expression_filter) >= (min_perc_sample_filter * ncol(dds))
  }else if(method=='counts'){
    keep <- rowSums(dds >= min_expression_filter) >= (min_perc_sample_filter * ncol(dds))
  }
  
  dds_keep <- dds[keep,]
  return(dds_keep)
}

#function to build model matrix and deseq()
deseq_studies=function(study_obj,
                       protect_col,
                       fix_col=NULL){
  results=c()
  if(is.null(fix_col)&!is.null(protect_col)){
    temp_formula=as.formula(paste0('~',protect_col))  
  }else if(!is.null(fix_col)&!is.null(protect_col)){
    temp_formula=as.formula(paste0('~',paste0(fix_col,sep='+'),protect_col))  
  }else if(!is.null(fix_col)&is.null(protect_col)){
    temp_str_1=paste0(fix_col,sep='+')
    temp_str_all=substr(temp_str_1,1,nchar(temp_str_1)-1)
    temp_formula=as.formula(paste0('~',temp_str_all,protect_col))
  }else if(is.null(fix_col)&is.null(protect_col)){
    temp_formula=as.formula('~1')
  }
  message('Using formula: ',temp_formula)
  dds_temp=DESeqDataSet(study_obj,design=temp_formula)
  
  temp_df=data.frame(colData(dds_temp))
  mm_return=model.matrix(data=temp_df,object=temp_formula)
  dds_temp=DESeq(dds_temp)
  results[['deseq_obj']]=dds_temp
  results[['mm']]=mm_return
  return(results)
}

#Function to normalise and batch-correct counts for plotting
normalise_counts=function(dds_obj_use,
                          batch_col=NULL,
                          protect_col=NULL,
                          method=c('combat','limma','wgcna',NULL),
                          normalisation=c('vst_full','vst','norm'),
                          blind=TRUE,
                          fitToSamples=NULL,
                          vst_n=1000){
  dds_obj_temp=dds_obj_use
  dds_transformed=NULL
  if(!is.null(method)){
    if(method=='combat'){
      message('Applying combat-seq() to raw counts')
      dds_transformed=combat_transform(dds_obj_use=dds_obj_temp,
                                       batch_col=batch_col,
                                       protect_col=protect_col)
    }
  }
  
  if(is.null(dds_transformed)){
    dds_transformed=dds_obj_temp
  }
  
  if(normalisation=='vst_full'){
    message('VST normalising all genes')
    if(!is.integer(assay(dds_transformed)[1])){
      dds_transformed=counts_to_integer(dds_transformed)
    }
    
    normalised_counts=varianceStabilizingTransformation(dds_transformed,
                                                        blind = blind)
  }else if(normalisation=='vst'){
    message(paste0('VST normalising ', vst_n, ' genes'))
    if(!is.integer(assay(dds_transformed)[1])){
      dds_transformed=counts_to_integer(dds_transformed)
    }
    normalised_counts=vst(dds_transformed,
                          blind = blind,
                          nsub=vst_n)
  }else if(normalisation=='norm'){
    message('Applying normTransform')
    normalised_counts <- normTransform(dds_transformed)
  }
  
  if(!is.null(method)&!is.null(batch_col)){
    if(method=='limma'){
      message('Applying removeBatchEffect to normalised counts')
      normalised_counts=limma_transform(dds_obj_use=normalised_counts,
                                        batch_col=batch_col,
                                        protect_col=protect_col)
    }else if(method=='wgcna'){
      message('Applying empiricalBayesLM() to normalised counts')
      if(!is.null(protect_col)){
        normalised_counts=wgcna_transform(dds_obj_use=normalised_counts,
                                          batch_col=batch_col,
                                          protect_col=protect_col,
                                          fitToSamples=fitToSamples)
      }
      
    }
  }
  return(normalised_counts)
}

#Function to PCA --> Can remove batch effect for plotting
pca_rds=function(deseq_obj,
                 pca_plot=c('study','cell_line',
                            'cell_type','tissue','cell_state','cell_substate',
                            'sra.platform_model'),
                 ntop=1000,
                 plot_title=NULL,
                 null_pca='cell_state',
                 remove_pca=NULL,
                 nrow=2,
                 legend_nrow=4){
  if(is.null(ntop)){
    ntop=nrow(deseq_obj)
  }
  if(!is.null(remove_pca)){
    pca_plot=pca_plot[!pca_plot%in%remove_pca]
  }
  relevant_pca=do.call(c,lapply(pca_plot,function(per_pca){
    if(length(unique(deseq_obj[[per_pca]]))>1){
      return(per_pca)
    }else{
      return(NULL)
    }
  }))
  
  if(is.null(relevant_pca)){
    relevant_pca=null_pca
  }
  
  temp_test=lapply(relevant_pca,function(per_pca){
    p_temp=custom_pca(count_data=assay(deseq_obj),
                      ntop=ntop,
                      pheno_data=colData(deseq_obj),
                      intgroup=per_pca)
    p_temp=p_temp+
      theme(legend.position='bottom',
            legend.title=element_blank(),
            axis.text=element_text(size=15),
            axis.title=element_text(size=15),
            legend.text = element_text(size=15))+
      guides(color=guide_legend(nrow=legend_nrow, byrow=TRUE))
    return(p_temp)
  })
  p1=cowplot::plot_grid(plotlist=temp_test,nrow = nrow)
}

#Function to plot PCA, derived from https://github.com/mikelove/DESeq2/blob/master/R/plots.R
custom_pca=function(count_data,
                    ntop=500,
                    pheno_data,
                    intgroup,
                    hide_legend=FALSE){
  # calculate the variance for each gene
  rv <- rowVars(count_data)
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(count_data[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(pheno_data))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(pheno_data[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    pheno_data[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=pheno_data[['external_id']])
  
  # if (returnData) {
  #   attr(d, "percentVar") <- percentVar[1:2]
  #   return(d)
  # }
  # 
  p=ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
  
  if(hide_legend){
    p=p+theme(legend.position = "none")
  }
  return(p)
}

#wgcna
wgcna_transform=function(dds_obj_use,
                         batch_col,
                         protect_col,
                         fitToSamples=NULL){
  temp_coldata=colData(dds_obj_use)
  transposed_temp=t(assay(dds_obj_use))
  temp_coldata=temp_coldata[match(temp_coldata$external_id,rownames(transposed_temp)),]
  if(!is.null(fitToSamples)){
    fit_samples_temp=temp_coldata[[protect_col]]==fitToSamples
  }else{
    fit_samples_temp=NULL
  }
  if(!is.null(protect_col)){
    temp_counts=empiricalBayesLM(data=transposed_temp,
                                 removedCovariates=temp_coldata[[batch_col]],
                                 retainedCovariates=temp_coldata[[protect_col]],
                                 fitToSamples=fit_samples_temp)[['adjustedData']]
  }else{
    temp_counts=empiricalBayesLM(data=transposed_temp,
                                 removedCovariates=temp_coldata[[batch_col]],
                                 retainedCovariates=NULL,
                                 fitToSamples=fit_samples_temp)[['adjustedData']]
  }
  
  temp_counts_transposed=t(temp_counts)
  if(max(temp_counts_transposed)>.Machine$integer.max){
    message('Dividing counts by 2')
    temp_counts_transposed=temp_counts_transposed/2
  }
  
  assay(dds_obj_use)=temp_counts_transposed
  return(dds_obj_use)
}

find_degs_new=function(rse_obj,
                       group_col,
                       group_1,
                       group_2='Proliferating',
                       pval_cutoff=0.05,
                       log2fc_cutoff=log2(1.5),
                       independentFiltering=TRUE,
                       batch_col=NULL,
                       heatmap_cols = c('tissue','study','cell_substate','cell_state'),
                       find_degs=TRUE,
                       sample_col='external_id'){
  rse_obj_temp=rse_obj
  rse_obj_temp=rse_obj_temp[,rse_obj_temp[[group_col]]%in%c(group_1,group_2)]
  results=c()
  temp_coldata=colData(rse_obj_temp)
  temp_coldata_df=data.frame(temp_coldata)
  #table temp
  if(!is.null(batch_col)){
    temp_coldata_df_relevant=temp_coldata_df[,colnames(temp_coldata_df)%in%
                                               c(group_col,batch_col,sample_col)]
  }else{
    temp_coldata_df_relevant=temp_coldata_df[,colnames(temp_coldata_df)%in%
                                               c(group_col,sample_col)]
  }
  temp_coldata_df_relevant=temp_coldata_df_relevant[temp_coldata_df_relevant[[group_col]]%in%
                                                      c(group_1,group_2),]
  sample_relevant=temp_coldata_df_relevant[[sample_col]]
  
  rse_obj_temp=rse_obj_temp[,colnames(rse_obj_temp)%in%sample_relevant]
  
  dds_obj=deseq_studies(study_obj=rse_obj_temp,
                        protect_col=group_col,
                        fix_col=batch_col)
  
  # dds_deseq=tryCatch(DESeq(dds_obj[['deseq_obj']]),
  #                    error=function(x)return(NULL))
  if(!find_degs){
    return(dds_obj)
  }
  if(is.null(dds_obj)){
    return(NULL)
  }
  
  contrast_matrix=c(group_col,as.character(group_1),as.character(group_2))
  
  all_degs=results(dds_obj$deseq_obj,
                   contrast_matrix,
                   independentFiltering = independentFiltering,
                   alpha = pval_cutoff)
  
  all_degs <- as.data.frame(all_degs)
  all_degs$sig=ifelse(abs(all_degs$log2FoldChange)>log2fc_cutoff&
                        all_degs$padj<pval_cutoff,'y','n')
  all_degs$direction=ifelse(all_degs$log2FoldChange<0,'down','up')
  all_degs=all_degs%>%rownames_to_column('gene')
  
  #remove outliers
  all_degs=all_degs[!is.na(all_degs$padj),]
  
  all_degs$group_col=group_col
  all_degs$group_1=group_1
  all_degs$group_2=group_2
  
  all_degs=clean_degs(degs_obj=all_degs)
  all_degs$dir_accession=clean_name(all_degs$dir_accession)
  
  results[['degs']]=all_degs
  results[['dds_deseq']]=dds_obj
  return(results)
}

#Wrapper for cleaning function
clean_degs=function(degs_obj){
  degs_obj_temp_test=clean_col_names(degs_obj)
  degs_obj_temp_test=degs_accession(degs_obj_temp_test)
  degs_obj_temp_test=degs_direction(degs_obj_temp_test)
  degs_obj_temp_test=degs_accession(degs_obj_temp_test,cols_accession=c('group_1_dir','group_2_dir'),
                                    accession_col_name='dir_accession')
  return(degs_obj_temp_test)
}

#Function to merge in gene nmaes
merge_gene_symbol_from_ensembl=function(df,
                                        ensembl_col){
  gene_order=unique(as.character(df[[ensembl_col]]))
  colnames(df)[colnames(df)==ensembl_col]='ensembl'
  df_merged=merge(df,human_pc,all.x=TRUE,by.x='ensembl',by.y='ensembl_gene_id')
  df_merged[["ensembl"]]=factor(df_merged[["ensembl"]],levels=gene_order)
  df_merged=df_merged%>%dplyr::arrange(ensembl)
  df_merged[['ensembl']]=as.character(df_merged[['ensembl']])
  return(df_merged)
}

#Function to clean group names in degs
clean_col_names=function(degs,
                         cols_clean=c('group_1','group_2')){
  for(i in 1:length(cols_clean)){
    degs[[cols_clean[i]]]=clean_name(degs[[cols_clean[i]]])
  }
  return(degs)
}

#Function to add accession to degs
degs_accession=function(degs,
                        cols_accession=c('group_1','group_2'),
                        accession_col_name='accession'){
  degs_temp=degs
  degs_temp[[accession_col_name]]=apply(degs_temp,MARGIN=1,function(per_row){
    all_cols=do.call(c,lapply(cols_accession,function(per_col){
      per_row[[per_col]]
    }))
    all_cols=paste0(all_cols,collapse=' vs\n')
  })
  return(degs_temp)
}

#Function to add direction to degs
degs_direction=function(degs){
  degs_temp=degs
  direction_vector=c('up','down')
  #Dir 1 group 1
  degs_temp[['group_1_dir']]=apply(degs_temp,MARGIN=1,function(dir_1){
    paste0(dir_1[['group_1']],'_',dir_1[['direction']])
  })
  degs_temp[['direction_2']]=apply(degs_temp,MARGIN=1,function(dir_2){
    direction_vector[direction_vector!=dir_2[['direction']]]
  })
  degs_temp[['group_2_dir']]=apply(degs_temp,MARGIN=1,function(dir_2){
    paste0(dir_2[['group_2']],'_',dir_2[['direction_2']])
  })
  colnames(degs_temp)[colnames(degs_temp)=='direction']='direction_1'
  return(degs_temp)
}

#Function to remove '-' and ' ' from a character vector
clean_name=function(name){
  name_temp=as.character(name)
  name_temp=gsub(name_temp,pattern='-',replacement='_')
  name_temp=gsub(name_temp,pattern=' ',replacement='_')
  name_temp=gsub(name_temp,pattern='\n',replacement='_')
  return(name_temp)
}

#Function to build heatmap
build_heatmap=function(
    degs,
    dds_obj,
    batch_col,
    protect_col,
    groups,
    group_col,
    method='wgcna',
    normalisation,
    group_col_heatmap=c('cell_state','cell_substate',
                        'study','tissue'),
    show_rownames=TRUE,
    show_colnames=FALSE,
    scale='row',
    colour=greenred(75),
    top_x_perc=NULL,
    top_x=50,
    vst_n=1000,
    gene_col='ensembl_gene_id',
    rownames_to_symbol=TRUE
){
  dds_obj=dds_obj[,dds_obj[[group_col]]%in%groups]
  normalised_corrected_counts_arrest=normalise_counts(
    dds_obj_use = dds_obj,
    batch_col=batch_col,
    protect_col=protect_col,
    method=method,
    normalisation=normalisation,
    vst_n = vst_n,
    blind = TRUE)
  conditions=unique(degs[['dir_accession']])
  relevant_genes_heatmap=do.call('rbind',lapply(conditions,function(per_condition){
    per_condition_temp=degs[degs$dir_accession==per_condition,]
    sig_degs_use_heatmap=pi_score_degs(sig_degs = per_condition_temp,
                                       top_x_perc=top_x_perc,
                                       top_x=top_x)[['filter_pi']]
    return(sig_degs_use_heatmap)
  }))
  relevant_genes_heatmap_genes=unique(relevant_genes_heatmap$gene)
  
  pheat_return=pheat_counts(dds_obj_use=normalised_corrected_counts_arrest,
                            group_col=group_col_heatmap,
                            sig_degs=relevant_genes_heatmap_genes,
                            show_rownames=show_rownames,
                            show_colnames=show_colnames,
                            scale=scale,
                            gene_col=gene_col,
                            rownames_to_symbol = rownames_to_symbol,
                            colour = colour)
  return(pheat_return)
}

#Function to find 'pi-score' by multiplying log10(adj-pval)*log2(fold change)
pi_score_degs=function(sig_degs,
                       top_x_perc=NULL,
                       top_x=NULL){
  results=c()
  min_adj=min(sig_degs[['padj']][sig_degs[['padj']]!=0])*10^-1
  sig_degs[['padj']]=sig_degs[['padj']]+min_adj
  sig_degs[['pi_score']]=abs(sig_degs[['log2FoldChange']])*-log(sig_degs[['padj']],
                                                                base=10)
  results[['all_pi']]=sig_degs
  
  if(!is.null(top_x_perc)){
    sig_degs=sig_degs%>%
      dplyr::filter(sig_degs[['pi_score']] >= quantile(sig_degs[['pi_score']], top_x_perc))
  }else if(!is.null(top_x)){
    sig_degs=sig_degs%>%
      dplyr::arrange(desc(pi_score))
    sig_degs=sig_degs[1:top_x,]
  }
  results[['filter_pi']]=sig_degs
  return(results)
}

#function to pheatmap normalised counts
pheat_counts=function(dds_obj_use,
                      sig_degs,
                      group_col,
                      scale='none',
                      rownames_to_symbol=TRUE,
                      show_rownames=TRUE,
                      show_colnames=FALSE,
                      gene_col='ensembl_gene_id',
                      colour=colorRampPalette(rev(brewer.pal(n = 7, 
                                                             name = "RdYlBu")))(100)){
  degs_unique=unique(sig_degs)
  gene_subset=assay(dds_obj_use)[rownames(dds_obj_use)%in%degs_unique,]
  
  #Need to do this in case group_col is just one column
  rownames_temp=rownames(colData(dds_obj_use))
  df <- as.data.frame(colData(dds_obj_use)[,c(group_col)])
  rownames(df)=rownames_temp
  
  df=df[match(rownames(df),colnames(dds_obj_use)),]
  df=data.frame(df)
  rownames(df)=colnames(dds_obj_use)
  
  colnames(df)=group_col
  
  list_temp=c()
  for(i in 1:length(group_col)){
    per_col=group_col[i]
    unique_group=unique(df[[per_col]])
    newCols <- colorRampPalette(grDevices::rainbow(length(unique_group)))
    mycolors <- newCols(length(unique_group))
    names(mycolors) <- unique(unique_group)
    mycolors <- list(mycolors)
    names(mycolors)=per_col
    list_temp=c(list_temp,mycolors)
  }
  
  if(rownames_to_symbol){
    human_pc_temp=human_pc[human_pc[[gene_col]]%in%rownames(gene_subset),]
    human_pc_temp=human_pc_temp[order(match(human_pc_temp[[gene_col]],
                                            rownames(gene_subset))), ]
    rownames(gene_subset)=human_pc_temp[['external_gene_name']]
  }
  
  heatmap_p=pheatmap(gene_subset, cluster_rows=TRUE, show_rownames=show_rownames,
                     show_colnames=show_colnames,cluster_cols=TRUE,
                     annotation_col=df,scale=scale,
                     annotation_colors=list_temp,
                     border_color =NA,color=colour,fontsize = 12)
  return(heatmap_p)
}

#Function to save pheatmap object
save_pheatmap=function(pheatmap_to_save,
                       file_name,
                       save_dir,
                       p_width=400,
                       p_height=400,
                       p_dpi=300){
  file_name_temp=paste0(save_dir,file_name,'.png')
  png(filename=file_name_temp, width = p_width, height = p_height,res=p_dpi)
  print(pheatmap_to_save)
  dev.off()
  message('File saved as ',file_name_temp)
}

#Upset plot
upsetR_plot_intersection = function(arrest_degs,
                                    dir_col = 'group_1_dir',
                                    gene_col = 'gene',
                                    dir = NULL,
                                    min_intersect_size = 100,
                                    mode='intersect') {
  # Filter if dir is not null
  if(!is.null(dir)){
    arrest_degs <- arrest_degs[grepl(dir, arrest_degs[[dir_col]]), ]
  }
  
  arrest_degs=arrest_degs[,colnames(arrest_degs)%in%
                            c(gene_col,dir_col)]
  
  arrest_degs=unique(arrest_degs)
  
  # Pivot the data to a wide format
  arrest_degs_wide <- arrest_degs %>%
    mutate(value = 1) %>% # add a column of 1s to spread
    spread(key = all_of(dir_col), value = 'value', fill = 0) # spread to wide format
  
  arrest_degs_wide=arrest_degs_wide[,colnames(arrest_degs_wide)!=gene_col]
  
  # Generate upset plot
  upset_plot <- upset(
    arrest_degs_wide,
    intersect = colnames(arrest_degs_wide),
    mode = mode)
  
  return(upset_plot)
}

# 
# upset_plot_all_conditions <- function(df,
#                                       dir_col = 'group_1_dir',
#                                       gene_col = 'external_gene_name') {
#   
#   df = df[, colnames(df) %in% c(gene_col, dir_col)]
#   
#   df = unique(df)
#   
#   # Pivot the data to a wide format
#   df_wide <- df %>%
#     mutate(value = 1) %>% # add a column of 1s to spread
#     spread(key = all_of(dir_col), value = 'value', fill = 0) # spread to wide format
#   
#   df_wide=df_wide[,2:ncol(df_wide)]
#   # Filter for the genes that change in all conditions
#   df_wide <- df_wide[rowSums(df_wide) == length(unique(df$group_1))/2, ]
#   
#   upset_plot <- upset(
#     df_wide,
#     intersect = colnames(df_wide),
#     min_size = 0,
#     mode='intersect'
#   )
#   
#   return(upset_plot)
# }

#Function to scramble degs
scramble_list=function(list_to_scramble,
                       scramble_name='gene',
                       n_scramble,
                       list_name='db',
                       tag){
  temp_frame=data.frame(sample(list_to_scramble,size=n_scramble))
  colnames(temp_frame)[1]=scramble_name
  temp_frame[[list_name]]=tag
  return(temp_frame)
}

simulate_overlaps = function(relevant_degs,
                             facet_col,
                             ingroup_col = 'group_1_dir',
                             outgroup_col = 'group_2_dir',
                             outgroup_relevant = c('Proliferating_up', 'Proliferating_down'),
                             gene_col,
                             seed = 1,
                             simulation_n = 10000,
                             sig_col='sig',
                             sig_label='y') {
  set.seed(seed)
  
  relevant_degs_to_use=relevant_degs[,colnames(relevant_degs)%in%
                                       c(gene_col,facet_col,ingroup_col,outgroup_col,sig_col)]
  
  relevant_degs_to_use=unique(relevant_degs_to_use)
  
  # Convert data.frame to data.table for efficient operations
  relevant_degs_func = data.table(relevant_degs_to_use)
  relevant_degs_func = relevant_degs_func[relevant_degs_func[[outgroup_col]] %in% outgroup_relevant,]
  
  # Run simulations for each outgroup
  sim_return = rbindlist(lapply(outgroup_relevant, function(per_outgroup) {
    message(per_outgroup)
    
    # Subset data for this outgroup
    relevant_degs_func_out = relevant_degs_func[get(outgroup_col) == per_outgroup,]
    unique_ingroup = unique(relevant_degs_func_out[[ingroup_col]])
    
    # Calculate counts for each ingroup
    ingroup_n = rbindlist(lapply(unique_ingroup, function(per_ingroup) {
      relevant_degs_func_in = relevant_degs_func_out[get(ingroup_col) == per_ingroup]
      facet_return = unique(relevant_degs_func_in[[facet_col]])
      return(list(db = per_ingroup, facet = facet_return, n = sum(relevant_degs_func_in[[sig_col]] == sig_label)))
    }))
    
    ingroup_n$outgroup = per_outgroup
    
    # Skip this outgroup if not enough genes
    if(nrow(ingroup_n) < 2){
      message('Not enough genes to overlap for ', per_outgroup)
      return(NULL)
    }
    
    all_sim = rbindlist(lapply(seq_len(simulation_n), function(sim_x) {
      # Scramble and count overlaps
      all_scramble = rbindlist(
        lapply(ingroup_n[['db']], function(per_ingroup) {
          n_scramble = as.numeric(ingroup_n[['n']][ingroup_n[['db']] == per_ingroup])
          list_to_scramble = relevant_degs_func[get(ingroup_col) == per_ingroup, get(gene_col)]
          temp_scramble=scramble_list(list_to_scramble,
                                      scramble_name=gene_col,
                                      n_scramble=n_scramble,
                                      list_name='db',
                                      tag=per_ingroup)
        }))
      
      sim_overlap = data.frame(table(all_scramble[[gene_col]]))
      sim_overlap_n = data.frame(table(sim_overlap$Freq),stringsAsFactors = FALSE)
      sim_overlap_n[['Var1']]=as.numeric(as.character(sim_overlap_n[['Var1']]))
      
      # Add rows with zero frequency for missing overlaps
      max_overlap = length(unique(ingroup_n[['db']]))
      missing_overlaps = setdiff(seq_len(max_overlap), sim_overlap_n$Var1)
      if(length(missing_overlaps) > 0){
        sim_overlap_n = rbind(sim_overlap_n, 
                              data.table(Var1 = missing_overlaps, Freq = rep(0, length(missing_overlaps))))
      }
      
      sim_overlap_n$sim = sim_x
      return(sim_overlap_n)
    }), fill = TRUE)
    
    all_sim$outgroup = per_outgroup
    return(all_sim)
  }))
  
  return(sim_return)
}

plot_top_simulation_cumsum <- function(df,
                                       Var1_col = "Var1",
                                       Freq_col = "Freq", 
                                       sim_col = "sim",
                                       outgroup_col = "outgroup",
                                       cutoff=0.05,
                                       plot_cutoff=TRUE,
                                       xlab='Overlap n',
                                       ylab='Reverse Overlap Cumulative Frequency',
                                       scale=1) {
  
  setDT(df)  # Convert df to a data.table if it isn't already
  
  df=df[df$Var1==max(df$Var1),]
  
  unique_group=unique(df[[outgroup_col]])
  
  to_plot=c()
  for(i in unique_group){
    df_temp=df[df[[outgroup_col]]==i,]
    table_freq=data.frame(table(df_temp$Freq))
    temp_cumsum=data.frame(rev_cumsum=(sum(table_freq$Freq)-cumsum(table_freq$Freq))/sum(table_freq$Freq))
    temp_cumsum$overlap_n=as.numeric(rownames(temp_cumsum))-1
    temp_cumsum$group=i
    temp_cumsum$accession=min(temp_cumsum$overlap_n[temp_cumsum$rev_cumsum<cutoff])
    to_plot=rbind(to_plot,
                  temp_cumsum)
  }
  
  max_overlap = ceiling(max(to_plot$overlap_n)/scale)*scale
  break_points = seq(0, max_overlap, by = scale)
  
  p=to_plot%>%ggplot(aes(overlap_n,y=rev_cumsum))+
    geom_line()+
    facet_grid(~group)+
    scale_x_continuous(breaks=break_points)+
    xlab(xlab)+
    ylab(ylab)
  if(plot_cutoff){
    p=p+
      geom_vline(aes(xintercept = accession),colour='red')
  }
  results=c()
  results[['p']]=p
  results[['cumsum']]=to_plot
  return(results)
}

#Overlap degs with db
overlap_function <- function(df_1, df_2, gene_col_1, gene_col_2,
                             carry_col_1=NULL,carry_col_2=NULL,
                             group_col_1=NULL, group_col_2=NULL,background) {
  background=unique(background)
  df_1_temp=df_1[df_1[[gene_col_1]]%in%background,]
  df_2_temp=df_2[df_2[[gene_col_2]]%in%background,]
  if(!is.null(group_col_1)){
    df_1_temp[[group_col_1]]=as.character(df_1_temp[[group_col_1]])
  }
  if(!is.null(carry_col_1)){
    for(i in carry_col_1){
      df_1_temp[[i]]=as.character(df_1_temp[[i]])
    }
  }
  if(!is.null(group_col_2)){
    df_2_temp[[group_col_2]]=as.character(df_2_temp[[group_col_2]])
  }
  if(!is.null(carry_col_2)){
    for(i in carry_col_2){
      df_2_temp[[i]]=as.character(df_2_temp[[i]])
    }
  }
  genome.size=length(background)
  # Function to calculate overlaps
  calculate_overlap <- function(df_1_sub, df_2_sub, gene_col_1, gene_col_2,
                                genome.size){
    gene_set_1 <- unique(df_1_sub[[gene_col_1]])
    gene_set_2 <- unique(df_2_sub[[gene_col_2]])
    
    overlap_obj <- testGeneOverlap(newGeneOverlap(gene_set_1, gene_set_2,
                                  genome.size = genome.size))
    intersect_genes <- paste0(getIntersection(overlap_obj),collapse='/')
    
    temp_cont=getContbl(overlap_obj)
    actual=temp_cont[4]
    expected=chisq.test(temp_cont)$exp[4]
    diff=actual-expected
    pval=fisher.test(temp_cont,alternative = 'two.sided')$p.val
    odds=fisher.test(temp_cont)$est
    result <- c(actual, expected, 
                diff,
                pval, odds, intersect_genes)
    names(result) <- c('actual','expected','diff','pval','odds','int')
    
    return(result)
  }
  
  # If group_cols are not null, calculate overlaps by group combinations
  # Getting unique combinations for both dataframes
  df_1_combinations <- df_1_temp %>% dplyr::select(all_of(group_col_1)) %>% distinct() %>% drop_na()
  df_2_combinations <- df_2_temp %>% dplyr::select(all_of(group_col_2)) %>% distinct() %>% drop_na()
  results=c()
  
  for (i in 1:nrow(df_1_combinations)){
    df_1_sub <- df_1_temp
    condition_1=c()
    carry_1=df_1_sub
    carry_1_carry=c()
    for(k in group_col_1){
      df_1_sub <- df_1_sub[df_1_sub[[k]] == df_1_combinations[[k]][i], ]
      condition_1=rbind(condition_1,df_1_combinations[[k]][i])
      condition_filter_1=df_1_sub[df_1_sub[[k]]%in%unique(df_1_sub[[k]]),]
      carry_1_carry=unique(do.call(c,lapply(carry_col_1,function(per_carry){
        temp_carry_1=unique(condition_filter_1[,colnames(condition_filter_1)%in%
                                    per_carry])
        if(length(temp_carry_1)>1){
          stop('carry cols are not unique for df_1')
        }
        return(temp_carry_1)
      })))
    }
    for (j in 1:nrow(df_2_combinations)){
      df_2_sub <- df_2_temp
      condition_2=c()
      carry_2=df_2_sub
      carry_2_carry=c()
      for(l in group_col_2){
        df_2_sub <- df_2_sub[df_2_sub[[l]] == df_2_combinations[[l]][j], ]
        condition_2=rbind(condition_2,df_2_combinations[[l]][j])
        condition_filter_2=df_2_sub[df_2_sub[[l]]%in%unique(df_2_sub[[l]]),]
        carry_2_carry=unique(do.call(c,lapply(carry_col_2,function(per_carry){
          temp_carry_2=unique(condition_filter_2[,colnames(condition_filter_2)%in%
                                                   per_carry])
          if(length(temp_carry_2)>1){
            stop('carry cols are not unique for df_2')
          }
          return(temp_carry_2)
        })))
      }
      result <- calculate_overlap(df_1_sub, df_2_sub,
                                  gene_col_1, gene_col_2,
                                  genome.size=genome.size)
      result <- c(result, condition_1, condition_2,
                  carry_1_carry,carry_2_carry)
      names(result)[names(result)=='']=c(group_col_1,group_col_2,carry_col_1,
                                         carry_col_2)
      
      result[['n1']]=nrow(df_1_sub)
      result[['n2']]=nrow(df_2_sub)
      results <- rbind(results,
                       result)
    }
  }
  results=data.frame(results)
    
  # Adjust p-values using BH correction
  results$actual=as.numeric(results$actual)
  results$expected=as.numeric(results$expected)
  results$diff=as.numeric(results$diff)
  results$pval=as.numeric(results$pval)
  results$odds=as.numeric(results$odds)
  results$adj=p.adjust(results$pval)
  results$n1=as.numeric(results$n1)
  results$n2=as.numeric(results$n2)
  rownames(results)=NULL
  
  if(!is.null(group_col_1)){
    if(is.factor(df_1[[group_col_1]])){
    results[[group_col_1]]=factor(results[[group_col_1]],
                                  levels=levels(df_1[[group_col_1]]))
    }
  }
  
  if(!is.null(group_col_2)){
    if(is.factor(df_2[[group_col_2]])){
    results[[group_col_2]]=factor(results[[group_col_2]],
                                  levels=levels(df_2[[group_col_2]]))
    }
  }
  
  if(!is.null(carry_col_1)){
    for(i in carry_col_1){
      if(!is.factor(df_1[[i]])){
        next
      }
      results[[i]]=factor(results[[i]],
                                       levels=levels(df_1[[i]]))
    }
  }
  
  if(!is.null(carry_col_2)){
    for(i in carry_col_2){
      if(!is.factor(df_2[[i]])){
        next
      }
      results[[i]]=factor(results[[i]],
                                       levels=levels(df_2[[i]]))
    }
  }
  
  results$background_n=genome.size
  
  return(results)
}

#wrapper for overlap_function() that p-value corrects for self-overlaps
self_overlaps_full_grid=function(df,
                                 gene_col,
                                 background,
                                 group_col,
                                 carry_col=NULL){
  self_overlaps=overlap_function(df_1 = df,
                                 df_2 = df,
                                 gene_col_1 = gene_col,
                                 gene_col_2 = gene_col,
                                 background = background,
                                 group_col_1 = group_col,
                                 group_col_2 = group_col,
                                 carry_col_1=carry_col,
                                 carry_col_2=carry_col)
  
  #adj
  n_temp=nrow(self_overlaps)/2
  
  adj_temp=do.call(c,lapply(self_overlaps[['pval']],function(per_pval){
    p.adjust(p = per_pval,n = n_temp)
  }))
  self_overlaps$adj=adj_temp
  self_col=do.call(c,as.list(apply(self_overlaps,MARGIN = 1,function(x){
    ifelse(x[[group_col]]==x[[paste0(group_col,'.1')]],
           TRUE,FALSE)
  })))
  self_overlaps[['self']]=self_col
  
  return(self_overlaps)
}

#Function to remove self overlaps and repeat overlaps
clean_self_overlaps=function(df,
                             group_1,
                             group_2){
  accessions=do.call(c,list(apply(df,MARGIN = 1,function(get_accession){
    x=get_accession[[group_1]]
    y=get_accession[[group_2]]
    accession=paste0(sort(c(x,y)),collapse='_')
  })))
  df$accession=accessions
  df=df[!duplicated(df$accession),]
  df=df[,colnames(df)!='accession']
  
  if('self'%in%colnames(df)){
    df=df[!df$self,]
  }
  return(df)
}

#create simple overlap plot
simple_overlap_plot=function(df,
                             x,
                             y,
                             actual_col='actual',
                             odds_col='odds',
                             pval_col='adj',
                             self=FALSE,
                             x_factor=NULL,
                             y_factor=NULL,
                             rev_y=FALSE,
                             xlab=NULL,
                             ylab=NULL,
                             x_tilt=0,
                             text_size=10,
                             add_almost_sig=FALSE,
                             remove_nonsig=FALSE){
  if(!is.null(x_factor)){
    df[[x]]=factor(df[[x]],
                   levels=x_factor)
  }
  
  if(!is.null(y_factor)){
    df[[y]]=factor(df[[y]],
                   levels=y_factor)
  }
  
  if(rev_y){
    df[[y]]=factor(df[[y]],
                   levels=rev(y_factor))
  }
  
  df[['log2odds']]=log2(df[[odds_col]])
  graph_max = max(ceiling(df$log2odds/2)[!is.infinite(ceiling(df$log2odds))]) * 2
  graph_min = min(floor(df$log2odds/2)[!is.infinite(ceiling(df$log2odds))])*2
  if(graph_min > 0){
    graph_min = 0
  }
  
  li <- c(graph_min, graph_max)
  la <- c(seq(graph_min,graph_max,2))
  br <- c(seq(graph_min,graph_max,2))
  
  df=df%>%rstatix::add_significance(pval_col,output.col = 'adj.signif')
  
  if(self){
    df[['adj.signif']][df[['self']]]=''
  }
  
  if(remove_nonsig){
    df$adj.signif[df$adj.signif=='ns']=''
  }
  if(add_almost_sig){
    df[[paste0(pval_col,'.signif')]][df[[pval_col]]>0.05&
                                       df[[pval_col]]<0.1]='^'
  }
  
  p_temp=df%>%ggplot(aes_string(x=x,y=y,fill='log2odds'))+
    geom_tile(colour='black')+
    geom_text(aes_string(label=actual_col),nudge_y=-0.25,size = text_size / 3)+
    geom_text(aes(label=adj.signif),nudge_y=0.25,size = text_size / 3)+
    scale_fill_gradient2(expression('log'[2]*'(Odds)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks = TRUE, 
                                                ticks.colour = 'black', title.position = 'top', title.hjust=0.5),
                         breaks=br,
                         labels=la,
                         limits=li)+
    xlab(xlab)+
    ylab(ylab)+
    scale_x_discrete(expand=expansion(c(0,0)))+
    scale_y_discrete(expand=expansion(c(0,0)))
  
  if(x_tilt==45) {
    p_temp = p_temp + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = text_size))
  } else if(x_tilt==0) {
    p_temp = p_temp + theme(axis.text.x = element_text(size = text_size))
  }else if(x_tilt==90){
    p_temp = p_temp + theme(axis.text.x = element_text(angle = 90,
                                                       hjust = 1,vjust=0.5, size = text_size))
  }
  
  p_temp + theme(axis.title.y = element_text(size = text_size),
                 axis.title.x = element_text(size = text_size),
                 strip.text = element_text(size = text_size),
                 legend.title = element_text(size = text_size),
                 legend.text = element_text(size = text_size),
                 axis.text.y=element_text(size=text_size))
  
  return(p_temp)
}

create_overlap_plot <- function(deg_db_overlap,
                                odds_column='odds',
                                pval_col='adj',
                                facet_col,
                                facet_2=NULL,
                                x,
                                y,
                                xlab=NULL,
                                ylab=NULL,
                                ggtitle=NULL,
                                remove_y=FALSE,
                                text_size=10,
                                x_tilt=0,
                                add_almost_sig=FALSE,
                                remove_nonsig=FALSE,
                                cols=NULL,
                                self=FALSE) {
  
  deg_db_overlap$log2odds = log2(deg_db_overlap[[odds_column]])
  
  # Determine graph_max and graph_min
  graph_max = max(ceiling(deg_db_overlap$log2odds/2)[!is.infinite(ceiling(deg_db_overlap$log2odds))]) * 2
  graph_min = min(floor(deg_db_overlap$log2odds/2)[!is.infinite(ceiling(deg_db_overlap$log2odds))])*2
  if(graph_min > 0){
    graph_min = 0
  }
  
  li <- c(graph_min, graph_max)
  la <- c(seq(graph_min,graph_max,2))
  br <- c(seq(graph_min,graph_max,2))
  
  # Add significance
  deg_db_overlap = deg_db_overlap %>% rstatix::add_significance(pval_col,
                                                                output.col = 'adj.signif')
  if(remove_nonsig){
    deg_db_overlap$adj.signif[deg_db_overlap$adj.signif=='ns']=''
  }
  if(add_almost_sig){
    deg_db_overlap[[paste0(pval_col,'.signif')]][deg_db_overlap[[pval_col]]>0.05&
                                deg_db_overlap[[pval_col]]<0.1]='^'
  }
  
  # Create plot
  deg_db_overlap[['facet']]=deg_db_overlap[[facet_col]]
  
  if (!is.null(facet_2)) {
    deg_db_overlap[['facet_2']]=deg_db_overlap[[facet_2]]
  }
  
  if(self){
    deg_db_overlap[[paste0(pval_col,'.signif')]][deg_db_overlap[[self]]]=''
  }
  
  p_overlap = deg_db_overlap %>%
    ggplot(aes_string(x=x, y=y)) +
    geom_tile(aes(fill=log2odds),colour='black') +
    scale_fill_gradient2(expression('log'[2]*'(Odds)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks = TRUE, 
                                                ticks.colour = 'black', title.position = 'top', title.hjust=0.5),
                         breaks=br,
                         labels=la,
                         limits=li) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggtitle(ggtitle)+
    geom_text(aes(label=adj.signif), nudge_y=0.25, size = text_size / 3) +
    geom_text(aes(label=actual), nudge_y=-0.25, size = text_size / 3) +
    xlab(xlab)+
    ylab(ylab)
  
  if(x_tilt==45) {
    p_overlap = p_overlap + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = text_size))
  } else if(x_tilt==0) {
    p_overlap = p_overlap + theme(axis.text.x = element_text(size = text_size))
  }else if(x_tilt==90){
    p_overlap = p_overlap + theme(axis.text.x = element_text(angle = 90,
                                                             hjust = 1,vjust=0.5, size = text_size))
  }
  
  # Use facet_grid if facet_2 is provided, else use facet_wrap
  if (!is.null(facet_2)) {
    deg_db_overlap[['facet_2']]=deg_db_overlap[[facet_2]]
    p_overlap = p_overlap + facet_grid(facet ~ facet_2,scales = 'free',
                                       space='free')
  } else {
    if(is.null(cols)){
      p_overlap = p_overlap + facet_wrap(.~facet)
    }else{
      p_overlap = p_overlap + facet_wrap(.~facet,ncol = cols)
    }
    
  }
  
  # Condition to remove y axis labels and tick marks
  p_overlap = p_overlap + theme(axis.title.y = element_text(size = text_size),
                                axis.title.x = element_text(size = text_size),
                                strip.text = element_text(size = text_size),
                                legend.title = element_text(size = text_size),
                                legend.text = element_text(size = text_size))
  
  if(remove_y){
    p_overlap = p_overlap + theme(axis.text.y = element_blank(), 
                                  axis.ticks.y = element_blank())
  }else{
    p_overlap=p_overlap+
      theme(axis.text.y = element_text(size = text_size))
  }
  
  return(p_overlap)
}

#Overlap within a df
overlap_within_df <- function(dataframe,
                              group_col,
                              other_cols, gene_col, background,
                              remove_self=FALSE) {
  background=unique(background)
  dataframe=dataframe[dataframe[[gene_col]]%in%background,]
  genome.size=length(background)
  
  # Function to calculate overlaps
  calculate_overlap <- function(df_1_sub, df_2_sub, gene_col){
    gene_set_1 <- unique(df_1_sub[[gene_col]])
    gene_set_2 <- unique(df_2_sub[[gene_col]])
    
    overlap_obj <- testGeneOverlap(newGeneOverlap(gene_set_1, gene_set_2,
                                                  genome.size = genome.size))
    intersect_genes <- paste0(getIntersection(overlap_obj),collapse='/')
    
    temp_cont=getContbl(overlap_obj)
    actual=temp_cont[4]
    expected=chisq.test(temp_cont)$exp[4]
    diff=actual-expected
    pval=fisher.test(temp_cont,alternative = 'two.sided')$p.val
    odds=fisher.test(temp_cont)$est
    result <- c(actual, expected, 
                diff,
                pval, odds, intersect_genes)
    names(result) <- c('actual','expected','diff','pval','odds','int')
    
    return(result)
  }
  
  # Getting unique combinations for dataframe
  if(is.null(other_cols)){
    df_combinations <- dataframe%>%dplyr::select(group_col) %>% distinct() %>% drop_na()
  }else{
    df_combinations <- dataframe[,colnames(dataframe)%in%
                                   c(group_col,other_cols)] %>% distinct() %>% drop_na()
  }
  
  
  # Cross-join df_combinations with itself
  df_combinations_cross <- merge(df_combinations, df_combinations, by = NULL)
  
  # Remove rows where group_col from both dataframes are the same
  df_combinations_cross <- df_combinations_cross[df_combinations_cross[[paste0(group_col, ".x")]] != 
                                                   df_combinations_cross[[paste0(group_col, ".y")]], ]
  colnames(df_combinations_cross)=gsub(colnames(df_combinations_cross),pattern='.x',replacement='')
  colnames(df_combinations_cross)=gsub(colnames(df_combinations_cross),pattern='.y',replacement='')
  temp_ncol_1=1:(ncol(df_combinations_cross)/2)
  colnames(df_combinations_cross)[temp_ncol_1]=paste0(colnames(df_combinations_cross)[temp_ncol_1],'_1')
  temp_ncol_2=(ncol(df_combinations_cross)/2+1):ncol(df_combinations_cross)
  colnames(df_combinations_cross)[temp_ncol_2]=paste0(colnames(df_combinations_cross)[temp_ncol_2],'_2')
  accession_filter=apply(df_combinations_cross,MARGIN = 1,function(get_accession){
    get_accession=matrix(get_accession,ncol=max(temp_ncol_2))
    access_1=paste0(get_accession[,temp_ncol_1],collapse='_')
    access_2=paste0(get_accession[,temp_ncol_2],collapse='_')
    return(paste0(sort(c(access_1,access_2)),collapse='_'))
  })
  df_combinations_cross$accession=accession_filter
  if(remove_self){
    df_combinations_cross=df_combinations_cross[!duplicated(df_combinations_cross$accession),]
  }
  
  results <- c()
  
  # Iterate over rows of df_combinations_cross
  for (i in 1:nrow(df_combinations_cross)){
    temp_accession=df_combinations_cross$accession[[i]]
    df_1_sub <- dataframe
    df_2_sub <- dataframe
    condition_1 = c()
    condition_2 = c()
    for(k in c(group_col, other_cols)){
      df_1_sub <- df_1_sub[df_1_sub[[k]] == df_combinations_cross[[paste0(k, "_1")]][i], ]
      df_2_sub <- df_2_sub[df_2_sub[[k]] == df_combinations_cross[[paste0(k, "_2")]][i], ]
      condition_1 = rbind(condition_1, df_combinations_cross[[paste0(k, "_1")]][i])
      condition_2 = rbind(condition_2, df_combinations_cross[[paste0(k, "_2")]][i])
    }
    result <- calculate_overlap(df_1_sub, df_2_sub, gene_col)
    result <- c(result, condition_1, condition_2,temp_accession)
    names(result)[names(result)==''] = c(paste0(c(group_col, other_cols), "_1"), paste0(c(group_col, other_cols), "_2"),
                                         'accession')
    
    results <- rbind(results, result)
  }
  
  results = data.frame(results)
  
  # Adjust p-values using BH correction
  results$actual=as.numeric(results$actual)
  results$expected=as.numeric(results$expected)
  results$diff=as.numeric(results$diff)
  results$pval=as.numeric(results$pval)
  results$odds=as.numeric(results$odds)
  results$adj=p.adjust(results$pval)
  
  rownames(results)=NULL
  
  return(results)
}

# plot overlap_within_df
plot_self_overlaps <- function(self_overlaps_recount_to_plot,
                               odds_column = 'odds',
                               x = 'direction_1_1',
                               y = 'outgroup_dir',
                               facet_col = 'group_1_1',
                               xlab = NULL,
                               ylab = NULL,
                               text_size = 10,  # added a text_size parameter with a default value of 10
                               angle_x = FALSE,
                               scale=2) {  # added an angle_x parameter with a default value of FALSE
  
  # Calculate log2odds
  self_overlaps_recount_to_plot$log2odds = log2(self_overlaps_recount_to_plot[[odds_column]])
  
  # Determine graph_max and graph_min
  graph_max = max(ceiling(self_overlaps_recount_to_plot$log2odds/scale)[!is.infinite(ceiling(self_overlaps_recount_to_plot$log2odds))]) * scale
  graph_min = min(floor(self_overlaps_recount_to_plot$log2odds/scale)[!is.infinite(ceiling(self_overlaps_recount_to_plot$log2odds))]) * scale
  if(graph_min > 0){
    graph_min = 0
  }
  
  li <- c(graph_min, graph_max)
  la <- c(seq(graph_min,graph_max,scale))
  br <- c(seq(graph_min,graph_max,scale))
  
  self_overlaps_recount_to_plot = self_overlaps_recount_to_plot %>% rstatix::add_significance('adj')
  if(!is.null(facet_col)){
    self_overlaps_recount_to_plot[['facet_col']]=self_overlaps_recount_to_plot[[facet_col]]
  }
  
  # Create plot
  p_overlap = self_overlaps_recount_to_plot %>%
    ggplot(aes_string(x=x, y=y)) +
    geom_tile(aes(fill=log2odds),colour='black') +
    scale_fill_gradient2(expression('log'[2]*'(Odds)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks = TRUE, 
                                                ticks.colour = 'black', title.position = 'top', title.hjust=0.5),
                         breaks=br,
                         labels=la,
                         limits=li) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    geom_text(aes(label=adj.signif), nudge_y=0.25, size = text_size / 3) +  # adjusted size for geom_text
    geom_text(aes(label=actual), nudge_y=-0.25, size = text_size / 3) +  # adjusted size for geom_text
    xlab(xlab) +
    ylab(ylab)
  
  if(!is.null(facet_col)){
    p_overlap=p_overlap+
      facet_wrap(.~facet_col, scales ='free')
  }
  
  # Adjust x-axis text angle if angle_x is TRUE
  if (angle_x) {
    p_overlap = p_overlap + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = text_size))
  } else {
    p_overlap = p_overlap + theme(axis.text.x = element_text(size = text_size))
  }
  
  # Rest of the text elements
  p_overlap = p_overlap + theme(axis.text.y = element_text(size = text_size),
                                axis.title.y = element_text(size = text_size),
                                axis.title.x = element_text(size = text_size),
                                strip.text = element_text(size = text_size),
                                legend.title = element_text(size = text_size),
                                legend.text = element_text(size = text_size))
  
  return(p_overlap)
}

plot_self_overlaps2=function(test_temp,
                             accession_col='accession',
                             accession_cols_1=c('group_1_1','direction_1_1'),
                             accession_cols_2=c('group_1_2','direction_1_2'),
                             odds_column='odds',
                             actual_n='actual',
                             pval_col='adj',
                             max_odds=6,
                             min_odds=-6){
  test_temp=test_temp%>%rstatix::add_significance(pval_col,output.col = 'pval_label')
  accessions_use=do.call('rbind',apply(test_temp,MARGIN = 1,function(get_accession){
    accession_1_temp=paste0(get_accession[[accession_cols_1[1]]],' ',stringr::str_to_sentence(get_accession[[accession_cols_1[2]]]))
    accession_2_temp=paste0(get_accession[[accession_cols_2[1]]],' ',stringr::str_to_sentence(get_accession[[accession_cols_2[2]]]))
    return(data.frame(matrix(c(accession_1_temp,accession_2_temp),ncol=2)))
  }))
  
  test_temp$accession_1=accessions_use$X1
  test_temp$accession_2=accessions_use$X2
  test_temp$log2odds = log2(test_temp[[odds_column]])
  
  # Determine graph_max and graph_min
  graph_max = max(ceiling(test_temp$log2odds/2)[!is.infinite(ceiling(test_temp$log2odds))]) * 2
  if(!is.null(max_odds)&max_odds<graph_max){
    graph_max=max_odds
    test_temp$log2odds[test_temp$log2odds>max_odds]=max_odds
    cut_max=TRUE
  }else{
    cut_max=FALSE
  }
  
  graph_min = min(floor(test_temp$log2odds/2)[!is.infinite(ceiling(test_temp$log2odds))])*2
  if(!is.null(min_odds)&min_odds>graph_min){
    graph_min=min_odds
    test_temp$log2odds[test_temp$log2odds<min_odds]=min_odds
    cut_min=TRUE
  }else{
    cut_min=FALSE
  }
  if(graph_min > 0){
    graph_min = 0
  }
  
  li <- c(graph_min, graph_max)
  la <- c(seq(graph_min,graph_max,2))
  br <- c(seq(graph_min,graph_max,2))
  
  if(cut_max){
    la[br==max(br)]=paste0('>',graph_max)
  }
  
  if(cut_min){
    la[br==min(br)]=paste0('<',graph_min)
  }
  
  test_temp=test_temp[!duplicated(test_temp[[accession_col]]),]
  
  order_x=rev(sort(table(test_temp$accession_1)))
  order_y=rev(sort(table(test_temp$accession_2)))
  test_temp$accession_1=factor(test_temp$accession_1,
                               levels=names(order_x))
  test_temp$accession_2=factor(test_temp$accession_2,
                               levels=names(order_y))
  test_temp%>%ggplot(aes(x=accession_1,y=accession_2))+
    geom_tile(aes(fill=log2odds),colour='black')+
    scale_fill_gradient2(expression('log'[2]*'(Odds)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks = TRUE, 
                                                ticks.colour = 'black', title.position = 'top', title.hjust=0.5),
                         breaks=br,
                         labels=la,
                         limits=li)+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    geom_text(aes_string(label=actual_n),nudge_y=-0.2)+
    geom_text(aes(label=pval_label),nudge_y = 0.2)+
    scale_x_discrete(expand=expansion(c(0,0)))+
    scale_y_discrete(expand=expansion(c(0,0)))+
    xlab('Condition 1')+
    ylab('Condition 2')
}

extract_up_down <- function(input_vector, pattern="down|up") {
  
  extract_pattern <- function(string) {
    match <- regexpr(pattern, string)
    if (match > 0) {
      extracted <- substr(string, match, match + attr(match, "match.length") - 1)
      return(extracted)
    } else {
      return(NA)
    }
  }
  
  result <- sapply(input_vector, extract_pattern)
  return(result)
}

#Function to find 'pi-score' and extract top x% of genes by group
pi_score_filter_by_group=function(sig_degs,
                                  group_col,
                                  top_x_perc=NULL,
                                  top_x=NULL){
  all_group=unique(sig_degs[[group_col]])
  pi_filtered_by_group=do.call('rbind',lapply(all_group,function(per_group){
    temp_degs=sig_degs[sig_degs[[group_col]]==per_group,]
    pi_return=pi_score_degs(sig_degs=temp_degs,
                            top_x_perc=top_x_perc,
                            top_x=top_x)[['filter_pi']]
    return(pi_return)
  }))
}

build_search=function(column,
                      search_terms){
  list(c(column,search_terms))
}

#Function to grab accessions
grab_accessions=function(cs_metadata,
                         search_terms,
                         search_terms2=NULL,
                         grab_all=FALSE,
                         search_within_search=NULL,
                         exclude_search=NULL,
                         just_search_within_search=FALSE,
                         sra_col='sra'){
  results=c()
  #Some checks
  if(grab_all&!is.null(search_within_search)){
    message('Must be either grab_all=TRUE OR search_within_search=TRUE, not both')
    return(NULL)
  }
  temp_meta=cs_metadata
  if(!is.null(search_terms2)){
    search_terms_all=c(search_terms,search_terms2)
  }else{
    search_terms_all=search_terms
  }
  for(i in 1:length(search_terms_all)){
    relevant_search=search_terms_all[i]
    all_search=unlist(relevant_search)
    search_col=all_search[1]
    search_temp=all_search[2:length(all_search)]
    temp_meta=temp_meta[temp_meta[[search_col]]%in%search_temp,]
  }
  sra_temp=unique(temp_meta[[sra_col]])
  if(grab_all){
    temp_meta=cs_metadata[cs_metadata[[sra_col]]%in%sra_temp,]
  }else if (!is.null(search_within_search)){
    cs_metadata_temp=cs_metadata[cs_metadata[[sra_col]]%in%sra_temp,]
    for(i in 1:length(search_terms)){
      relevant_search=search_terms_all[i]
      all_search=unlist(relevant_search)
      search_col=all_search[1]
      search_temp=all_search[2:length(all_search)]
      cs_metadata_temp=cs_metadata_temp[cs_metadata_temp[[search_col]]%in%search_temp,]
    }
    
    for(i in 1:length(search_within_search)){
      relevant_search=search_within_search[i]
      all_search=unlist(relevant_search)
      search_col=all_search[1]
      search_temp=all_search[2:length(all_search)]
      cs_metadata_temp=cs_metadata_temp[cs_metadata_temp[[search_col]]%in%search_temp,]
    }
    if(!just_search_within_search){
      temp_meta=rbind(temp_meta,cs_metadata_temp)
    }else if(just_search_within_search){
      temp_meta=cs_metadata_temp
    }
    
  }
  if(!is.null(exclude_search)){
    for(i in 1:length(exclude_search)){
      relevant_search=exclude_search[i]
      all_search=unlist(relevant_search)
      search_col=all_search[1]
      search_temp=all_search[2:length(all_search)]
      temp_meta=temp_meta[temp_meta[[search_col]]%!in%search_temp,]
    }
  }
  results[['study']]=unique(temp_meta[[sra_col]])
  results[['sample']]=unique(temp_meta$sample_ID)
  return(results)
}

#Function to label custom batch to remove with combat-seq
label_batch=function(deseq2_obj,
                     batch_accession_1){
  colData(deseq2_obj)[['batch']]=ifelse(colData(deseq2_obj)[['external_id']]%in%batch_accession_1,'batch_1','batch_2')
  return(deseq2_obj)
}

find_degs_multi=function(rse_obj,
                         group_col,
                         independentFiltering=TRUE,
                         pval_cutoff=0.05,
                         log2fc_cutoff=log2(1.5),
                         outgroup=NULL,
                         batch_col=NULL,
                         find_degs=TRUE,
                         sample_col='external_id',
                         only_outgroup=TRUE){
  
  results=c()
  all_groups=unique(rse_obj[[group_col]])
  full_merge=vector_to_df(all_groups)
  
  if(!is.null(outgroup)){
    full_merge=do.call('rbind',apply(full_merge,MARGIN=1,function(find_outgroup){
      if(find_outgroup[['x']]==outgroup|
         find_outgroup[['y']]==outgroup){
        both_group=c(as.character(find_outgroup[['x']]),
                     as.character(find_outgroup[['y']]))
        x=both_group[both_group!=outgroup]
        y=outgroup
        data.frame(cbind(x,y))
      }else{
        x=as.character(find_outgroup[['x']])
        y=as.character(find_outgroup[['y']])
        data.frame(cbind(x,y))
      }
    }))
    if(only_outgroup){
      full_merge=full_merge[full_merge[['y']]==outgroup,]
    }
  }
  
  all_degs=c()
  for(i in 1:nrow(full_merge)){
    temp_name_1=full_merge[['x']][i]
    temp_name_1=clean_name(temp_name_1)
    temp_name_2=full_merge[['y']][i]
    temp_name_2=clean_name(temp_name_2)
    message(temp_name_1,' vs ',temp_name_2)
    temp_degs=find_degs_new(rse_obj=rse_obj,
                            group_col=group_col,
                            group_1=full_merge[['x']][i],
                            group_2=full_merge[['y']][i],
                            independentFiltering=independentFiltering,
                            pval_cutoff=pval_cutoff,
                            log2fc_cutoff=log2fc_cutoff,
                            batch_col=batch_col,
                            find_degs=find_degs,
                            sample_col=sample_col)
    all_degs=rbind(all_degs,temp_degs$degs)
    
    temp_name=paste0(temp_name_1,'_vs_',temp_name_2,'_p')
    results[[temp_name]]=temp_degs$dds_deseq
  }
  results[['all_degs']]=all_degs
  return(results)
}


vector_to_df=function(vector_to_use){
  vector_df=merge(vector_to_use,vector_to_use)
  vector_df=vector_df[vector_df[['x']]!=vector_df[['y']],]
  vector_df[['accession']]=apply(vector_df,MARGIN=1,function(per_row){
    temp_sort=paste0(sort(c(per_row[['x']],per_row[['y']])),collapse='_')
    return(temp_sort)
  })
  vector_df=vector_df[!duplicated(vector_df[['accession']]),]
  vector_df=vector_df[,colnames(vector_df)!='accession']
  return(vector_df)
}

#Function to enrich DEGs
enrich_degs=function(DEGs,
                     by_col='dir_accession',
                     return_cols=c('accession','group_1','direction_1','group_1_dir',
                                   'group_2','direction_2','group_2_dir'),
                     gene_col='gene',
                     entrez_dictionary,
                     entrez_dictionary_ensembl_col='ensembl_gene_id',
                     entrez_dictionary_entrez_col='entrezgene_id',
                     entrez_dictionary_gene_symbol_col = 'external_gene_name',
                     simplify=FALSE,
                     pval_cutoff=0.05,
                     use_ensembl=TRUE){
  unique_accession=unique(DEGs[[by_col]])
  do.call('rbind',lapply(unique_accession,function(per_accession){
    message(per_accession)
    per_accession_title=gsub(per_accession,pattern='\n',replacement='_')
    per_accession_title=gsub(per_accession_title,pattern=' ',replacement='_')
    DEGs_relevant=DEGs[DEGs[[by_col]]==per_accession,]
    DEGs_relevant_sig=DEGs_relevant[DEGs_relevant[['sig']]=='y',]
    testing_p=enrich_genes(gene_list=unique(DEGs_relevant_sig[[gene_col]]),
                           background=unique(DEGs[[gene_col]]),
                           gene_dictionary=entrez_dictionary,
                           ensembl_symbol=entrez_dictionary_ensembl_col,
                           gene_entrez=entrez_dictionary_entrez_col,
                           gene_symbol = entrez_dictionary_gene_symbol_col,
                           graph_enrichment=FALSE,
                           simplify=simplify,
                           enrichment_pval=pval_cutoff,
                           use_ensembl=use_ensembl)
    all_enriched=testing_p$enrichment
    #Get relevant cols back
    DEGs_relevant_sig_temp=unique(DEGs_relevant_sig[,colnames(DEGs_relevant_sig)%in%c(return_cols,by_col)])
    DEGs_relevant_sig_temp=DEGs_relevant_sig_temp%>%dplyr::select(c(return_cols,by_col))
    
    all_enriched=cbind(all_enriched,DEGs_relevant_sig_temp)
    return(all_enriched)
  }))
}

enrich_genes=function(gene_list,
                      background,
                      gene_dictionary,
                      ensembl_symbol='ensembl_gene_id',
                      gene_symbol='external_gene_name',
                      gene_entrez='entrezgene_id',
                      enrichment_correction='BH',
                      enrichment_pval=0.05,
                      simplify=FALSE,
                      space_indent=60,
                      axis_indent=10,
                      graph_enrichment=TRUE,
                      # max_n=30,
                      simplify_cutoff=0.7,
                      use_ensembl=TRUE){
  results=c()
  #Get proper dictionary
  gene_dictionary=gene_dictionary[!is.na(gene_dictionary[[gene_entrez]]),]
  if(use_ensembl){
    #Get background entrez
    background_entrez=as.character(gene_dictionary[gene_dictionary[[ensembl_symbol]]%in%background,][[gene_entrez]])
    #get gene list entrez
    gene_list_entrez=as.character(gene_dictionary[gene_dictionary[[ensembl_symbol]]%in%gene_list,][[gene_entrez]])
  }else{
    #Get background entrez
    background_entrez=as.character(gene_dictionary[gene_dictionary[[gene_symbol]]%in%background,][[gene_entrez]])
    #get gene list entrez
    gene_list_entrez=as.character(gene_dictionary[gene_dictionary[[gene_symbol]]%in%gene_list,][[gene_entrez]])
  }
  
  gene_list_entrez=gene_list_entrez[gene_list_entrez%in%background_entrez]
  #Kegg
  kegg_temp=data.frame(enrichKEGG(gene = unique(gene_list_entrez),
                                  organism = 'human',
                                  pAdjustMethod = enrichment_correction,
                                  pvalueCutoff = enrichment_pval,
                                  universe = unique(background_entrez)))
  
  if(nrow(kegg_temp)>0){
    kegg_temp$enrichment='KEGG'
    kegg_temp=enrichment_gene_symbol(enrichment_object=kegg_temp,
                                     entrez_ID='geneID',
                                     entrez_dictionary=gene_dictionary,
                                     entrez_dictionary_entrez_col=gene_entrez,
                                     entrez_dictionary_gene_symbol_col=gene_symbol)
    kegg_temp$simplify=NA
    kegg_temp$simplify_cutoff=NA
    if(graph_enrichment){
      KEGG_dotplot=enrichment_dotplot(enrichment_table=kegg_temp,
                                      enrichment_type = 'KEGG',
                                      space_indent=space_indent,
                                      axis_indent=axis_indent)
      results[['kegg_plot']]=KEGG_dotplot
    }
  }
  
  GO_temp=enrichGO(gene = unique(gene_list_entrez),
                   keyType = 'ENTREZID',
                   OrgDb = 'org.Hs.eg.db',
                   ont='BP',
                   pvalueCutoff = enrichment_pval,
                   pAdjustMethod = enrichment_correction,
                   universe = unique(background_entrez),
                   readable = FALSE)
  
  if(simplify){
    message('Simplifying GO terms')
    GO_temp=clusterProfiler::simplify(GO_temp,
                                      cutoff=simplify_cutoff)
    
  }
  
  GO_temp=data.frame(GO_temp)
  
  if(nrow(GO_temp)>0){
    GO_temp$enrichment='GO'
    GO_temp=enrichment_gene_symbol(enrichment_object=GO_temp,
                                   entrez_ID='geneID',
                                   entrez_dictionary=gene_dictionary,
                                   entrez_dictionary_entrez_col=gene_entrez,
                                   entrez_dictionary_gene_symbol_col=gene_symbol)
    GO_temp$simplify=simplify
    if(simplify){
      GO_temp$simplify_cutoff=simplify_cutoff
    }else{
      GO_temp$simplify_cutoff=NA
    }
    if(graph_enrichment){
      GO_dotplot=enrichment_dotplot(GO_temp,enrichment_type = 'GO',
                                    space_indent=space_indent,
                                    axis_indent=axis_indent)
      results[['go_plot']]=GO_dotplot
    }
    GO_temp$category=NA
    GO_temp$subcategory=NA
    GO_temp <- GO_temp[, colnames(kegg_temp)]
  }
  
  all_enrichment=rbind(kegg_temp,GO_temp)
  if(nrow(all_enrichment)==0){
    return(NULL)
  }else{
    #Stop microsoft being dumb
    all_enrichment[['GeneRatio']]=paste0(' ',all_enrichment[['GeneRatio']])
    all_enrichment[['BgRatio']]=paste0(' ',all_enrichment[['BgRatio']])
    all_enrichment$pval_cutoff=enrichment_pval
    results[['enrichment']]=all_enrichment 
  }
  
  return(results)
}

#Function to merge in gene symbols to entrez enrichment results
enrichment_gene_symbol=function(enrichment_object,
                                entrez_ID,
                                entrez_dictionary,
                                entrez_dictionary_entrez_col,
                                entrez_dictionary_gene_symbol_col){
  enrichment_temp=enrichment_object
  gene_symbol=apply(enrichment_temp,MARGIN = 1,function(bind_genes){
    temp_genes=data.frame(gene=unlist(strsplit(bind_genes[[entrez_ID]],split = '/')))
    gene_entrez=merge(temp_genes,
                      entrez_dictionary,
                      all.x=TRUE,
                      by.x='gene',
                      by.y=entrez_dictionary_entrez_col)
    
    gene_entrez$gene=factor(gene_entrez$gene,levels=temp_genes$gene)
    gene_entrez=gene_entrez%>%dplyr::arrange(gene)
    genes_return=paste0(unique(gene_entrez[[entrez_dictionary_gene_symbol_col]]),collapse='/')
    return(genes_return)
  })
  gene_row=cbind(enrichment_temp,gene_symbol)
  return(gene_row)
}

#Function to graph DEGs by various time points
time_graph_deg_db=function(degs,
                           group_1_deg_col='deg_group_1',
                           group_1_dir_col='group_1_dir',
                           group_2_deg_col='deg_group_2',
                           pval_col='adj',
                           odds_col='odds',
                           n_col='actual',
                           order=NULL,
                           facet_col='db',
                           exclude_db_facet=NULL,
                           nudge_y=35,
                           label_ns=FALSE,
                           xlab='Days',
                           ylab='Number of overlapping DEGs',
                           cell_type=NULL,
                           ggtitle=NULL,
                           facet_order=NULL,
                           tilt_x=FALSE){
  if(!is.null(cell_type)){
    ggtitle_temp=paste(ggtitle, cell_type)
  }else{
    ggtitle_temp=ggtitle
  }
  degs_temp=degs
  if(length(exclude_db_facet)>0){
    degs_temp=degs_temp[degs_temp[[facet_col]]%!in%exclude_db_facet,]
  }
  degs_temp=degs_temp[,colnames(degs_temp)!='intersect']
  
  if(!is.null(order)){
    degs_temp[[group_1_deg_col]]=factor(degs_temp[[group_1_deg_col]],
                                        levels=order)
  }
  degs_temp$log2odds=log2(degs_temp[[odds_col]])
  degs_temp=degs_temp%>%rstatix::add_significance(pval_col)
  if(!label_ns){
    degs_temp[['adj.signif']][degs_temp[['adj.signif']]=='ns']=''
  }
  degs_temp[['x']]=degs_temp[[group_1_deg_col]]
  degs_temp[['y']]=degs_temp[[n_col]]
  degs_temp[['facet_temp']]=degs_temp[[facet_col]]
  degs_temp[['dir_temp']]=degs_temp[[group_1_dir_col]]
  
  graph_max=max(ceiling(degs_temp$log2odds/2)[!is.infinite(ceiling(degs_temp$log2odds))])*2
  graph_min=min(floor(degs_temp$log2odds)[!is.infinite(ceiling(degs_temp$log2odds))])
  if(graph_min>0){
    graph_min=0
  }
  
  li <- c(graph_min, graph_max)
  la <- c(seq(graph_min,graph_max,2))
  br <- c(seq(graph_min,graph_max,2))
  if(!is.null(facet_order)){
    degs_temp[['facet_temp']]=factor(degs_temp[['facet_temp']],
                                     levels=facet_order)
  }
  graph = degs_temp%>%ggplot(aes(x=x,y=y,fill=log2odds))+
    geom_point(shape = 21,size=4)+
    facet_grid(facet_temp~dir_temp)+
    geom_text(aes(label=adj.signif),nudge_y = nudge_y)+
    scale_fill_gradient2(expression('log'[2]*'(Odds)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks = TRUE, 
                                                ticks.colour = 'black',title.position = 'top',title.hjust=0.5),
                         breaks=br,
                         labels=la,
                         limits=li)+
    xlab(xlab)+
    ylab(ylab)+
    ggtitle(ggtitle_temp)
  
  if (tilt_x) {
    graph = graph + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(graph)
}

compare_logfc_outgroup=function(degs,
                                fc_col='log2FoldChange',
                                degs_gene_col='gene',
                                degs_dir_ingroup='direction_1',
                                outgroup='none',
                                outgroup_col='group_2',
                                ingroup_col='group_1',
                                database,
                                database_gene_col='ensembl_gene_id',
                                database_group_col,
                                xlab=NULL,
                                ylab='log2FC vs. Proliferating',
                                facet_order=NULL,
                                plot_title=TRUE){
  database_unique=unique(database[[database_group_col]])
  degs_temp=degs[degs[[outgroup_col]]==outgroup,]
  wilcox_return=c()
  kw_return=c()
  for(i in 1:length(database_unique)){
    per_db=database_unique[i]
    db_temp=database[database[[database_group_col]]==per_db,][[database_gene_col]]
    db_fc=degs_temp[degs_temp[[degs_gene_col]]%in%db_temp,]
    # unique_dir=unique(db_fc[[degs_dir_ingroup]])
    # lapply(unique_dir,function(per_dir){
    # db_fc_dir=db_fc[db_fc[[degs_dir_ingroup]]==per_dir,]
    db_fc[['logfc_col']]=db_fc[[fc_col]]
    db_fc[['comparison_col']]=db_fc[[ingroup_col]]
    temp_merge=merge(unique(db_fc[['comparison_col']]),unique(db_fc[['comparison_col']]))
    temp_merge=temp_merge[temp_merge[['x']]!=temp_merge[['y']],]
    temp_merge$accession=apply(temp_merge,MARGIN=1,function(accession_row){
      paste0(sort(c(accession_row[['x']],accession_row[['y']])),collapse='_')
    })
    temp_merge=temp_merge[!duplicated(temp_merge$accession),]
    wilcox_temp=do.call('rbind',apply(temp_merge,MARGIN=1,function(per_row){
      db_fc_dir_group=db_fc[db_fc[[ingroup_col]]%in%c(per_row[['x']],per_row[['y']]),]
      db_fc_dir_group[['comparison_col']]=droplevels(db_fc_dir_group[['comparison_col']])
      temp_return=db_fc_dir_group%>%rstatix::wilcox_test(formula=logfc_col~comparison_col)
      temp_return$db=per_db
      return(temp_return)
    }))
    wilcox_return=rbind(wilcox_return,wilcox_temp)
    
    kw_temp=db_fc%>%rstatix::kruskal_test(formula=logfc_col~comparison_col)
    kw_temp$db=per_db
    kw_return=rbind(kw_return,kw_temp)
    # })
  }
  #correct
  wilcox_return$adj=p.adjust(wilcox_return$p,method='BH')
  kw_return$adj=p.adjust(kw_return$p,method='BH')
  
  wilcox_return=wilcox_return%>%rstatix::add_significance('adj')
  kw_return=kw_return%>%rstatix::add_significance('adj')
  
  plot_temp=lapply(database_unique,function(per_db){
    db_temp=database[database[[database_group_col]]==per_db,][[database_gene_col]]
    db_fc=degs_temp[degs_temp[[degs_gene_col]]%in%db_temp,]
    db_fc[['logfc_col']]=db_fc[[fc_col]]
    db_fc[['comparison_col']]=db_fc[[ingroup_col]]
    temp_wilcox=wilcox_return[wilcox_return[['db']]==per_db,]
    temp_kw=kw_return[kw_return[['db']]==per_db,]
    label_kw=paste0('KW: ',temp_kw[['adj.signif']])
    temp_max=max(db_fc[['logfc_col']])+0.5
    # temp_max_kw=temp_max+nrow(temp_wilcox)
    temp_wilcox$y.position=ceiling(temp_max)
    if(!is.null(facet_order)){
      db_fc[['comparison_col']]=factor(db_fc[['comparison_col']],
                                       levels=facet_order)
    }
    
    scale_temp_min=floor(min(db_fc[['logfc_col']]))
    # if(scale_temp_min>0){
    #   scale_temp_min=0
    # }
    # scale_temp_max=floor(min(db_fc[['logfc_col']]))
    # if(scale_temp_max<0){
    #   scale_temp_max=0
    # }
    db_fc$facet=label_kw
    p_temp=db_fc%>%ggplot(aes(y=logfc_col,x=comparison_col))+
      geom_boxplot(aes(fill=comparison_col),show.legend=FALSE)+
      stat_pvalue_manual(temp_wilcox,label='adj.signif',
                         step.increase=0.1,tip.length=0.01,size = 9)+
      # ylim(c(scale_temp_min,scale_temp_max))+
      # scale_y_continuous(breaks=seq(scale_temp_min,scale_temp_max,by=1))+
      # geom_text(aes(x=1,y=temp_max_kw,label=label_kw),inherit.aes=FALSE,size=9)+
      xlab(xlab)+
      ylab(ylab)+
      theme(axis.text = element_text(size=20),
            axis.title=element_text(size=20),
            strip.text=element_text(size=20))+
      facet_wrap(~facet)
    if(plot_title){
      p_temp=p_temp+
        ggtitle(per_db)+
        theme(plot.title = element_text(size=20))
    }
    return(p_temp)
  })
  names(plot_temp)=paste0(clean_name(database_unique),'_p')
  plot_temp$wilcox=wilcox_return
  plot_temp$kw=kw_return
  return(plot_temp)
}

#Function to replace underscore with spaces while maintaining factor order
replace_underscore_with_space=function(df, old_col, new_col) {
  # Replace underscores with spaces in the old column
  new_col_values = gsub("_", " ", df[[old_col]])
  
  # Get the levels of the old column
  old_levels = levels(df[[old_col]])
  
  # Replace underscores with spaces in the old levels
  new_levels = gsub("_", " ", old_levels)
  
  # Convert the new column to a factor with the same order as the old column
  df[[new_col]] = factor(new_col_values, levels = new_levels)
  
  return(df)
}

#Run whole analysis for time point data
time_analysis=function(cell_type,
                       recount_pheno,
                       search_term_study,
                       ensembl_dictionary,
                       enrich_degs_test=FALSE,
                       entrez_dictionary=NULL,
                       db,
                       db_col,
                       db_gene_col='gene_name',
                       results_dir,
                       exclude_db=NULL,
                       grab_all=FALSE,
                       exclude_search=NULL,
                       ensembl_gene_col='ensembl_gene_id',
                       treatment_time_col='time_after_treatment',
                       time_levels=c('none', '4_days', '10_days', '20_days'),
                       pval_cutoff=0.05,
                       log2fc_cutoff=log2(1.5),
                       independentFiltering=TRUE,
                       pca_plot_group=c('cell_state','time_after_treatment'),
                       sep='/',
                       fix_batch=TRUE,
                       batch_1=NULL,
                       run_simulations=FALSE,
                       simulation_n=10000,
                       vst_n=500){
  if(fix_batch&is.null(batch_1)){
    message('Must provide accessions for custom batch')
    return(NULL)
  }
  if(enrich_degs_test&is.null(entrez_dictionary)){
    message('Must provide entrez dictionary')
    return(NULL)
  }
  study_temp=search_term_study[[1]][2]
  temp_dir=paste0(results_dir,study_temp,sep)
  if(!dir.exists(temp_dir)){
    dir.create(temp_dir)
  }
  temp_dir_cell=paste0(temp_dir,cell_type,sep)
  if(!dir.exists(temp_dir_cell)){
    dir.create(temp_dir_cell)
  }
  ##
  time_cell=build_search(column='cell_type',
                         search_terms=c(cell_type))
  
  accession_time=grab_accessions(cs_metadata=recount_pheno,
                                 search_terms=search_term_study,
                                 search_terms2=time_cell,
                                 grab_all=grab_all,
                                 search_within_search=NULL,
                                 exclude_search=exclude_search)
  
  message('Finding samples')
  accession_time_study=download_studies(studies=study_temp,
                                        sra_organism='human')
  
  accession_time_study_processed=process_rse(rse=accession_time_study,
                                             meta_add=recount_pheno,
                                             meta_id_col = 'sample_ID',
                                             filter_pc=TRUE,
                                             condition_col = 'cell_state',
                                             ensembl_dictionary=ensembl_dictionary,
                                             filter_duplicates=TRUE,
                                             sample_filter=accession_time[['sample']],
                                             low_expression_filter='per',
                                             scale_counts=FALSE,
                                             remove_cq = TRUE)
  
  time_coldata=data.frame(colData(accession_time_study_processed))
  
  save_csv(data = time_coldata,file_name = 'sample_pheno',
           path=temp_dir_cell)
  # nrow(time_coldata)
  # 
  # #Number of samples by CS type
  # table(time_coldata$cell_state)
  # table(time_coldata$time_after_treatment)
  # 
  # #Number of unique studies
  # length(unique(time_coldata$study))
  
  accession_time_study_processed[[treatment_time_col]][is.na(accession_time_study_processed[[treatment_time_col]])]=time_levels[1]
  accession_time_study_processed[[treatment_time_col]]=gsub(accession_time_study_processed[[treatment_time_col]],pattern = ' ',replacement = '_')
  accession_time_study_processed[[treatment_time_col]]=factor(accession_time_study_processed[[treatment_time_col]],
                                                              levels=time_levels)
  ##no correction
  if(fix_batch){
    accession_time_study_processed=label_batch(deseq2_obj=accession_time_study_processed,
                                               batch_accession_1=batch_1)
    pca_plot_group_temp=unique(c(pca_plot_group,'batch'))
  }else{
    pca_plot_group_temp=pca_plot_group
  }
  
  time_no_correction=deseq_studies(study_obj=accession_time_study_processed,
                                   protect_col=treatment_time_col,
                                   fix_col=NULL)
  
  saveRDS(time_no_correction,file = paste0(temp_dir_cell,
                                           'deseq.rds'))
  
  time_no_correction_keep_vst=normalise_counts(
    dds_obj_use=time_no_correction[['deseq_obj']],
    batch_col=NULL,
    protect_col=NULL,
    method=NULL,
    normalisation='vst',
    vst_n=vst_n,
    blind=TRUE)
  
  if(fix_batch){
    pca_temp=c(treatment_time_col,'batch')
  }else{
    pca_temp=c(treatment_time_col)
  }
  
  time_no_correction_keep_vst_pca=pca_rds(deseq_obj = time_no_correction_keep_vst,
                                          ntop=vst_n,
                                          pca_plot=pca_temp)
  
  save_p(time_no_correction_keep_vst_pca,
         file_name='pca_no_batch_correction',save_dir = temp_dir_cell,p_width=12,p_height=10)
  ##correction
  if(fix_batch){
    message('Batch correcting')
    time_correction=deseq_studies(study_obj=accession_time_study_processed,
                                  protect_col=treatment_time_col,
                                  fix_col='batch')
    
    saveRDS(time_correction,file = paste0(temp_dir_cell,
                                             'deseq.rds'))
    
    time_correction_keep_vst=normalise_counts(
      dds_obj_use=time_correction[['deseq_obj']],
      batch_col='batch',
      protect_col=treatment_time_col,
      method='wgcna',
      normalisation='vst',
      vst_n=vst_n,
      blind=TRUE)
    
    time_correction_keep_vst_pca=pca_rds(
      time_correction_keep_vst,
      ntop=vst_n,
      pca_plot=pca_temp)
    
    dds_use=time_correction
    
    save_p(time_correction_keep_vst_pca,
           file_name = 'pca_batch_correction',save_dir=temp_dir_cell,p_height = 9,
           p_width=9)
    batch_col='batch'
  }else{
    dds_use=time_no_correction
    batch_col=NULL
  }
  
  accession_time_study_DEGs=find_degs_multi(rse_obj=dds_use[['deseq_obj']],
                                            group_col=treatment_time_col,
                                            independentFiltering=independentFiltering,
                                            pval_cutoff=pval_cutoff,
                                            log2fc_cutoff=log2fc_cutoff,
                                            outgroup='none',
                                            batch_col=batch_col,
                                            find_degs=TRUE,
                                            only_outgroup = TRUE,
                                            sample_col='external_id')
  
  time_degs=accession_time_study_DEGs[['all_degs']]
  
  if(fix_batch){
    time_degs$batch_corrected='y'
  }else{
    time_degs$batch_corrected='n'
  }
  
  save_csv(data = time_degs,
           file_name ='degs',
           path=temp_dir_cell,
           merge_ensembl_symbol=TRUE,gene_col='gene')
  
  #####
  #heatmap
  message('building heatmaps and volcano plots')
  time_degs$accession_find=paste0(time_degs$group_1,'_vs_',
                                  time_degs$group_2)
  # time_degs=merge(time_degs,human_pc,by.x='gene',by.y='ensembl_gene_id')
  temp_deseq=names(accession_time_study_DEGs[names(accession_time_study_DEGs)!='all_degs'])
  lapply(temp_deseq,function(per_deseq){
    find_accession=gsub(per_deseq,pattern='_p',replacement='')
    time_degs_volcano=time_degs[time_degs$accession_find==find_accession,]
    plot_title_temp=gsub(find_accession,pattern='_',replacement=' ')
    plot_title_temp=gsub(plot_title_temp,pattern='none',replacement='proliferating')
    plot_title_temp=paste0(cell_type,' ',plot_title_temp)
    volcano_plot_temp=volcano_function(degs=time_degs_volcano,
                                       gene_col = 'gene',
                                       # custom_labs='group_1',
                                       title = plot_title_temp)
    save_p(volcano_plot_temp,file_name=paste0(find_accession,'_volcano'),save_dir=temp_dir_cell)
    temp_degs_heatmap=time_degs_volcano[time_degs_volcano$sig=='y',]
    group_1_heatmap=unique(temp_degs_heatmap$group_1)
    group_2_heatmap=unique(temp_degs_heatmap$group_2)
    group_col_heatmap=c('time_after_treatment',batch_col)
    temp_heat=build_heatmap(degs = temp_degs_heatmap,
                            dds_obj = accession_time_study_DEGs[[per_deseq]]$deseq_obj,
                            group_col = 'time_after_treatment',
                            groups = c(group_1_heatmap,group_2_heatmap),
                            protect_col = 'time_after_treatment',
                            batch_col = batch_col,
                            normalisation = 'vst_full',
                            group_col_heatmap=group_col_heatmap,
                            gene_col = 'gene',
                            rownames_to_symbol = FALSE)
    temp_pi=pi_score_degs(sig_degs =temp_degs_heatmap,
                          top_x = 50)
    save_csv(data = temp_pi$filter_pi,
             file_name =paste0(find_accession,'_pi'),path = temp_dir_cell)
    
    save_pheatmap(pheatmap_to_save=temp_heat,
                  file_name=paste0(find_accession,'_heatmap'),
                  save_dir=temp_dir_cell,
                  p_width=3000,
                  p_height=6000)
  })
  message('making whole heatmap')
  time_degs_sig=time_degs[time_degs$sig=='y',]
  temp_heat=build_heatmap(degs = time_degs_sig,
                            dds_obj = dds_use[['deseq_obj']],
                            group_col = 'time_after_treatment',
                            groups = unique(c(time_degs_sig$group_1,
                                       time_degs_sig$group_2)),
                            protect_col = 'time_after_treatment',
                            batch_col = batch_col,
                            normalisation = 'vst_full',
                          group_col_heatmap=c('time_after_treatment',batch_col),
                            gene_col = 'gene',
                            rownames_to_symbol = FALSE,
                          top_x = 25)
  conditions_pi=unique(time_degs_sig[['dir_accession']])
  relevant_genes_heatmap=do.call('rbind',lapply(conditions_pi,function(per_condition){
    per_condition_temp=time_degs_sig[time_degs_sig$dir_accession==per_condition,]
    sig_degs_use_heatmap=pi_score_degs(sig_degs = per_condition_temp,
                                       top_x=25)[['filter_pi']]
    return(sig_degs_use_heatmap)
  }))
  save_pheatmap(pheatmap_to_save=temp_heat,
                file_name=paste0(cell_type,'_heat'),
                save_dir=temp_dir_cell,
                p_width=3000,
                p_height=6000)
  save_csv(relevant_genes_heatmap,
           file_name = paste0(cell_type,'_heat_pi'),
           path = temp_dir_cell)
  #####
  
  volcano_plot=volcano_function(degs=time_degs,
                                custom_labs='group_1',
                                gene_col='gene',
                                title = cell_type)
  
  save_p(volcano_plot,file_name='volcano_all',save_dir=temp_dir_cell)
  
  time_degs_sig=time_degs[time_degs$sig=='y',]
  time_degs_sig_pheat=time_degs_sig
  time_degs_sig_pheat=pi_score_filter_by_group(sig_degs=time_degs_sig_pheat,
                                               group_col='group_1',
                                               top_x=50)
  dds_obj_use=dds_use[['deseq_obj']]
  if(!fix_batch){
    normalised_corrected_counts=normalise_counts(
      dds_obj_use,
      batch_col=NULL,
      protect_col=NULL,
      method='wgcna',
      normalisation='vst_full')
    pheat_all_p=pheat_counts(dds_obj_use=normalised_corrected_counts,
                             group_col=treatment_time_col,
                             sig_degs=time_degs_sig_pheat[['gene']],
                             scale='row',
                             rownames_to_symbol = FALSE,
                             gene_col = 'gene')
  }else if(fix_batch){
    normalised_corrected_counts=normalise_counts(
      dds_obj_use,
      batch_col='batch',
      protect_col=treatment_time_col,
      method='wgcna',
      normalisation='vst_full')
    pheat_all_p=pheat_counts(dds_obj_use=normalised_corrected_counts,
                             group_col=c(treatment_time_col,'batch'),
                             sig_degs=time_degs_sig_pheat[['gene']],
                             scale='row',
                             rownames_to_symbol = FALSE,
                             gene_col = 'gene')
  }
  
  save_pheatmap(pheatmap_to_save=pheat_all_p,
                file_name=paste0(cell_type,'_all_condition'),
                save_dir=temp_dir_cell,
                p_width=2000,
                p_height=4000)
  
  if(enrich_degs_test){
    ##Enrich DEGs
    time_degs_enriched=enrich_degs(DEGs=time_degs,
                                   by_col='dir_accession',
                                   gene_col='gene',
                                   entrez_dictionary=entrez_dictionary,
                                   use_ensembl = FALSE)
    
    save_csv(data =time_degs_enriched,
             file_name = paste0('time_deg_enrichment'),
             path =temp_dir_cell)
  }
  
  time_degs_sig_count=data.frame(table(time_degs_sig$dir_accession))
  time_degs_sig_count$Var1=clean_name(time_degs_sig_count$Var1)
  colnames(time_degs_sig_count)[1]='DEG'
  
  save_csv(data=time_degs_sig_count,
           file_name ='deg_count',path=temp_dir_cell)
  
  message('Overlapping DEGs with each other')
  all_deg_overlaps_time=overlap_within_df(dataframe = time_degs_sig,
                                          gene_col = 'gene',
                                          group_col = 'group_1',
                                          other_cols = 'group_1_dir',
                                          background = time_degs$gene)
  
  save_csv(data = all_deg_overlaps_time,
           file_name = 'deg_self_overlaps',path = temp_dir_cell)
  
  get_pval=all_deg_overlaps_time[,colnames(all_deg_overlaps_time)%in%
                                   c('accession','adj')]
  
  self_overlaps_to_plot=overlap_within_df(dataframe = time_degs_sig,
                                          gene_col = 'gene',
                                          group_col = 'group_1',
                                          other_cols = 'group_1_dir',
                                          background = unique(time_degs$gene),
                                          remove_self = FALSE
  )
  
  self_overlaps_to_plot=self_overlaps_to_plot[,colnames(self_overlaps_to_plot)!='adj']
  self_overlaps_to_plot=merge(self_overlaps_to_plot,get_pval)
  
  self_overlaps_to_plot$outgroup_dir=paste0(self_overlaps_to_plot$group_1_2,' ',self_overlaps_to_plot$direction_1_2)
  self_overlaps_to_plot$outgroup_dir=gsub(self_overlaps_to_plot$outgroup_dir,pattern='_',replacement=' ')
  self_overlaps_to_plot$group_1_dir_1=
    extract_up_down(self_overlaps_to_plot$group_1_dir_1)
  self_overlaps_to_plot$group_1_dir_1=tools::toTitleCase(self_overlaps_to_plot$group_1_dir_1)
  self_overlaps_to_plot$group_1_dir_1=paste0(self_overlaps_to_plot$group_1_dir_1,'\nin Arrest')
  
  self_overlaps_to_plot[['group_1_1']]=factor(self_overlaps_to_plot[['group_1_1']],
                                              levels=time_levels[2:length(time_levels)])
  self_overlaps_to_plot=factor_column_and_modify(df = self_overlaps_to_plot,
                           column = 'group_1_dir_2',
                           old_list = time_levels[2:length(time_levels)],
                           keyword = '\nin Arrest')
  
  self_overlaps_to_plot=factor_column_and_modify(df = self_overlaps_to_plot,
                                                 column = 'group_1_1',
                                                 old_list = levels(self_overlaps_to_plot$group_1_1),
                                                 keyword = NULL)
  
  self_overlap_p=plot_self_overlaps(self_overlaps_recount_to_plot = self_overlaps_to_plot,
                                    odds_column = 'odds',
                                    x = 'group_1_dir_1',
                                    y = 'group_1_dir_2',
                                    facet_col = 'group_1_1',
                                    xlab = 'Temporal DEG Facet vs Proliferating',
                                    ylab = 'Temporal DEG Y-Axis vs Proliferating')
  save_p(self_overlap_p,file_name = 'deg_self_overlaps',save_dir=temp_dir_cell)
  
  # time_degs_sig_temp=time_degs_sig[time_degs_sig$group_2!='none',]
  message('Overlapping DEGs with database')
  if(!is.null(exclude_db)){
    db=db[db[[db_col]]!=exclude_db,]
  }
  ##all overlaps
  # db_deg_overlaps_time=overlap_function(df_1 = time_degs_sig,
  #                                       gene_col_1 = 'gene',
  #                                       group_col_1 = 'group_1_dir',
  #                                       df_2 = db,
  #                                       group_col_2 = db_col,
  #                                       gene_col_2 = db_gene_col,
  #                                       background = time_degs$gene,
  #                                       carry_col_1 = c('group_1','direction_1'))
  # 
  # save_csv(data=db_deg_overlaps_time,file_name = 'deg_db_overlaps_all',
  #          path=temp_dir_cell)
  # db_deg_overlaps_time=factor_column_and_modify(df = db_deg_overlaps_time,
  #                                               column = 'group_1',
  #                                               old_list = time_levels[2:length(time_levels)],
  #                                               keyword = NULL)
  # 
  # db_deg_overlaps_time$direction_1[db_deg_overlaps_time$direction_1=='up']=
  #   'Up in\nArrest'
  # db_deg_overlaps_time$direction_1[db_deg_overlaps_time$direction_1=='down']=
  #   'Down in\nArrest'
  # 
  # overlap_plot=create_overlap_plot(deg_db_overlap = db_deg_overlaps_time,
  #                                  odds_column = 'odds',
  #                                  pval_col = 'adj',
  #                                  facet_col = 'group_1',
  #                                  x = 'direction_1',
  #                                  y = 'Senescence',
  #                                  xlab = paste0(cell_type,' DEGs'),
  #                                  ylab = 'CellAge')
  # save_p(overlap_plot,file_name = 'deg_db_overlaps_all',save_dir =temp_dir_cell)
  # 
  if(run_simulations){
    message('Running simulations')
    time_simulation_none_outgroup=simulate_overlaps(relevant_degs=time_degs,
                                                    facet_col='accession',
                                                    ingroup_col='group_1_dir',
                                                    outgroup_col='group_2_dir',
                                                    outgroup_relevant=c('none_up',
                                                                        'none_down'),
                                                    gene_col='gene',
                                                    seed=1,
                                                    simulation_n=simulation_n)
    
    save_csv(time_simulation_none_outgroup,
             'simulation_overlap.csv',temp_dir_cell)
    
    time_simulation_none_outgroup$outgroup[time_simulation_none_outgroup$outgroup=='none_up']='Down vs Proliferating'
    time_simulation_none_outgroup$outgroup[time_simulation_none_outgroup$outgroup=='none_down']='Up vs Proliferating'
    sim_cumsum=plot_top_simulation_cumsum(time_simulation_none_outgroup,
                                          plot_cutoff = TRUE,scale = 3)
    save_csv(sim_cumsum$cumsum,file_name = 'reverse_cumsum_recount',
             path = temp_dir_cell)
    
    save_p(sim_cumsum$p,file_name = 'reverse_cumsum_recount',
           save_dir  = temp_dir_cell)
  }
  
  ##non-zero overlaps
  # db_deg_overlaps_time_no_zero=deg_db_fisher(
  #   degs=time_degs_sig_temp,
  #   degs_facet='accession',
  #   degs_groups='dir_accession',
  #   db=db,
  #   db_col='db')
  # 
  # save_ggplot(db_deg_overlaps_time_no_zero$p,obj_name = 'deg_db_overlaps_no_proliferative',
  #             dir=temp_dir_cell)
  # p_over_time=time_graph_deg_db(
  #   degs=db_deg_overlaps_time,
  #   group_1_deg_col='group_1',
  #   group_1_dir_col='direction_1',
  #   group_2_deg_col='Senescence',
  #   pval_col='adj',
  #   odds_col='odds',
  #   n_col='actual',
  #   order= NULL,
  #   facet_col='Senescence',
  #   exclude_db_facet=exclude_db,
  #   nudge_y=45,
  #   label_ns=FALSE,
  #   xlab='Days',
  #   ylab='Number of overlapping DEGs',
  #   cell_type=cell_type,
  #   ggtitle='Differentially expressed genes compared to proliferating controls')
  # save_p(plot = p_over_time,file_name = 'deg_db_overlaps_vs_proliferative',
  #        save_dir = temp_dir_cell)
  
  #Compare expression amounts
  message('Comparing logFC between groups')
  time_degs=factor_column_and_modify(df = time_degs,
                                     column = 'group_1',
                                     old_list = time_levels[2:length(time_levels)],
                                     keyword = NULL)
  # testing_fc=compare_logfc_outgroup(ingroup_col = 'group_1',
  #                                   outgroup_col = 'group_2',
  #                                   degs_gene_col = 'gene',
  #                                   degs=time_degs,
  #                                   xlab='Day',
  #                                   facet_order=levels(time_degs$group_1),
  #                                   database=db,
  #                                   database_gene_col = db_gene_col,
  #                                   database_group_col=db_col)
  # testing_fc_p=testing_fc[1:(length(testing_fc)-2)]
  # for(i in 1:length(testing_fc_p)){
  #   temp_name=names(testing_fc_p[i])
  #   save_p(plot=testing_fc_p[[temp_name]],
  #          file_name = paste0(temp_name,'_log2FC'),save_dir = temp_dir_cell)
  # }
  # save_csv(data=testing_fc$wilcox,file_name='wilcox_fc_db',path=temp_dir_cell)
  # save_csv(data=testing_fc$kw,file_name='kw_fc_db',path=temp_dir_cell)
  # 
  #Compare expression amounts of triple intersect
  time_degs_sig_up=time_degs_sig[time_degs_sig$direction_1=='up',]
  time_degs_sig_up=data.frame(table(time_degs_sig_up[['gene']]))
  time_degs_sig_up=time_degs_sig_up[time_degs_sig_up[['Freq']]==3,]
  time_degs_sig_down=time_degs_sig[time_degs_sig$direction_1=='down',]
  time_degs_sig_down=data.frame(table(time_degs_sig_down[['gene']]))
  time_degs_sig_down=time_degs_sig_down[time_degs_sig_down[['Freq']]==3,]
  time_degs_sig_graph=unique(c(as.character(time_degs_sig_up[['Var1']]),
                               as.character(time_degs_sig_down[['Var1']])))
  time_degs_triple=time_degs_sig[time_degs_sig[['gene']]%in%time_degs_sig_graph,]
  
  time_degs$label=paste0(time_degs$direction_1,' vs Proliferating')
  testing_fc_triple=compare_logfc_outgroup(degs=time_degs,
                                           xlab='Day',
                                           facet_order=NULL,
                                           database=time_degs_triple,
                                           database_group_col='direction_1',
                                           database_gene_col='gene')
  graphs_temp=list(testing_fc_triple[['up_p']],testing_fc_triple[['down_p']])
  plot_both=cowplot::plot_grid(plotlist=graphs_temp)
  save_p(plot = plot_both,file_name = 'triple_intersect_range',
         save_dir = temp_dir_cell,p_width = 10)
  save_csv(data = testing_fc_triple$wilcox,file_name='wilcox_fc_triple',
           path=temp_dir_cell)
  save_csv(data = testing_fc_triple$kw,file_name='kw_fc_triple',
           path=temp_dir_cell)
}

#Function to compile time degs
compile_time_files=function(
    cell_types=c('fibroblast','melanocyte','keratinocyte'),
    dir_use,
    sep_slash='/',
    deg_file
){
  all_degs=do.call('rbind',lapply(cell_types,function(read_degs){
    dir_use_temp=paste0(dir_use,read_degs,sep_slash)
    temp_file=read.csv(paste0(dir_use_temp,deg_file))
    temp_file$cell_type=read_degs
    return(temp_file)
  }))
}

factor_list_based_on_other_list=function(old_list, new_list) {
  
  # Initialize a list to store the new levels for the new_list
  new_levels = list()
  
  # For each level in old_list
  for(i in seq_along(old_list)) {
    # Find the corresponding level in new_list that contains the level in old_list
    matching_new_levels = grep(old_list[i], new_list, value=TRUE)
    
    # If matching levels are found
    if(length(matching_new_levels) > 0) {
      # Add the matching levels to the new_levels
      new_levels[[i]] = matching_new_levels
    }
  }
  
  # Unlist the new_levels to create a character vector
  new_levels = unlist(new_levels)
  
  # Remove duplicates from new_levels
  new_levels = unique(new_levels)
  
  # Create a new factor from new_list, with levels ordered by new_levels
  new_factor = factor(new_list, levels = new_levels)
  
  # If some levels didn't match, make them NA
  new_factor[!(new_list %in% new_levels)] <- NA
  
  return(new_factor)
}

factor_column_and_modify=function(df, column, old_list,
                                  keyword = "\nin Arrest",
                                  caps=TRUE,
                                  remove_underscores=TRUE) {
  df[[column]]=as.character(df[[column]])
  # Convert underscores to spaces
  if(remove_underscores){
    df[[column]] = gsub("_", " ", df[[column]])
    # Convert the old list underscore to spaces
    old_list = gsub("_", " ", old_list)
  }
  
  # Convert to title case
  if(caps){
    df[[column]] = tools::toTitleCase(df[[column]])
    # Convert the old list to title case
    old_list = tools::toTitleCase(old_list)
  }
  
  # Add the keyword
  if(!is.null(keyword)){
    df[[column]] = paste(df[[column]], keyword)
  }
  
  # Now factor the column based on the old list
  df[[column]] = factor_list_based_on_other_list(old_list, df[[column]])
  
  return(df)
}

#Function to find common DEGs across groups
find_common_genes = function(df, group_col, gene_col) {
  
  # Split the dataframe into subsets by group
  group_lists = split(df, df[[group_col]])
  
  # Extract the genes from each group
  gene_lists = lapply(group_lists, function(group_df) {
    unique(group_df[[gene_col]])
  })
  
  # Find the intersection of all gene sets
  common_genes = Reduce(intersect, gene_lists)
  
  return(common_genes)
}

#Function to make dotplot of enrichment when markers aren't divided into groups
enrichment_dotplot=function(enrichment_table,
                            enrichment_type=NULL,
                            log2pval_filter=30,
                            pval_scale=5,
                            text_size=15,
                            legend_position='right',
                            space_indent=60,
                            axis_indent=10,
                            plot_max=50){
  results=c()
  #copy of markers table
  markers_table_copy=enrichment_table
  
  #add in proper gene ratio
  ##total genes first
  total_genes_relevant=as.numeric(unlist(strsplit(as.character(markers_table_copy$GeneRatio[1]),split = '/'))[2])
  gene_ratio_temp=apply(markers_table_copy,MARGIN = 1,function(gene_ratio_func){
    as.numeric(unlist(strsplit(gene_ratio_func[['GeneRatio']],split = '/'))[1])/total_genes_relevant
  })
  markers_table_copy$gene_ratio=gene_ratio_temp
  
  #Max pval
  markers_table_copy$log2pval=log2(markers_table_copy$p.adjust)*-1
  markers_table_copy$log2pval[markers_table_copy$log2pva>log2pval_filter]=log2pval_filter
  
  #Fix descriptions that are very long
  if (!is.null(space_indent)){
    markers_table_copy$Description=as.character(markers_table_copy$Description)
    for (j in 1:nrow(markers_table_copy)){
      markers_table_copy$Description[j]=paste(stri_wrap(markers_table_copy$Description[j], space_indent, 0.0), collapse="\n")
      
    }
  }
  markers_table_copy=markers_table_copy%>%dplyr::arrange(desc(gene_ratio))
  markers_table_copy$Description=factor(markers_table_copy$Description,levels=rev(markers_table_copy$Description))
  
  legend_title=paste('Enriched',enrichment_type,'terms')
  
  li <- c(0, log2pval_filter)
  la <- c(seq(0,log2pval_filter-pval_scale,by=pval_scale),paste0('>', log2pval_filter))
  br <- seq(0, log2pval_filter, pval_scale)
  
  if(!is.null(plot_max)&nrow(markers_table_copy)>plot_max){
    markers_table_copy=markers_table_copy[1:plot_max,]
  }
  
  p=markers_table_copy%>%ggplot(aes(x=gene_ratio,y=Description))+
    geom_point(aes(fill=log2pval,size=Count),colour='black',pch=21)+
    xlab('Gene Ratio')+
    ggtitle(legend_title)+
    scale_fill_gradient(expression('-log'[2]*'(p-val)'),low='white', high = "red2",
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks = TRUE, 
                                               ticks.colour = 'black',
                                               title.position = 'top',
                                               title.hjust=0.5),
                        breaks=br,
                        labels=la,
                        limits=li)+
    guides(size = guide_legend(title.position="top", title.hjust = 0.5))+
    theme(axis.text=element_text(size=text_size),
          axis.title=element_text(size=text_size),
          legend.key.width=unit(1,"cm"),
          plot.title = element_text(size=text_size),
          legend.position="bottom",
          plot.subtitle = element_text(size=text_size-3),
          panel.background = element_rect(fill = "transparent",
                                          colour = "white",
                                          size = 0.5),
          axis.text.y=element_text(size=text_size-2),
          legend.title = element_text(size=text_size),
          legend.text = element_text(size=text_size-2),
          # axis.ticks.x=element_blank(),
          panel.grid.major=element_line(colour='lightgrey',size=0.2),
          axis.line = element_line(size = 0.5),
          plot.margin = margin(r = 1,unit =  "cm")
    )+
    coord_cartesian(clip = "off")
  return(p)
}

#Function to make revigo tree map
rrvgo_treemap=function(ont_sem_sim,
                       go_terms,
                       file_name,
                       save_dir){
  go_terms=go_terms[go_terms[['enrichment']]=='GO',]
  sim_matrix=calculateSimMatrix(go_terms$ID,
                                semdata = ont_sem_sim,
                                orgdb = "org.Hs.eg.db",
                                method="Rel")
  
  scores <- setNames(-log10(go_terms$qvalue), go_terms$ID)
  
  reduce <- reduceSimMatrix(sim_matrix,
                            scores,
                            threshold=0.7,
                            orgdb="org.Hs.eg.db")
  
  pdf(file = paste0(save_dir,file_name,'.pdf'),
      width = 15, 
      height = 15)
  treemapPlot(reduce,
              inflate.labels = TRUE)
  dev.off()
  message('File saved as ',paste0(save_dir,file_name,'.pdf'))
}

get_sem_sim=function(orgdb="org.Hs.eg.db",
                     ont="BP"
){
  temp_sem_sim=GOSemSim::godata(orgdb, ont=ont)
}

label_genes <- function(df, new_column, gene_vector, gene_column, gene_label) {
  # Check if new_column exists
  if(!new_column %in% colnames(df)) {
    # If not, create it with NA values
    df[[new_column]] <- NA
  }
  
  # Replace NA values in new_column with gene_label if gene_vector is in gene_column
  df[[new_column]][df[[gene_column]] %in% gene_vector] <- gene_label
  
  return(df)
}

# extract genes before labelling them
extract_label_genes=function(df,
                             new_df=NULL,
                             new_column,
                             gene_vector,
                             gene_column,
                             gene_label){
  temp_genes=df[df[[gene_column]]%in%gene_vector,]
  temp_genes[[new_column]]=gene_label
  if(!is.null(new_df)){
    new_df=rbind(new_df,temp_genes)
    return(new_df)
  }else{
    return(temp_genes)
  }
}

#Function to plot gene expression by group
compare_expression=function(degs,
                            genes,
                            gene_col='ensembl',
                            group,
                            facet=NULL,
                            facet2=NULL,
                            fill_col='log2FoldChange',
                            pval_col='padj',
                            scale_temp=1,
                            log2_cutoff=log2(1.5),
                            flip_axis=FALSE,
                            fill_cap=NULL,
                            tilt_x=FALSE,
                            hide_ns=TRUE,
                            text_size=10, 
                            facet_text_size=10,
                            hide_legend=FALSE,
                            legend_bottom=FALSE,
                            gene_lab='Gene',
                            condition_lab='Condition',
                            graph_min_new=NULL,
                            graph_max_new=NULL,
                            barwidth=10){
  
  degs_relevant=degs[degs[[gene_col]]%in%genes,]
  
  # Cap the fill values
  if (!is.null(fill_cap)) {
    degs_relevant[[fill_col]] = pmin(degs_relevant[[fill_col]], fill_cap)
  }
  
  degs_relevant[['x']]=degs_relevant[[gene_col]]
  degs_relevant[['y']]=degs_relevant[[group]]
  degs_relevant[['fill']]=degs_relevant[[fill_col]]
  degs_relevant=degs_relevant%>%rstatix::add_significance(pval_col)
  
  if(!is.null(facet)){
    degs_relevant[['facet_temp']]=degs_relevant[[facet]]
  }
  if(!is.null(facet2)){
    degs_relevant[['facet_temp2']]=degs_relevant[[facet2]]
  }
  
  graph_max=max(ceiling(degs_relevant[[fill_col]]/scale_temp)[!is.infinite(ceiling(degs_relevant[[fill_col]]))])*scale_temp
  graph_min=min(floor(degs_relevant[[fill_col]]/scale_temp)[!is.infinite(ceiling(degs_relevant[[fill_col]]))])*scale_temp
  if(graph_min>0){
    graph_min=0
  }
  
  if(!is.null(graph_min_new)){
    graph_min=graph_min_new
  }
  
  if(!is.null(graph_max_new)){
    graph_max=graph_max_new
  }
  
  li <- c(graph_min, graph_max)
  la <- c(seq(graph_min,graph_max,scale_temp))
  br <- c(seq(graph_min,graph_max,scale_temp))
  
  # Modify labels if fill_cap is not NULL
  if (!is.null(fill_cap)) {
    la[la == fill_cap] = paste0(">", fill_cap)
  }
  
  if(!is.null(log2_cutoff)){
    degs_relevant$padj.signif[abs(degs_relevant[[fill_col]])<log2_cutoff]='ns'
  }
  
  if(hide_ns){
    degs_relevant$padj.signif[degs_relevant$padj.signif=='ns']=''
    # degs_relevant$padj.signif=paste0('\n',degs_relevant$padj.signif)
    nudge_y=-0.2
  }else{
    nudge_y=0
  }
  
  if(!flip_axis){
    p_temp=degs_relevant%>%ggplot(aes(x=x,y=y,fill=fill))+
      geom_tile(colour='black')+
      geom_text(aes(label=padj.signif),nudge_y=nudge_y)+
      scale_fill_gradient2(expression('log'[2]*'(FC)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                           guide = guide_colorbar(frame.colour = "black",
                                                  ticks = TRUE, 
                                                  ticks.colour = 'black',title.position = 'top',title.hjust=0.5,
                                                  barwidth = barwidth),
                           breaks=br,
                           labels=la,
                           limits=li)+
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0))+
      xlab(gene_lab)+
      ylab(condition_lab)
  }else{
    p_temp=degs_relevant%>%ggplot(aes(x=y,y=x,fill=fill))+
      geom_tile(colour='black')+
      geom_text(aes(label=padj.signif),nudge_y=nudge_y)+
      scale_fill_gradient2(expression('log'[2]*'(FC)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                           guide = guide_colorbar(frame.colour = "black",
                                                  ticks = TRUE, 
                                                  ticks.colour = 'black',title.position = 'top',title.hjust=0.5),
                           breaks=br,
                           labels=la,
                           limits=li)+
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0))+
      xlab(condition_lab)+
      ylab(gene_lab)
  }
  if (tilt_x) {
    p_temp = p_temp + theme(axis.text.x=element_text(angle=45,hjust=1, size=text_size))
  }
  
  # add size control for axis labels, axis titles, legend and geom_text
  p_temp = p_temp + theme(
    text = element_text(size=text_size),
    axis.title = element_text(size=text_size),
    axis.text = element_text(size=text_size),
    legend.text = element_text(size=text_size),
    strip.text = element_text(size=facet_text_size),
    # modify geom_text size
    plot.subtitle = element_text(size=text_size),
    plot.caption = element_text(size=text_size)
  )
  
  # move the legend to the bottom if legend_bottom is TRUE
  if (legend_bottom) {
    p_temp = p_temp + theme(legend.position = "bottom")
  }
  
  if(!is.null(facet)&is.null(facet2)){
    if(!flip_axis){
      p_temp=p_temp+
        facet_wrap(facet_temp~.,scales = 'free_x')
    }else{
      p_temp=p_temp+
        facet_wrap(facet_temp~.,scales = 'free_y')
    }
  }else if(!is.null(facet)&!is.null(facet2)){
    p_temp=p_temp+
      facet_grid(facet_temp~facet_temp2,scales = 'free',space='free')
  }
  
  if(hide_legend){
    p_temp=
      p_temp + theme(legend.position="none")
  }
  return(p_temp)
}

compare_degs_between_groups=function(all_degs,
                                     group_comparison_col_1,
                                     group_comparison_col_2,
                                     outgroup_col,
                                     group_col,
                                     gene_col='ensembl',
                                     accession_col='dir_accession',
                                     return_overlap=TRUE,
                                     outgroup,
                                     direction_col_1='direction_1',
                                     direction_col_2='direction_2',
                                     xlab=NULL,
                                     ylab=NULL,
                                     temp_order=NULL,
                                     graph_scale=1,
                                     text_size=12){
  all_comparison=unique(all_degs[[group_col]])
  all_degs[[group_comparison_col_1]]=as.character(all_degs[[group_comparison_col_1]])
  all_degs[[group_comparison_col_2]]=as.character(all_degs[[group_comparison_col_2]])
  all_degs[[group_col]]=as.character(all_degs[[group_col]])
  comparison_df=vector_to_df(all_comparison)
  all_overlaps=do.call('rbind',apply(comparison_df,MARGIN=1,function(per_row){
    group_1=per_row[['x']]
    group_2=per_row[['y']]
    
    degs_1=all_degs[all_degs[[group_col]]==group_1,]
    degs_1_relevant=degs_1[degs_1[[group_comparison_col_2]]==outgroup&degs_1[['sig']]=='y',]
    degs_2=all_degs[all_degs[[group_col]]==group_2,]
    degs_2_relevant=degs_2[degs_2[[group_comparison_col_2]]==outgroup&degs_2[['sig']]=='y',]
    
    background_temp=unique(c(degs_1[[gene_col]],
                             degs_2[[gene_col]]))
    
    accession_1=unique(degs_1_relevant[[accession_col]])
    accession_2=unique(degs_2_relevant[[accession_col]])
    overlap_temp=do.call('rbind',lapply(accession_1,function(per_accession_1){
      degs_1_relevant_temp=degs_1_relevant[degs_1_relevant[[accession_col]]==per_accession_1,]
      degs_1_relevant_temp_gene=degs_1_relevant_temp[[gene_col]]
      
      n_1=length(degs_1_relevant_temp_gene)
      n_1_name=degs_1_relevant_temp[[group_col]][1]
      n_1_dir_ingroup=degs_1_relevant_temp[[direction_col_1]][1]
      n_1_dir_time=degs_1_relevant_temp[[group_comparison_col_1]][1]
      do.call('rbind',lapply(accession_2,function(per_accession_2){
        degs_2_relevant_temp=degs_2_relevant[degs_2_relevant[[accession_col]]==per_accession_2,]
        degs_2_relevant_temp_gene=degs_2_relevant_temp[[gene_col]]
        overlap=testGeneOverlap(newGeneOverlap(listA=degs_1_relevant_temp_gene,
                                               listB=degs_2_relevant_temp_gene,
                                               genome.size=length(background_temp)))
        con_tbl=getContbl(overlap)
        actual=con_tbl[4]
        expected=chisq.test(con_tbl)$exp[4]
        n_diff=actual-expected
        odds=fisher.test(con_tbl)$est
        pval=fisher.test(con_tbl)$p.val
        overlap=paste0(getIntersection(overlap),collapse='/')
        
        n_2=length(degs_2_relevant_temp_gene)
        
        n_2_name=degs_2_relevant_temp[[group_col]][1]
        n_2_dir_ingroup=degs_2_relevant_temp[[direction_col_1]][1]
        n_2_dir_time=degs_2_relevant_temp[[group_comparison_col_1]][1]
        
        if(return_overlap){
          full_return=matrix(c(actual,expected,n_diff,odds,pval,
                               n_1_name,n_1_dir_ingroup,n_1_dir_time,n_1,
                               n_2_name,n_2_dir_ingroup,n_2_dir_time,n_2,
                               length(background_temp),overlap),ncol=15)
          colnames(full_return)=c('actual','expected','diff','odds','pval',
                                  'group_1','dir_ingroup_1','time_1','n_1',
                                  'group_2','dir_ingroup_2','time_2','n_2',
                                  'background_n','overlap')
        }else{
          full_return=matrix(c(actual,expected,n_diff,odds,pval,
                               n_1_name,n_1_dir_ingroup,n_1_dir_time,n_1,
                               n_2_name,n_2_dir_ingroup,n_2_dir_time,n_2,
                               length(background_temp)),ncol=14)
          colnames(full_return)=c('actual','expected','diff','odds','pval',
                                  'group_1','dir_ingroup_1','time_1','n_1',
                                  'group_2','dir_ingroup_2','time_2','n_2',
                                  'background_n')
        }
        
        return(full_return)
      }))
    }))
    overlap_temp=data.frame(overlap_temp)
    return(overlap_temp)
  }))
  all_overlaps[['actual']]=as.numeric(as.character(all_overlaps[['actual']]))
  all_overlaps[['expected']]=as.numeric(as.character(all_overlaps[['expected']]))
  all_overlaps[['diff']]=as.numeric(as.character(all_overlaps[['diff']]))
  all_overlaps[['pval']]=as.numeric(as.character(all_overlaps[['pval']]))
  all_overlaps[['odds']]=as.numeric(as.character(all_overlaps[['odds']]))
  all_overlaps[['n_1']]=as.numeric(as.character(all_overlaps[['n_1']]))
  all_overlaps[['n_2']]=as.numeric(as.character(all_overlaps[['n_2']]))
  all_overlaps[['background_n']]=as.numeric(as.character(all_overlaps[['background_n']]))
  
  all_overlaps$p.adj=p.adjust(all_overlaps[['pval']])
  all_overlaps$log2odds=log2(all_overlaps[['odds']])
  all_overlaps=all_overlaps%>%rstatix::add_significance('p.adj')
  
  all_overlaps[['accession']]=paste0(all_overlaps[['group_1']],' ',
                                     all_overlaps[['dir_ingroup_1']],
                                     '\n',all_overlaps[['group_2']],' ',
                                     all_overlaps[['dir_ingroup_2']]
  )
  all_overlaps[['x']]=paste0(all_overlaps[['group_1']],' ',gsub(all_overlaps[['time_1']],pattern='_',replacement=' '))
  all_overlaps[['y']]=paste0(all_overlaps[['group_2']],' ',gsub(all_overlaps[['time_2']],pattern='_',replacement=' '))
  
  #Fix some legend parameters
  graph_max=max(ceiling(all_overlaps$log2odds/2)[!is.infinite(ceiling(all_overlaps$log2odds))])*2
  graph_min=min(floor(all_overlaps$log2odds)[!is.infinite(ceiling(all_overlaps$log2odds))])
  if(graph_min>0){
    graph_min=0
  }
  li <- c(graph_min, graph_max)
  la <- c(seq(graph_min,graph_max,graph_scale))
  br <- c(seq(graph_min,graph_max,graph_scale))
  
  if(!is.null(temp_order)){
    x_relevant=unique(all_overlaps[['x']])
    x_order=do.call(c,lapply(temp_order,function(per_order){
      x_relevant[grep(x_relevant,pattern=per_order)]
    }))
    all_overlaps[['x']]=factor(all_overlaps[['x']],
                               levels=x_order)
    
    y_relevant=unique(all_overlaps[['y']])
    y_order=do.call(c,lapply(temp_order,function(per_order){
      y_relevant[grep(y_relevant,pattern=per_order)]
    }))
    all_overlaps[['y']]=factor(all_overlaps[['y']],
                               levels=rev(y_order))
  }
  
  p_temp=all_overlaps%>%ggplot(aes(x=x,y=y,fill=log2odds))+
    geom_tile(colour='black')+
    geom_text(aes(label=p.adj.signif),nudge_y = 0.25, size = text_size / 2)+
    geom_text(aes(label=actual),nudge_y=-0.25, size = text_size / 2)+
    scale_fill_gradient2(expression('log'[2]*'(Odds)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks = TRUE, 
                                                barwidth = 8,
                                                ticks.colour = 'black',
                                                title.position = 'top',title.hjust=0.5),
                         breaks=br,
                         labels=la,
                         limits=li)+
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0))+
    facet_wrap(~accession,scales='free')+
    theme(axis.text.x=element_text(angle=45,hjust=1, size = text_size),
          axis.text.y = element_text(size = text_size),
          legend.title = element_text(size = text_size),
          legend.text = element_text(size = text_size),
          strip.text = element_text(size = text_size),
          legend.position='bottom')+
    xlab(xlab)+
    ylab(ylab)
  
  return_me=c()
  return_me[['df']]=all_overlaps
  return_me[['p']]=p_temp
  return(return_me)
}

# function to make summary plot of overlaps
summarise_overlaps = function(db,
                              x,
                              y,
                              odds_col = 'odds',
                              pval_col = 'adj',
                              xlab = NULL,
                              ylab = NULL,
                              add_line = TRUE,
                              legend_title = NULL,
                              label_cq = TRUE,
                              remove_self = FALSE,
                              self_col = c('group_1','group_1.1'),
                              legend_position_bottom = FALSE,
                              rotate_x_text = TRUE) {
  
  db$log2odds = log2(db[[odds_col]])
  db$fill = ifelse(db[[pval_col]] > 0.05, 'Not Sig', ifelse(db$log2odds > 0, 'Sig Over-represented', 'Sig Under-represented'))
  
  if(!label_cq) {
    db$fill = factor(db$fill,
                     levels = c('Sig Over-represented', 'Sig Under-represented', 'Not Sig'))
    fill_manual = c("Not Sig" = "#f1f1f1", "Sig Over-represented" = "#f8766d",
                    "Sig Under-represented" = "#00bfc4")
  } else {
    new_legend = apply(db, 1, function(per_row) {
      temp_fill = per_row[['fill']]
      if(temp_fill == 'Not Sig') {
        return(temp_fill)
      } else {
        paste0(ifelse(grepl(per_row[[y]], pattern = 'CQ'), 'CQ', 'CS'), ' ', per_row[['fill']])
      }
    })
    
    db$fill = new_legend
    db$fill = factor(db$fill,
                     levels = c("CS Sig Over-represented",
                                "CQ Sig Over-represented",
                                "CS Sig Under-represented",
                                "CQ Sig Under-represented",
                                'Not Sig'))
    fill_manual = c("Not Sig" = "#f1f1f1", "CS Sig Over-represented" = "#f8766d",
                    'CQ Sig Over-represented' = '#ffb3b3',
                    'CQ Sig Under-represented' = '#b3ffff',
                    "CS Sig Under-represented" = "#00bfc4")
  }
  
  db[['x']] = db[[x]]
  db[['y']] = db[[y]]
  
  if(remove_self) {
    db$fill[db[[self_col[1]]] == db[[self_col[2]]]] = NA
  }
  
  temp_p = db %>% ggplot(aes(x = x, y = y, fill = fill)) +
    geom_tile(colour = 'black') +
    scale_x_discrete(expand = expansion(c(0, 0))) +
    scale_y_discrete(expand = expansion(c(0, 0))) +
    ylab(ylab) +
    xlab(xlab) +
    scale_fill_manual(name = legend_title,
                      values = fill_manual)
  
  if(add_line) {
    yint = (length(unique(db[['y']])) / 2) + 0.5
    temp_p = temp_p + geom_hline(yintercept = yint, size = 1)
  }
  
  # Adjust x-axis text rotation if rotate_x_text is TRUE
  if(rotate_x_text) {
    temp_p = temp_p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # Adjust legend position if legend_position_bottom is TRUE
  if(legend_position_bottom) {
    temp_p = temp_p + theme(legend.position = "bottom")
  }
  
  return(temp_p)
}

summarise_overlaps_facet=function(db,
                            x,
                            y,
                            odds_col='odds',
                            pval_col='adj',
                            xlab=NULL,
                            ylab=NULL,
                            add_line=TRUE,
                            legend_title=NULL,
                            label_cq=TRUE,
                            facet_col,
                            angle=90){
  db$log2odds=log2(db[[odds_col]])
  db$fill=ifelse(db[[pval_col]]>0.05,'Not Sig',ifelse(db$log2odds>0,'Sig Over-represented','Sig Under-represented'))
  
  db$fill=factor(db$fill,
                   levels=c('Sig Over-represented','Sig Under-represented','Not Sig'))
  fill_manual=c("Not Sig" = "#f1f1f1", "Sig Over-represented" = "#f8766d",
                  "Sig Under-represented" = "#00bfc4")
  
  if(angle==90){
    angle_x=element_text(angle=90,vjust=0.5,hjust=1)
  }else if(angle==45){
    angle_x=element_text(angle=45,hjust=1)
  }
  db$facet_temp=db[[facet_col]]
  db[['x']]=db[[x]]
  db[['y']]=db[[y]]
  temp_p=
    db%>%ggplot(aes(x=x,y=y,fill=fill))+
    geom_tile(colour='black')+
    scale_x_discrete(expand=expansion(c(0,0)))+
    scale_y_discrete(expand=expansion(c(0,0)))+
    theme(axis.text.x=angle_x)+
    ylab(ylab)+
    xlab(xlab)+
    scale_fill_manual(name=legend_title,
                      values = fill_manual)+
  facet_grid(facet_temp~.)
    
  if(add_line){
    yint=(length(unique(db[['y']]))/2)+0.5
    
    temp_p=temp_p+
      geom_hline(yintercept = yint,size=1)
  }
  return(temp_p)
}

# ---------------------------------------------------------------------------
# filter_dataframe()
#
# Carried over verbatim from Final/Scripts/for_github/study_degs.R (origin/main),
# which was never merged into restructure-bulk-rna-seq. It is needed to
# reproduce SI Figure 16 (the SIPS/OIS/CQ PCA): it keeps only those studies
# that contain at least one sample matching each of two label patterns, so the
# CQ-versus-senescence comparison is never made across studies that contribute
# only one side of it. Replicative senescence drops out of that figure as a
# consequence, not by an explicit exclusion -- no study in the pool contributes
# both CQ and RS samples.
# ---------------------------------------------------------------------------
filter_dataframe <- function(dataframe,
                             iterate_col,
                             filter_col,
                             col_label_1,
                             col_label_2,
                             col_1_n_filter = 0,
                             col_2_n_filter = 0) {
  dataframe = dataframe[grepl(dataframe[[filter_col]], pattern = col_label_1) |
                          grepl(dataframe[[filter_col]], pattern = col_label_2), ]
  filter_val = do.call(c, lapply(unique(dataframe[[iterate_col]]), function(value) {
    subset_df <- dataframe[dataframe[[iterate_col]] == value, ]
    counts_col_label_1 <- sum(grepl(subset_df[[filter_col]], pattern = col_label_1)) > col_1_n_filter
    counts_col_label_2 <- sum(grepl(subset_df[[filter_col]], pattern = col_label_2)) > col_2_n_filter
    if (counts_col_label_1 & counts_col_label_2) {
      return(value)
    } else {
      return(NULL)
    }
  }))
  filtered_df <- dataframe[dataframe[[iterate_col]] %in% filter_val, ]
  return(filtered_df)
}

library(magrittr)
library(stringi)
#Save ggplot
save_p <- function(plot, file_name, save_dir, p_width = 7, p_height = 7,
                   add_date=FALSE,
                   add_time=FALSE) {
  if(!dir.exists(save_dir)){
    stop('Directory does not exist')
  }
  
  if(add_date){
    date_add=paste0(Sys.Date(),'_')
  }else{
    date_add=''
  }
  
  if(add_time){
    time_add=paste0(Sys.time(),'_')
  }else{
    time_add=''
  }
  file_name=paste0(date_add,time_add,file_name)
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
                     gene_col='gene',
                     add_date=FALSE,
                     add_time=FALSE) {
  
  if(merge_ensembl_symbol){
    data <- merge_gene_symbol_from_ensembl(df = data,
                                           ensembl_col = gene_col)
  }
  
  if(!dir.exists(path)){
    stop('Directory does not exist')
  }
  
  # Add date and time if requested
  date_add <- if (add_date) paste0(Sys.Date(), "_") else ""
  time_add <- if (add_time) paste0(format(Sys.time(), "%H-%M-%S"), "_") else ""
  
  # Modify file name
  file_name <- paste0(date_add, time_add, file_name)
  
  # Ensure .csv extension
  if (!grepl("\\.csv$", file_name, ignore.case = TRUE)) {
    file_name <- paste0(file_name, ".csv")
  }
  
  # Construct full file path
  file_path <- file.path(path, file_name)
  
  # Write the data frame to a CSV file
  write.csv(data, file = file_path, row.names = row.names)
  
  print(paste("Data saved to:", file_path))
}

#
replace_ensembl_with_gene_names <- function(rse, biomart_df, 
                                            ensembl_col_in_biomart="ensembl_gene_id", 
                                            gene_name_col_in_biomart="external_gene_names") {
  rse_use <- rse
  temp_rownames <- rownames(rse_use)
  
  biomart_df_use <- biomart_df[biomart_df[[ensembl_col_in_biomart]] %in% temp_rownames, ]
  
  # Reorder biomart_df_use to match the order of temp_rownames
  biomart_df_use <- biomart_df_use[match(temp_rownames, biomart_df_use[[ensembl_col_in_biomart]]), ]
  
  # Update row names of rse_use
  rownames(rse_use) <- biomart_df_use[[gene_name_col_in_biomart]]
  
  # Update row names of the assays within rse_use
  for (assay_name in assayNames(rse_use)) {
    assay_data <- assay(rse_use, assay_name)
    rownames(assay_data) <- rownames(rse_use)
    assay(rse_use, assay_name) <- assay_data
  }
  
  return(rse_use)
}

#Function to filter counts by min counts and samples
process_rse=function(rse,
                     group_col='group',
                     ensembl_dictionary=NULL,
                     ensembl_gene_col='ensembl_gene_id',
                     ensembl_symbol_col='external_gene_name',
                     filter_duplicates=TRUE,
                     sample_filter=NULL,
                     low_expression_filter=c('all','per'),
                     min_expression_filter=1,
                     min_perc_sample_filter=0.3,
                     method='cpm',
                     replace_ensembl=TRUE){
  meta_condition=data.frame(colData(rse))
  # filter for samples
  if(!is.null(sample_filter)){
    rse=rse[,colnames(rse)%in%sample_filter]
  }
  
  if(replace_ensembl&is.null(ensembl_dictionary)){
    message('must give ensembl_dictionary')
    return(NULL)
  }
  #clean gene names
  # rownames(rse)=clean_id(rownames(rse))
  
  if(replace_ensembl){
    rse=replace_ensembl_with_gene_names(rse,
                                        biomart_df = ensembl_dictionary,
                                        ensembl_col_in_biomart = ensembl_gene_col,
                                        gene_name_col_in_biomart = ensembl_symbol_col)
  }
  
  if(!is.null(low_expression_filter)){
    message('before removing lowly expressed genes: ',length(rownames(rse)))
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
      message('after removing lowly expressed genes: ',length(rownames(rse)))
    }else if(low_expression_filter=='per'){
      all_condition=unique(rse[[group_col]])
      counts_temp=assay(rse)
      genes_keep=do.call(c,lapply(all_condition,function(per_cond){
        temp_samples=rownames(meta_condition[meta_condition[[group_col]]==per_cond,,drop=FALSE])
        
        counts_temp_assay=counts_temp[,colnames(counts_temp)%in%temp_samples]
        
        counts_temp_assay=filter_counts(dds=counts_temp_assay,
                                        min_expression_filter=min_expression_filter,
                                        min_perc_sample_filter=min_perc_sample_filter,
                                        method=method)
        genes_temp=rownames(counts_temp_assay)
      }))
      genes_keep=unique(genes_keep)
      
      rse=rse[rownames(rse)%in%genes_keep,]
      message('after removing lowly expressed genes: ',length(rownames(rse)))
    }
  }
  
  rse=rse[rowSums(assay(rse)) != 0, ]
  #remove duplicates
  if(filter_duplicates){
    message('before duplicates: ',length(rownames(rse)))
    gene_duplicates=unique(rownames(rse)[duplicated(rownames(rse))])
    rse=rse[!rownames(rse)%in%gene_duplicates,]
    message('after duplicates: ',length(rownames(rse)))
  }
  
  return(rse)
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
                 null_pca='group',
                 remove_pca=NULL,
                 nrow=2,
                 legend_nrow=4,
                 id_col='sample'){
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
                      intgroup=per_pca,
                      id_col = id_col)
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
                    hide_legend=FALSE,
                    id_col='sample'){
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
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=pheno_data[[id_col]])
  
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

find_degs_new=function(rse_obj,
                       group_col,
                       group_1,
                       group_2,
                       pval_cutoff=0.05,
                       log2fc_cutoff=log2(1.5),
                       independentFiltering=TRUE,
                       batch_col=NULL,
                       find_degs=TRUE,
                       sample_col='sample'){
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

#Wrapper for cleaning function
clean_degs=function(degs_obj){
  degs_obj_temp_test=clean_col_names(degs_obj)
  degs_obj_temp_test=degs_accession(degs_obj_temp_test)
  degs_obj_temp_test=degs_direction(degs_obj_temp_test)
  degs_obj_temp_test=degs_accession(degs_obj_temp_test,cols_accession=c('group_1_dir','group_2_dir'),
                                    accession_col_name='dir_accession')
  return(degs_obj_temp_test)
}

#Function to remove '-' and ' ' from a character vector
clean_name=function(name){
  name_temp=as.character(name)
  name_temp=gsub(name_temp,pattern='-',replacement='_')
  name_temp=gsub(name_temp,pattern=' ',replacement='_')
  name_temp=gsub(name_temp,pattern='\n',replacement='_')
  return(name_temp)
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
  results$adj=p.adjust(results$pval,method = 'BH')
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
                              remove_self = FALSE,
                              self_col = c('group_1','group_1.1'),
                              legend_position_bottom = FALSE,
                              rotate_x_text = TRUE,
                              text_size = 12) {  # Added text_size argument
  
  db$log2odds = log2(db[[odds_col]])
  db$fill = ifelse(db[[pval_col]] > 0.05, 'Not Sig', ifelse(db$log2odds > 0, 'Sig Over-represented', 'Sig Under-represented'))
  
  db$fill = factor(db$fill,
                   levels = c('Sig Over-represented', 'Sig Under-represented', 'Not Sig'))
  fill_manual = c("Not Sig" = "#f1f1f1", "Sig Over-represented" = "#f8766d",
                  "Sig Under-represented" = "#00bfc4")
  
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
                      values = fill_manual) +
    theme(
      axis.text = element_text(size = text_size),  # Apply text_size
      axis.title = element_text(size = text_size - 2),  # Slightly smaller title text
      legend.text = element_text(size = text_size),
      legend.title = element_text(size = text_size + 2),
      strip.text = element_text(size = text_size + 2)  # Larger facet text
    )
  
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


#function to filter lowly expressed counts
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

#function to plot pca from rse object
plot_pca_rse <- function(rse, assay = "counts",
                         group_by = NULL,
                         top_n_genes = NULL,
                         label_text=FALSE) {
  
  # Extract assay data
  assay_data <- assay(rse, assay)
  
  # Option to filter top genes with the highest variance (if provided)
  if (!is.null(top_n_genes)) {
    variances <- apply(assay_data, 1, var)
    top_genes <- order(variances, decreasing = TRUE)[1:top_n_genes]
    assay_data <- assay_data[top_genes, ]
  }
  
  # Perform PCA
  pca <- prcomp(t(assay_data), scale. = TRUE)
  
  # Calculate variance explained by each PC
  percent_variance <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  pc1_var <- round(percent_variance[1], 1)
  pc2_var <- round(percent_variance[2], 1)
  
  # Create a data frame with PCA results
  pca_data <- as.data.frame(pca$x)
  pca_data$Sample <- rownames(pca_data)
  
  # Extract grouping information from colData
  coldata <- as.data.frame(colData(rse))
  if (!is.null(group_by)) {
    if (!all(group_by %in% colnames(coldata))) {
      stop("One or more specified group_by columns are not present in colData.")
    }
    pca_data <- cbind(pca_data, coldata[, group_by, drop = FALSE])
  }
  
  # Prepare ggplot2 PCA plot
  p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes_string(color = group_by)) +
    theme_minimal() +
    labs(
      title = "PCA Plot",
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    )
  
  if(label_text){
    
    p <- p + geom_text(aes(label = Sample), size = 3)
    
  }
  return(p)
}

#get edger degs for one outgroup by ingroup
edger_degs=function(rse,
                    group_col,
                    ingroup,
                    outgroup,
                    assay_slot,
                    log_transform=FALSE,
                    fc_cutoff=log2(1.5)){
  if(log_transform){
    assay_slot_old=assay_slot
    assay_slot=paste0('log2',assay_slot)
    rse@assays@data[[assay_slot]] <-
      log2(rse@assays@data[[assay_slot_old]]+1)
    
  }
  get_degs=do.call('rbind',lapply(ingroup,function(per_group){
    rse_temp=rse[,rse[[group_col]]%in%c(per_group,outgroup)]
    
    temp_counts=rse_temp@assays@data[[assay_slot]]
    
    # Extract group information
    group_labels <- colData(rse_temp)[[group_col]]
    
    # Ensure group labels are a factor and set the reference level to `outgroup`
    group_labels <- factor(group_labels, levels = c(outgroup, per_group))
    
    # Create a design matrix for the model
    design <- model.matrix(~ group_labels)
    
    counts <- rse_temp@assays@data[[assay_slot]]
    
    # Run the limma-trend pipeline
    fit <- lmFit(counts, design)
    
    if(log_transform){
      fit <- eBayes(fit, trend = TRUE)
      
      results <- topTable(fit, coef = 2, adjust = "BH", number = Inf)
      results$group_1=per_group
      results$group_2=outgroup
      
      results$sig=ifelse(results$adj.P.Val<0.05&
                           abs(results$logFC)>fc_cutoff,
                         TRUE,FALSE)
      results=results%>%rownames_to_column('gene')
      return(results)
    }else{
      message('need to code')
      return(NULL)
    }
  }))
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
  gene_list=gene_list%>%dplyr::select(dplyr::all_of(c('group','gene_col')))
  
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
                           minGSSize=0,
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

#function to only plot significant enrichment IDs
plot_gsea_enrichment_sig=function(enrichment_df,
                                  condition_col,
                                  pathway_col,
                                  pval_col='p.adjust',
                                  factor_condition_levels_x=NULL,
                                  factor_condition_levels_facet=NULL,
                                  tilt_x=90,
                                  pval_filter=0.05){
  sig_id=unique(enrichment_df[[pathway_col]][enrichment_df[[pval_col]]<pval_filter])
  plot_gsea_enrichment(enrichment_df = enrichment_df[enrichment_df[[pathway_col]]%in%sig_id,],
                       pathway_col = pathway_col,condition_col = condition_col,
                       pval_col=pval_col,factor_condition_levels_x=factor_condition_levels_x,
                       factor_condition_levels_facet=factor_condition_levels_facet,
                       tilt_x=tilt_x)
}

#Function to build the msigdb db using a file path
build_msig_db=function(file_path){
  # Read the file into a dataframe
  data <- readLines(file_path)
  
  # Create an empty dataframe to store the results
  result <- data.frame(pathway = character(), genes = character(), stringsAsFactors = FALSE)
  
  # Loop through each line of the file and process it
  for (line in data) {
    # Split the line by tab
    parts <- strsplit(line, "\t")[[1]]
    
    # Extract the pathway and the genes
    pathway <- parts[1]
    genes <- parts[-c(1,2)]
    
    # Create a temporary dataframe for the current line
    temp_df <- data.frame(pathway = rep(pathway, length(genes)), genes = genes, stringsAsFactors = FALSE)
    
    # Bind the temporary dataframe to the result dataframe
    result <- rbind(result,temp_df)
  }
  result$class='MSigDB'
  return(result)
}

# Function to transform a dataframe into a list of pathways with their genes
create_pathway_gene_list <- function(df, pathway_col, gene_col) {
  # Split the dataframe into a list based on unique pathway names
  pathway_list <- split(df[[gene_col]], df[[pathway_col]])
  
  # Return the resulting list
  return(pathway_list)
}

#Function to label deg direction and significance
label_degs=function(df,
                    fc_col='log2FoldChange',
                    pval_col='padj',
                    log_fc_cutoff=1.5,
                    pval_cutoff=0.05){
  df$sig=ifelse(abs(df[[fc_col]])>log2(log_fc_cutoff)&
                  df[[pval_col]]<pval_cutoff,TRUE,FALSE)
  df$dir=ifelse(df[[fc_col]]>0,'up','down')
  df[['pval_filter']]=pval_cutoff
  df[['fc_filter']]=log_fc_cutoff
  return(df)
}

#function to make upset plots
create_upset_plot <- function(deg_df,
                              group,
                              dir,
                              gene_col,
                              sig_col='sig') {
  
  deg_df=deg_df[deg_df[[sig_col]],]
  
  # Transform the input data frame and create the upset plot
  df_wide <- deg_df %>%
    mutate(set_id = paste(.data[[group]], .data[[dir]], sep = "_")) %>%
    select(.data[[gene_col]], set_id) %>%
    mutate(value = 1) %>%
    tidyr::pivot_wider(names_from = set_id, values_from = value, values_fill = 0)
  
  # Create the upset plot
  ComplexUpset::upset(df_wide, intersect = colnames(df_wide)[!colnames(df_wide)==gene_col])
}

#function to compare gene expression across samples
compare_sample_expression=function(expression_matrix,
                                   goi,
                                   gene_col='external_gene_name',
                                   col_compare){
  expression_matrix_filtered=expression_matrix[expression_matrix[[gene_col]]%in%
                                                 goi,]
  
  expression_matrix_filtered=expression_matrix_filtered[,colnames(expression_matrix_filtered)%in%
                                                          c(gene_col,col_compare)]
  
  expression_matrix_filtered=reshape2::melt(expression_matrix_filtered,id=gene_col)
  
  expression_matrix_filtered$x=expression_matrix_filtered[[gene_col]]
  
  expression_matrix_filtered%>%ggplot(aes(x=variable,y=value))+
    geom_point()+
    facet_wrap(~x)+
    theme(axis.text.x=element_text(angle=45,hjust=1))
}

#function to perform GO gsea
per_group_gsea_go=function(df,
                           gene_col='ensembl',
                           group_col,
                           pval_col=NULL,
                           log2fc_col=NULL,
                           rank_type = c("pval", "fc", "both"),
                           keyType='ENSEMBL',
                           minGSSize = 0,
                           maxGSSize = Inf,
                           pvalueCutoff = 1,
                           OrgDb,
                           pAdjustMethod = "BH"){
  ont_temp=c("BP", "MF", "CC")
  
  all_group=unique(df[[group_col]])
  all_enrich=do.call('rbind',lapply(all_group,function(per_group){
    message(per_group)
    df_temp=df[df[[group_col]]==per_group,]
    
    df_temp=rank_genelist(genes = df_temp,
                          gene_col = gene_col,
                          pval_col = pval_col,
                          log2fc_col = log2fc_col,
                          rank_type = rank_type)
    
    df_temp=df_temp%>%dplyr::arrange(desc('ranking_score'))
    
    temp_list=df_temp[['ranking_score']]
    names(temp_list)=df_temp[['gene']]
    
    temp_gse=do.call('rbind',lapply(ont_temp,function(per_ont){
      message(per_ont)
      gse <- data.frame(gseGO(geneList=temp_list, 
                              ont = per_ont, 
                              keyType = keyType, 
                              minGSSize = minGSSize, 
                              maxGSSize = maxGSSize, 
                              pvalueCutoff = pvalueCutoff, 
                              verbose = TRUE, 
                              OrgDb = OrgDb, 
                              pAdjustMethod = pAdjustMethod)
      )
      gse$ont=per_ont
      return(gse)
    }))
    temp_gse$group=per_group
    return(temp_gse)
  }))
  all_enrich$dir=ifelse(all_enrich$enrichmentScore>0,'up','down')
  return(all_enrich)
}

# Function to perform KEGG GSEA for mice
per_group_gsea_kegg <- function(df,
                                gene_col = 'entrezgene_id',
                                group_col,
                                pval_col = NULL,
                                log2fc_col = NULL,
                                rank_type = c("pval", "fc", "both"),
                                organism = 'mmu',  # Use 'mmu' for mouse
                                keyType = 'kegg',
                                minGSSize = 0,
                                maxGSSize = Inf,
                                pvalueCutoff = 1,
                                pAdjustMethod = "BH") {
  
  # Get unique groups
  all_group <- unique(df[[group_col]])
  
  # Perform enrichment for each group
  all_enrich <- do.call('rbind', lapply(all_group, function(per_group) {
    message(per_group)
    df_temp <- df[df[[group_col]] == per_group,]
    
    gene_dup=df_temp[[gene_col]][duplicated(df_temp[[gene_col]])]
    
    df_temp=df_temp[!df_temp[[gene_col]]%in%gene_dup,]
    
    # Rank genes for enrichment analysis
    df_temp <- rank_genelist(genes = df_temp,
                             gene_col = gene_col,
                             pval_col = pval_col,
                             log2fc_col = log2fc_col,
                             rank_type = rank_type)
    
    df_temp <- df_temp %>% dplyr::arrange(desc(ranking_score))
    
    # Prepare ranked list
    temp_list <- df_temp[['ranking_score']]
    names(temp_list) <- df_temp[['gene']]
    
    # Perform KEGG enrichment
    temp_gse <- data.frame(gseKEGG(geneList = temp_list,
                                   organism = organism,
                                   keyType = 'ncbi-geneid',
                                   minGSSize = minGSSize,
                                   maxGSSize = maxGSSize,
                                   pvalueCutoff = pvalueCutoff,
                                   pAdjustMethod = pAdjustMethod))
    temp_gse$group <- per_group
    return(temp_gse)
  }))
  
  # Add direction of enrichment
  all_enrich$dir <- ifelse(all_enrich$enrichmentScore > 0, 'up', 'down')
  
  all_enrich$ont='kegg'
  
  return(all_enrich)
}

#function to plot per_group_gsea_go results
##TODO finish find_opp
generate_enrichment_plots <- function(data,
                                      top_n = 20,
                                      p_adjust_threshold = 0.01,
                                      find_opp=FALSE,
                                      remove_dash=FALSE) {
  results <- list()  # Initialize an empty list to store the plots
  
  lapply(unique(data$ont), function(per_ont) {
    deg_all_merged_go_temp <- data[data$ont == per_ont, ]
    
    lapply(c('up', 'down'), function(per_dir) {
      deg_all_merged_go_temp_dir <- deg_all_merged_go_temp[deg_all_merged_go_temp$dir == per_dir, ]
      deg_all_merged_go_temp_dir <- deg_all_merged_go_temp_dir %>% dplyr::arrange(enrichmentScore)
      
      deg_all_merged_go_temp_dir_sig <- unique(deg_all_merged_go_temp_dir[deg_all_merged_go_temp_dir$p.adjust < p_adjust_threshold, ][['ID']])
      if(length(deg_all_merged_go_temp_dir_sig)==0){
        return(NULL)
      }
      if (length(deg_all_merged_go_temp_dir_sig) > top_n) {
        deg_all_merged_go_temp_dir_sig <- deg_all_merged_go_temp_dir_sig[1:top_n]
        p_title_temp <- top_n
      } else {
        p_title_temp <- length(deg_all_merged_go_temp_dir_sig)
      }
      
      enrich_use <- deg_all_merged_go_temp_dir[deg_all_merged_go_temp_dir$ID %in% deg_all_merged_go_temp_dir_sig, ]
      
      if(remove_dash){
        cleaned_strings <- sub(" - [^-]+$", "", enrich_use$Description)
        enrich_use$Description=cleaned_strings
      }
      
      gsea_msig_p_filtered <- plot_gsea_enrichment_sig(enrichment_df = enrich_use,
                                                       pathway_col = 'Description',
                                                       condition_col = 'group') +
        ggtitle(paste0('Top ', p_title_temp, ' ', per_ont, ' ', per_dir))
      
      # Save the plot with a name in the results list
      results[[paste0(per_ont, '_', per_dir)]] <<- gsea_msig_p_filtered
    })
    
    # if(find_opp){
    #   #look at opposites
    #   deg_all_merged_go_temp_opp=deg_all_merged_go_temp
    #   
    #   deg_all_merged_go_temp_opp=filter_opposing_directions(deg_all_merged_go_temp_opp)
    # }
  })
  
  return(results)
}

#function to filter for processes that change in opposing directions between two conditions
##TODO finish
filter_opposing_directions <- function(data, p_adjust_threshold = 0.01) {
  # Filter for significant processes based on p.adjust threshold
  significant_data <- data[data$p.adjust < p_adjust_threshold, ]
  
  ##filter for 2
  temp_table=table(significant_data$ID)
  temp_table=names(temp_table[temp_table==2])
  
  significant_data_plot_temp=significant_data_plot[significant_data_plot$ID%in%
                                                     temp_table,]
  ##filter for 1
  temp_table_dir=data.frame(table(significant_data_plot_temp))
  
  # Filter the original data for these opposing IDs
  filtered_data <- data[data$ID %in% opposing_ids, ]
  
  return(filtered_data)
}

#Function to filter for just significant groups to plot in enrichment dfs
filter_sig_enrich=function(enrichment_df,
                           pval_col='p.adjust',
                           id='ID',
                           pval_cutoff=0.05,
                           keep_nonsig=TRUE){
  if(keep_nonsig){
    sig_terms=enrichment_df[[id]][enrichment_df[[pval_col]]<pval_cutoff]
    return(enrichment_df[enrichment_df[[id]]%in%sig_terms,])
  }else{
    enrichment_df=enrichment_df[enrichment_df[[pval_col]]<pval_cutoff,]
    return(enrichment_df)
  }
}

#dotplot GSEA
plot_dotplot_enrichment <- function(enrichment_df,
                                    x_axis_col = "group",
                                    y_axis_col = "ID",
                                    size_col = NULL, # This will be calculated from pvalue_col if NULL
                                    pvalue_col = "p.adjust", # Specify the column for p-values
                                    colour_col = "enrichmentScore",
                                    label_sig = FALSE,
                                    tilt_x = 90,
                                    space_indent=60,
                                    axis_indent=10,
                                    filter_sig = TRUE,
                                    pval_cutoff = 0.05,
                                    colour_sig = FALSE,
                                    nudge_sig = 0,
                                    text_size = 10,
                                    str_remove=NULL,
                                    top_x=NULL,
                                    size_legend_title=NULL,
                                    colour_legend_title=NULL,
                                    plot_title=NULL,
                                    title_add=NULL,
                                    title_sep='_',
                                    xlab=NULL,
                                    ylab=NULL,
                                    facet=NULL,
                                    col_legend=TRUE,
                                    order_pval=TRUE) {
  if(!is.null(facet)){
    enrichment_df[['facet']]=enrichment_df[[facet]]
  }
  if(!is.null(str_remove)){
    enrichment_df[[y_axis_col]]=gsub(pattern=str_remove,
                                     replacement='',
                                     x=enrichment_df[[y_axis_col]],
                                     fixed = TRUE)
  }
  
  if (filter_sig) {
    enrichment_df_use <- filter_sig_enrich(enrichment_df = enrichment_df,
                                           pval_col = pvalue_col,
                                           id = y_axis_col,
                                           pval_cutoff = pval_cutoff,
                                           keep_nonsig = FALSE)
  } else {
    enrichment_df_use <- enrichment_df
  }
  
  # Calculate log2(pvalue) if size_col is not provided
  if (is.null(size_col)) {
    enrichment_df_use[['log2_pvalue']] <- -log2(enrichment_df_use[[pvalue_col]])
    size_col <- "log2_pvalue"
  }
  
  # Adjust colour based on p-value if colour_sig is TRUE
  if (colour_sig) {
    enrichment_df_use[['fill_colour']] <- ifelse(enrichment_df_use[[pvalue_col]] > pval_cutoff, NA, enrichment_df_use[[colour_col]])
    enrichment_df_use[['fill_colour']] <- as.numeric(enrichment_df_use[['fill_colour']])
    fill_aes <- "fill_colour"
  } else {
    fill_aes <- colour_col
  }
  if(is.null(colour_legend_title)){
    colour_legend_use=colour_col
  }else{
    colour_legend_use=colour_legend_title
  }
  
  #check whether label is discrete or continuous
  class_col=ifelse(class(enrichment_df_use[[colour_col]])=='character','disc','cont')
  if(class_col=='disc'){
    enrichment_df_use[[colour_col]]=tools::toTitleCase(enrichment_df_use[[colour_col]])
    class_col_use= scale_fill_manual(name=colour_legend_title,values = c("Down" = "#1E90FF", "Up" = "#FF4500"))
  }else{
    class_col_use=scale_fill_gradient2(name=colour_legend_title,low = "#1E90FF", mid = "white", high = '#FF4500', midpoint = 0, na.value = "transparent")
  }
  
  # Add significance labels if label_sig is TRUE
  if (label_sig) {
    enrichment_df_use[['sig']] <- ifelse(enrichment_df_use[[pvalue_col]] < pval_cutoff, '*', '')
  }
  
  if(!is.null(top_x)){
    enrichment_df_use=enrichment_df_use%>%dplyr::arrange(dplyr::across(dplyr::all_of(pvalue_col)))
    if(top_x<length(unique(enrichment_df_use[[y_axis_col]]))){
      enrichment_df_use_id=unique(enrichment_df_use[[y_axis_col]])[1:top_x]
    }else{
      enrichment_df_use_id=unique(enrichment_df_use[[y_axis_col]])
    }
    
    enrichment_df_use=enrichment_df_use[enrichment_df_use[[y_axis_col]]%in%
                                          enrichment_df_use_id,]
  }
  
  if(is.null(size_col)){
    legend_name=expression('-log'[2]*'(p-value)')
  }else{
    if(is.null(size_legend_title)){
      legend_name=size_col
    }else{
      legend_name=size_legend_title
    }
  }
  
  if(!is.null(title_add)&!is.null(plot_title)){
    plot_title=paste0(plot_title,title_sep,title_add)
  }
  
  if (!is.null(space_indent)){
    enrichment_df_use[[y_axis_col]]=as.character(enrichment_df_use[[y_axis_col]])
    for (j in 1:nrow(enrichment_df_use)){
      enrichment_df_use[[y_axis_col]][j]=paste(stringi::stri_wrap(enrichment_df_use[[y_axis_col]][j], space_indent, 0.0), collapse="\n")
    }
  }
  
  if(order_pval){
    enrichment_df_use=enrichment_df_use[order(enrichment_df_use[[pvalue_col]]),]
    enrichment_df_use[[y_axis_col]]=factor(enrichment_df_use[[y_axis_col]],
                                           levels=rev(unique(enrichment_df_use[[y_axis_col]])))
  }
  
  ## ---- NEW PART: angle-dependent hjust & vjust ----
  if (tilt_x == 0) {
    x_hjust <- 0.5
    x_vjust <- 0.5
  } else if (tilt_x == 45) {
    x_hjust <- 1
    x_vjust <- 1
  } else if (tilt_x == 90) {
    x_hjust <- 1
    x_vjust <- 0.5
  } else {
    # sensible default for other angles
    x_hjust <- 1
    x_vjust <- 1
  }
  ## -----------------------------------------------  
  
  # Create the ggplot object
  dotplot <- ggplot(enrichment_df_use, aes_string(x = x_axis_col, 
                                                  y = y_axis_col, 
                                                  size = size_col, 
                                                  fill = fill_aes)) +
    geom_point(alpha = 0.6, colour = 'black', pch = 21) +
    class_col_use +
    scale_size_continuous(range = c(2, 10), name = legend_name) +
    theme(axis.title.x = element_text(size = text_size),
          axis.title.y = element_text(size = text_size),
          plot.title = element_text(size=text_size+2),
          axis.text.x = element_text(size = text_size, 
                                     angle = tilt_x, 
                                     hjust = x_hjust, 
                                     vjust = x_vjust),
          axis.text.y = element_text(size = text_size),
          legend.title = element_text(size = text_size),
          legend.text = element_text(size = text_size),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey90"),
          panel.grid.minor = element_blank()) +
    labs(x = xlab, y = ylab)+
    ggtitle(plot_title)
  
  # Add significance labels if label_sig is TRUE
  if (label_sig) {
    dotplot <- dotplot +
      geom_text(aes(label = sig), size = text_size-6, nudge_y = nudge_sig)
  }
  
  if(!is.null(facet)){
    dotplot <- dotplot +
      facet_wrap(~facet)
  }
  
  if(!(col_legend)){
    dotplot=dotplot+
      guides(fill = "none")
  }
  
  return(dotplot)
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

#internal for enrich_things to calculate fractions
calc_fractions=function(fraction){
  temp_frac=unlist(strsplit(fraction,split='/'))
  temp_frac=as.numeric(temp_frac)
  temp_frac[1]/temp_frac[2]
}

#calc_fractions wrapper
calc_fractions_wrapper=function(fractions){
  unlist(lapply(fractions,calc_fractions))
}

#function to enrich a dataframe based on a group using enricher()
enricher_degs=function(deg_df,
                     deg_gene_col,
                     deg_group_col=NULL,
                     sig_col='sig',
                     enrichment_df,
                     enrichment_gene_col,
                     enrichment_term_col,
                     enrichment_group_col=NULL,
                     background=NULL,
                     pAdjustMethod='BH',
                     minGSSize=0,
                     maxGSSize=Inf,
                     enrichment_sig=0.01,
                     gene_dictionary = NULL,
                     to_replace_gene_col='ensembl_gene_id',
                     gene_symbol_col = 'external_gene_name',
                     carry_col=NULL,
                     verbose=TRUE){
  deg_df[[deg_gene_col]]=as.character(deg_df[[deg_gene_col]])
  if(is.null(deg_group_col)){
    deg_df[['group']]='sig'
    deg_group_col='group'
  }
  if(is.null(enrichment_group_col)){
    enrichment_df[['enrichment_group']]='enrichment'
    enrichment_group_col='enrichment_group'
  }
  
  if(is.null(background)){
    background=unique(deg_df[[deg_gene_col]])
  }else{
    background=unique(background)
  }
  background=as.character(background)
  
  deg_df_sig=deg_df[deg_df[[sig_col]],]
  all_group_deg=unique(deg_df[[deg_group_col]])
  all_enrich_group=unique(enrichment_df[[enrichment_group_col]])
  all_enrich_group=all_enrich_group[all_enrich_group!='']
  all_enriched=do.call('rbind',lapply(all_group_deg,function(per_group){
    if(verbose){
      message(per_group)
    }
    temp_deg=deg_df_sig[deg_df_sig[[deg_group_col]]==per_group,]
    
    enriched_groups=do.call('rbind',lapply(all_enrich_group,function(per_enrich){
      if(verbose){
        message(per_enrich)
      }
      temp_enrich_df=enrichment_df[enrichment_df[[enrichment_group_col]]==per_enrich,]
      
      temp_enrich_df=temp_enrich_df%>%dplyr::select(!!sym(enrichment_term_col), !!sym(enrichment_gene_col))
      test_enrich=enricher(gene = unique(temp_deg[[deg_gene_col]]),
                           universe = background,
                           TERM2GENE = temp_enrich_df,
                           pvalueCutoff = 1,
                           pAdjustMethod = pAdjustMethod,
                           minGSSize = minGSSize,
                           maxGSSize = maxGSSize,
                           qvalueCutoff = 1)
      if(is.null(test_enrich)){
        return(NULL)
      }
      test_enrich=data.frame(test_enrich)
      test_enrich$enrich_class=per_enrich
      rownames(test_enrich)=NULL
      return(test_enrich)
    }))
    if(is.null(enriched_groups)){
      return(NULL)
    }
    enriched_groups$group=per_group
    
    if(!is.null(carry_col)){
      for(i in carry_col){
        enriched_groups[[i]]=temp_deg[[i]][1]
      }
    }
    return(enriched_groups)
  }))
  
  all_enriched$sig=ifelse(all_enriched$p.adjust<enrichment_sig,TRUE,FALSE)
  
  all_enriched$gene_ratio=calc_fractions_wrapper(all_enriched$GeneRatio)
  
  if(!is.null(gene_dictionary)){
    message('Adding gene IDs')
    all_enriched=enrichment_gene_symbol(enrichment_object = all_enriched,
                                        entrez_ID = 'geneID',
                                        entrez_dictionary = gene_dictionary,
                                        entrez_dictionary_entrez_col = to_replace_gene_col,
                                        entrez_dictionary_gene_symbol_col = gene_symbol_col)
  }
  
  return(all_enriched)
}

plot_dotplot_enrichment_wrapper=function(enrichment_df,
                                         x_axis_col = "group",
                                         y_axis_col = "ID",
                                         size_col = NULL, # This will be calculated from pvalue_col if NULL
                                         pvalue_col = "p.adjust", # Specify the column for p-values
                                         colour_col = "enrichmentScore",
                                         label_sig = FALSE,
                                         tilt_x = 90,
                                         filter_sig = TRUE,
                                         pval_cutoff = 0.01,
                                         colour_sig = FALSE,
                                         nudge_sig = 0,
                                         text_size = 10,
                                         str_remove=NULL,
                                         top_x=30,
                                         size_legend_title=NULL,
                                         colour_legend_title=NULL,
                                         enrich_class_col='namespace_1003',
                                         title_add=NULL,
                                         title_sep='_',
                                         xlab=NULL,
                                         ylab=NULL,
                                         facet=NULL,
                                         col_legend=TRUE,
                                         order_pval=TRUE){
  require(stringi)
  all_class=unique(enrichment_df[[enrich_class_col]])
  all_plots=lapply(all_class,function(per_class){
    message(per_class)
    enrichment_df_temp=enrichment_df[enrichment_df[[enrich_class_col]]==per_class,]
    p=plot_dotplot_enrichment(enrichment_df = enrichment_df_temp,
                            x_axis_col = x_axis_col,
                            y_axis_col = y_axis_col,
                            size_col = size_col, 
                            pvalue_col = pvalue_col,
                            colour_col = colour_col,
                            label_sig = label_sig,
                            tilt_x = tilt_x,
                            filter_sig = filter_sig,
                            pval_cutoff = pval_cutoff,
                            colour_sig = colour_sig,
                            nudge_sig = nudge_sig,
                            text_size = text_size,
                            str_remove=str_remove,
                            top_x=top_x,
                            size_legend_title=size_legend_title,
                            colour_legend_title=colour_legend_title,
                            plot_title = per_class,
                            title_add=title_add,
                            title_sep=title_sep,
                            xlab = xlab,
                            ylab = ylab,
                            facet=facet,
                            col_legend=col_legend,
                            order_pval=order_pval)
  })
  names(all_plots)=all_class
  return(all_plots)
}

#Function to get lengths for fpkm normalisation
summarise_gene_lengths <- function(data, method = c("mean",
                                                    'median',
                                                    'max')) {
  # Validate method input
  if (!method %in% c("mean", "median", "max")) {
    stop("Invalid method. Choose one of 'mean', 'median', or 'max'.")
  }
  
  # Group by gene and summarise based on the chosen method
  summarised_data <- data %>%
    group_by(ensembl_gene_id) %>%
    summarise(
      transcript_length = case_when(
        method == "mean" ~ mean(transcript_length, na.rm = TRUE),
        method == "median" ~ median(transcript_length, na.rm = TRUE),
        method == "max" ~ max(transcript_length, na.rm = TRUE)
      )
    ) %>%
    ungroup()
  
  return(summarised_data)
}

#function to add fpkm assay to a deseq2 object given a gene_length_df
add_fpkm_assay <- function(deseq2_obj, gene_length_df, gene_col = "gene", length_col = "length") {
  # Check if the DESeq2 object has the required counts
  if (!"counts" %in% assayNames(deseq2_obj)) {
    stop("The DESeq2 object must have a 'counts' assay.")
  }
  
  # Extract counts matrix
  counts <- counts(deseq2_obj)
  
  # Ensure gene_length_df has the required columns
  if (!all(c(gene_col, length_col) %in% colnames(gene_length_df))) {
    stop(paste("gene_length_df must contain", gene_col, "and", length_col, "columns."))
  }
  
  # Ensure genes in gene_length_df and counts match
  gene_lengths <- gene_length_df[[length_col]]
  names(gene_lengths) <- gene_length_df[[gene_col]]
  
  if (!all(rownames(counts) %in% names(gene_lengths))) {
    stop("All genes in the DESeq2 object must have corresponding lengths in gene_length_df.")
  }
  
  # Subset gene lengths to match counts rownames
  gene_lengths <- gene_lengths[rownames(counts)]
  
  # Compute FPKM values
  # Formula: FPKM = (counts / (gene length in kb * total reads)) * 10^6
  gene_lengths_kb <- gene_lengths / 1000
  total_counts <- colSums(counts)
  
  fpkm <- sweep(counts, 2, total_counts, FUN = "/") # Normalize by library size
  fpkm <- sweep(fpkm, 1, gene_lengths_kb, FUN = "/") # Normalize by gene length
  fpkm <- fpkm * 1e6
  
  # Add FPKM as a new assay
  assays(deseq2_obj)[["fpkm"]] <- fpkm
  
  return(deseq2_obj)
}

#wrapper for add_fpkm_assay and add_tpm_assay to add fpkm/tpm to a deseq obj
calc_normalised_assay <- function(deseq2_obj,
                                  gene_length_df,
                                  gene_col = "ensembl_gene_id",
                                  length_col = "transcript_length",
                                  method = 'tpm', # Default to TPM
                                  summarisation_method = 'mean') {
  # Ensure gene lengths are available in the gene_length_df
  gene_length_df <- gene_length_df[gene_length_df[[gene_col]] %in% rownames(deseq2_obj), ]
  
  # Summarise gene lengths (if needed)
  gene_length_df_mean <- summarise_gene_lengths(data = gene_length_df,
                                                method = summarisation_method)
  
  # Filter DESeq2 object for relevant genes
  message('Before filtering DESeq2 gene n: ', nrow(deseq2_obj))
  deseq2_obj_match_lengths <- deseq2_obj[rownames(deseq2_obj) %in% 
                                           gene_length_df_mean[[gene_col]]]
  message('After filtering DESeq2 gene n: ', nrow(deseq2_obj_match_lengths))
  
  # Arrange gene length dataframe to match DESeq2 object
  gene_length_df_mean <- gene_length_df_mean %>% dplyr::arrange(rownames(deseq2_obj_match_lengths))
  
  # Call the appropriate function based on the method argument
  if (method == 'tpm') {
    result_obj <- add_tpm_assay(deseq2_obj = deseq2_obj_match_lengths,
                                gene_length_df = gene_length_df_mean,
                                gene_col = gene_col,
                                length_col = length_col)
  } else if (method == 'fpkm') {
    result_obj <- add_fpkm_assay(deseq2_obj = deseq2_obj_match_lengths,
                                 gene_length_df = gene_length_df_mean,
                                 gene_col = gene_col,
                                 length_col = length_col)
  } else {
    stop("Invalid method specified. Use 'tpm' or 'fpkm'.")
  }
  
  return(result_obj)
}

add_tpm_assay <- function(deseq2_obj,
                          gene_length_df,
                          gene_col = "gene",
                          length_col = "length") {
  # Check if the DESeq2 object has the required counts
  if (!"counts" %in% assayNames(deseq2_obj)) {
    stop("The DESeq2 object must have a 'counts' assay.")
  }
  
  # Extract counts matrix
  counts <- counts(deseq2_obj)
  
  # Ensure gene_length_df has the required columns
  if (!all(c(gene_col, length_col) %in% colnames(gene_length_df))) {
    stop(paste("gene_length_df must contain", gene_col, "and", length_col, "columns."))
  }
  
  # Ensure genes in gene_length_df and counts match
  gene_lengths <- gene_length_df[[length_col]]
  names(gene_lengths) <- gene_length_df[[gene_col]]
  
  if (!all(rownames(counts) %in% names(gene_lengths))) {
    stop("All genes in the DESeq2 object must have corresponding lengths in gene_length_df.")
  }
  
  # Subset gene lengths to match counts rownames
  gene_lengths <- gene_lengths[rownames(counts)]
  
  # Compute TPM values
  # Formula: TPM = (counts / gene length in kb) / scaling factor
  gene_lengths_kb <- gene_lengths / 1000
  norm_counts <- sweep(counts, 1, gene_lengths_kb, FUN = "/") # Normalize by gene length
  
  scaling_factors <- colSums(norm_counts) # Sum normalised counts across genes
  tpm <- sweep(norm_counts, 2, scaling_factors, FUN = "/") # Normalize by scaling factor
  tpm <- tpm * 1e6 # Scale to 1,000,000
  
  # Add TPM as a new assay
  assays(deseq2_obj)[["tpm"]] <- tpm
  
  return(deseq2_obj)
}

# Function to make stacked barcharts out of different conditions
stacked_barchart <- function(
    df,
    group,
    fill,
    cutoff = 100,
    percent_cutoff = 0,
    fill_order = NULL,
    group_order = NULL,
    label_type = "both",
    legend_position = "right",
    text_size = 12,
    ylab = NULL,           # controls VERTICAL axis label (post-flip)
    legend_title = NULL,
    n_split=' ',
    show_n = TRUE          # <---- NEW ARGUMENT
){
  group_sym <- rlang::sym(group)
  fill_sym  <- rlang::sym(fill)
  
  df[[fill]]  <- as.character(df[[fill]])
  df[[group]] <- as.character(df[[group]])
  
  # total N per group (for axis labels)
  group_sizes <- df %>%
    dplyr::count(!!group_sym, name = "total_n")
  
  df_percent <- df %>%
    dplyr::count(!!group_sym, !!fill_sym) %>%
    dplyr::group_by(!!group_sym) %>%
    dplyr::mutate(percent = n / sum(n) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      label = dplyr::case_when(
        n > .env$cutoff & percent > .env$percent_cutoff & .env$label_type == "both"    ~ paste0(n, n_split,"(", round(percent, 1), "%)"),
        n > .env$cutoff & percent > .env$percent_cutoff & .env$label_type == "number"  ~ as.character(n),
        n > .env$cutoff & percent > .env$percent_cutoff & .env$label_type == "percent" ~ paste0(round(percent, 1), "%"),
        TRUE ~ NA_character_
      )
    )
  
  if (!is.null(fill_order)) {
    df_percent[[fill]] <- factor(df_percent[[fill]], levels = fill_order)
  }
  
  if (!is.null(group_order)) {
    df_percent[[group]] <- factor(df_percent[[group]], levels = rev(group_order))
  }
  
  # base plot
  p <- ggplot(df_percent, aes(x = !!group_sym, y = percent, fill = !!fill_sym)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5),
              size = text_size * 0.4, na.rm = TRUE) +
    coord_flip() +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(c(0, 0)))
  
  # axis labels with N if requested
  if (show_n) {
    # named vector: names = original group values, values = "group (N)"
    label_map <- setNames(
      paste0(group_sizes[[group]], " (", group_sizes$total_n, ")"),
      group_sizes[[group]]
    )
    
    p <- p + scale_x_discrete(labels = label_map)
  } else {
    p <- p + scale_x_discrete()
  }
  
  p +
    labs(
      x    = if (is.null(ylab)) group else ylab,  # vertical axis label (after flip)
      y    = "Percentage",                        # horizontal axis label (after flip)
      fill = if (is.null(legend_title)) fill else legend_title
    ) +
    theme(
      legend.position = legend_position,
      legend.text  = element_text(size = text_size * 0.9),
      legend.title = element_text(size = text_size),
      axis.title.x = element_text(size = text_size),   # horizontal title (Percentage)
      axis.title.y = element_text(size = text_size),   # vertical title (ylab/group)
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(size = text_size * 0.9),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

#build enrichment rds files
build_enrichment_rds=function(df,
                              gene_col,
                              group_col,
                              output_file,
                              output_dir,
                              extra_col=TRUE){
  enrich_me=df%>%dplyr::select(dplyr::all_of(c(gene_col,group_col)))
  colnames(enrich_me)=c('gene','group')
  if(extra_col&ncol(df)>2){
    enrich_me=cbind(enrich_me,
                    df[,!colnames(df)%in%c(gene_col,
                                           group_col)])
  }
  saveRDS(object = enrich_me,
          file = paste0(output_dir,'/',output_file,'.rds'))
}

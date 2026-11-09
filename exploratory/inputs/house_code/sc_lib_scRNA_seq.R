#scRNAseq functions
#functions
customFeaturePlot <- function(seur.obj, alpha.use = 0.2, gradient.use = c("#adadad",viridis(100, direction = 1, end=0.99, option = "C")),
                              expression.threshold = 0, feature = NULL, pt.size = 0.5, reduction="umap", dims = c(1,2),
                              xlimset = NULL,ylimset=NULL,name=NULL) {
  require(ggplot2)
  require(viridis)
  
  coords <- data.frame(seur.obj[[reduction]]@cell.embeddings)[,dims] ## in case reduction is pca select only first 2
  colnames(coords) = c("X","Y")
  coords$Cell = rownames(coords)
  if (feature %in% rownames(seur.obj@assays$SCT@data)){
    coords$Feature = as.numeric(seur.obj@assays$SCT@data[feature,])
  } else{
    stop("Feature not found!")
  }
  coords$Cell = factor(coords$Cell,levels = coords$Cell[order(coords$Feature)]) ## try to plot first cells without expression
  coords$Feature = coords$Feature*10/max(coords$Feature)
  coords$Transparent = factor(ifelse(coords$Feature <= expression.threshold,"transparent","full"),levels=c("transparent","full"))
  
  if (is.null(xlimset)){
    xlimset = c(min(coords$X),max(coords$X))
  }
  if (is.null(ylimset)){
    ylimset = c(min(coords$Y),max(coords$Y))
  }
  
  if(is.null(name)){
    name=element_blank()
  }
  p <- ggplot(coords,aes(x=X,y=Y)) +
    geom_point(aes(colour=Feature,alpha = Transparent), size=pt.size, shape=16) +
    
    
    scale_colour_gradientn(colours = gradient.use,  
                           guide = guide_colourbar(direction = "vertical",
                                                   barwidth = unit(0.2, "in"),
                                                   title = name)) +
    #label.theme = element_text(size=16,family='Helvetica',hjust=0,vjust=0.5))) +
    scale_alpha_manual(values=c("transparent"=alpha.use,"full"=1),guide="none") +
    xlim(xlimset) + ylim(ylimset) +
    xlab(paste(reduction,1,sep="_")) + ylab(paste(reduction,2,sep="_")) +
    theme_classic() +
    ggtitle(feature) +
    theme(text=element_text(size=20,family="Helvetica"),
          plot.title = element_text(size=20,family="Helvetica"),
          legend.position="right",
          legend.text.align = 0,
          legend.text = element_text(size=14,family='Helvetica'),
          legend.justification = "left")
  
  
  ####
  return(p)
  
}

customSignaturePlot <- function(seur.obj,
                                alpha.use = 0.2,
                                grey_zero = FALSE,
                                expression.threshold = 0,
                                feature = NULL,
                                pt.size = 0.5,
                                reduction = "umap",
                                dims = c(1,2),
                                xlimset = NULL,
                                ylimset = NULL,
                                name = NULL,
                                remove_axes = FALSE) {
  
  require(ggplot2)
  require(viridis)
  
  if (grey_zero) {
    gradient.use <- c("#adadad",
                      viridis(100, direction = 1, end = 0.99, option = "C"))
  } else {
    gradient.use <- viridis(100, direction = 1, end = 0.99, option = "C")
  }
  
  coords <- data.frame(seur.obj[[reduction]]@cell.embeddings)[, dims]
  colnames(coords) <- c("X", "Y")
  coords$Cell <- rownames(coords)
  
  if (feature %in% colnames(seur.obj@meta.data)) {
    coords$Feature <- as.numeric(seur.obj@meta.data[, feature])
  } else {
    stop("Feature not found!")
  }
  
  coords$Cell <- factor(coords$Cell,
                        levels = coords$Cell[order(coords$Feature)])
  
  coords$Feature <- coords$Feature - min(coords$Feature)
  coords$Feature <- coords$Feature * 10 / max(coords$Feature)
  coords$Transparent <- factor(
    ifelse(coords$Feature <= expression.threshold,
           "transparent", "full"),
    levels = c("transparent", "full")
  )
  
  if (is.null(xlimset)) xlimset <- range(coords$X)
  if (is.null(ylimset)) ylimset <- range(coords$Y)
  if (is.null(name)) name <- element_blank()
  
  p <- ggplot(coords, aes(x = X, y = Y)) +
    geom_point(aes(colour = Feature, alpha = Transparent),
               size = pt.size, shape = 16) +
    scale_colour_gradientn(
      colours = gradient.use,
      guide = guide_colourbar(
        direction = "vertical",
        barwidth = unit(0.2, "in"),
        title = name
      )
    ) +
    scale_alpha_manual(
      values = c("transparent" = alpha.use, "full" = 1),
      guide = "none"
    ) +
    xlim(xlimset) +
    ylim(ylimset) +
    xlab(paste(reduction, 1, sep = "_")) +
    ylab(paste(reduction, 2, sep = "_")) +
    ggtitle(feature) +
    theme_classic() +
    theme(
      text = element_text(size = 20, family = "Helvetica"),
      plot.title = element_text(size = 20, family = "Helvetica"),
      legend.position = "right",
      legend.text.align = 0,
      legend.text = element_text(size = 14, family = "Helvetica"),
      legend.justification = "left"
    )
  
  ## ---- Remove all axes if requested ----
  if (remove_axes) {
    p <- p + theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    )
  }
  
  return(p)
}

#wrapper for customfeatureplot and customsignatureplot
custom_plot_wrapper <- function(seurat_obj,
                                wrap = c('customFeaturePlot',
                                         'customSignaturePlot'),
                                goi,
                                ncol = NULL,  # Number of plots per row
                                ...) {
  
  wrap <- match.arg(wrap)  # Ensure valid selection of wrap function
  
  if (wrap == 'customFeaturePlot') {
    all_plots <- lapply(goi, function(per_goi) {
      customFeaturePlot(seur.obj = seurat_obj,
                        feature = per_goi,
                        name = 'Relative Expression',
                        ...)
    })
  } else if (wrap == 'customSignaturePlot') {
    all_plots <- lapply(goi, function(per_goi) {
      customSignaturePlot(seur.obj = seurat_obj,
                          feature = per_goi,
                          name = 'Relative Expression',
                          ...)
    })
  }
  
  patchwork::wrap_plots(all_plots, ncol = ncol) + 
    patchwork::plot_layout(guides = "collect")
}

# function to wilcox scores in metadata
wilcox_scores=function(df,
                       group_col,
                       score_col){
  df[['group_col']]=df[[group_col]]
  df[['score_col']]=df[[score_col]]
  wilcox_group=df%>%wilcox_test(score_col~group_col)%>%add_xy_position()
  wilcox_group[['group']]=group_col
  wilcox_group[['score']]=score_col
  
  if(length(unique(df[[group_col]]))>2){
    kw=df%>%rstatix::kruskal_test(score_col~group_col)
    kw=kw%>%add_significance()
    wilcox_group$kw=kw$p
    wilcox_group$kw_signif=kw$p.signif
  }
  return(wilcox_group)
}

#wrapper for wilcox_scores
wilcox_scores_by_group=function(df,
                                group_col = 'cluster',
                                score_col = 'entropy',
                                subset_col='pathway'){
  all_subset=unique(df[[subset_col]])
  wilcox_all=do.call('rbind',lapply(all_subset,function(per_subset){
    df_temp=df[df[[subset_col]]==per_subset,]
    wilcox_temp=wilcox_scores(df = df_temp,
                              group_col = group_col,
                              score_col = score_col)
    wilcox_temp[['pathway']]=per_subset
    return(wilcox_temp)
  }))
  wilcox_all[['p.adj']]=p.adjust(wilcox_all[['p']],method='BH')
  wilcox_all=wilcox_all%>%add_significance(p.col = 'p.adj',
                                           output.col = 'p.adj.signif')
  return(wilcox_all)
}

#wrapper for wilcox_scores
wilcox_scores_wrapper=function(df,
                               group_col,
                               score_cols){
  all_scores=do.call('rbind',lapply(score_cols,function(per_score){
    wilcox_scores(df = df,
                  group_col = group_col,
                  score_col = per_score)
  }))
}

#function to make dotplot of scores
dotplot_scores=function(df,
                        group_col,
                        score_cols,
                        xlab=NULL,
                        ylab=NULL,
                        kw_scores=NULL){
  df[['group_col']]=df[[group_col]]
  df=df[,colnames(df)%in%c('group_col',score_cols)]
  df_melted=reshape2::melt(df)
  
  if(is.null(kw_scores)){
    p=df_melted%>%ggplot(aes(x=group_col,y=value))+
      geom_point(aes(colour=group_col),alpha=0.3,position = 'jitter',show.legend = FALSE)+
      geom_violin(aes(fill=group_col),show.legend = FALSE)+
      xlab(xlab)+
      ylab(ylab)+
      facet_wrap(~ variable)
  }else{
    df_melted=
      merge(df_melted,kw_scores,by.x='variable',by.y='score')
    df_melted$label=paste0(df_melted$variable,' (kw: ',df_melted$kw_signif,')')
    
    p=
      df_melted%>%ggplot(aes(x=group_col,y=value))+
      geom_point(aes(colour=group_col),alpha=0.3,position = 'jitter',show.legend = FALSE)+
      geom_violin(aes(fill=group_col),show.legend = FALSE)+
      xlab(xlab)+
      ylab(ylab)+
      facet_wrap(~ label)
  }
  
  
  return(p)
}

#wrapper combining wilcox_scores_wrapper and dotplot_scores
wilcox_dotplot_wrapper=function(df,
                                group_col,
                                score_cols,
                                xlab=NULL,
                                ylab=NULL){
  df_melt=reshape2::melt(df[,colnames(df)%in%c(group_col,score_cols)])
  colnames(df_melt)[1]='group'
  
  #get kw
  kw=do.call('rbind',lapply(score_cols,function(per_col){
    df_melt_filtered=df_melt[df_melt$variable==per_col,]
    df_melt_filtered=df_melt_filtered%>%
      kruskal_test(value~group)
    
    df_melt_filtered$variable=per_col
    return(df_melt_filtered)
  }))
  
  kw=kw%>%adjust_pvalue(p.col = 'p')%>%add_significance(p.col = 'p.adj')
  kw$var_new=paste0(kw$variable,' (kw: ',kw$p.adj.signif,')')
  kw_merge=kw[,colnames(kw)%in%c('variable','var_new')]
  df_melt=merge(df_melt,kw_merge)
  
  score_cols_new=unique(df_melt$var_new)
  all_wilcox=do.call('rbind',lapply(score_cols_new,function(per_col){
    df_melt_filtered=df_melt[df_melt$var_new==per_col,]
    df_melt_filtered=df_melt_filtered%>%
      wilcox_test(value~group)%>%
      add_xy_position(scales = 'free_y')
    
    df_melt_filtered$variable=per_col
    return(df_melt_filtered)
  }))
  
  all_wilcox=all_wilcox%>%adjust_pvalue(p.col = 'p',
                                        output.col = 'p.adj')%>%
    add_significance(p.col = 'p.adj',output.col = 'p.adj.signif')
  
  colnames(all_wilcox)[colnames(all_wilcox)=='variable']='var_new'
  
  p=df_melt%>%ggplot(aes(x=group,y=value))+
    geom_point(aes(colour=group),alpha=0.3,position = 'jitter',show.legend = FALSE)+
    geom_violin(aes(fill=group),show.legend = FALSE)+
    xlab(xlab)+
    ylab(ylab)+
    facet_wrap(~ var_new,scales='free_y')+
    stat_pvalue_manual(all_wilcox,hide.ns = TRUE,tip.length = 0,
                       step.group.by = 'var_new',step.increase = 0.03)
  
  results=c()
  results[['p']]=p
  results[['wilcox']]=all_wilcox
  results[['kw']]=kw
  return(results)
}

#function to plot the fc of set of genes
plot_goi_fc=function(markers,
                     goi,
                     x_col='cluster',
                     pval_col='p_val_adj',
                     legend_title=NULL,
                     xlab=NULL,
                     ylab=NULL){
  temp_markers=markers[markers[['gene']]%in%goi,]
  temp_markers$fill = ifelse(temp_markers[[pval_col]] > 0.05,
                             'Not Sig', ifelse(temp_markers$avg_log2FC > 0, 'Sig Overexpressed',
                                               'Sig Underexpressed'))
  
  temp_markers$fill = factor(temp_markers$fill,
                             levels = c('Sig Overexpressed', 'Sig Underexpressed', 'Not Sig'))
  fill_manual = c("Not Sig" = "#f1f1f1", "Sig Overexpressed" = "#f8766d",
                  "Sig Underexpressed" = "#00bfc4")
  temp_markers=temp_markers%>%rstatix::add_significance(p.col = 'p_val_adj')
  temp_markers[['x']]=temp_markers[[x_col]]
  temp_markers%>%ggplot(aes(x=x,y=gene))+
    geom_tile(aes(fill=fill),colour='black')+
    scale_x_discrete(expand = expansion(c(0, 0))) +
    scale_y_discrete(expand = expansion(c(0, 0))) +
    ylab(ylab) +
    xlab(xlab)+
    scale_fill_manual(name = legend_title,
                      values = fill_manual)
}

#function to subsample celltypes from a sc obj
subset_sce_by_celltype <-
  function(sce,
           cell_type_column,
           x = NULL,
           include_types = NULL,
           exclude_types = NULL,
           allow_all_cells = FALSE,
           seed = 1) {
    # Ensure the cell_type_column exists in colData
    if (!(cell_type_column %in% colnames(colData(sce)))) {
      stop("The specified cell_type_column does not exist in colData.")
    }
    
    # Set seed for reproducibility
    set.seed(seed)
    
    # Extract the cell type information
    cell_types <- colData(sce)[[cell_type_column]]
    
    # Filter cell types based on include and exclude arguments
    if (!is.null(include_types)) {
      cell_types <- cell_types[cell_types %in% include_types]
    }
    if (!is.null(exclude_types)) {
      cell_types <- cell_types[!(cell_types %in% exclude_types)]
    }
    
    # Determine the minimum cell count if x is not provided
    cell_type_counts <- table(cell_types)
    if (is.null(x)) {
      x <- min(cell_type_counts)
    }
    
    # Identify cell types with fewer cells than x
    insufficient_cells <- names(cell_type_counts[cell_type_counts < x])
    
    # Check if x exceeds the minimum cell type count
    if (x > min(cell_type_counts) && !allow_all_cells) {
      stop(paste("The value of 'x' exceeds the minimum cell count of a specific cell type.",
                 "Set 'allow_all_cells = TRUE' to include all cells for those types.",
                 "Cell types with fewer cells than 'x':", paste(insufficient_cells, collapse = ", ")))
    }
    
    # Inform about cell types with fewer cells if allow_all_cells is TRUE
    if (allow_all_cells && length(insufficient_cells) > 0) {
      message(paste("Cell types with fewer cells than 'x' (including all cells for these types):",
                    paste(insufficient_cells, collapse = ", ")))
    }
    
    # Create a list to store indices for selected cells
    selected_indices <- vector("list", length(unique(cell_types)))
    names(selected_indices) <- unique(cell_types)
    
    # Iterate through each cell type and randomly sample x cells
    for (type in unique(cell_types)) {
      type_indices <- which(colData(sce)[[cell_type_column]] == type)
      if (length(type_indices) <= x || (allow_all_cells && length(type_indices) < x)) {
        # If fewer cells are available than x or allow_all_cells is TRUE, take all
        selected_indices[[type]] <- type_indices
      } else {
        # Randomly sample x cells
        selected_indices[[type]] <- sample(type_indices, x)
      }
    }
    
    # Combine all selected indices
    final_indices <- unlist(selected_indices)
    
    # Subset the SCE object
    subset_sce <- sce[, final_indices]
    
    return(subset_sce)
  }

#function to scran normalize data
scran_by_celltypes=function(sce,
                            cell_type_col,
                            seed=123,
                            log=FALSE){
  require(scran)
  set.seed(seed)  # For reproducibility
  clusters <- sce[[cell_type_col]]  # Use cell types for clustering
  sce_temp <- computeSumFactors(sce, clusters = clusters)
  sce_temp <- normalizeCounts(sce_temp, log = log)
  heart_sc@assays@data[[paste0('sce_log_',log)]]=sce_temp
  return(heart_sc)
}

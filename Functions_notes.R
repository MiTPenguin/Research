##Source script for commonly used functions, presets, and other notes for general R programming.

install_if_missing <- function(packages) {
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
  }
}

Axis_themes <- theme(plot.title = element_text(size = 8),
                     axis.title = element_text(size = 8), 
                     axis.text = element_text(size = 6),
                     axis.text.x = element_text(size = 6),
                     legend.text = element_text(size =6),
                     legend.title = element_text(size = 8),
                     strip.text.x = element_text(size = 8))

#theme(axis.text.x = element_text(angle = 90, hjust = 1))

Read_10x_data <- function(Sample_list, data_dir){
#basically take a lits of samples, then use a baes data_dir to pull out the barocdes, genes, and matrix files
#remember to gunzip everything. the Read10X function is not really quite that smart.
    d10x.data <- sapply(Sample_list, function(i){
      d10x <- Read10X(file.path(data_dir,i))
      colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-") 
    #this line is a bit tougher to understand
      d10x
    })
return(d10x.data)
}

find_pcCut <- function(seurat){
  #pull out the pca deviation from each PCA components
  PCA_dev <- seurat@reductions$pca@stdev
  #This part is doing some geometry calculation. we should write this into a function.
  allCoor<-cbind(seq_along(PCA_dev),PCA_dev)
  line_Vec<-allCoor[length(PCA_dev),]-allCoor[1,]
  ## normalize the line vector
  lineVecN = line_Vec / sqrt(sum(line_Vec^2));
  ## find the distance from each point to the line:
  ## vector between all points and first point
  vecFromFirst<-allCoor-do.call("rbind", rep(list(allCoor[1,]), length(PCA_dev)))
  q<-do.call("rbind", rep(list(lineVecN), length(PCA_dev)))
  scalarProduct<-q[,1]
  for (i in 1:length(PCA_dev)) {
    scalarProduct[i]<- vecFromFirst[i,] %*% q[i,]
  }
  vecFromFirstParallel = scalarProduct * q
  vecToLine = vecFromFirst - vecFromFirstParallel
  distToLine = sqrt(rowSums(vecToLine^2))
  ##Point in elbow is point furthest from line connecting first and last point
  pcCut<-which.max(distToLine) 
  return(pcCut)
}


limit_vargenes <- function(seurat,n_var_genes = 1000, cf = 0.2){
  length_var <- length(seurat@var.genes)
  while (length_var>n_var_genes) {
    seurat <- FindVariableGenes(object = seurat, x.low.cutoff = 0.1, y.cutoff = cf, do.plot = FALSE)
    length_var <- length(seurat@var.genes)
    cf=cf+.05
  }
  return(seurat)
}

#process assumes that we have already made a seurat object to push in (in most cases, one where some cells
#were just removed, and just needs to be re-processed)

Calculate_mito <- function(Seurat, Species = 'mouse'){
    regex_script <- case_when(Species == 'mouse' ~ "^mt-",
                              Species == 'human' ~ "^MT-")
    mito.genes <- grep(regex_script, rownames(GetAssayData(Seurat)), value = T)
    percent.mito <- Matrix::colSums(GetAssayData(Seurat, slot = "counts")[mito.genes, ]) / 
                    Matrix::colSums(GetAssayData(Seurat, slot = "counts"))
    Seurat <- AddMetaData(object = Seurat,
                           metadata = percent.mito,
                           col.name = "percent.mito")
    return(Seurat)
}

#update: with the new Seurat, I think find Variable genes, and other processes would have to be re-written
# Seurat_process <- function(seurat, perp = 30, cutoff = NULL, limit = FALSE,cf = 0.2, mf = 0.1,umap.run = TRUE,...) {
#   seurat <- NormalizeData(object = seurat)
#   if(limit!=FALSE){
#     seurat <- limit_vargenes(seurat,limit,cf)
#   }else{
#     seurat <-  FindVariableGenes(seurat,x.low.cutoff = mf, y.cutoff = cf, do.plot = FALSE)
#   }
#   seurat <-  ScaleData(seurat, vars.to.regress = c('nUMI', 'percent.mito'), 
#                        model.use = 'poisson', genes.use = seurat@var.genes)
#   seurat <-  RunPCA(seurat, do.print = FALSE)
#   if(is.null(cutoff)){
#     cutoff <- find_pcCut(seurat)
#   }
#   print(paste("PCA cutoff used is",cutoff,"components"))
#   seurat <-  RunTSNE(seurat, dims.use = 1:cutoff, perplexity = perp)
#   #seurat = SetAllIdent(seurat, 'expt')
#   TSNEPlot(seurat)
#   if(umap.run){
#     seurat <-  RunUMAP(seurat, dims.use = 1:cutoff)
#     DimPlot(seurat, 'umap')
#   }
#   return(seurat)
# }
Seurat_preprocess <- function(seurat, species = "human", normalize = TRUE,
  vars.to.regress = c('nCount_RNA','percent.mito'), model.use = 'linear',...) {
  seurat <- Calculate_mito(seurat, species)
  if(normalize){
  seurat <- NormalizeData(object = seurat)
  }
  seurat <- FindVariableFeatures(seurat) #currently, this just does the default findvariable feature plot
  
  seurat <-  ScaleData(seurat, vars.to.regress = vars.to.regress, model.use = model.use)
  seurat <-  RunPCA(seurat)

  cutoff <- find_pcCut(seurat)
  print(paste("PCA cutoff used is",cutoff,"components"))
return(seurat)
}

Seurat_cluster<- function(seurat, cutoff = 10, perp = 30,
  tsne.run = FALSE, umap.run = TRUE, res = 0.8,...){
  seurat <- FindNeighbors(seurat,dims = 1:cutoff)
  seurat <- FindClusters(seurat, resolution = res)
  if(tsne.run){
  seurat <-  RunTSNE(seurat, dims = 1:cutoff, perplexity = perp)
  DimPlot(seurat, 'tsne', label = TRUE)
  }
  if(umap.run){
    seurat <-  RunUMAP(seurat, dims = 1:cutoff,umap.method = 'uwot-learn')
    DimPlot(seurat, reduction = 'umap', label = TRUE)
  }
  return(seurat)
}


# Seurat_preprocess <- function(seurat, species = "human", perp = 30, 
#   vars.to.regress = c('nCount_RNA','percent.mito'), model.use = 'linear',
#   tsne.run = FALSE, umap.run = TRUE, ...) {
#   seurat <- Calculate_mito(seurat, species)
#   seurat <- NormalizeData(object = seurat)
#   seurat <- FindVariableFeatures(seurat) #currently, this just does the default find
  
#   seurat <-  ScaleData(seurat, vars.to.regress = vars.to.regress, model.use = model.use)
#   seurat <-  RunPCA(seurat)

#   cutoff <- find_pcCut(seurat)
#   print(paste("PCA cutoff used is",cutoff,"components"))
#   if(tsne.run){
#   seurat <-  RunTSNE(seurat, dims.use = 1:cutoff, perplexity = perp)
#   DimPlot(seurat, 'tsne', label = TRUE)
#   }
#   if(umap.run){
#     seurat <-  RunUMAP(seurat, dims.use = 1:cutoff,umap.method = 'uwot-learn')
#     DimPlot(seurat, 'umap', label = TRUE)
#   }
#   return(seurat)
# }

#  seurat <- FindNeighbors(Petti_Seurat,dims = 1:cutoff)
#  seurat <- FindClusters(Petti_Seurat, resolution = res)
#  if(tsne.run){
#  seurat <-  RunTSNE(seurat, dims.use = 1:cutoff, perplexity = perp)
#  DimPlot(seurat, 'tsne', label = TRUE)
#  }
#  if(umap.run){
#    seurat <-  RunUMAP(seurat, dims.use = 1:cutoff,umap.method = 'uwot-learn')
#    DimPlot(seurat, 'umap', label = TRUE)
#  }
#  return(seurat)
#}
#This is the more defualt DE test that I usually use to wrap around 
DE_test <- function(seurat, group_1, group_2 = NULL, method = "bimod",max_p = 10e-30, min_p = 0.001, 
                    foldchangethresh = 2,logfc.threshold = 0,max_cells = Inf, ...){
  require(tidyverse)
  require(Seurat)
  
  DE_output <- FindMarkers(seurat,logfc.threshold = logfc.threshold, ident.1 = group_1, ident.2 = group_2,
                           test.use = method, verbose = TRUE, max.cells.per.ident = max_cells, ...)
  DE_output <- as_tibble(DE_output, rownames = "genes")
  DE_output <- DE_output %>% mutate(Graph_adj_p = ifelse(p_val_adj > max_p, p_val_adj, max_p),
                                    significance = ifelse(p_val_adj< min_p, "p-significant","not significant"))
  DE_output <- DE_output %>% mutate(significance = ifelse((abs(avg_logFC)> log(foldchangethresh)) & significance == "p-significant", "significant",significance)) %>%
    mutate(significance = factor(significance, levels = c("not significant","p-significant","significant")))
  
  return(DE_output)
}
#Asuuming we ran the DE_test function above.
#modified so that if no list of genes are supplied, automatically look for top differnetially expressing genes in terms of increasing
#and decrreasing fold-change genes. also include the nudge parameter to put in, so that we can 
plot_volcano <- function(DE_results,gene_list_1=NULL, gene_list_2 = NULL, max_p = 10e-30, 
                         min_p = 0.001, right_nudge = 0.8, left_nudge = 0, label_size = 2){
  require(ggrepel)
  require(ggplot2)

  volcano <- ggplot(DE_results, aes(avg_logFC,-log10(Graph_adj_p))) + 
    geom_point(aes(color = significance), size = 0.3) + geom_hline(yintercept = -log10(min_p), linetype = 2) +
  geom_hline(yintercept = -log10(max_p),linetype = 2) + scale_color_manual(values =c("grey87","magenta3","darkturquoise")) +
  labs(x = "ln(fold change)", y = "-log10(adjusterd p-value)")
  
  if(!is.null(gene_list_1)){
    right_labeled_genes <- DE_results %>% filter(genes %in% gene_list_1)
    volcano <- volcano + geom_text_repel(data = right_labeled_genes, aes(label = genes),
                                         segment.size = 0.2, segment.color = "grey50", size =label_size, nudge_y = -2,
                                         nudge_x = right_nudge - right_labeled_genes$avg_logFC, direction = "y", hjust = 1) 
  } else{
    right_labeled_genes <- DE_results %>% arrange(desc(p_val_adj),desc(avg_logFC)) %>% top_n(10,avg_logFC)
    volcano <- volcano + geom_text_repel(data = right_labeled_genes, aes(label = genes),
                                         segment.size = 0.2, segment.color = "grey50", size =label_size, nudge_y = -2,
                                         nudge_x = right_nudge - right_labeled_genes$avg_logFC, direction = "y", hjust = 1) 
  }
  if(!is.null(gene_list_2)){
    left_labeled_genes <- DE_results %>% filter(genes %in% gene_list_2)
    volcano <- volcano + geom_text_repel(data = left_labeled_genes, aes(label = genes),
                                         segment.size = 0.2, segment.color = "grey50", size =label_size, nudge_y = -2,
                                         nudge_x = left_nudge + left_labeled_genes$avg_logFC, direction = "y", hjust = 1)
  } else{
    left_labeled_genes <- DE_results %>% arrange(desc(p_val_adj),(avg_logFC)) %>% top_n(10,-avg_logFC)
    volcano <- volcano + geom_text_repel(data = left_labeled_genes, aes(label = genes),
                                         segment.size = 0.2, segment.color = "grey50", size =label_size, nudge_y = -2,
                                         nudge_x = left_nudge + left_labeled_genes$avg_logFC, direction = "y", hjust = 1)  
  }
  return(volcano)
}
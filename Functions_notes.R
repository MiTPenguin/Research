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

#update: with the new Seurat, I think find Variable genes, 
Seurat_process <- function(seurat, perp = 30, cutoff = NULL, limit = FALSE,cf = 0.2, mf = 0.1,umap.run = TRUE,...) {
  seurat <- NormalizeData(object = seurat)
  if(limit!=FALSE){
    seurat <- limit_vargenes(seurat,limit,cf)
  }else{
    seurat <-  FindVariableGenes(seurat,x.low.cutoff = mf, y.cutoff = cf, do.plot = FALSE)
  }
  seurat <-  ScaleData(seurat, vars.to.regress = c('nUMI', 'percent.mito'), 
                       model.use = 'poisson', genes.use = seurat@var.genes)
  seurat <-  RunPCA(seurat, do.print = FALSE)
  if(is.null(cutoff)){
    cutoff <- find_pcCut(seurat)
  }
  print(paste("PCA cutoff used is",cutoff,"components"))
  seurat <-  RunTSNE(seurat, dims.use = 1:cutoff, perplexity = perp)
  #seurat = SetAllIdent(seurat, 'expt')
  TSNEPlot(seurat)
  if(umap.run){
    seurat <-  RunUMAP(seurat, dims.use = 1:cutoff)
    DimPlot(seurat, 'umap')
  }
  return(seurat)
}

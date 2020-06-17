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



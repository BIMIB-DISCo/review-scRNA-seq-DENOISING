if (!require("devtools")) {
  install.packages("devtools")
}
if (!require("Seurat")) {
  install.packages("Seurat")
}

if (!require("SymSim")) {
  library(devtools)
  devtools::install_github("YosefLab/SymSim")
}

if (!require("NMF")) {
  install.packages("NMF")
}

if (!require("utils")) {
  install.packages("utils")
}
if (!require("stringi")) {
  install.packages("stringi")
}
if (!require("aricode")) {
  install.packages("aricode")
}
if (!require("stringr")) {
  install.packages("stringr")
}
library(ggplot2)
# Function to perform clustering and evaluating its performances
calculate_clustering <- function(tmp_so, res = "default", clustering_performances=TRUE) {
  # PREPROCESSING FOR CLUSTERING
  tmp_so <- NormalizeData(tmp_so, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
  all.genes <- rownames(tmp_so)
  tmp_so <- ScaleData(tmp_so, 
                      features = all.genes)
  tmp_so <- RunPCA(tmp_so, 
                   npcs = 50, 
                   features = all.genes, 
                   verbose = F)
  tmp_so <- RunTSNE(tmp_so, 
                    reduction = "pca", 
                    check_duplicates = FALSE, 
                    dims = 1:50)
  #tmp_so <- RunUMAP(tmp_so, reduction = "pca", dims = 1:50)
  tmp_so <- FindNeighbors(tmp_so, dims = 1:50)
  
  # Run Louvain clustering
  tmp_so <- FindClusters(tmp_so, 
                         algorithm = 2, 
                         verbose = F, 
                         resolution = if (res=="default") 0.8 else res)
  
  # Now plot the result
  p1 = DimPlot(tmp_so, 
               reduction = "tsne", 
               group.by = "pop", 
               label = T, 
               pt.size = 1.0, 
               do.return=T)
  p2 = DimPlot(tmp_so, 
               reduction = "tsne", 
               label = T, 
               pt.size = 1.0, 
               do.return=T)
  
  trueLab = unlist(tmp_so@meta.data$pop)
  if (clustering_performances) {
    
    
    # Calculate cluster entropy
    en = round(NMF::entropy(x = tmp_so$seurat_clusters, 
                            y = trueLab), 
               digits = 3)
    
    # Calculate Rand Index for clustering solution
    rand = round(ARI(tmp_so$seurat_clusters, 
                     trueLab), 
                 digits = 3)
    
    # Calculate Mutual Information for clustering solution
    mutual = round(AMI(tmp_so$seurat_clusters, 
                       trueLab), 
                   digits = 3)
    ##### COMMENTATO PER NON PORTARE VIA TROPPO TEMPO
    # Calculate the percentage of zeros in the count matrix
    #counts = as.matrix(tmp_so@assays$counts@counts@x)
    zeros = 0#round(sum(counts == 0) / prod(dim(counts)), 
    #digits = 3)
    
    # Calculate silhouette coefficient, using the known cell population as labels
    #distance_m_o = dist(t(as.matrix(tmp_so@assays$counts@counts)))
    #pop = as.numeric(tmp_so@meta.data$pop)
    #sil = silhouette(pop, distance_m_o)
    avg_sil = 0#round(summary(sil)$avg.width, digits = 3)
    
    newList <- list("entropy" = en, 
                    "percZero" = zeros, 
                    "plotPop" = p1, 
                    "plotCluster" = p2, 
                    "clusters" = tmp_so$seurat_clusters, 
                    "ari" = rand, 
                    "ami" = mutual, 
                    "sil" = avg_sil,
                    "seurat_object" = tmp_so)
  } else {
    newList <- list("plotPop" = p1, 
                    "plotCluster" = p2)
  }
  return(newList)
}

# Function to convert summarizedExperiment to SeuratObject
summarizedToSeurat = function(experiment) {
  counts = experiment$counts
  colnames(counts) = experiment$cell_meta$cellid
  rownames(counts) = c(1:nrow(counts))
  # Should be a data frame where the rows are cell names and the columns are additional metadata fields.
  meta = experiment$cell_meta$pop
  npops = length(unique(meta))
  meta = as.data.frame(meta)
  rownames(meta) = colnames(counts)
  colnames(meta) = c("pop")
  tmp_so <- CreateSeuratObject(counts, 
                               assay = "counts",
                               names.field = 1, 
                               names.delim = "-", 
                               meta.data = meta)
  return(tmp_so)
}

# Function to convert a list with three elements (count matrix, 
# list of populations and list of cell ids) to SeuratObject
listToSeurat = function(experiment) {
  counts = experiment$counts
  colnames(counts) = experiment$cell_id
  rownames(counts) = c(1:nrow(counts))
  # Should be a data frame where the rows are cell names and the columns are additional metadata fields.
  meta = experiment$pop
  npops = length(unique(meta))
  meta = as.data.frame(meta)
  rownames(meta) = colnames(counts)
  colnames(meta) = c("pop")
  tmp_so <- CreateSeuratObject(counts, 
                               assay = "counts",
                               names.field = 1, 
                               names.delim = "-", 
                               meta.data = meta)
  return(tmp_so)
}

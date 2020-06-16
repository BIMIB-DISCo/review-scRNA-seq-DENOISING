#Normalize and filter datasets


cleanFile <- function(symsimFile){
  cleanSS <- symsimFile
  cleanSS$nreads_perUMI <- NULL
  cleanSS$nUMI2seq <- NULL
  cleanSS$cell_meta <- cleanSS$cell_meta[, c("cellid","pop","alpha","depth","batch")]
  return(cleanSS)
}

#Filtering of genes with expression == 0 in more than minGenePerc cells, and select top nVaribale most variable genes
# Return a filtered dataset, with the information about remoed and kept genes
filterDataset <- function(counts, minCellPerc, minGenePerc){
  idRm <- which((rowSums(counts==0)/dim(counts)[2]) >= (1-minCellPerc))
  idKeep <- setdiff(1:dim(counts)[1], idRm)
  idCRm <- which((colSums(counts==0)/dim(counts)[1]) >= (1-minGenePerc))
  idCKeep <- setdiff(1:dim(counts)[2], idCRm)
  res <- list()
  res$counts <- counts[idKeep,idCKeep]
  res$preprocess_info$minCellPerc <- minCellPerc
  res$preprocess_info$idG_Keep <- idKeep
  res$preprocess_info$idG_Rm <- idRm
  res$preprocess_info$minGenePerc <- minCellPerc
  res$preprocess_info$idC_Keep <- idCKeep
  res$preprocess_info$idC_Rm <- idCRm
  print(paste0("Keept ",length(idKeep), " genes and ", length(idCKeep), " cells."))
  return(res)
}

#Normalize SMARTSEQ over gene length and library size
normalizeSmartseq <- function(counts, geneL, scaleFac, logSc){
  normCount <- counts / geneL
  if(logSc){
    normCount <- round(log1p(t(t(normCount) * scaleFac / colSums(normCount))), digits = 2)
  } else {
    normCount <- round(t(t(normCount) * scaleFac / colSums(normCount)), digits = 2)
  }
  return(normCount)
}

#Normalzie 10x only over library size
normalize10x <- function(counts, scaleFac, logSc){
  if(logSc){
    normCount <- round(log1p(t(t(counts) * scaleFac / colSums(counts))), digits = 2)
  } else {
    normCount <- round(t(t(counts) * scaleFac / colSums(counts)), digits = 2)
  }
  return(normCount)
}


#_________________________________________________________________#
#Script sui dati

working_dir = ""# Set this to the path of the folder where "Review denoising" was downloaded
setwd(working_dir)
truedir = "./final_GT/rdata/"
dirCSV_true = "./final_data/csv/true_filt/"
dirCSV_noise = "./final_data/csv/noise_filt/"
dirCSV_scaled = "./final_data/csv/noise_filt_normalized/"

if(!dir.exists(dirCSV_true)){dir.create(dirCSV_true, recursive = T)}
if(!dir.exists(dirCSV_noise)){dir.create(dirCSV_noise, recursive = T)}
if(!dir.exists(dirCSV_scaled)){dir.create(dirCSV_scaled, recursive = T)}


geneLengthFile = "tuning/gene_len.rdata"
load(geneLengthFile)

minCellPerc <- 0.05
minGenePerc <- 0.1

scaleFactor = 10^4

logScale <- T

trueNames <- dir(path = truedir, pattern = "UMI|nonUMI", ignore.case = T)

for(GT_f in trueNames){
  
  datasetName = gsub(".rds|GT_", "", GT_f)
  
  noiseNames <- dir(path = "noise_mixed/", 
                    pattern = paste0("^", datasetName), 
                    ignore.case = T)
  true_counts = readRDS(paste0(truedir, GT_f))
  
  true_counts = list(counts = true_counts$counts,
                     cell_id = as.character(true_counts$cell_meta$cellid),
                     pop = true_counts$cell_meta$pop)
  colnames(true_counts$counts) = true_counts$cell_id
  names(true_counts$pop) = true_counts$cell_id
  
  for(N_f in noiseNames){
    print(paste0("Processing dataset ", N_f))
    trueFile = paste0(dirCSV_true, 
                      gsub(".rds", "_true.csv", N_f))
    
    noiseFile = paste0(dirCSV_noise, 
                       gsub(".rds", "_noise.csv", N_f))
    logscaledFile = paste0(dirCSV_scaled, 
                           gsub(".rds", "_noise_logScaled.csv", N_f))
    if(
      file.exists(trueFile) &
      file.exists(noiseFile) &
      file.exists(logscaledFile)
    ) {next}
    
    observed_counts = readRDS(paste0("noise_mixed/", N_f))
    colnames(observed_counts$counts) = observed_counts$cell_id
    names(observed_counts$pop) = observed_counts$cell_id
    # Re-order the cells in the ground truth so that they are in the same order as in the 
    # observed datasets. First re-order and then fix the cell_id order
    true_counts$counts = true_counts$counts[,observed_counts$cell_id]
    true_counts$cell_id = colnames(true_counts$counts)
    true_counts$pop = true_counts$pop[observed_counts$cell_id]
    
    resFilt <- filterDataset(observed_counts$counts, 
                             minCellPerc, 
                             minGenePerc)
    idxGenes <- resFilt$preprocess_info$idG_Keep
    idxCells = resFilt$preprocess_info$idC_Keep
    
    #Filtering true
    countFilt_True <- as.data.frame(true_counts$counts[idxGenes, idxCells])
    
    #Concatanate cell number and population info
    colnames(countFilt_True) <- paste(as.character(true_counts$cell_id[idxCells]), 
                                      true_counts$pop[idxCells], 
                                      sep = "p")
    row.names(countFilt_True) <- paste0("G_", idxGenes)
    write.table(x = countFilt_True,
                file = trueFile,
                row.names = T, 
                col.names = T, 
                sep = ",")
    rm(countFilt_True)
    
    #Filtering noise
    countFilt_Noise <- as.data.frame(resFilt$counts)
    colnames(countFilt_Noise) <- paste0(observed_counts$cell_id[idxCells], 
                                        "p", 
                                        observed_counts$pop[idxCells])
    row.names(countFilt_Noise) <- paste0("G_", idxGenes)
    write.table(x = countFilt_Noise,
                file = noiseFile,
                row.names = T, 
                col.names = T, sep = ",")
    rm(countFilt_Noise)
    
    #Normalize noise
    if(grepl(pattern = "nonUMI", N_f, ignore.case = T)){
      NormCount <- as.data.frame(normalizeSmartseq(counts = resFilt$counts, 
                                                   geneL = gene_len[idxGenes],
                                                   scaleFac = scaleFactor, 
                                                   logSc = logScale))
    } else if(grepl(pattern = "UMI", N_f, ignore.case = T)) {
      NormCount <- as.data.frame(normalize10x(counts = resFilt$counts, 
                                              scaleFac = scaleFactor,
                                              logSc = logScale))
    } else {
      stop("Dataset not smartseq or 10x")
    }
    colnames(NormCount) <- paste0(observed_counts$cell_id[idxCells], 
                                  "p", 
                                  observed_counts$pop[idxCells])
    row.names(NormCount) <- paste0("G_", idxGenes)
    write.table(x = NormCount,
                file = logscaledFile,
                row.names = T, 
                col.names = T, 
                sep = ",")
    rm(NormCount)
  }
}
library(Seurat)
library(scImpute)

repo_dir = "" # Set this to the folder containing the cloned repository (i.e. local path to "review-scRNA-seq-DENOISING")
source(paste0(repo_dir, "/register_RAM_time.R"))
# Set this to the path of the folder containing the csv of datasets with noise
# e.g., to denoise synthetic datasets with the most variables genes, set this to "localPath/review-scRNA-seq-DENOISING/final_data/csv_mostvar"
work_dir = "" 
setwd(work_dir)

algorithm = "SCimpute"

inputDir = "noise_filt/"
outputDir = paste0("denoised_",algorithm,"/")

if(!dir.exists(outputDir)){dir.create(outputDir)}

filenames = dir(inputDir, pattern = "_noise.csv")
write("experiment,method,cpuTime,elapsedTime,RAMpeak(MiB)", file=paste0(outputDir,"runtime.csv"), append=T)

for (f in filenames){
  experiment_name = gsub("_noise.csv$","", f)
  print(experiment_name)
  denoised_name <- paste0(outputDir, experiment_name, "_denoised_",algorithm,".csv")
  if(file.exists(denoised_name)){next}
  # 
  
  if (grepl("3_pop", f)) {
    nclust = 3
  } else if (grepl("5_pop", f)) {
    nclust = 5
  } else {
    # Use Seurat to identify the number of clusters
    tmp = read.csv(file=paste0(inputDir, f), row.names=1)
    tmpS <- CreateSeuratObject(tmp)
    tmpS <- ScaleData(tmpS, verbose = F)
    tmpS <- FindVariableFeatures(tmpS, verbose = F)
    tmpS <- RunPCA(tmpS, verbose = F)
    tmpS <- FindNeighbors(tmpS, verbose = F)
    tmpS <- FindClusters(tmpS, verbose = F)
    nclust = length(unique(tmpS@meta.data$seurat_clusters))
  }
  

  print(paste0("Cluster: ", nclust))
  print(paste0("Running ", algorithm))
  
  start = Sys.time()
  peakMemory = peakRAM_custom(t <- try(scimpute(count_path = normalizePath(paste0(inputDir, f)), 
                                                infile = "csv",          # format of input file
                                                outfile = "csv",         # format of output file
                                                out_dir = paste0(getwd(), "/",outputDir, "SCimpTmp"),     # full path to output directory
                                                labeled = FALSE,         # cell type labels not available
                                                drop_thre = 0.5,         # threshold set on dropout probability
                                                Kcluster = nclust,         #
                                                ncores = 8),              # number of cores used in parallel computation
                                       silent = T))
  if(class(t) != "try-error") {
    end = Sys.time()
    duration = end-start
    duration_sec = as.numeric(end) - as.numeric(start)
    print("Saving results")

    file.rename(from = "denoised_SCimpute/SCimpTmpscimpute_count.csv", to = denoised_name)
    file.remove(paste0("denoised_SCimpute/",dir("denoised_SCimpute/", pattern = "SCimpTmp*")))
    write(paste(experiment_name, 
                algorithm, 
                peakMemory$user_Time_sec, 
                peakMemory$Elapsed_Time_sec,
                peakMemory$Peak_RAM_Used_MiB,  
                sep=','), 
          file=paste0(outputDir,"runtime.csv"), append=T)   
  }else {
    print("ERROR!")
    write(paste(experiment_name, algorithm, "error", "error","error", sep=','), file=paste0(outputDir,"runtime.csv"), append=T)
  }
}



            
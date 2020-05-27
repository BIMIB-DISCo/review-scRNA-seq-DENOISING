library(Seurat)
library(Matrix)

repo_dir = "" # Set this to the folder containing the cloned repository (i.e. local path to "review-scRNA-seq-DENOISING")
source(paste0(repo_dir, "/register_RAM_time.R"))
source(paste0(repo_dir, "/denoising_script/enhance.R"))
# Set this to the path of the folder containing the csv of datasets with noise
# e.g., to denoise synthetic datasets with the most variables genes, set this to "localPath/review-scRNA-seq-DENOISING/final_data/csv_mostvar"
work_dir = "" 
setwd(work_dir)

inputDir = "noise_filt/"
outputDir = "denoised_ENHANCE/"

if(!dir.exists(outputDir)){dir.create(outputDir, recursive = T)}

filenames = dir(inputDir, pattern = "_noise.csv")
write("experiment,method,cpuTime,elapsedTime,RAMpeak(MiB)", 
      file=paste0(outputDir,"runtime.csv"), 
      append=T)

for (f in filenames){
  print(paste0("Denoising ", f, " with ENHANCE"))
  experiment_name = gsub("_noise.csv$","", f)
  print(experiment_name)
  denoised_name <- paste0(outputDir, 
                          experiment_name, 
                          "_denoised_ENHANCE.csv")
  if(file.exists(denoised_name)){next}
  tmp = read.csv(file=paste0(inputDir, f), row.names=1)

  tmpA <- CreateSeuratObject(counts = tmp)
  print("Running ENHANCE")
  start = Sys.time()
  set.seed(2)
  peakMemory = peakRAM_custom(t <- try(tmpA <- enhance_seurat_wrapper(tmpA, 
                                                                     assay = 'RNA',
                                                                     #ratio_pcs = 2,
                                                                     percent_cells_max = 2,
                                                                     setDefaultAssay = TRUE),
                                       silent=T))

  if(class(t) != "try-error") {
    end = Sys.time()
    duration = end-start
    duration_sec = as.numeric(end) - as.numeric(start)
    print("Saving results")
    tmp <- as.data.frame(as.matrix(tmpA@assays$enhance@data))
    tmp <- round(tmp, digits = 2)
    write.table(tmp, 
                file = denoised_name, 
                sep = ",", 
                col.names = T, 
                row.names = T)
    write(paste(experiment_name, 
                "ENHANCE", 
                peakMemory$user_Time_sec, 
                peakMemory$Elapsed_Time_sec,
                peakMemory$Peak_RAM_Used_MiB, 
                sep=','), 
          file=paste0(outputDir,"runtime.csv"), append=T)   
    } else {
      print("ERROR!")
      write(paste(experiment_name, "ENHANCE", "error", "error","error", sep=','), 
            file=paste0(outputDir,"runtime.csv"), 
            append=T)
    }
  
}

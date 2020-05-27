library(SAVERX)
repo_dir = "" # Set this to the folder containing the cloned repository (i.e. local path to "review-scRNA-seq-DENOISING")
source(paste0(repo_dir, "/register_RAM_time.R"))
# Set this to the path of the folder containing the csv of datasets with noise
# e.g., to denoise synthetic datasets with the most variables genes, set this to "localPath/review-scRNA-seq-DENOISING/final_data/csv_mostvar"
work_dir = "" 
setwd(work_dir)

algorithm <- "saverx"

inputDir = "noise_filt/"

filenames = dir(inputDir, pattern = "_noise.csv")

pre_trained_paths = c("","human_Tcells.hdf5", "human_Immune.hdf5")

for (weights_file in pre_trained_paths) {
  outputDir = paste0("denoised_",
                     algorithm,
                     "-",
                     gsub(".hdf5", "", weights_file),
                     "/")
  
  if(!dir.exists(outputDir)){dir.create(outputDir, recursive = T)}
  write("experiment,method,cpuTime,elapsedTime,RAMpeak(MB)", 
        file=paste0(outputDir,"runtime.csv"), 
        append=T)

  
  for (f in filenames){
  
    experiment_name = gsub("_noise.csv$","", f)
    print(experiment_name)
  
    denoised_name <- paste0(outputDir, 
                            experiment_name, 
                            "_denoised_", 
                            algorithm,
                            "-",
                            gsub(".hdf5", "", weights_file),
                            ".csv")
  
  
    if(file.exists(denoised_name)){next}
    
    print("Running saverx")
    start = Sys.time()

    pretrain_flag = if (weights_file == '') F else T
    peakMemory = peakRAM_custom(e <- try(saverx(paste0(inputDir, f), 
                    verbose=T, 
                    data.species = "Human", 
                    is.large.data = F, 
                    clearup.python.session = F, 
                    use.pretrain = pretrain_flag,
                    pretrained.weights.file = weights_file,
                    ncores = 8)))
    if(class(e) %in% 'try-error'){
      write(paste(gsub("_noise.csv$","", f), "saverx", "error","error", "ERROR!", sep=','), 
            file="denoised_saverx/runtime.csv", append=T)
      print(e)
      next
    }
    end = Sys.time()
    duration = end-start
    duration_sec = as.numeric(end) - as.numeric(start)
    print("Saving results")
    
    tmp <- readRDS(file = e)
    tmp <- round(tmp$estimate, digits = 2)
    
    write.table(tmp, file = denoised_name, 
                sep = ",", 
                col.names = T, 
                row.names = T)
    unlink(gsub("denoised.rds", "", e), recursive = T)
    write(paste(experiment_name, 
                "saverx", 
                peakMemory$user_Time_sec, 
                peakMemory$Elapsed_Time_sec,
                peakMemory$Peak_RAM_Used_MiB * 1.04858, sep=','), 
                file=paste0(outputDir,"runtime.csv"), append=T)
  }
}

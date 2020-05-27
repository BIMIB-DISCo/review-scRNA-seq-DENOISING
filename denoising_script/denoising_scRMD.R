
library(scRMD)

repo_dir = "" # Set this to the folder containing the cloned repository (i.e. local path to "review-scRNA-seq-DENOISING")
source(paste0(repo_dir, "/register_RAM_time.R"))

algorithm = "scRMD"

# Set this to the path of the folder containing the csv of datasets with noise
# e.g., to denoise synthetic datasets with the most variables genes, set this to "localPath/review-scRNA-seq-DENOISING/final_data/csv_mostvar"
work_dir = "" 
setwd(work_dir)


inputDir = "noise_filt/"
outputDir = paste0("denoised_",algorithm,"/")


if(!dir.exists(outputDir)){dir.create(outputDir)}

filenames = dir(inputDir, pattern = "_noise.csv")
write("experiment,method,cpuTime,elapsedTime,RAMpeak(MiB)", 
      file=paste0(outputDir,"runtime.csv"), append=T)

for (f in filenames){
  experiment_name = gsub("_noise.csv$","", f)
  print(experiment_name)
  denoised_name <- paste0(outputDir, experiment_name, "_denoised_",algorithm,".csv")
  if(file.exists(denoised_name)){next}
  tmp = read.csv(file=paste0(inputDir, f), row.names=1)
  
  print(paste0("Running ", algorithm))
  start = Sys.time()
  gc1 <- gc(reset = TRUE)
  peakMemory = peakRAM_custom(t <- try(tmpD <- rmd(as.matrix(tmp)), silent = T))
  gc2 <- gc()
  gc2[,6] - gc1[,2]
  end = Sys.time()
  if(class(t) != "try-error") {
    end = Sys.time()
    duration = end-start
    duration_sec = as.numeric(end) - as.numeric(start)
    print("Saving results")
    
    denoised = tmpD$exprs
    denoised <- round(denoised, digits = 2)
    
    write.table(denoised, file = denoised_name, sep = ",", col.names = T, row.names = T)
    write(paste(experiment_name, 
                algorithm, 
                peakMemory$user_Time_sec, 
                peakMemory$Elapsed_Time_sec,
                peakMemory$Peak_RAM_Used_MiB,
                sep=','), 
          file=paste0(outputDir,"runtime.csv"), append=T)   
  }else {
    print("ERROR!")
    write(paste(experiment_name, algorithm, "error", "error", "error", sep=','), 
          file=paste0(outputDir,"runtime.csv"), append=T)
  }
}
            
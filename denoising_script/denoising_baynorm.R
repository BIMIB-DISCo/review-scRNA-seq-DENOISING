library(bayNorm)
library(Seurat)

repo_dir = "" # Set this to the folder containing the cloned repository (i.e. local path to "review-scRNA-seq-DENOISING")
source(paste0(repo_dir, "/register_RAM_time.R"))
# Set this to the path of the folder containing the csv of datasets with noise
# e.g., to denoise synthetic datasets with the most variables genes, set this to "localPath/review-scRNA-seq-DENOISING/final_data/csv_mostvar"
work_dir = "" 
setwd(work_dir)

inputDir = "noise_filt/"
outputDir = "denoised_baynorm/"

if(!dir.exists(outputDir)){dir.create(outputDir, recursive = T)}

filenames = dir(inputDir, pattern = "_noise.csv")
write("experiment,method,cpuTime,elapsedTime,RAMpeak(MiB)", 
      file=paste0(outputDir,"runtime.csv"), append=T)

#f = filenames[1]

for (f in filenames){
  experiment_name = gsub("_noise.csv$","", f)
  print(experiment_name)
  denoised_name <- paste0(outputDir, experiment_name, "_denoised_baynorm.csv")
  if(file.exists(denoised_name)){next}
  tmp = read.csv(file=paste0(inputDir, f), row.names=1)
  
  conditions = NULL
  Prior_type = NULL

  tmpA <- CreateSeuratObject(counts = tmp)
  print("Running bayNorm")
  start = Sys.time()
  peakMemory = peakRAM_custom(t <- try(tmpA <- bayNorm(Data=tmpA@assays$RNA@counts,
                                                       Conditions = conditions,
                                                       Prior_type = Prior_type,
                                                       BETA_vec = NULL,
                                                       mode_version=FALSE,
                                                       mean_version = TRUE,#S=20,
                                                       verbose =FALSE,
                                                       parallel = TRUE,
                                                       NCores = 8), 
                                       silent = T))
  if(class(t) != "try-error") {
    end = Sys.time()
    duration = end-start
    duration_sec = as.numeric(end) - as.numeric(start)
    print("Saving results")
      tmp <- as.data.frame(as.matrix(tmpA$Bay_out))
      tmp <- round(tmp, digits = 2)
      write.table(tmp, file = denoised_name, sep = ",", col.names = T, row.names = T)
      write(paste(experiment_name, 
                  "baynorm", 
                  peakMemory$user_Time_sec, 
                  peakMemory$Elapsed_Time_sec,
                  peakMemory$Peak_RAM_Used_MiB, 
                  sep=','), 
            file=paste0(outputDir, "runtime.csv"), append=T)      
  } else {
    print("ERROR!")
    write(paste(experiment_name, "baynorm", "error", "error", "error", sep=','), file=paste0(outputDir, "runtime.csv"), append=T)   
  }
  
}

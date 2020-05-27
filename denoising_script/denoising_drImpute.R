# library(devtools)
# install_github('gongx030/DrImpute')

library(DrImpute)
repo_dir = "" # Set this to the folder containing the cloned repository (i.e. local path to "review-scRNA-seq-DENOISING")
source(paste0(repo_dir, "/register_RAM_time.R"))

# Set this to the path of the folder containing the csv of datasets with noise
# e.g., to denoise synthetic datasets with the most variables genes, set this to "localPath/review-scRNA-seq-DENOISING/final_data/csv_mostvar"
work_dir = "" 
setwd(work_dir)
algorithm = "DRimpute"

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
  tmp = read.csv(file=paste0(inputDir, f), row.names=1)
  
  print(paste0("Running ", algorithm))
  
  LB <- colSums(tmp)
  SF <- 1
  tmpN <- log1p((sweep(tmp,2,LB,"/")*SF))
  
  start = Sys.time()
  peakMemory = peakRAM_custom(t <- try(tmpD <- DrImpute(as.matrix(tmpN),  
                                                        mc.cores = 8),
                                       silent = T))
  if(class(t) != "try-error") {
    end = Sys.time()
    duration = end-start
    duration_sec = as.numeric(end) - as.numeric(start)
    print("Saving results")
    
    tmpD <- log1p((exp(tmpD)-1)*10000) #Riscalo per avere dei valori più alti
    tmpD <- round(tmpD, digits = 2)
    
    write.table(tmpD, file = denoised_name, sep = ",", col.names = T, row.names = T)
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





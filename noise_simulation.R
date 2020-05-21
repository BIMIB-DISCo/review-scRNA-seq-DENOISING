# This script simulates the sequencing process. For each ground truth (5 and 3 populations) it simulates 50 UMI sequencing and 50 nonUMI sequencing experiments.
# For nonUMI it first samples 100 cells from the ground truth and then it simulates sequencing considering only the sampled cells.

working_dir = ""# Set this to the path of the folder where "Review denoising" was downloaded
setwd(working_dir)
outdir = paste0(working_dir, "noise_mixed/")


library(foreach)
library(doParallel)

rndCellSempler <- function(true_counts, idxCell){
  sampledCount <- true_counts
  sampledCount$counts <- sampledCount$counts[, idxCell]
  sampledCount$cell_meta <- sampledCount$cell_meta[idxCell,]
  sampledCount$cell_meta$cellid <- droplevels(sampledCount$cell_meta$cellid)
  sampledCount$kinetic_params[[1]] <- sampledCount$kinetic_params[[1]][,idxCell]
  sampledCount$kinetic_params[[2]] <- sampledCount$kinetic_params[[2]][,idxCell]
  sampledCount$kinetic_params[[3]] <- sampledCount$kinetic_params[[3]][,idxCell]
  return(sampledCount)
}

if(!dir.exists(outdir)){dir.create(outdir, recursive = T)}

if(!dir.exists(paste0(outdir, "final_noise/"))) {
  dir.create(paste0(outdir, "final_noise/"), recursive = T)
}
if(!dir.exists(paste0(outdir, "final_noise/figures/"))) {
  dir.create(paste0(outdir, "final_noise/figures/"), recursive = T)
}


ground_truth_list = list.files("final_GT/rdata/", pattern = "^GT")

for (ground_truth_file in ground_truth_list) {
  ground_truth = readRDS(paste0("./final_GT/rdata/", ground_truth_file))
  
  for (protocol in c("UMI", "nonUMI")) {
    sym_name = gsub("GT_|.rds", "", ground_truth_file)
    result_folder = paste0(outdir, "final_noise/", protocol, "_", sym_name, "/")
    if (!dir.exists(result_folder)) {
      dir.create(result_folder, recursive = T)
    }
    if (!dir.exists(paste0(result_folder, "figures/"))) {
      dir.create(paste0(result_folder, "figures/"), recursive = T)
    }
    # Set the values for the simulation parameters
    if (protocol=="UMI") {
      set.seed(1)
      m = 0.009
      alpha_mean = rnorm(50, mean = m, sd = 0.001)
      set.seed(2)
      d = 85000
      depht_mean = round(rnorm(50, mean = d, sd = d*0.20))
      ground_truth_new = ground_truth
    } else {
      set.seed(1)
      m = 0.1
      alpha_mean = rnorm(50, mean = m, sd = 0.01)
      set.seed(2)
      d = 200000
      depht_mean = round(rnorm(50, mean = d, sd = d*0.20))
      set.seed(3)
      # Sample 100 cells and save the nonUMI ground truth file
      id_sampled_cells = sample(1:ncol(ground_truth$counts), 100)
      ground_truth_new = rndCellSempler(ground_truth, 
                                        id_sampled_cells)
      saveRDS(object = ground_truth_new,
              file=paste0("final_GT/rdata/GT_", sym_name, "_nonUMI.rds"))
    }
    load("gene_len.rdata")
    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(15)
    registerDoParallel(cl)
    
    # Simulate the 50 sequencing experiments
    foreach(r=1:length(alpha_mean), .combine=cbind,
            .packages = c("SymSim", "Seurat")) %dopar% {
      alpha_mean_r = alpha_mean[r]
      depht_mean_r = depht_mean[r]
      depht_sd_r = depht_mean_r * 0.25
      alpha_sd_r = alpha_mean_r * 0.25
      # File name where the simulation result will be stored in .rds format
      fileN = paste0(protocol, 
                     "_alpha", 
                     alpha_mean_r, 
                     "_depht", 
                     depht_mean_r, 
                     "_", 
                     "noise.rds")
      
      if (!file.exists(paste0(result_folder, fileN))) {
        source("tmp_functions.R")
        
        to_study_n <- True2ObservedCounts(true_counts = ground_truth_new[[1]], 
                                          meta_cell = ground_truth_new[[3]],
                                          protocol = protocol, 
                                          alpha_mean= alpha_mean_r, 
                                          alpha_sd = alpha_sd_r, 
                                          gene_len = gene_len, 
                                          depth_mean = depht_mean_r, 
                                          nPCR1 = 14, 
                                          nbatch = 1, 
                                          depth_sd = depht_sd_r)
  
        saveRDS(object = to_study_n, 
                file=paste0(result_folder, fileN))
        tmp_so <- summarizedToSeurat(to_study_n)
        clustering_result = calculate_clustering(tmp_so)
        
        ggsave(paste0(result_folder, "figures/population_", gsub(".rds", ".png", fileN)), 
               plot = clustering_result$plotPop, 
               dpi = 500)
        ggsave(paste0(result_folder, "figures/clustering_", gsub(".rds", ".png", fileN)), 
               plot = clustering_result$plotCluster, 
               dpi = 500)
      }
    }
  } 
}

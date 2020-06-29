# This script simulates the ground truth datasets with 3 and 5 subpopulations
library(foreach)
library(doParallel)

working_dir = ""# Set this to the path of the folder where "Review denoising" was downloaded
setwd(working_dir)
outdir = working_dir

if(!dir.exists(outdir)){dir.create(outdir, recursive = T)}

if(!dir.exists(paste0(outdir, "ground_truth/figures/"))) {
  dir.create(paste0(outdir, "ground_truth/figures/"), recursive = T)
}

if(!dir.exists(paste0(outdir, "ground_truth/rdata/"))) {
  dir.create(paste0(outdir, "ground_truth/rdata/"), recursive = T)
}

source("tmp_functions.R")

SEED_SYM = 10

#load(paste0(outdir, "tuning/tuning_gt/g.rdata"))

load("gene_len.rdata")

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(2)
registerDoParallel(cl)

# Simulate Ground Truth
n_pop = c(3,5)
foreach(pop=n_pop, .combine=cbind,
  .packages = c("SymSim", "Seurat")) %dopar% {
  source("./tmp_functions.R")
  
  print("Simulating biological counts")
  to_study = SimulateTrueCounts(ncells_total = 3000,
                                ngenes = 10000,
                                min_popsize = 600,
                                nevf = 60,
                                n_de_evf = 9,
                                evf_type = 'discrete',
                                vary = 's',
                                gene_effect_prob = 0.1,
                                gene_effects_sd = 0.5,
                                Sigma = 0.5,
                                randseed = SEED_SYM,
                                phyla = if (pop == 5) Phyla5() else Phyla3())
  
  fileN = paste0("GT_", pop, "_pop")
  saveRDS(object=to_study, 
          file=paste0(outdir, 
                      "ground_truth/rdata/", 
                      fileN,
                      ".rds"))
  tmp_so <- summarizedToSeurat(to_study)
  clustering_result = calculate_clustering(tmp_so)
  
  ggsave(paste0(outdir, "ground_truth/figures/", fileN, "_TPOP.png"), 
         plot = clustering_result$plotPop, 
         dpi = 500)
  ggsave(paste0(outdir, "ground_truth/figures/", fileN, "_TCLUS.png"), 
         plot = clustering_result$plotCluster, 
         dpi = 500)
  }

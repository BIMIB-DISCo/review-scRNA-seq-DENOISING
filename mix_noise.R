# This script mixes cells from different sequencing experiments in order to produce the final noisy datasets.

working_dir = ""# Set this to the path of the folder where "Review denoising" was downloaded
setwd(working_dir)
outdir = paste0(working_dir, "noise_mixed/")

if(!dir.exists(outdir)){dir.create(outdir, recursive = T)}

source("tmp_functions.R")
load_symsim_experiment = function(experiment_path) {
  curr_obj = readRDS(experiment_path)
  curr_obj$perc_zeros = round(sum(curr_obj$counts == 0) / prod(dim(curr_obj$counts)),
                              digits=3)
  # Convert cell_id to list, this is necessary later when the ground truths will be mixed 
  curr_obj$cell_meta$cellid = as.character(curr_obj$cell_meta$cellid)
  zeros_cells = c()
  for (c in 1:ncol(curr_obj$counts)) {
    curr_zeros = round(sum(curr_obj$counts[,c] == 0) / nrow(curr_obj$counts),
                       digits=3)
    zeros_cells = c(zeros_cells, curr_zeros)
  }
  curr_obj$zeros_per_cell = zeros_cells
  return(curr_obj)
}

experiments_noise = list.files("./final_noise/", pattern='^UMI|^nonUMI')

for (experiment in experiments_noise) {
  if (!dir.exists(paste0(outdir, "figures/"))) {
    dir.create(paste0(outdir, 
                      "figures/"), 
               recursive = T)
  }

  experiment_dir = paste0("./final_noise/", 
                          experiment, "/")
  sequencing_list = list.files(experiment_dir, 
                               pattern = '^UMI|^nonUMI')
  
  sequencing_objects = list()
  # Read all the sequencing experiments over the considered ground truth
  for (i in seq(sequencing_list)) {
    obj_name = gsub(".rds|UMI_|nonUMI", 
                    "", 
                    sequencing_list[i])
    sequencing_objects[[obj_name]] = load_symsim_experiment(paste0(experiment_dir, 
                                                                   sequencing_list[i]))
  }
  # Check that the cells are all in the same order in all objects
  for (i in 2 : length(sequencing_objects)) {
    cell_id_i = sequencing_objects[[i]]$cell_meta$cellid 
    cell_id_i1 = sequencing_objects[[i-1]]$cell_meta$cellid
    equal_cells = cell_id_i == cell_id_i1
    if (sum(equal_cells) < ncol(sequencing_objects[[i]]$counts)) {
      stop("Careful, cells are not in the same order")
    }
  }

  cell_ids = 1:length(sequencing_objects[[1]]$cell_meta$cellid)
  
  
  nonUMI_flag = grepl("nonUMI",experiment)
  ncells = length(cell_ids)
  seed_i = 1
  n_experiments = if (nonUMI_flag) c(2,4,6,8,10) else c(10,20,30,40,50)
  
  for (n in n_experiments) {
    set.seed(seed_i)
    permutation = sample(x = 1:length(cell_ids), 
                         size = ncells)
    seed_i = seed_i + 1
    new_object = list()
    set.seed(seed_i)
    cuts = sample(1:n*2, n)
    seed_i = seed_i + 1
    cuts = cuts/sum(cuts)
    # obj_considered is the list of indeces of the objects from which cells are samples. Thus, 
    # I randomly sample from them in order to have 10, 20, 30, 40 and finally 50 objects to sample from
  
    set.seed(seed_i)
    objs_considered = sample(1:length(sequencing_objects), n)
    seed_i = seed_i + 1
    start = 1
    for (i in 1:length(cuts)) {
      # Calculated the width of the window of cells that need to be considered in this step
      window_width = floor(ncells * cuts[i]) 
      if (i == length(cuts)) {
        end = ncells
      } else {
        end = start + window_width -1
      }
      
      # Indices of the window of cells that were permuted
      current_cut = start:end 
      # Take the ids of cells that need to be considered in this step
      current_cells = permutation[current_cut] 
      start = end + 1

      j = objs_considered[i]
      new_object$counts = cbind(new_object$counts, 
                                sequencing_objects[[j]]$counts[,current_cells])
      new_object$cell_id[current_cut] = sequencing_objects[[j]]$cell_meta$cellid[current_cells]
      new_object$pop[current_cut] = sequencing_objects[[j]]$cell_meta$pop[current_cells]
      new_object$source_obj[current_cut] = j

    }
    
    tmp_so = listToSeurat(new_object)
    clustering_result = calculate_clustering(tmp_so, clustering_performances = FALSE)
    
    fileN = paste0(experiment, "_", n, "simulations")
    saveRDS(file=paste0(outdir, fileN, ".rds"), new_object)
    
    ggsave(paste0(outdir, "figures/clustering_", fileN, ".png"), 
           plot = clustering_result$plotCluster, 
           dpi = 500)
    ggsave(paste0(outdir, "figures/populations_", fileN, ".png"), 
           plot = clustering_result$plotPop, 
           dpi = 500)
  }
  
}












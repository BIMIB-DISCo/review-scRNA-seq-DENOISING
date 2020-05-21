# review-scRNA-seq-DENOISING

This repository contains the source code that was used for the paper *A review of computational strategies for denoisingand imputation of single-cell transcriptomic data*

### Folder organization

**final_GT**: this folder contains the rds files with the datasets used as ground truth

**noise_mixed**: this folder contains the rds files with the full datasets with noise

After simulating the counts, for each dataset we applied a basic preprocessing, removing
genes expressed in less than 5% of cells and cells with more than 90% of zero counts. 
Then we also selected for each dataset the top 500,1000,2000 and 5000 most variable genes.

**final_data/csv/noise_filt** contains the csv of the pre-processed noise files

**final_data/csv_mostvar/noise_filt** contains the csv files of the pre-processed noise files with the selection
				  of the most variable genes

In Both csv_mostvar and csv there is a folder named noise_filt_normalized: this folder contains
the noise datsets that were also library size normalized (with a scaling factor of 1e4) and then
log-transformed (log(X+1)).

Both in csv_mostvar and csv there is a folder named "true_filt". This file contains the ground truth
for each dataset with noise, without those genes and cells that were removed during the pre-processing
of the noise files.

IMPORTANT: in each csv, the information about each cell's subpopulation is stored in the cell name.
In particular, each cell has a name cellxxxxpx, the number after p indicates the population.

### REPRODUCIBILITY:
the scripts need to be run in the following order:
GT_simulation.R (produces the files located in the folder "final_GT"
noise_simulation.R (it produces the files located in the folder "final_noise")
mix_noise.R (it produces the files located in folder "noise_mixed")
preprocess_dataset.R (it produces the files located in folder "final_data/csv")
preprocess_dataset_highlyvar.R (it produces the files located in folder "final_data/csv_mostvar")





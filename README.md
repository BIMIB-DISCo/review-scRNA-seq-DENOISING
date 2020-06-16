# review-scRNA-seq-DENOISING

This repository contains the source code that was used for the paper
*A review of computational strategies for denoising and imputation of
single-cell transcriptomic data.*


### Folder organization

In order to download the datasets used in the analyses, after cloning
this repository in your local machine, download the zip files located
[here](https://drive.google.com/drive/folders/1IGREoIP3jGzWah5mnS13xdnZlSTo29wZ?usp=sharing)
and extract them in the local folder that contains this
repository. After downloading all files, you will have the following
folders:

- **`final_GT`**: this folder contains the rds files with the datasets
  used as ground truth.

- **`noise_mixed`**: this folder contains the rds files with the full
  datasets with noise.

After simulating the counts, for each dataset we applied a basic
preprocessing, removing genes expressed in less than 5% of cells and
cells with more than 90% of zero counts.  Then we also selected for
each dataset the top 500, 1000, 2000 and 5000 most variable genes.

- **`final_data/csv/noise_filt`** contains the csv of the
  pre-processed noise files.

- **`final_data/csv_mostvar/noise_filt`** contains the csv files of
  the pre-processed noise files with the selection of the most
  variable genes.

In Both `csv_mostvar` and `csv` there is a folder named
`noise_filt_normalized`: this folder contains the noise datsets that
were also library size normalized (with a scaling factor of 1e4) and
then log-transformed (log(X+1)).

Both in `csv_mostvar` and `csv` there is a folder named `true_filt`. This
folder contains the ground truth for each dataset with noise, without
those genes and cells that were removed during the pre-processing of
the noise files.

IMPORTANT: in each csv, the information about each cell's
subpopulation is stored in the cell name. In particular, each cell has
a name `cellxxxxpx`, the number after `p` indicates the population.


### Reproducibility

To reproduce the **generation of synthetic datasets**, the scripts
need to be run in the following order:

* `GT_simulation.R`, which produces the files located in the folder `final_GT`.

* `noise_simulation.R`, which produces the files located in the folder
  `final_noise`.
  
* `mix_noise.R`, which produces the files located in folder `noise_mixed`.

* `preprocess_dataset.R`, which produces the files located in folder
  `final_data/csv`.
  
* `preprocess_dataset_highlyvar.R`, which produces the files located in
  folder `final_data/csv_mostvar`.
  

To reproduce the results of the different **denoising** methods,
follow these steps.

* For all denoisers written in **R**, execute the corresponding script
  located in folder `denoising_script`. Note that, before running the
  script, you need to change to working directory to the folder that
  contains the folder `noise_filt`
  (e.g. `review-scRNA-seq-DENOISING/final_data/csv` to denoise synthetic
  datasets).
  
* For all denoisers wirtten in **Python**, you need to execute the script
  `all_python.py`, executing it inside the directory
  `denoising_script`.







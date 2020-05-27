"""

This script executes the code to run DCA, MAGIC, DeepImpute and kNN-smoothing. Before running it, set the value of variable data_directory.

"""

import pandas as pd
import os
import re
import time
import datetime
import dca.api as d
import scanpy as sc
import anndata
from anndata import AnnData
from convert_adata import anndata_to_csv
import numpy as np
from memory_profiler import memory_usage
from resource import getrusage as resource_usage, RUSAGE_SELF
from time import time as timestamp


# Set this to the path of the folder containing the csv of datasets with noise
# e.g., to denoise synthetic datasets with the most variables genes, set this to "localPath/review-scRNA-seq-DENOISING/final_data/csv_mostvar"
data_directory = "" 

algorithm = "DCA"
inputDir = f"{data_directory}/noise_filt"
outputDir = f"{data_directory}/denoised_{algorithm}"

if not os.path.isdir(outputDir):
    os.mkdir(outputDir)

filenames = os.listdir(inputDir)
filenames = [x for x in filenames if "noise.csv" in x]
f = open(f"{outputDir}/{algorithm}_runtime.csv", "a+")
f.write("experiment,method,cpuTime,elapsedTime,RAMpeak(MB)")
f.close()


# In[3]:


for f in filenames:
    experiment_name = re.sub("_noise.csv$","", f)
    denoised_name = f"{outputDir}/{experiment_name}_denoised_{algorithm}.csv"
    print(denoised_name)
    if os.path.isfile(denoised_name):
        continue
        
    adata = sc.read(f"{inputDir}/{f}")
    adata = adata.transpose()
    loss = "zinb-conddisp" if ("nonUMI" in experiment_name) else "nb-conddisp"

    adata.X = np.ceil(adata.X)

    start_time, start_resources = timestamp(), resource_usage(RUSAGE_SELF)
    (mem_registered, adata_denoised) = memory_usage((d.dca, (adata,),
                                                    {"copy":True, 
                                                    "log1p":False, 
                                                    "return_info":True, 
                                                    "verbose":True, 
                                                    "reduce_lr":5, 
                                                    "scale":False, 
                                                    "ae_type":loss}),
                                                    retval = True, 
                                                    max_usage = True, 
                                                    include_children =True)
    end_resources, end_time = resource_usage(RUSAGE_SELF), timestamp()

    real = end_time - start_time
    systime =  end_resources.ru_stime - start_resources.ru_stime
    usertime = end_resources.ru_utime - start_resources.ru_utime
    cpu_time = systime + usertime
    
    adata_denoised.X = np.round(adata_denoised.X)
    print("Saving denoised_data")
    anndata_to_csv(adata_denoised.transpose(), 
                  denoised_name)
    print("Saving scale factors")
    os.makedirs(f"{outputDir}/size/", 
                exist_ok=True)
    adata_denoised.obs['size_factors'].to_csv(f"{outputDir}/size/{experiment_name}_dca_size.csv")

    if "nonUMI" in f:
        print("Saving dropout probabilities")
        tmp = AnnData(adata_denoised.obsm["X_dca_dropout"])
        tmp.obs.index = adata_denoised.obs.index
        tmp.var.index = adata_denoised.var.index
        os.makedirs(f"{outputDir}/dropout/", exist_ok=True)
        anndata_to_csv(tmp.transpose(), 
                         f"{outputDir}/dropout/{experiment_name}_dca_dropout.csv")
    file = open(f"{outputDir}/{algorithm}_runtime.csv", "a+")
    file.write(f"\n{experiment_name},{algorithm},{str(cpu_time)},{str(real)},{str(mem_registered)}")
    file.close()


from deepimpute.multinet import MultiNet


algorithm = "DeepImpute"
inputDir = f"{data_directory}/noise_filt"
outputDir = f"{data_directory}/denoised_{algorithm}"

if not os.path.isdir(outputDir):
    os.mkdir(outputDir)

filenames = os.listdir(inputDir)
filenames = [x for x in filenames if "noise.csv" in x]
f = open(f"{outputDir}/{algorithm}_runtime.csv", "a+")
f.write("experiment,method,cpuTime,elapsedTime,RAMpeak(MB)")
f.close()

for f in filenames:
    experiment_name = re.sub("_noise.csv$","", f)
    denoised_name = f"{outputDir}/{experiment_name}_denoised_{algorithm}.csv"
    print(denoised_name)
    if os.path.isfile(denoised_name):
        continue
    tmp = pd.read_csv(f"{inputDir}/{f}", index_col=0)
    tmp = tmp.transpose()
    
    # Using custom parameters
    outputdim = 512 if tmp.shape[1] > 500 else 256
    intermediate = 200 if tmp.shape[1] > 500 else 128
    NN_params = {
            'learning_rate': 1e-4,
            'batch_size': 64,
            'max_epochs': 500,
            'ncores': 8,
            'sub_outputdim': outputdim,
            'architecture': [
                {"type": "dense", "activation": "relu", "neurons": intermediate},
                {"type": "dropout", "activation": "dropout", "rate": 0.3}]
        }

    multinet = MultiNet(**NN_params)

    start_time, start_resources = timestamp(), resource_usage(RUSAGE_SELF)
    mem_registered = memory_usage((multinet.fit,
                                (tmp,),
                                {'cell_subset': 1,
                                 'minVMR':0.5}),
                                retval = False, 
                                max_usage = True, 
                                include_children =True)
    end_resources, end_time = resource_usage(RUSAGE_SELF), timestamp()
    imputedData = multinet.predict(tmp)

    real = end_time - start_time
    systime =  end_resources.ru_stime - start_resources.ru_stime
    usertime = end_resources.ru_utime - start_resources.ru_utime
    cpu_time = systime + usertime

    imputedData = np.round(imputedData)
    imputedData.transpose().to_csv(denoised_name, index=True, header=True) 
    file = open(f"{outputDir}/{algorithm}_runtime.csv", "a+")
    file.write(f"\n{experiment_name},{algorithm},{str(cpu_time)},{str(real)},{str(mem_registered)}")
    file.close()


import magic


def find_pca_comp(adata, figName, figTitle):
    adata_scaled = adata.copy()
    sc.pp.scale(adata_scaled)
    pca = PCA().fit(adata_scaled.X)
    #Plotting the Cumulative Summation of the Explained Variance
    plt.figure()
    n = next(i for i,x in enumerate(np.cumsum(pca.explained_variance_ratio_)) if x >= 0.7)
    plt.plot(np.cumsum(pca.explained_variance_ratio_))
    plt.vlines(n, 0, 1, linestyles="dashed")
    plt.xlabel('Number of Components')
    plt.ylabel('Variance (%)') #for each component
    plt.title(figTitle)
    plt.savefig(figName, dpi=200)
    plt.close()
    return n



algorithm = "MAGIC"
inputDir = f"{data_directory}/noise_filt_normalized"
outputDir = f"{data_directory}/denoised_{algorithm}"

if not os.path.isdir(outputDir):
    os.mkdir(outputDir)

filenames = os.listdir(inputDir)
filenames = [x for x in filenames if "noise_logScaled.csv" in x]
f = open(f"{outputDir}/{algorithm}_runtime.csv", "a+")
f.write("experiment,method,cpuTime,elapsedTime,RAMpeak(MB)")
f.close()

os.makedirs(f"{outputDir}/figures/", exist_ok=True)

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
for f in filenames:
    experiment_name = re.sub("_noise_logScaled.csv$","", f)
    denoised_name = f"{outputDir}/{experiment_name}_denoised_{algorithm}.csv"
    print(denoised_name)
    if os.path.isfile(denoised_name):
        continue

    adata = sc.read(f"{inputDir}/{f}")
    adata = adata.transpose()
    adata.X = np.expm1(adata.X)
    sc.pp.sqrt(adata)

    n = find_pca_comp(adata, 
                      figName = f"{outputDir}/figures/{experiment_name}_variance.png",
                      figTitle = f'{experiment_name} Explained Variance')

    magic_op = magic.MAGIC(t=6, n_pca = n)

    start_time, start_resources = timestamp(), resource_usage(RUSAGE_SELF)
    (mem_registered, adata_denoised) = memory_usage((magic_op.fit_transform, 
                                                    (adata,), 
                                                    {'genes': 'all_genes'}), 
                                                    retval = True, 
                                                    max_usage = True, 
                                                    include_children =True)

    end_resources, end_time = resource_usage(RUSAGE_SELF), timestamp()
    real = end_time - start_time
    systime =  end_resources.ru_stime - start_resources.ru_stime
    usertime = end_resources.ru_utime - start_resources.ru_utime
    cpu_time = systime + usertime


    adata_denoised.X = adata_denoised.X**2
    print("Saving denoised_data")
    anndata_to_csv(adata_denoised.transpose(), 
                  denoised_name)
    
    file = open(f"{outputDir}/{algorithm}_runtime.csv", "a+")
    file.write(f"\n{experiment_name},{algorithm},{str(cpu_time)},{str(real)},{str(mem_registered)}")
    file.close()


from knn_smoothing.knn_smooth import knn_smoothing

algorithm = "knn_smoothing"
inputDir = f"{data_directory}/noise_filt"
outputDir = f"{data_directory}/denoised_{algorithm}"

if not os.path.isdir(outputDir):
    os.mkdir(outputDir)

filenames = os.listdir(inputDir)
filenames = [x for x in filenames if "noise.csv" in x]
f = open(f"{outputDir}/{algorithm}_runtime.csv", "a+")
f.write("experiment,method,cpuTime,elapsedTime,RAMpeak(MB)")
f.close()

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
for f in filenames:
    experiment_name = re.sub("_noise.csv$","", f)
    denoised_name = f"{outputDir}/{experiment_name}_denoised_{algorithm}.csv"
    print(denoised_name)
    if os.path.isfile(denoised_name):
        continue

    adata = sc.read(f"{inputDir}/{f}")

    start_time, start_resources = timestamp(), resource_usage(RUSAGE_SELF)


    (mem_registered, S) = memory_usage((knn_smoothing,
                                    (adata.X,),
                                    {'k': 16,
                                    'd':10,
                                    'dither':0.1,
                                    'seed':0}), 
                                retval = True, 
                                max_usage = True, 
                                include_children =True)

    end_resources, end_time = resource_usage(RUSAGE_SELF), timestamp()
    real = end_time - start_time
    systime =  end_resources.ru_stime - start_resources.ru_stime
    usertime = end_resources.ru_utime - start_resources.ru_utime
    cpu_time = systime + usertime

    denoised = pd.DataFrame(S, index=adata.obs_names, columns=adata.var_names)

    print("Saving denoised_data")
    denoised.to_csv(denoised_name, header=True, index=True)

    file = open(f"{outputDir}/{algorithm}_runtime.csv", "a+")
    file.write(f"\n{experiment_name},{algorithm},{str(cpu_time)},{str(real)},{str(mem_registered[0])}")
    file.close()
import pandas as pd
def anndata_to_csv(adata, path):
    """
    This function saves an anndata object to csv, with header and rownames. 
    Input:
        - adata: AnnData object to save
        - path: path to the csv file to save
    """
    df = pd.DataFrame(data=adata.X, index=adata.obs.index, columns=adata.var.index)
    df.to_csv(path)
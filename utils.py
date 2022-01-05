from pandas.io.pytables import ClosedFileError
import scvi
import scanpy as sc
import torch
import numpy as np
import pandas as pd
def gene2ct(adata,cell_type_df_):
    """\
    Read 10x-Genomics-formatted mtx directory.
    Parameters
    ----------
    adata
        The AnnData object to obtain scores for the genes
    cell_type_df_
        The dataframe with columns [official gene symbol,	cell type,	organ] 
        for data points belonging to a cluster
    Returns
    -------
    A tuple containing organ with maximum ocurrences for a given gene, cell type
    with highest cumulative score for all genes (organ corresponds to that cell type)
    """
    genestoscore={}
    x=list(adata.uns['rank_genes_groups']['names'])
    y=list(adata.uns['rank_genes_groups']['scores'])
    for i in range(len(x)):
        x1=list(x[i])
        y1=list(y[i])
    for j in range(len(x1)):
        if(not(x1[j] in genestoscore.keys())):
            genestoscore[x1[j]]=y1[j]
    genestoorgan={}
    for i in range(len(cell_type_df_)):
        if(cell_type_df_["official gene symbol"][i] in genestoorgan.keys()):
            genestoorgan[cell_type_df_["official gene symbol"][i]].append([cell_type_df_["organ"][i],cell_type_df_["cell type"][i],(genestoscore[cell_type_df_["official gene symbol"][i]] if cell_type_df_["official gene symbol"][i] in genestoscore.keys() else 0)])
        else:
            genestoorgan[cell_type_df_["official gene symbol"][i]]=[]
            genestoorgan[cell_type_df_["official gene symbol"][i]].append([cell_type_df_["organ"][i],cell_type_df_["cell type"][i],(genestoscore[cell_type_df_["official gene symbol"][i]] if cell_type_df_["official gene symbol"][i] in genestoscore.keys() else 0)])
    genestoorgan1={}
    for i in genestoorgan.keys():
        genestoorgan1[i]=sorted(genestoorgan[i], key = list(np.array(genestoorgan[i])[:,0]).count,
                                    reverse = True)  
    genestocelltype={}
    for i in genestoorgan1.keys():
        genestocelltype[i]=[genestoorgan1[i][0]]
        for j in range(len(genestoorgan1[i])-1):
            if(genestoorgan1[i][0][0]==genestoorgan1[i][j+1][0]):
                genestocelltype[i].append(genestoorgan1[i][j+1])
            else:
                break
    celltypetoscore={}
    for i in genestocelltype.keys():
        for j in genestocelltype[i]:
            if(j[1] in celltypetoscore.keys()):
                celltypetoscore[j[1]][0]+=j[2]
            else:
                celltypetoscore[j[1]]=[0,j[0]]
                celltypetoscore[j[1]][0]+=j[2]
    mc=""
    m=-1000
    mo=""
    for i in celltypetoscore.keys():
        if(celltypetoscore[i][0]>m):
            mc=i
            m=celltypetoscore[i][0]
            mo=celltypetoscore[i][1]
    return (mo,mc)

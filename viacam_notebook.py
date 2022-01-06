from pandas.io.pytables import ClosedFileError
import scvi
import scanpy as sc
import torch
import numpy as np
import pandas as pd
import preprocessing as prep
import models
import utils
import argparse


def viacam(dir_name="./input/",
    dataset_name = "brain",
    n_models = 2,
    train_iters = 20,
    n_latent=10,
    use_pca = 'true',
    n_pca = 20):
    results_file = 'write/dat.h5ad'

    
    
    adata=prep.preprocess(dir_name)
        
    adata,latents=models.train(adata, n_models,train_iters)
    
    adata.obsm["X_aggr"] = np.concatenate((latents), axis=1)
    


    sc.pp.neighbors(adata, use_rep="X_aggr")
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, key_added="leiden_aggr", resolution=0.5)

    sc.pl.umap(
        adata,
        color=["leiden_aggr"],
        frameon=True,
        title = dataset_name, 
        save = "_vae_only.pdf"
    )

    if use_pca:
        sc.tl.pca(adata, svd_solver='arpack')
        x_pca = adata.obsm['X_pca'][:,:n_pca]
        latents.append(x_pca)
        adata.obsm["X_aggr"] = np.concatenate((latents), axis=1)
        sc.pp.neighbors(adata, use_rep="X_aggr")
        sc.tl.umap(adata, min_dist=0.3)
        sc.tl.leiden(adata, key_added="leiden_aggr", resolution=0.5)

        sc.pl.umap(
            adata,
            color=["leiden_aggr"],
            frameon=True,
            title = dataset_name, 
            save = "_vae+pca.pdf"
        )
    sc.tl.rank_genes_groups(adata, 'leiden_aggr', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

    clusters = list(set(adata.obs['leiden_aggr']))
    dat = [list(ele) for ele in list(adata.uns['rank_genes_groups']['names'])] 
    dat_score = [list(ele) for ele in list(adata.uns['rank_genes_groups']['scores'])] 
    ranked_genes_df = pd.DataFrame(data=dat, columns=clusters)
    ranked_genes_df['scores'] = dat_score

    cell_type_df = pd.read_csv("./data/markers.tsv", sep="\t")
    cell_type_df_ = cell_type_df[["official gene symbol", "cell type","organ"]].fillna('')

    gene2celltype = {}

    for i in clusters:
        i = int(i)
        cluster_i = pd.DataFrame(list(ranked_genes_df.iloc[:, i]), columns=["official gene symbol"])
        cluster_i['scores'] = [s[i] for s in list(ranked_genes_df['scores'])]
        cell_type_i = pd.merge(cluster_i,cell_type_df_, how ='inner', on =["official gene symbol"])
        gene2celltype[i] = utils.gene2ct(adata, cell_type_i)
        

        
        # save the file for prototyping ranking function
        filename = "cluster_" + str(i) + ".csv"
        cluster_i = pd.merge(cluster_i,cell_type_i, how ='left', on =["official gene symbol"])
        cluster_i.to_csv(filename, index=False)
    markers = list(ranked_genes_df.iloc[0])[:-1]
    sc.pl.heatmap(
        adata,
        markers,
        groupby='leiden_aggr',
        standard_scale="var",
        dendrogram=True,
        figsize=(6,8),
        save = ".pdf"
    )
    labels = {}
    for k,v in gene2celltype.items():
        labels[k+1] = str(k) +" " + str(v[1])
        
    l=[]
    m=labels
    for i in range(len(adata.obs)):
      l.append(m[int(adata.obs["leiden_aggr"][i])+1])
    adata.obs["celltype"]=l

    sc.pl.umap(
        adata, 
        color=["celltype"],
        frameon=False,
        title = dataset_name,
        save = '_annot.pdf'
    )
    sc.tl.pca(adata, svd_solver='arpack')   
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, key_added="leiden_aggr", resolution=0.5)

    sc.pl.umap(
        adata,
        color=["leiden_aggr"],
        frameon=True,
        title = dataset_name, 
        save = "_pca_only.pdf"
    )

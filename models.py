from pandas.io.pytables import ClosedFileError
import scvi
import scanpy as sc
import torch
import numpy as np
import pandas as pd
def train(adata, n_models=2,train_iters=1):
    latents = []
    for i in range(1,n_models+1):
        layer_name = "counts_" + str(i)
        z_name = "X_latent_"+ str(i)
        scvi.model.SCVI.setup_anndata(adata, layer=layer_name)
        model = scvi.model.SCVI(adata, n_latent=10, n_layers=i)
        model.train(max_epochs=train_iters, plan_kwargs={'lr':5e-3}, check_val_every_n_epoch=5)
        z = model.get_latent_representation()
        adata.obsm[z_name] = z
        adata.layers["scvi_normalized" + str(i)] = model.get_normalized_expression(library_size=10e4)
        latents.append(z)
    return (adata,latents)
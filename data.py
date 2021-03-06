import os
import torch
import numpy as np
import pandas as pd
import scanpy as sc

def get_data(dir_name='./input/',verbosity=3):
    #os.mkdir("data")
    #os.chdir("data")
    os.mkdir("write")
    torch.random.manual_seed(42)
    sc.settings.verbosity = verbosity             # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=120, facecolor='white', dpi_save=300, vector_friendly=True)
    sc._settings.ScanpyConfig(autosave=True)

    adata = sc.read_10x_mtx(
        dir_name,  # the directory with the `.mtx` file
        var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
        cache=False)  
    
    adata.var_names_make_unique()
    return adata
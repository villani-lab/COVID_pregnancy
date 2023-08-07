import anndata as an
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(dpi=300, dpi_save=300, vector_friendly=True)
from matplotlib.backends.backend_pdf import PdfPages

exec(open("functions.py").read())


# Supp Figure 8h - T+NK interim object - UMAP
data = an.read_h5ad("T_NK_interim.h5ad")
sc.set_figure_params(vector_friendly=True, dpi_save=300, figsize=(5,4))
g = sc.pl.umap(data, color = "leiden_res_1.8", return_fig=True, show=False, size=5)
g.figure.tight_layout()
g.figure.savefig("CD4_CD8_split_umap.pdf")

# Supp Figure 8f - CD8T Feature plots for doublets
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad")
feature_plots_2(cd8_data, ["CD4","CD8A","CD8B","S100A9","S100A8","VCAN"], title="doublets", dot_size=10, figsize=(5,6))
plt.clf()
sc.set_figure_params(figsize=(3,3))
g= sc.pl.umap(cd8_data, color= 'percent_mito', cmap="YlOrRd", size=10, return_fig=True, show=False)
g.figure.tight_layout()
g.figure.savefig("percent_mito.pdf")





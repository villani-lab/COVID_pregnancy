import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc
import pandas as pd
import anndata as an

exec(open("functions.py").read())


# Reading all anndata objects
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
cd4_dat = an.read_h5ad("CD4T_GEX_TCR.h5ad") # From GSE239452
myeloid_data = an.read_h5ad("Myeloid.h5ad") # From GSE239452
B_data = an.read_h5ad("B_Plasma_GEX_BCR.h5ad") # From GSE239452

# UMAP per lineage - Figure 6A, C, E, G
sc.pl.umap(cd8_data, color = 'cluster_title')
sc.pl.umap(cd4_dat, color = 'cluster_title')
sc.pl.umap(myeloid_data, color = 'cluster_title')
sc.pl.umap(B_data, color = 'cluster_title')


# ISG signature per lineage - Figure 6B, D, F, Supp Figure 11C
# For every lineage we first compute Z-scores across all cells in the anndata object and then compute mean z-score for ISG gene set.
ISGs = pd.read_excel("Supp Table 2 - ISGs.xlsx", header=None)[0].to_list() # Supplementary table 2

# CD8 T cells
signature_score_per_cell_2(data=cd8_data, gene_set=ISGs, score_title="ISG_score")
plt.clf()
with PdfPages("ISG_score_CD8T.pdf") as pdf:
    plot_ISG_score_boxplot(data=cd8_data, cluster_col="cluster_title", title="CD8T", pdf=pdf, doublet_clusters=['Doublets (CD4/CD8)','High-MT','Doublets (Myeloid/CD8)'], shorten=True)

# Myeloid cells
signature_score_per_cell_2(data=myeloid_data, gene_set=ISGs, score_title="ISG_score")
plt.clf()
with PdfPages("ISG_score_Myeloid.pdf") as pdf:
    plot_ISG_score_boxplot(data=myeloid_data, cluster_col="cluster_title", title="Myeloid", pdf=pdf, doublet_clusters=['Doublets (B/myeloid)','Doublets (T/myeloid)',"Doublets (platelets/ Myeloid)"])

# CD4 T cells
signature_score_per_cell_2(data=cd4_dat, gene_set=ISGs, score_title="ISG_score")
plt.clf()
with PdfPages("ISG_score_CD4T.pdf") as pdf:
    plot_ISG_score_boxplot(data=cd4_dat, cluster_col="cluster_title", title="CD4T", pdf=pdf, doublet_clusters=["high MT/ribo"])

# B cells
signature_score_per_cell_2(data=B_data, gene_set=ISGs, score_title="ISG_score")
plt.clf()
with PdfPages("ISG_score_B.pdf") as pdf:
    plot_ISG_score_boxplot(data=B_data, cluster_col="cluster_title", title="B_Plasma", pdf=pdf, doublet_clusters = ["Doublets (platelets/B)","high_ribo","high_mito","Doublets (T/B)"])


# Figure 6H - see SuppFile_11.py file
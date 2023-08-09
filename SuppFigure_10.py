import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc
import seaborn as sns
import anndata as an

exec(open("functions.py").read())



# Supp Figure 10A - Myeloid UMAP
myeloid_data = an.read_h5ad("Myeloid.h5ad") # From GSE239452
plt.clf()
with PdfPages("Myeloid_UMAP.pdf") as pdf:
    sc.set_figure_params(figsize=(8,6))
    g = sc.pl.umap(myeloid_data, color='cluster_title', show=False, return_fig=True, size=5)
    plt.legend(ncol=1, loc='upper left', bbox_to_anchor=(1, 1))
    g.figure.tight_layout()
    pdf.savefig(g)

# Supp Figure 10A - Myeloid Heatmap
genes1 = ["CD14","S100A8",	"S100A9","LYZ","VCAN", # MNP1
          "BHLHE40","METRNL","ID1","THBD", # MNP2
          "TMEM176B", "TMEM176A", "FOS", # MNP3
          "CXCL10","IFIT2","TNFSF10", # MNP4
          "FCGR3A","ITGAL","CDKN1C","LST1", # MNP5
          "FCGR3A","C1QA","C1QB","C1QC", # MNP6
          "CLEC10A","CD1C","CD74","PLD3" ,# cDC
          "CLEC4C","TCF4","IRF8"] # pDC
# Removing doublets for heatmap
myeloid_data = myeloid_data[[x not in ["Doublets (B/myeloid)","Doublets (T/myeloid)","Doublets (platelets/ Myeloid)"] for x in myeloid_data.obs["cluster_title"]]].copy()
df1= myeloid_data[:,genes1].to_df().groupby(myeloid_data.obs['cluster_title']).mean()

sns.set(font_scale = 0.5)
plt.clf()
g=sns.clustermap(df1, standard_scale=1, cmap="YlOrRd", col_cluster=False, row_cluster=False,xticklabels=True, yticklabels=True)
plt.xticks(rotation=90)
plt.yticks(rotation=90)
g.fig.set_size_inches(5,3)
plt.tight_layout()
g.figure.savefig("myeloid_markers_heatmap.pdf")


# Supp Figure 10B - marker genes for doublets and platelets
myeloid_data = an.read_h5ad("Myeloid.h5ad") # From GSE239452
sc.set_figure_params(transparent =True)
feature_plots_2(myeloid_data, gene_list = ["CD3E","CD3D","MS4A1","CD79A","PPBP","PF4"], title="myeloid", dot_size=3)


# Supp Figure 10C,D - cell subset abundances
myeloid_data = an.read_h5ad("Myeloid.h5ad") # From GSE239452

with PdfPages("Myeloid_abundances.pdf") as pdf:
    compute_DA_and_plot_per_category(data=myeloid_data, category="Category_2", groupby='cluster_title', pdf=pdf,figsize=(7, 4), title="Myeloid",omit_clusters=["Doublets (T/myeloid)", "Doublets (B/myeloid)","Doublets (platelets/ Myeloid)"])
    compute_DA_and_plot_per_category(data=myeloid_data, category="Category", groupby='cluster_title', pdf=pdf,figsize=(5, 4), title="Myeloid",omit_clusters=["Doublets (T/myeloid)", "Doublets (B/myeloid)","Doublets (platelets/ Myeloid)"])
    plot_per_category_separate_scale(abundances_file="Myeloid_DA.csv", groupby='cluster_title', title="Myeloid", pdf=pdf, figsize=(5, 4))


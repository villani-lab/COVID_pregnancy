import anndata as an
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(dpi=300, dpi_save=300, vector_friendly=True)
from matplotlib.backends.backend_pdf import PdfPages



# Figure 5a - CD8T UMAP
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad")
colors = list(matplotlib.colormaps.get("tab10").colors + matplotlib.colormaps.get("Accent").colors)
colors =[colors[x] for x in [0,1,2,3,4,5,6,15,8,9,10,11,16,13,12,7]]# rearraning colors to match patterns in UMAP
ordered_clusters = cd8_data.obs["cluster_title"].unique().categories.to_list()
cd8_colors = dict(zip(ordered_clusters, list(colors)))

with PdfPages("CD8T_UMAP.pdf") as pdf:
    sc.set_figure_params(figsize=(11,6))
    g = sc.pl.umap(cd8_data, color='cluster_title', show=False, return_fig=True, palette = cd8_colors, size=5)
    plt.legend(ncol=1, loc='upper left', bbox_to_anchor=(1, 1))
    g.figure.tight_layout()
    pdf.savefig(g)

# Figure 5a - CD8T heatmap
genes2 = ['CD4','CD40LG','ITGB1', #c4
          'GZMH','NKG7',"CST7", #c3
          'KLRD1','FCGR3A','GZMB', #c2
          "CXCR4","TGFB1","NR4A2", #c10
          "KLRF1","FCER1G","STMN1", # c15
          "CXCR3", "HMGB2", "PRDX3", "TALDO1", #c9
          "CMC1", "CD74", "GZMK", #c1
          "KLRB1", "IL7R", "CXCR6", # c12
          "TNFAIP3", #c13
          "TRAV1-2", "IL23R", "AQP3", #c14
          "LTB", "TCF7",	"TIGIT", #c6
          "SELL","CCR7", # c5
          "MX1", "ISG15" ] # c11
df= cd8_data[:,genes2].to_df().groupby(cd8_data.obs['final_res_1.8']).mean()

sns.set(font_scale = 0.5)
plt.clf()
g=sns.clustermap(df, standard_scale=1, cmap="YlOrRd", col_cluster=False, row_cluster=False,xticklabels=True, yticklabels=True)
plt.xticks(rotation=90)
g.fig.set_size_inches(5,3)
g.figure.tight_layout()
plt.tight_layout()
g.figure.savefig("CD8_markers_heatmap.pdf")


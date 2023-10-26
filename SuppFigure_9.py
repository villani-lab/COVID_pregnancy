import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc
import seaborn as sns
import anndata as an
sc.set_figure_params(dpi=300, dpi_save=300, vector_friendly=True)

exec(open("functions.py").read())

B_data = an.read_h5ad("B_Plasma_GEX_BCR.h5ad") # From GSE239452

# Supp Figure 9A - UMAP
plt.clf()
with PdfPages("B_Plasma_UMAP.pdf") as pdf:
    sc.set_figure_params(figsize=(6,4))
    g = sc.pl.umap(B_data, color='cluster_title', show=False, return_fig=True, size=5)
    plt.legend(ncol=1, loc='upper left', bbox_to_anchor=(1, 1))
    g.figure.tight_layout()
    pdf.savefig(g)


# Supp Figure 9A - heatmap
genes =["CD27","IGHG1","FCER2","TCL1A", "IGHD","IGHM", # all
        "IL4R", # B_1 - naive
        "IFITM1", "IFI44L", "IRF7", # B_3 - naive ISG
        "PLD4","SOX4","CD9", # B_6 - immature
        "IGLC7", "IGLC6", # B_8
        "CD83", "SNX9", # B_2 - early activated
        "CD83","GRASP","PELI1","COTL1", # B_7
        "CRIP1", "ITGB1", # B_4
        "FGR", "CIB1", "ITGAX", # B_7
        "JCHAIN", "MZB1", "CD38"] # Plasmablasts
        #"NRGN", "PPBP", "SPARC" #c10   ]

B_data = B_data[[x not in ["Doublets (platelets/B)","high_ribo","high_mito","Doublets (T/B)"] for x in B_data.obs["cluster_title"]],:].copy() # removing doublets and high-MT for heatmap
df= B_data[:,genes].to_df().groupby(B_data.obs['cluster_title']).mean()

sns.set(font_scale = 0.5)
plt.clf()
g=sns.clustermap(df, standard_scale=1, cmap="YlOrRd", col_cluster=False, row_cluster=False,xticklabels=True, yticklabels=True)
plt.xticks(rotation=90)
plt.yticks(rotation=90)
g.fig.set_size_inches(5,3)
plt.tight_layout()
g.figure.savefig("B_markers_heatmap.pdf")



# Supp Figure 9B, Figure 4H - B cell subset abundances
with PdfPages("B_plasma_abundances.pdf") as pdf:
    compute_DA_and_plot_per_category(B_data, "Category_2", groupby='cluster_title', pdf=pdf, figsize=(7, 4), title="B+Plasma cells", omit_clusters = ["Doublets (platelets/B)","Doublets (T/B)","high_mito","high_ribo"])
    compute_DA_and_plot_per_category(B_data, "Category", groupby='cluster_title', pdf=pdf, figsize=(5, 4), title="B+Plasma cells", omit_clusters = ["Doublets (platelets/B)","Doublets (T/B)","high_mito","high_ribo"])


# Supp Figure 9D - BCR data projected on a UMAP
B_data.obs.loc[B_data.obs["has_kap"]=="","has_kap"]='False' # correcting 1 cell
B_data.obs["has_kap"] =B_data.obs["has_kap"].cat.remove_unused_categories()
B_data.uns["has_lam_colors"]=["lightgrey","red"]
B_data.uns["has_kap_colors"]=["lightgrey","red"]
B_data.uns["has_heavy_colors"]=["lightgrey","red"]
B_data.uns["has_bcr_colors"]=["lightgrey","red"]
fig = sc.pl.umap(B_data, color=['has_lam', 'has_kap', 'has_heavy', 'has_bcr'], return_fig=True, show=False, size=4)
fig.savefig("has_BCR_umap.pdf")


# Supp Figure 9E - BCR data per cell subset
cluster_col = 'cluster_title'
df_lam = get_perc_out_of_cluster(B_data, cluster_col=cluster_col, feature="has_lam")
df_kap = get_perc_out_of_cluster(B_data, cluster_col=cluster_col, feature="has_kap")
df_heavy = get_perc_out_of_cluster(B_data, cluster_col=cluster_col, feature="has_heavy")
df_bcr = get_perc_out_of_cluster(B_data, cluster_col=cluster_col, feature="has_bcr")

df = df_lam[[cluster_col, "total_has_lam", "prop_has_lam"]].merge(df_kap[[cluster_col, "total_has_kap", "prop_has_kap"]], on=cluster_col).merge(df_heavy[[cluster_col, "total_has_heavy", "prop_has_heavy"]], on=cluster_col).merge(df_bcr[[cluster_col, "total_has_bcr", "prop_has_bcr"]], on=cluster_col)
df = df.loc[[x not in ["3","5","10","12"] for x in df[cluster_col]],:].copy()
df = df[[cluster_col,"prop_has_kap","prop_has_lam","prop_has_heavy","prop_has_bcr"]].melt(id_vars=cluster_col,value_name="Percent",var_name="Feature")
df= df.replace({"prop_has_kap":"Kappa", "prop_has_lam":"Lambda","prop_has_heavy":"Heavy","prop_has_bcr":"BCR"})
df[cluster_col] = df[cluster_col].astype(str)

plt.clf()
plt.figure(figsize=(5,5))
sns.set_style("white")
g = sns.barplot(data =df, y=cluster_col, x="Percent", hue="Feature", orient='h')
plt.tight_layout()
g.figure.savefig("has_BCR_barplot.pdf")



# Supp Figure 9F - BCR clone frequency projected over a UMAP
clone_prop_umap_per_category(B_data, pdf_prefix="BCR_clone_prop_umap")



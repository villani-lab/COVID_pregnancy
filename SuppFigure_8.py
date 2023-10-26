import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc
import seaborn as sns
import anndata as an

exec(open("functions.py").read())
cd4_data = an.read_h5ad("CD4T_GEX_TCR.h5ad") # From GSE239452

# Supp Figure 8A UMAP
plt.clf()
with PdfPages("CD4T_UMAP.pdf") as pdf:
    sc.set_figure_params(figsize=(8,6))
    g = sc.pl.umap(cd4_data, color='cluster_title', show=False, return_fig=True, size=5)
    plt.legend(ncol=1, loc='upper left', bbox_to_anchor=(1, 1))
    g.figure.tight_layout()
    pdf.savefig(g)

# Supp Figure 8A - heatmaqp
genes =["CCR7", "FHIT","SELL", #c2
"KRT1",  # c3
"FOXP3", "TIGIT", "IKZF2", "CTLA4","IL2RA", #c9
"CD8A","CD8B", "GATA3", #c8
"AIRE", "CCR10", "LMNA", #c5
"IL7R", "TNFSF13B","AQP3", "KLRB1", #c4
"GZMK" ,"TNFAIP3","CCL5",  #c1
"TRDV2", "TRGV9", "PRF1"] #c6

cd4_data = cd4_data[[x not in ["high MT/ribo"] for x in cd4_data.obs["cluster_title"]],:].copy() # removing high-MT for heatmap
df= cd4_data[:,genes].to_df().groupby(cd4_data.obs['cluster_title']).mean()
sns.set(font_scale = 0.5)
plt.clf()
g=sns.clustermap(df, standard_scale=1, cmap="YlOrRd", col_cluster=False, row_cluster=False,xticklabels=True, yticklabels=True)
plt.xticks(rotation=90)
plt.yticks(rotation=90)
g.fig.set_size_inches(5,3)
plt.tight_layout()
g.figure.savefig("CD4T_markers_heatmap.pdf")


# Supp Figure 8B,C - cell subset abundances
with PdfPages("CD4T_abundances.pdf") as pdf:
    compute_DA_and_plot_per_category(cd4_data, "Category_2", groupby='cluster_title', pdf=pdf, figsize=(7,4), title="CD4 T cells",omit_clusters = ["high MT/ribo"])
    compute_DA_and_plot_per_category(cd4_data, "Category", groupby='cluster_title', pdf=pdf, figsize=(5,4), title="CD4 T cells",omit_clusters = ["high MT/ribo"])


# Supp Fig 8D - TCR data projected on UMAP
cd4_data.uns["has_alpha_colors"]=["lightgrey","red"]
cd4_data.uns["has_beta_colors"]=["lightgrey","red"]
cd4_data.uns["has_tcr_colors"]=["lightgrey","red"]
fig = sc.pl.umap(cd4_data, color=['has_alpha', 'has_beta', 'has_tcr'], return_fig=True, show=False, size=4)
fig.savefig("has_TCR_CD4_umap.pdf")



# Supp Fig 8E - has_tcr in barplot
cluster_col = "cluster_title"
df_alpha = get_perc_out_of_cluster(cd4_data, cluster_col=cluster_col, feature="has_alpha")
df_beta = get_perc_out_of_cluster(cd4_data, cluster_col=cluster_col, feature="has_beta")
df_tcr = get_perc_out_of_cluster(cd4_data, cluster_col=cluster_col, feature="has_tcr")

df = df_alpha[[cluster_col, "total_has_alpha", "prop_has_alpha"]].merge(df_beta[[cluster_col, "total_has_beta", "prop_has_beta"]], on=cluster_col).merge(df_tcr[[cluster_col, "total_has_tcr", "prop_has_tcr"]], on=cluster_col)
df = df.loc[[x not in ['high MT/ribo'] for x in df[cluster_col]],:].copy()
df = df[[cluster_col,"prop_has_alpha","prop_has_beta","prop_has_tcr"]].melt(id_vars=cluster_col,value_name="Percent",var_name="Feature")
#df= df.replace({"prop_has_kap":"Kappa", "prop_has_lam":"Lambda","prop_has_heavy":"Heavy","prop_has_bcr":"BCR"})
df[cluster_col] = df[cluster_col].astype(str)

plt.clf()
plt.figure(figsize=(5,5))
sns.set_style("white")
g = sns.barplot(data =df, y=cluster_col, x="Percent", hue="Feature", orient='h')
plt.tight_layout()
g.figure.savefig("has_CD4_TCR_barplot.pdf")


# Supp Fig 8F - CD4 TCR clonal expansion on UMAP
clone_prop_umap_per_category(data, pdf_prefix="CD4_TCR_clone_prop_umap")
clone_prop_umap_per_sample(data, pdf_prefix="CD4_TCR_clone_prop_umap_sample")
clone_prop_umap_per_sample(data, pdf_prefix="CD4_TCR_clone_prop_umap_sample_no_vmax", vmax=False)


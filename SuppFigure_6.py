import anndata as an
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import matplotlib
import seaborn as sns
matplotlib.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(dpi=300, dpi_save=300, vector_friendly=True)
from matplotlib.backends.backend_pdf import PdfPages

exec(open("functions.py").read())


# Supp Figure 6A - all cells UMAP
data = an.read_h5ad("AllCells_GEX.h5ad") # From GSE239452
sc.pl.umap(data, color="Lineage")


# Supp Figure 6B - all cells lineage composition per patient
data = an.read_h5ad("AllCells_GEX.h5ad") # From GSE239452
df = pd.crosstab(data.obs.loc[:, "sample"], data.obs.loc[:, "Lineage"])
df = df.div(df.sum(axis=1), axis=0) * 100.0
# Sample order - by Category_2, then by T cekk frequency
order_df = data.obs[["sample","Category_2"]].drop_duplicates().reset_index(drop=True)
order_df["Category_2"] = pd.Categorical(order_df["Category_2"], ordered=True, categories=["NonPreg_ASX","NonPreg_SEV","Preg_CTRL","Preg_ASX","Preg_SEV"])
df = df.merge(order_df, right_on="sample",left_index=True, how="left")
df = df.sort_values(by=["Category_2","T-NK"]).set_index("sample")
plt.clf()
fig, axes = plt.subplots(1,2,sharey=True,gridspec_kw={'width_ratios': [8, 2]}, figsize=(5, 5))
df.plot(kind="barh", stacked=True, legend=True, grid=False, fontsize=10, width=1,style={'border': None}, ax=axes[0])#, color=color_dict)
axes[0].set_xticks([0,25,50,75,100],fontsize=8)
axes[0].legend(loc='lower center', bbox_to_anchor=(1, 1), fontsize=8)
counts = data.obs["sample"].value_counts().reset_index().rename(columns={"index":"sample", "sample":"cell_no"})
counts["sample"] = pd.Categorical(counts["sample"], ordered=True, categories=df.index.to_list())
counts = counts.sort_values(by="sample")
counts.plot(kind="barh", y="cell_no", x="sample", stacked=True, legend=True, grid=False, fontsize=10, width=1, color="grey", ax=axes[1], )
plt.tight_layout()
plt.show()


# Supp Figure 6C - lineage marker feature plots
data = an.read_h5ad("AllCells_GEX.h5ad") # From GSE239452
feature_plots_2(data, ["CD3D","CD8A","CD4","NCAM1","MS4A1","LYZ"], title="doublets", dot_size=10, figsize=(5,6))


# Supp Figure 6D,E - all cells differential analysis
data = an.read_h5ad("AllCells_GEX.h5ad") # From GSE239452
with PdfPages("GlobalUMAP_DA.pdf") as pdf:
    compute_DA_and_plot_per_category(data, category="Category", groupby="Lineage", pdf=pdf, title="Global UMAP", figsize=(5,4))
    compute_DA_and_plot_per_category(data, category="Category_2", groupby="Lineage", pdf=pdf, title="Global UMAP")


# Supp Figure 6F - CD8T Feature plots for doublets
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
feature_plots_2(cd8_data, ["CD4","CD8A","CD8B","S100A9","S100A8","VCAN"], title="doublets", dot_size=10, figsize=(5,6))
plt.clf()
sc.set_figure_params(figsize=(3,3))
g= sc.pl.umap(cd8_data, color= 'percent_mito', cmap="YlOrRd", size=10, return_fig=True, show=False)
g.figure.tight_layout()
g.figure.savefig("percent_mito.pdf")

# Supp Figure 6G - Expression of CD4 and CD8 in various cell subsets
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
cd4_data = an.read_h5ad("CD4T_GEX_TCR.h5ad") # From GSE239452
cluster_CD8T_3 = cd8_data[cd8_data.obs["cluster_title"]=='CD8T_3: CD4, CD40LG, ITGB1, GZMH','CD8T_7: CXCR4, TGFB1, NR4A2',:].copy()
#cluster_CD8T_7 = cd8_data[cd8_data.obs["cluster_title"]=='CD8T_7: CXCR4, TGFB1, NR4A2',:].copy()
cd8_pos = cd8_data[cd8_data.obs["cluster_title"]=='CD8T_2: NKG7, GZMH, CST7',:].copy()
cd4_pos = cd4_data[cd4_data.obs['leiden_res_1']=='5',:].copy()
clear_data_2(cluster_CD8T_3, ["sample"])
#clear_data_2(cluster_CD8T_7, ["sample"])
clear_data_2(cd8_pos, ["sample"])
clear_data_2(cd4_pos, ["sample"])
#x = cluster_CD8T_3.concatenate(cluster_CD8T_7, cd8_pos, cd4_pos, batch_categories = ["CD8T_3","CD8T_7","CD8T_2","CD4_5"], join="inner")
x = cluster_CD8T_3.concatenate(cd8_pos, cd4_pos, batch_categories = ["CD8T_3","CD8T_7","CD8T_2","CD4_5"], join="inner")
x.obs["group"] = x.obs["sample"].astype(str)+"$"+x.obs["batch"].astype(str)
x2 = x.to_df().groupby(x.obs["group"]).mean()
x3 = x2[["CD4","CD8A","CD8B"]].reset_index().melt(id_vars="group")
x3["sample"] = [x.split("$")[0] for x in x3["group"]]
x3["batch"] = [x.split("$")[1] for x in x3["group"]]
plt.clf()
fig,ax = plt.subplots(figsize=(5,5))
sns.boxplot(data = x3, x="featurekey",y="value", hue="batch",ax=ax)
sns.stripplot(data = x3, x="featurekey",y="value", hue="batch", dodge=True, palette=["black","black","black"], size=2)
plt.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.show()
fig.savefig("CD8T7_CD8T3_boxplots.pdf")


# Supp Figure 6H - T+NK interim object - UMAP
data = an.read_h5ad("T_NK_interim.h5ad") # From GSE239452
sc.set_figure_params(vector_friendly=True, dpi_save=300, figsize=(5,4))
g = sc.pl.umap(data, color = "leiden_res_1.8", return_fig=True, show=False, size=5)
g.figure.tight_layout()
g.figure.savefig("CD4_CD8_split_umap.pdf")


# Supp Figure 6I - T+NK interim object - marker genes feature plots
data = an.read_h5ad("T_NK_interim.h5ad") # From GSE239452
feature_plots_2(data, ["CD3D","CD3G","CD8A","CD8B","CD4","NCAM1"], title="T_NK", dot_size=10, figsize=(5,6))


# Supp Figure 6J,K - differential abundance CD8 T cells
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
with PdfPages("CD8T_DA.pdf") as pdf:
    compute_DA_and_plot_per_category(cd8_data, "Category", groupby="cluster_title", pdf=pdf, figsize=(5,4), title="CD8 T cells", omit_clusters = ["CD8T_7","CD8T_8","CD8T_16"])
    compute_DA_and_plot_per_category(cd8_data, "Category_2", groupby="cluster_title", pdf=pdf, figsize=(7,4), title="CD8 T cells", omit_clusters = ["CD8T_7","CD8T_8","CD8T_16"])


import anndata as an
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(dpi=300, dpi_save=300, vector_friendly=True)
from matplotlib.backends.backend_pdf import PdfPages

exec(open("/projects/COVID_pregnancy/Code_publication/functions.py").read())
##############

exec(open("functions.py").read())

# Supp Figure 9A - where we have CD8 TCR data
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
cd8_data.uns["has_alpha_colors"]=["lightgrey","red"]
cd8_data.uns["has_beta_colors"]=["lightgrey","red"]
cd8_data.uns["has_tcr_colors"]=["lightgrey","red"]
fig = sc.pl.umap(cd8_data, color=["cluster_title","has_alpha","has_beta","has_tcr"], return_fig=True, show=False, size=4)
fig.savefig("has_tcr_umap.pdf")

# Supp Figure 9B - TCR data per cluster
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
df_alpha = get_perc_out_of_cluster(cd8_data, cluster_col="cluster_title", feature="has_alpha")
df_beta = get_perc_out_of_cluster(cd8_data, cluster_col="cluster_title", feature="has_beta")
df_tcr = get_perc_out_of_cluster(cd8_data, cluster_col="cluster_title", feature="has_tcr")
df = df_alpha.merge(df_beta, on="cluster_title").merge(df_tcr, on="cluster_title")
df = df.loc[[x not in ["CD8T_7","CD8T_8","CD8T_16"] for x in df["cluster_title"]],:].copy()
df = df[["cluster_title","prop_has_alpha","prop_has_beta","prop_has_tcr"]].melt(id_vars="cluster_title",value_name="Percent",var_name="Feature")
df= df.replace({"prop_has_beta":"Beta", "prop_has_alpha":"Alpha","prop_has_tcr":"TCR"})
df["cluster_title"] = df["cluster_title"].astype(str)

plt.clf()
plt.figure(figsize=(5,5))
sns.set_style("white")
g = sns.barplot(data =df, y="cluster_title", x="Percent", hue="Feature", orient='h')
plt.tight_layout()
g.figure.savefig("has_tcr_barplot.pdf")


# Supp figure 9 C,D - Where are the expanded clones?
#---------------------------------------------
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
cd8_data.obs["expanded_t"] = cd8_data.obs["expanded_t"].map({False:"False", True:"True"})
cluster_col="cluster_title"
with PdfPages("expnaded_clones_umap_bar_plots.pdf") as pdf:
    g=sc.pl.umap(cd8_data, color="expanded_t", palette={"False":"lightgrey", "True":"red"}, return_fig=True, show=False)
    pdf.savefig(g)

    # Hoa many exoanded clones do we have per cluster?
    expanded_data = cd8_data[cd8_data.obs["expanded_t"]=="True",:].copy()
    df = expanded_data.obs[[cluster_col, "expanded_t"]].value_counts().reset_index().rename(columns = {0:"expanded_count"})
    has_tcr_data = cd8_data[cd8_data.obs["has_tcr"]=="True",:].copy()
    df2 = has_tcr_data.obs[[cluster_col,'has_tcr']].value_counts().reset_index().rename(columns = {0:"total_count"})
    df = df[[cluster_col,"expanded_count"]].merge(df2[[cluster_col,"total_count"]], on=cluster_col)
    df["perc_expanded"] = (df["expanded_count"]*100)/df["total_count"]
    #df[cluster_col] = pd.Categorical(df[cluster_col], ordered=True, categories=sorted(list(set(df[cluster_col])), key=int))
    df= df.loc[[x not in ['CD8T_7','CD8T_8','CD8T_16'] for x in df[cluster_col]],:].copy()
    df[cluster_col] = df[cluster_col].cat.remove_unused_categories()

    plt.clf()
    ax = sns.barplot(data =df, x=cluster_col, y="perc_expanded", color="grey")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=8)
    pdf.savefig(ax.figure)

# Supp Figure 9E - Clonal expansion projected on a UMAP per sample
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
clone_prop_umap_per_sample(cd8_data, pdf_prefix="TCR_clone_prop_umap_sample")


# Supp Figure 9F,G - percent of expanded clones per cluster - 5 categories and 2 categories
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
cd8_data.obs["cluster_title"] = [x.split(":")[0] for x in cd8_data.obs["cluster_title"]]
cd8_data.obs["cluster_title"]  = pd.Categorical(cd8_data.obs["cluster_title"] )
cd8_data_wo_doublets = cd8_data[[x not in ["CD8T_7","CD8T_8","CD8T_16"] for x in cd8_data.obs["cluster_title"]],:].copy()

with PdfPages("supp_barplot_expanded_TCRs_5categories.pdf") as pdf:
    plot_expanded_clones_per_clust(data = cd8_data_wo_doublets, pdf=pdf, figsize=(15,4), hue_order=["NonPreg_ASX","NonPreg_SEV", "Preg_CTRL", "Preg_ASX" ,"Preg_SEV"], category="Category_2", pairs=["Preg_SEV", "NonPreg_SEV"], colors = Category_2_colors, cluster_col='cluster_title')

samp_category2 = cd8_data_wo_doublets.obs[["sample","Category_2"]].drop_duplicates().reset_index(drop=True)
covid_samples =  list(samp_category2.loc[[x != "Preg_CTRL" for x in samp_category2["Category_2"]],"sample"])
cd8_data_wo_doublets = cd8_data_wo_doublets[[x in covid_samples for x in cd8_data_wo_doublets.obs["sample"]],:].copy()
cd8_data_wo_doublets.obs["Category"] = cd8_data_wo_doublets.obs["Category_2"].map({"Preg_SEV":"Preg_COVID","Preg_ASX":"Preg_COVID","NonPreg_SEV":"NonPreg_COVID","NonPreg_ASX":"NonPreg_COVID"})
with PdfPages("barplot_expanded_TCRs_2categories.pdf") as pdf:
    plot_expanded_clones_per_clust(data = cd8_data_wo_doublets, pdf=pdf, figsize=(15,4), hue_order=["NonPreg_COVID", "Preg_COVID"], category="Category", pairs=["Preg_COVID", "NonPreg_COVID"], colors ={"Preg_COVID":"#EE6677", "NonPreg_COVID":"#4477AA"}, cluster_col='cluster_title')





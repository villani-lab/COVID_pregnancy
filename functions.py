import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc
import seaborn as sns
import pandas as pd
import numpy as np
from statannotations.Annotator import Annotator


Category_2_colors = {"NonPreg_ASX":"lightsteelblue","NonPreg_SEV":"royalblue", "Preg_CTRL":"mistyrose", "Preg_ASX":"lightcoral" ,"Preg_SEV":"firebrick"}
Category_colors = {"Preg_COVID":"#EE6677", "NonPreg_COVID":"#4477AA"}


# A function to plot feature plots on a UMAP - expression of a list of genes/ proteins projected on a UMAP
# adata - anndata object with all the data
# gene_list - a list of genes/ protein that appear in adata.var_names to plot on the UMAP
# title - the title of the main figure and the prefix of the PDF file that will be created in case pdf==True
# uae_raw - a boolean parameter indicating whether to use raw values
# hex - a boolean parameter indicating wither to bin the expression values rather than plot individual cells
# gridsize - number of hexagons on the x axis - the larger this number the smaller the hexagons will be
# umap - X_<umap> is the attribute name in adata.obsm that will be the basis for the UMAP
# pdf - a boolean parameter indicating whether to plot the figure or to save it in a PDF file.
def feature_plots_2(adata, gene_list, title, use_raw=False, hex=False, gridsize=100, umap="umap", pdf=True, figsize=(10, 7), dot_size=1, axes_titlesize= 35):
    genes = [g for g in gene_list if g in adata.var_names]
    print(str(len(genes)) + " out of " + str(len(gene_list)) + " were measured")
    page_no = math.floor(len(genes)/16)
    if(len(genes)%16 != 0):
        page_no+=1
    rows=4

    if (umap != "umap"):
        adata = adata.copy()
        adata.obsm["X_umap"] = adata.obsm[
            "X_" + umap]  # sc.pl.scatter doesn't have all the feature os sc.pl.umap to plot correctly

    def plot(dot_size):
        fig_list = []
        count = 0
        for page in range(page_no):
            print("page: " + str(page))
            if (hex):
                x = adata.obsm["X_umap"][:, 0]
                y = adata.obsm["X_umap"][:, 1]
                fig, axs = plt.subplots(ncols=4, nrows=4, figsize=figsize)
                for i in range(4):
                    for j in range(4):
                        if (count > (len(genes) - 1)):
                            break
                        gene = genes[count]
                        count += 1
                        ax = axs[i][j]
                        hb = ax.hexbin(x, y, C=adata[:, gene].to_df()[gene], cmap='YlOrRd', gridsize=gridsize)
                        ax.set_title(gene, fontsize=11)
                        ax.set_yticklabels([])
                        ax.set_xticklabels([])
                        ax.set_ylabel([])
                        ax.set_xlabel([])
                        ax.set_rasterization_zorder(10)

                        ax.tick_params(labelsize=2)
                        fig.colorbar(hb, ax=ax)
            else:
                parameters = {'axes.labelsize': 0, 'axes.titlesize': 55, 'legend.fontsize': 40, 'axes.grid': False,'legend.fontsize': "medium"}
                plt.rcParams.update(parameters)
                fig = sc.pl.umap(adata, color=genes[16 * page:16 * (page + 1)], size=dot_size, return_fig=True, show=False,cmap="YlOrRd", use_raw=use_raw)
            fig_list.append(fig)
        return(fig_list)

    if (pdf):
        if (dot_size==""):
            dot_size=25

        with PdfPages(title.replace(" ","_")+"_feature_plots.pdf") as pdf:
            fig_list = plot(dot_size=dot_size)
            for fig in fig_list:
                fig.suptitle(title, fontsize=axes_titlesize)
                pdf.savefig(fig, bbox_inches="tight")
    else:
        if (dot_size==""):
            dot_size=20
        fig_list = plot(dot_size=dot_size)
        for fig in fig_list:
            plt.suptitle(title, fontsize=rows*15, y=1)
            plt.subplots_adjust(hspace=0.8)
            fig.show()



# data - anndata object
# category - the category in obs to compute the differential analysis in - either "Category" which consists of the categories 'Preg_COVID'/ 'NonPreg_COVID', or "Category_2 that includes the categories 'NonPreg_ASX', 'NonPreg_SEV', 'Preg_ASX', 'Preg_CTRL', 'Preg_SEV'
# pdf - pdf object to print the plots into
# figsize - the figure size
# omit_clusters - clusters to omit from plot and multiple hypothesis testing
# focus_clusters - focusing on specific clusters to plot
def compute_DA_and_plot_per_category(data, category, groupby, pdf, title, figsize=(8,4), omit_clusters=[], focus_clusters=[]):
    if ((category != "Category") & (category != "Category_2")):
        print("Error: 'category' should either be 'Category' or 'Category_2'")
        return ()

    data_filtered = data[[x not in omit_clusters for x in data.obs[groupby]],:].copy()
    samples = list(set(data_filtered.obs["sample"]))

    # df.groupby(['Symbol', 'Year']).count().unstack(fill_value=0).stack() - another way to count without a function (from Stack Overflow)
    # This function counts cells per cluster for a piece of the dataframe fromd ata.obs that includes only that cluster, fills with 0s values of sampels with no counts
    def count_df(x):
        return (pd.DataFrame(x["sample"].value_counts().reindex(samples, fill_value=0)))
    df_init = data_filtered.obs.groupby(by =groupby).apply(count_df).reset_index()
    df_init = df_init.rename(columns={"level_1":"sample", "sample":"Count"})
    df_init = df_init.merge(data_filtered.obs[["sample","Category","Category_2"]].drop_duplicates().reset_index(drop=True), on="sample", how="left")

    if (category == "Category"):
        df_init = df_init.loc[df_init["Category_2"] != "Preg_CTRL", :].copy()
        df = df_init[["sample", groupby, "Category", "Category_2", "Count"]].copy()
        df["Category"] = df["Category"].map({'NonPregnant': "NonPreg_COVID", 'Pregnant': "Preg_COVID",'NonPreg': "NonPreg_COVID"})
        df["Category"] = pd.Categorical(df["Category"], ordered=True, categories=["NonPreg_COVID", "Preg_COVID"])
        curr_colors = Category_colors
    else:
        df = df_init[["sample", groupby, "Category", "Category_2", "Count"]].copy()
        df["Category_2"] = pd.Categorical(df["Category_2"], ordered=True,categories=["NonPreg_ASX", "NonPreg_SEV", "Preg_CTRL", "Preg_ASX","Preg_SEV"])
        curr_colors = Category_2_colors

    total_per_donor = data_filtered.obs["sample"].value_counts().reset_index().rename(columns={"index": "sample", "sample": "total_cells"})
    df = df.merge(total_per_donor, on="sample", how="left")
    df["freq_group"] = df["Count"] * 100 / df["total_cells"]
    if (category == "Category_2"):
        df.to_csv(title + "_DA.csv", index=False)

    if (len(focus_clusters)!=0): # If we want to plot only some of the clusters we should remove them after computing their frequency out of all the cells
        df = df.loc[[x in focus_clusters for x in df[groupby]],:].copy()
        df[groupby] =pd.Categorical(df[groupby]).remove_unused_categories()

    plt.clf()
    sns.set_style("ticks")
    plt.figure(figsize=figsize)
    ax = sns.boxplot(data=df, x=groupby, y="freq_group", hue=category, palette=curr_colors, linewidth=0.25)
    sns.stripplot(data=df, x=groupby, y="freq_group", hue=category, size=2.5, dodge=True, palette=["black"] * 5)
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.ylabel("Frequency")
    plt.xlabel(title)
    ax.figure.suptitle(title, fontsize=12)

    groups = sorted(list(set(df[groupby])))
    if (category == "Category"):
        pairs = [((x, "Preg_COVID"), (x, "NonPreg_COVID")) for x in groups]
    else:
        pairs = [((x, "Preg_SEV"), (x, "Preg_ASX")) for x in groups]
        pairs = pairs + [((x, "Preg_SEV"), (x, "Preg_CTRL")) for x in groups]
        pairs = pairs + [((x, "NonPreg_SEV"), (x, "NonPreg_ASX")) for x in groups]
        pairs = pairs + [((x, "NonPreg_SEV"), (x, "Preg_SEV")) for x in groups]
    annotator = Annotator(data=df, x=groupby, y="freq_group", hue=category, pairs=pairs, ax=ax)
    annotator.configure(test='Mann-Whitney', comparisons_correction="BH", line_width=0.25, fontsize='x-small').apply_and_annotate()
    handles, labels = ax.get_legend_handles_labels()

    no_categories = len(set(df[category]))
    plt.legend(handles[0:no_categories], labels[0:no_categories], bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., fontsize=7)
    plt.tight_layout()

    pdf.savefig(ax.figure)


def clear_data_2(adata, obs_keep = []) :
    try :
        data = adata.to_anndata()
    except :
        data = adata
    data.obsm.clear()
    data.obsm.clear()
    data.varm.clear()
    data.uns.clear()
    data.var = data.var[[]]
    data.obs = data.obs[obs_keep]
    return data


# This functions plots TCR/BCR clone proportions on a given UMAP split by categories of conditions
# adata - the anndaat object that contains the UMAP coordinates and a "clone_prop" attribute in obs
# pdf_prefix - the prefix of the PDF file to which the plot will be saved
def clone_prop_umap_per_category(adata, pdf_prefix=None):
    fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(15, 3))
    ax = ax.ravel()
    clone_prop_max =  np.nanmax(adata.obs["clone_prop"],)
    categories = ['NonPreg_ASX', 'NonPreg_SEV', 'Preg_CTRL',  'Preg_ASX', 'Preg_SEV']
    sc.set_figure_params(dpi=300, dpi_save=300, vector_friendly=True)
    for i in range(len(categories)):
        cat = categories[i]
        sc.pl.umap(adata[(adata.obs["Category_2"]==cat),:], groups=cat, color="clone_prop", cmap="Reds",vmin=0, vmax=clone_prop_max, ax=ax[i], show=False, size=20)
        ax[i].set_title(cat, size=15)
        ax[i].set(xlabel=None)
        ax[i].set(ylabel=None)
    plt.tight_layout()
    if (pdf_prefix==None):
        plt.show()
    else:
        with PdfPages(pdf_prefix+".pdf") as pdf:
            pdf.savefig(fig, bbox_inches='tight')



# This functions plots TCR/BCR clone proportions on a given UMAP split by sample
# adata - the anndaat object that contains the UMAP coordinates and a "clone_prop" attribute in obs
# pdf_prefix - the prefix of the PDF file to which the plot will be saved
# vmax - a boolean vairable depicting whether to have a different color scale per sample based on its maximum value
def clone_prop_umap_per_sample(adata, pdf_prefix=None, vmax=True):
    clone_prop_max = np.nanmax(adata.obs["clone_prop"], )
    categories = ['NonPreg_ASX', 'NonPreg_SEV', 'Preg_CTRL', 'Preg_ASX', 'Preg_SEV']
    sc.set_figure_params(dpi=300, dpi_save=300, vector_friendly=True)
    with PdfPages(pdf_prefix + ".pdf") as pdf:
        for cat in set(categories):
            curr_data = adata[(adata.obs["Category_2"] == cat), :]
            samples = list(set(curr_data.obs["sample"]))
            if (len(samples)>6):
                rows=2
                cols=np.ceil(len(samples)/2).astype(int)
            else:
                rows = 1
                cols = len(samples)
            fig, ax = plt.subplots(nrows=rows, ncols=cols, figsize=(cols*2.5, rows*2.5))
            ax = ax.ravel()
            for i in range(len(samples)):
                samp=samples[i]
                samp_data = adata[(adata.obs["sample"] == samp), :]
                if vmax:
                    sc.pl.umap(samp_data, groups=samp, color="clone_prop", cmap="Reds", vmin=0, vmax=clone_prop_max, ax=ax[i], show=False, size=20)
                else:
                    sc.pl.umap(samp_data, groups=samp, color="clone_prop", cmap="Reds", vmin=0,ax=ax[i], show=False, size=20)
                ax[i].set_title(samp, size=10)
                ax[i].set(xlabel=None)
                ax[i].set(ylabel=None)
                cbar = ax[i].collections[0].colorbar
                cbar.ax.tick_params(labelsize=8)

            fig.suptitle(cat)
            #plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight', transparent=True)

# Computes the percentage of cells with a feature per cluster
# data - the anndata object
# cluster_col - the name of the clustering column in obs
# feature - the feature in obs to enumerate
def get_perc_out_of_cluster(data, cluster_col, feature):
    df_feature = data.obs[[cluster_col, feature]].value_counts().reset_index()
    df_feature = df_feature.loc[df_feature[feature]=="True",[cluster_col,0]].copy().rename(columns={0:"total_"+feature})
    df_cluster = data.obs[cluster_col].value_counts().reset_index().rename(columns={cluster_col:"total", "index":cluster_col})

    df = df_feature.merge(df_cluster, on = cluster_col)
    df["prop_"+feature] = (df["total_"+feature]*100)/df["total"]
    return(df)


# Plots a barplot of the percent of cells that are part of expanded clones out of the total number of cells that have TCR/BCR data, per cluster
# data - the anndata object
# pdf - the pdf object to plot into
# figsize - the figure size
# hue_order - category order to plot
# category - the category to split each cluster by
# pairs - what pairs to run statistical tests on
# colors - category colors to fill the barplots in
# cluster_col - the cluster column names in obs
# cluster_focus - a list of clusters to plot
def plot_expanded_clones_per_clust(data, pdf, figsize, hue_order, category, pairs, colors, cluster_col, cluster_focus=[]):
    # Computing number of expanded cells per cluster per sample out of total number of expanded clones per sample
    data = data[data.obs["expanded_t"], :].copy()
    summary_per_expanded = pd.DataFrame()
    clusters = list(set(data.obs[cluster_col]))
    for samp in list(set(data.obs["sample"])):
        curr = data[data.obs["sample"] == samp, :].copy()
        tot = curr.shape[0]

        curr_df = curr.obs[cluster_col].value_counts().reindex(clusters, fill_value=0).reset_index()
        curr_df = curr_df.rename(columns={"index": cluster_col, cluster_col:"no_expanded"})
        curr_df["perc_expanded"] = (curr_df["no_expanded"]*100) / tot
        curr_df["sample"] = samp
        curr_df["total_cells_per_sample"] = tot
        summary_per_expanded = summary_per_expanded.append(curr_df[["sample", "total_cells_per_sample", cluster_col, "no_expanded", "perc_expanded"]])

    df = summary_per_expanded.pivot_table(index=cluster_col, columns="sample", values="perc_expanded")
    df[df.isna()] = 0

    # sample order
    samps = data.obs[["sample", category]].drop_duplicates().reset_index(drop=True)
    samps = samps.sort_values(by=category)
    df = df[samps["sample"]]
    df = df.loc[data.obs[cluster_col].cat.categories.to_list(),:]

    if (len(cluster_focus)!=0):
        df = df.loc[cluster_focus,:].copy()

    plt.clf()
    plt.figure(figsize=figsize)
    sns.set_style('ticks')
    df_box3 = df.reset_index().melt(id_vars=[cluster_col])
    df_box3 = df_box3.merge(samps, on="sample")
    ax= sns.boxplot(data=df_box3, hue=category, y="value", x=cluster_col, palette=colors, hue_order=hue_order, linewidth=0.5)
    sns.stripplot(data=df_box3,y="value", x=cluster_col, hue=category, dodge=True, size=2.5, linewidth=0.25, ax=ax,hue_order=hue_order,palette=["black"]*len(hue_order))
    plt.ylabel("Percent of expanded clones per cluster out of total expanded clones")
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:len(hue_order)], labels[0:len(hue_order)], bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0.,fontsize=10)
    plt.show()

    pairs = [((x, pairs[0]),(x,pairs[1])) for x in list(set(df_box3[cluster_col]))]
    annotator = Annotator(data=df_box3, x=cluster_col, y="value", hue=category, ax=ax, pairs = pairs, hue_order=hue_order)
    annotator.configure(test='Mann-Whitney', comparisons_correction="BH", line_width=0.25).apply_and_annotate()

    pdf.savefig(ax.figure, bbox_inches='tight')



# Plot feature plots using hex-bin summary
# data - the anndata object
# genes - a list of genes or protein to plot the expression of
# ncols, nrows - number of columns and rows of figure panels
# figsize - the figure size
# grid size - the grid size for hex-bin summary
# pdf - the pdf object to plot to figures into
def plot_hex_feature_plots(data, genes, ncols, nrows, figsize, gridsize, pdf):
    x = data.obsm["X_umap"][:, 0]
    y = data.obsm["X_umap"][:, 1]

    plt.clf()
    count = 0
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=figsize, dpi=500)
    for i in range(nrows):
        for j in range(ncols):
            # After ploting all genes - remove the rest of the unused axes
            if (count > (len(genes) - 1)):
                fig.delaxes(axs[i][j])
            else:
                print(count)
                gene = genes[count]
                count += 1
                ax = axs[i][j]
                hb = ax.hexbin(x, y, C=data[:, gene].to_df()[gene], cmap='YlOrRd', gridsize=gridsize, linewidths=0.2)
                ax.set_title(gene, fontsize=8, fontweight='bold')
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                ax.set_yticks([])
                ax.set_xticks([])
                ax.set_rasterization_zorder(10)
                cbar = fig.colorbar(hb, ax=ax)
                cbar.ax.tick_params(labelsize=6)

    plt.tight_layout()
    pdf.savefig(fig)

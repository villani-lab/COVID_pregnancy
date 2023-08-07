import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc


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
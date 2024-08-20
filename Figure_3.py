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

# Figure 3A - CD8T UMAP
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
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


# Figure 3A - CD8T heatmap
genes2 = ["CD8A","CD8B","CD6",'CD4','CD40LG','ITGB1',"PRF1","GNLY","FGFBP2", #c4
          'GZMH','NKG7',"CST7", #c3
          'KLRD1','FCGR3A','GZMB', #c2
          "CXCR4","TGFB1","NR4A2", #c10
          "KLRF1","FCER1G","STMN1", # c15
          "CXCR3", "HMGB2", "PRDX3", "TALDO1", #c9
          "CMC1", "CD74", "GZMK", #c1
          "KLRB1", "IL7R", "CXCR6","RORC", # c12
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
g.figure.savefig("CD8_markers_heatmap_v3.pdf")


# Figure 3B - differential abundance
compute_DA_and_plot_per_category(data=cd8_data, category="Category", groupby="cluster_title", pdf=pdf, figsize=(5, 4), title="CD8 T cells", omit_clusters=["CD8T_7", "CD8T_8", "CD8T_16"], focus_clusters=["CD8T_2", "CD8T_3", "CD8T_4", "CD8T_12", "CD8T_13", "CD8T_14"])

# Figure 3C - Clonal expansion projected on a UMAP per category
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
clone_prop_umap_per_category(cd8_data, pdf_prefix="TCR_clone_prop_umap")


# Figure 3D - bar plot of clonal expansion per cluster
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
cd8_data.obs["cluster_title"] = [x.split(":")[0] for x in cd8_data.obs["cluster_title"]]
cd8_data.obs["cluster_title"]  = pd.Categorical(cd8_data.obs["cluster_title"] )
cd8_data_wo_doublets = cd8_data[[x not in ["CD8T_7","CD8T_8","CD8T_16"] for x in cd8_data.obs["cluster_title"]],:].copy()
samp_category2 = cd8_data_wo_doublets.obs[["sample","Category_2"]].drop_duplicates().reset_index(drop=True)
covid_samples =  list(samp_category2.loc[[x != "Preg_CTRL" for x in samp_category2["Category_2"]],"sample"])
cd8_data_wo_doublets = cd8_data_wo_doublets[[x in covid_samples for x in cd8_data_wo_doublets.obs["sample"]],:].copy()
cd8_data_wo_doublets.obs["Category"] = cd8_data_wo_doublets.obs["Category_2"].map({"Preg_SEV":"Preg_COVID","Preg_ASX":"Preg_COVID","NonPreg_SEV":"NonPreg_COVID","NonPreg_ASX":"NonPreg_COVID"})
with PdfPages("barplot_expanded_TCRs_2categories.pdf") as pdf:
    plot_expanded_clones_per_clust(cd8_data_wo_doublets, pdf=pdf, figsize=(4,4), hue_order=["NonPreg_COVID", "Preg_COVID"], category="Category", pairs=["Preg_COVID", "NonPreg_COVID"], colors ={"Preg_COVID":"#EE6677", "NonPreg_COVID":"#4477AA"}, cluster_col='cluster_title', cluster_focus = ["CD8T_2","CD8T_3","CD8T_4"])


# Figure 3E - Gene expression and protein expression characterizing CD8T_3
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
cd8_adt = an.read_h5ad("CD8T_ADT.h5ad") # From GSE239452

with PdfPages("GEX_features.pdf") as pdf:
    plot_hex_feature_plots(data=cd8_data, genes=["CD40LG", "GZMH", "ITGB1", "IL7R"], ncols=2, nrows=2, figsize=(5,4), gridsize=(200,200), pdf=pdf)
with PdfPages("ADT_features.pdf") as pdf:
    plot_hex_feature_plots(data=cd8_adt, genes=["prot_CD8", "prot_CD4","prot_CD45RO","prot_CD150","prot_CD57","prot_KLRG1"], ncols=3, nrows=3, figsize=(5, 4),gridsize=(200, 200), pdf=pdf)


# Figure 3F, G - COVID epitope hits projected on a UMAP
hits_A = pd.read_excel("VDJdb_hits_TRA_wAlleles.xlsx")
hits_B = pd.read_excel("VDJdb_hits_TRB_wAlleles.xlsx")

# We're missing J gene in the hits file - to get them back let's merge it with the "expanded_clones" file
hits_A = hits_A.rename(columns={"CDR3":"TRA_cdr3","J":'TRA_j_gene'})
hits_B = hits_B.rename(columns={"CDR3":"TRB_cdr3","J":'TRB_j_gene'})
hits = hits_A[['sample', 'TRB_cdr3', 'TRB_v_gene', 'TRA_j_gene', 'TRA_cdr3','TRA_v_gene', 'Category_2',"Epitope", "Epitope.gene", "Epitope.species"]].append(hits_B[['sample', 'TRB_cdr3', 'TRB_v_gene', 'TRB_j_gene', 'TRA_cdr3','TRA_v_gene', 'Category_2', "Epitope", "Epitope.gene", "Epitope.species"]])
expanded_clones = pd.read_csv("expanded_clones.csv")
hits_final = hits.merge(expanded_clones, on=['sample', 'TRB_cdr3', 'TRB_v_gene', 'TRA_cdr3','TRA_v_gene', 'Category_2'], how="left")
hits_final.pop("TRB_j_gene_x")
hits_final.pop("TRA_j_gene_x")
hits_final = hits_final.rename(columns={"TRA_j_gene_y":"TRA_j_gene", "TRB_j_gene_y":"TRB_j_gene"})
clone_comps = ["TRA_cdr3", "TRB_cdr3", "TRA_v_gene", "TRB_v_gene", "TRA_j_gene", "TRB_j_gene"]
hits_final["clone"] = hits_final[ clone_comps].apply(lambda row: "".join(row.values.astype(str)), axis=1)
hits_final = hits_final.drop_duplicates()
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
with PdfPages("COVID_hits.pdf") as pdf:
    covid_hits = hits_final.loc[hits_final["Epitope.species"]=="SARS-CoV-2",:].copy() # 33 lines
    covid_hits = covid_hits[["clone","Epitope.species"]].drop_duplicates() # 13 clones
    cd8_data.obs = cd8_data.obs.reset_index().merge(covid_hits, on=["clone"], how="left").set_index("index")
    g = sc.pl.umap(cd8_data, color="Epitope.species", size=7, palette=["red"], title="clones with hits to COVID", return_fig=True, show=False)
    pdf.savefig(g)

    # Quantify COVID hits per cluster
    cd8_data = cd8_data[[x not in ['Doublets (CD4/CD8)','High-MT', 'Doublets (Myeloid/CD8)'] for x in cd8_data.obs["cluster_title"]],:].copy()
    df_hits = cd8_data.obs[["Epitope.species","cluster_title"]].value_counts().reset_index().rename(columns={0:"COVID_hits"})
    df_has_tcr = cd8_data.obs.loc[cd8_data.obs["has_tcr"]=="True",["cluster_title"]].value_counts().reset_index().rename(columns={0:"has_tcr_counts"})
    df = df_hits.merge(df_has_tcr, on="cluster_title", how="left")
    df["per_COVID_hit"] = (df["COVID_hits"]*100) / df["has_tcr_counts"]
    df = df.sort_values(by="per_COVID_hit", ascending=False)
    df["CD8T_subsets"] = [x.split(":")[0] for x in df["cluster_title"]]
    df["CD8T_subsets"] = pd.Categorical(df["CD8T_subsets"], ordered=True, categories=df["CD8T_subsets"])
    df = df.reset_index(drop=True)

    plt.clf()
    sns.set(rc={'figure.figsize': (4,4)})
    sns.set_style("ticks")
    ax = sns.barplot(data = df, y="CD8T_subsets", x ="per_COVID_hit", palette=["grey"], orient='h')
    for i in range(len(ax.patches)):
        p = ax.patches[i]
        ax.annotate(str(df.loc[i, "COVID_hits"]), xy=(p.get_width(), p.get_y() + p.get_height() / 2), xytext=(5, 0), textcoords='offset points', ha="left", va="center", size=11)
    ax.figure.tight_layout()
    pdf.savefig(ax.figure)




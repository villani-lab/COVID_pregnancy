import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc
import pandas as pd
import anndata as an
import numpy as np
import statsmodels.stats.multitest as stats

exec(open("functions.py").read())


# Reading all anndata objects
cd8_data = an.read_h5ad("CD8T_GEX_TCR.h5ad") # From GSE239452
cd4_dat = an.read_h5ad("CD4T_GEX_TCR.h5ad") # From GSE239452
myeloid_data = an.read_h5ad("Myeloid.h5ad") # From GSE239452
B_data = an.read_h5ad("B_Plasma_GEX_BCR.h5ad") # From GSE239452

# UMAP per lineage - Figure 4A, C, E, G
sc.pl.umap(cd8_data, color = 'cluster_title')
sc.pl.umap(cd4_dat, color = 'cluster_title')
sc.pl.umap(myeloid_data, color = 'cluster_title')
sc.pl.umap(B_data, color = 'cluster_title')





# ISG figures
#=======================================================================================
# 1) Create BP for all objects without doublet clusters and withou control
#-------------------------------------------------------------------------
data_B_init = an.read_h5ad("B_Plasma_GEX_BCR.h5ad")
data_B = data_B_init[~data_B_init.obs['cluster_title'].isin(['high_ribo','Doublets (T/B)','Doublets (platelets/B)','high_mito']),:].copy()
data_B = data_B[data_B.obs["Category_2"].isin(['Preg_SEV', 'Preg_ASX', 'NonPreg_SEV', 'NonPreg_ASX']),:].copy()
calc_pseudobulk(data_B, save_name="B_woCont", cluster_label='cluster_title', sample_label="sample", condition_label="Category", thresh_cells=20, save=True)

data_CD8T_init = an.read_h5ad("CD8T_GEX_TCR.h5ad")
data_CD8T_init.obs["cluster_name"] = [x.split(":")[0] for x in data_CD8T_init.obs["cluster_title"]]
data_CD8T = data_CD8T_init[~data_CD8T_init.obs['cluster_name'].isin(['Doublets (Myeloid/CD8)','Doublets (CD4/CD8)','High-MT']),:].copy()
data_CD8 = data_CD8T[data_CD8T.obs["Category_2"].isin(['Preg_SEV', 'Preg_ASX', 'NonPreg_SEV', 'NonPreg_ASX']),:].copy()
calc_pseudobulk(data_CD8, save_name="CD8T_woCont", cluster_label='cluster_name', sample_label="sample", condition_label="Category", thresh_cells=20, save=True)

data_CD4T_init = an.read_h5ad("CD4T_GEX_TCR.h5ad")
data_CD4T = data_CD4T_init[~data_CD4T_init.obs['cluster_title'].isin(['high MT/ribo']),:].copy()
data_CD4T = data_CD4T[data_CD4T.obs["Category_2"].isin(['Preg_SEV', 'Preg_ASX', 'NonPreg_SEV', 'NonPreg_ASX']),:].copy()
calc_pseudobulk(data_CD4T, save_name="CD4T_woCont", cluster_label='cluster_title', sample_label="sample", condition_label="Category", thresh_cells=20, save=True)

data_myeloid_init = an.read_h5ad("Myeloid.h5ad")
data_myeloid = data_myeloid_init[~data_myeloid_init.obs['cluster_title'].isin(['Doublets (platelets/ Myeloid)', 'Doublets (B/myeloid)', 'Doublets (T/myeloid)']),:].copy()
data_myeloid = data_myeloid[data_myeloid.obs["Category_2"].isin(['Preg_SEV', 'Preg_ASX', 'NonPreg_SEV', 'NonPreg_ASX']),:].copy()
calc_pseudobulk(data_myeloid, save_name="myeloid_woCont", cluster_label='cluster_title', sample_label="sample", condition_label="Category", thresh_cells=20, save=True)

# Read PB files and concatenate
myeloid_PB = an.read_h5ad("myeloid_woCont_PB_norm_counts.h5ad")
CD4T_PB = an.read_h5ad("CD4T_woCont_PB_norm_counts.h5ad")
CD8T_PB = an.read_h5ad("CD8T_woCont_PB_norm_counts.h5ad")
B_PB = an.read_h5ad("B_woCont_PB_norm_counts.h5ad")
myeloid_PB.obs["Category"] = myeloid_PB.obs["Category"].replace({'SEV':'Preg', 'ASX':'Preg'})
CD4T_PB.obs["Category"] = CD4T_PB.obs["Category"].replace({'SEV':'Preg', 'ASX':'Preg'})
CD8T_PB.obs["Category"] = CD8T_PB.obs["Category"].replace({'SEV':'Preg', 'ASX':'Preg'})
B_PB.obs["Category"] = B_PB.obs["Category"].replace({'SEV':'Preg', 'ASX':'Preg'})

all_PB = an.concat([myeloid_PB, CD4T_PB, CD8T_PB, B_PB])
all_PB.write_h5ad("all_PB.h5ad")


# 2) Reading ISGs from MsigDB
#-------------------------------------------------------------------------
g_list = pd.read_csv("gamma_msigDB.csv", index_col=0)["x"].to_list() # from R code Figure4.R
a_list = pd.read_csv("alpha_msigDB.csv", index_col=0)["x"].to_list() # from R code Figure4.R

# how many genes do these two lists share?
len(set(g_list).intersection(set(a_list))) # no genes are overlapping
msigDB = g_list + a_list
ISGs = pd.read_csv("/projects/COVID_pregnancy/Results/exhaustion/ISGs.txt", header=None)[0].to_list()
# how many genes for MSigDB are shared with our 99 genes?
len(set(ISGs).intersection(set(msigDB))) # 51 genes are overlapping
pd.DataFrame(list(set(ISGs).union(set(msigDB)))).to_csv("ISG_long.csv")


# 3) Out of the extended list of ISGs - which are different between preg and non-preg
#---------------------------------
all_PB=an.read_h5ad("all_PB.h5ad")
ISGs=pd.read_csv("ISG_long.csv",index_col=0)["0"].to_list()
ISGs = list(set(ISGs).intersection(set(all_PB.var_names)))
all_PB = all_PB[:, ISGs].copy()
all_PB.obs["Lineage"] = [x.split("_")[0] for x in all_PB.obs["CellPopulation"]]
all_PB.obs["Lineage"] = all_PB.obs["Lineage"].replace({"pDC": "MNP", "MNP3": "MNP", "MNP2": "MNP", "MNP1": "MNP", "MNP4": "MNP", "MNP5": "MNP", "NK": "CD8T","gdT": "CD4T", "cDC": "MNP", "cDC": "MNP", "MAIT": "CD8T", "NK/NKT": "CD8T", "MNP6": "MNP",'Plasmablasts': "B"})

# test difference between pregnant and non-preg expression
# Making sure that there are at least 3 samples in each
df_pvals = pd.DataFrame(columns=["CellSubset","Pval"])
df_FC = pd.DataFrame(columns=["CellSubset","FC"])
for isg in ISGs:
    for subset in list(set(all_PB.obs["CellPopulation"])):
        #print(isg+" "+subset)
        curr = all_PB[all_PB.obs["CellPopulation"].isin([subset]), isg].copy()
        if ((sum(curr.obs["Category"]=="Preg")>2) & (sum(curr.obs["Category"]=="NonPreg")>2)):
            preg_values = curr[curr.obs["Category"] == "Preg", :].to_df()
            nonpreg_values = curr[curr.obs["Category"] == "NonPreg", :].to_df()
            ststistic, pvalue = sp.stats.mannwhitneyu(preg_values, nonpreg_values)
            fc = preg_values.mean()/ (nonpreg_values.mean()+0.001)
            df_pvals = df_pvals.append(pd.DataFrame([[subset+":"+isg, pvalue[0]]],columns=["CellSubset","Pval"]))
            df_FC = df_FC.append(pd.DataFrame([[subset+":"+isg, fc.loc[isg,]]] ,columns=["CellSubset","FC"]))
        else:
            df_pvals = df_pvals.append(pd.DataFrame([[subset + ":" + isg, np.nan]], columns=["CellSubset", "Pval"]))
            df_FC = df_FC.append(pd.DataFrame([[subset + ":" + isg, np.nan]], columns=["CellSubset", "FC"]))
df_pvals.to_csv("ISG_per_clust_pvals.csv", index=False)
df_FC.to_csv("ISG_per_clust_FC.csv", index=False)
# These 2 files above will comprise supp table 9

df_pvals = pd.read_csv("ISG_per_clust_pvals.csv")
df_FC = pd.read_csv("ISG_per_clust_FC.csv")

df_pvals["cluster"] = [x.split(":")[0] for x in df_pvals["CellSubset"]]
df_pvals["gene"] = [x.split(":")[1] for x in df_pvals["CellSubset"]]
df_FC["cluster"] = [x.split(":")[0] for x in df_FC["CellSubset"]]
df_FC["gene"] = [x.split(":")[1] for x in df_FC["CellSubset"]]
# What clusters don't have enough samples for a robust comparison?
df_na = df_pvals.loc[df_pvals["Pval"].isna(),:].copy()
df_na["cluster"].value_counts() # these are the clusters
na_clusters = df_na["cluster"].value_counts().index.to_list()

df_pvals_final = df_pvals.loc[~df_pvals["Pval"].isna(),:].copy()
df_FC_final = df_FC.loc[~df_FC["FC"].isna(),:].copy()
#df_FC_final["FC"].plot.hist(bins=2000, xlim=(0,4), ylim=(0,1000))

pvals_mat = pd.pivot_table(df_pvals_final, index="cluster", columns="gene", values="Pval")
FC_mat = pd.pivot_table(df_FC_final, index="cluster", columns="gene", values="FC")

pvals_adj_mat = pd.DataFrame(columns = pvals_mat.columns)
for clust in pvals_mat.index:
    x = pd.DataFrame(data=stats.multipletests(pvals_mat.loc[clust,:], method="fdr_bh")[1], columns=[clust], index=pvals_mat.columns).T
    pvals_adj_mat = pvals_adj_mat.append(x)

# 4) Filtering genes - adj_pval<0.05 AND FC>1.3
#-----------------------------------------------------------------------
FC_mat = np.log10(FC_mat+0.0001)
FC_mat.melt()["value"].plot.hist(bins=200)
FC_mat.to_csv("FC_mat.csv")
pvals_adj_mat.to_csv("Pval_mat.csv")
# These files will be used to plot the heatmap in SuppFigure 10I


bool_FC = abs(FC_mat)>0.2 # FC>1.58
bool_pval = pvals_adj_mat<0.05
relevant_comps = bool_FC & bool_pval

gene_sums = relevant_comps.sum()
sum(gene_sums==0) # 172 genes are not significant anywhere
sig_indx = ~(gene_sums==0)
df_diff = pvals_adj_mat.loc[:,sig_indx.to_list()] # We're left with 90 genes that have at least one significant comparison between preg and non-preg



#=======================================================================================
# ISG signature per lineage - Figure 4B, D, F, Supp Figure 9C
# For every lineage we first compute Z-scores across all cells in the anndata object and then compute mean z-score for ISG gene set.

# Lets now focus on the ISGs that show some kind of difference between preg and non-preg and re-plot the main figures
FC_mat = pd.read_csv("../Expand_ISG_list/FC_mat.csv", index_col=0)
pvals_adj_mat = pd.read_csv("../Expand_ISG_list/Pval_mat.csv", index_col=0)
bool_FC = abs(FC_mat)>0.2 # FC>1.58
bool_pval = pvals_adj_mat<0.05
relevant_comps = bool_FC & bool_pval
gene_sums = relevant_comps.sum()
sig_indx = ~(gene_sums==0)
sig_ISGs = sig_indx[sig_indx].index.to_list()


ISGs = pd.read_excel("Supp Table 9 - ISGs.xlsx", header=None)[0].to_list() # Supplementary table 9

# CD8 T cells
signature_score_per_cell_2(data=cd8_data, gene_set=sig_ISGs, score_title="ISG_score")
plt.clf()
with PdfPages("ISG_score_CD8T.pdf") as pdf:
    plot_ISG_score_boxplot(data=cd8_data, cluster_col="cluster_title", title="CD8T", pdf=pdf, doublet_clusters=['Doublets (CD4/CD8)','High-MT','Doublets (Myeloid/CD8)'], shorten=True)

# Myeloid cells
signature_score_per_cell_2(data=myeloid_data, gene_set=sig_ISGs, score_title="ISG_score")
plt.clf()
with PdfPages("ISG_score_Myeloid.pdf") as pdf:
    plot_ISG_score_boxplot(data=myeloid_data, cluster_col="cluster_title", title="Myeloid", pdf=pdf, doublet_clusters=['Doublets (B/myeloid)','Doublets (T/myeloid)',"Doublets (platelets/ Myeloid)"])

# CD4 T cells
signature_score_per_cell_2(data=cd4_dat, gene_set=sig_ISGs, score_title="ISG_score")
plt.clf()
with PdfPages("ISG_score_CD4T.pdf") as pdf:
    plot_ISG_score_boxplot(data=cd4_dat, cluster_col="cluster_title", title="CD4T", pdf=pdf, doublet_clusters=["high MT/ribo"])

# B cells
signature_score_per_cell_2(data=B_data, gene_set=sig_ISGs, score_title="ISG_score")
plt.clf()
with PdfPages("ISG_score_B.pdf") as pdf:
    plot_ISG_score_boxplot(data=B_data, cluster_col="cluster_title", title="B_Plasma", pdf=pdf, doublet_clusters = ["Doublets (platelets/B)","high_ribo","high_mito","Doublets (T/B)"])



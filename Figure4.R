library(fgsea)

# Getting full list of ISGs
gene_sets <- gmtPathways("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/msigdb_symbols.gmt")

gamma_list = gene_sets$HALLMARK_INTERFERON_GAMMA_RESPONSE
alpha_list = gene_sets$HALLMARK_INTERFERON_ALPHA_RESPONSE

write.csv(gamma_list, "/projects/COVID_pregnancy/Results/16_Paper/Revision/ISGs/Expand_ISG_list/gamma_msigDB.csv")
write.csv(alpha_list, "/projects/COVID_pregnancy/Results/16_Paper/Revision/ISGs/Expand_ISG_list/alpha_msigDB.csv")

#**************************************
library(ggplot2)
library(ComplexHeatmap)
require(circlize)

# Files below are based on code in Figure4.py
FC_mat = read.csv("FC_mat.csv", row.names = 1)
pval_mat = read.csv("Pval_mat.csv", row.names = 1)

bool_FC = abs(FC_mat)>0.2 # FC>1.58
bool_pval = pval_mat<0.05
relevant_comps = bool_FC & bool_pval


plot_heatmap = function(FC_mat_final, pval_mat_final, tit, width=40, height=10)
{
  FC_mat_final[FC_mat_final>1]=1
  FC_mat_final[FC_mat_final<(-1)]=-1

  col_FC <- colorRamp2(c(min(FC_mat_final),0, max(FC_mat_final)), c("blue","white", "red"))
  cell_fun = function(j, i, x, y, width, height, fill)
                    {
                         grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = col_FC(as.numeric(FC_mat_final[i,j])), col='black'))
                         if ((pval_mat_final[i, j] < 0.05) & (abs(FC_mat_final[i,j])>0.2))
                            {grid.circle(x, y, r = unit(1, "mm"), gp = gpar(fill = "light grey"))}
                    }

  ht <- Heatmap(FC_mat_final,
          cell_fun = cell_fun,
          rect_gp = gpar(type = "none"),
          border_gp = gpar(col = "black"),
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          show_row_names = TRUE,
          row_title = NULL,
          column_names_side = "top",
          row_names_side = "left",
          show_column_names = TRUE,
          width = ncol(FC_mat_final)*unit(6, "mm"),
          height = nrow(FC_mat_final)*unit(6, "mm"))

  pdf(paste0(tit,".pdf"), width = width, height=height)
  draw(ht, heatmap_legend_side = "left")
  dev.off()
}


# Sorting for genes with at least 1 comparison (Figure S10I)
#---------------------------------------------------------------------------
sig_genes = !(apply(relevant_comps,2, sum)==0) # only 91 genes have at least one significant comparison
FC_mat_final = FC_mat[,sig_genes]
pval_mat_final = pval_mat[,sig_genes]
plot_heatmap(FC_mat_final, pval_mat_final, tit="Filtered_1comp_heatmap", width=25, height=10)


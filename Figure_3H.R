library(ggplot2)
library("alakazam")
setwd("/projects/COVID_pregnancy/Results/16_Paper/CD8T_clonal_diversity/")

df <- read.csv("clones_per_patient.csv")
df <- na.omit(df)
df <- df[df["clone"]!="",]
res <- data.frame()
samples_category = unique(df[c("sample","Category_2")])
rownames(samples_category) <- 1:nrow(samples_category)

data <- estimateAbundance(df, group="sample", clone="clone", nboot=100, min_n=200)
div  <- alphaDiversity(data)
saveRDS(div, "diversity_object_min200.rds")

ggplot(div@diversity, aes(x=q, y=!!rlang::sym("d"), group=!!rlang::sym(data@group_by)))+
  geom_line(aes(color=!!rlang::sym(data@group_by)), size=1)+
  scale_color_manual(name="sampels",  values=sample_color) +
  baseTheme()

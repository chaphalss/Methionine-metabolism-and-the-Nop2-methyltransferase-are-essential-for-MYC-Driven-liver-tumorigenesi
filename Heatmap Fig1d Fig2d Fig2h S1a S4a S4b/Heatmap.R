setwd("D:\\Sensen\\Methionine\\Goga\\EC4 tumor microarray\\methionine related genes in EC4 cells")
res <- read.table("methionine related genes MYC high vs low.txt", header=TRUE, sep = "\t")
View(res)
row.names(res) <- res$GeneName
res <- res[,2:7]
res_matrix <- data.matrix(res)
library(stats)
library(gplots)
library(pheatmap)
pheatmap(res_matrix, cluster_rows=TRUE, cluster_cols = TRUE, show_rownames=TRUE, scale = "row", fontsize_row=6, cellwidth=20, color=colorRampPalette(c("navy", "white", "red"))(100))
library(made4)
heatplot(res,margins=c(8,20),)
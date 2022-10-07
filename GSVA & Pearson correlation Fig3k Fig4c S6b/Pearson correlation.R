# Pearson correlation
library(ggpubr)
my_data <- read.table('NOP2 VS cell cycle.txt', header=TRUE, row.names=1, sep='\t')
my_data$name <- rownames(my_data)
ggscatter(my_data, x = "NOP2", y = "CELL_CYCLE", title=" NOP2 vs cell cycle", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", label = "name", repel = TRUE, xlab = "Normalized NOP2 expression", ylab = "Pathway score")

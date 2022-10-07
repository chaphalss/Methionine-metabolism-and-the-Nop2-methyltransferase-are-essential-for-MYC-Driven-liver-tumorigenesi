rData<-read.table('I:\\RNA seq datasets\\star\\EC4\\self analysis data\\htseq_counts.txt', header=TRUE, row.names=1, sep='\t')
dim(rData)
condition <- c('Myc high Methionine high', 'Myc high Methionine high', 'Myc low Methionine high', 'Myc low Methionine high', 'Myc low Methionine high', 'Myc high Methionine low', 'Myc high Methionine low', 'Myc high Methionine low', 'Myc low Methionine low', 'Myc low Methionine low', 'Myc low Methionine low')
colData <- data.frame (row.names=colnames (rData), condition=factor(condition, levels=c('Myc high Methionine high', 'Myc low Methionine high', 'Myc high Methionine low', 'Myc low Methionine low')))
library(DESeq2)
dataset<-DESeqDataSetFromMatrix(countData=rData,colData=colData,design=~condition)
dds<-DESeq(dataset)
result<-results(dds, contrast=c('condition', 'Myc high Methionine high', 'Myc low Methionine high'))
result <- result[complete.cases(result),] 
head(result)
summary(result)
plotMA(result, main='Myc high vs. low', ylim=c(-5,5))
result_sig<-subset(result, padj<.05&abs(log2FoldChange)>0.58&(!is.na(result$padj)))
View(result_sig)
write.csv(as.data.frame(result_sig), file="MYC high vs low DEGs.csv")
library(clusterProfiler)
library(org.Mm.eg.db)
df2 = read.csv("MYC high vs low DEGs.csv", header=TRUE)
original_gene_list <- df2$log2FoldChange
names(original_gene_list) <- df2$Gene
gene_list<-na.omit(original_gene_list) 
gene_list = sort(gene_list, decreasing = TRUE)
gse1 <- gseGO(geneList=gene_list,
              ont ="BP",                         
              keyType = "ENSEMBL",     
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              OrgDb ="org.Mm.eg.db",     
              pAdjustMethod = "BH")
require(DOSE)
dotplot(gse1, showCategory=10, font.size=10, split=".sign") + facet_grid(.~.sign)
library(enrichplot)
gseaplot2(gse1, geneSetID = 1, title = gse1$Description[1])


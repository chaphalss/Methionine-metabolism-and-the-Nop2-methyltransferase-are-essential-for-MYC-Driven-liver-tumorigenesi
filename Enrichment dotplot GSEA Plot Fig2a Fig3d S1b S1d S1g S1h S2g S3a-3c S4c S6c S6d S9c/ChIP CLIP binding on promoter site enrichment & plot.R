SRR16936850Peaks <- "LT2-MYC MYC on_peaks.csv"
SRR16936850Peaks_DF <- read.delim(SRR16936850Peaks, header=TRUE, sep=",")
SRR16936850Peaks_DF[1:2, ]
library(GenomicRanges)
SRR16936850Peaks_GR <- GRanges(seqnames = SRR16936850Peaks_DF[, "chr"], IRanges(SRR16936850Peaks_DF[, "start"], SRR16936850Peaks_DF[, "end"]))      
mcols(SRR16936850Peaks_GR) <- SRR16936850Peaks_DF[, c("abs_summit", "fold_enrichment")] 
SRR16936850Peaks_GR
# download gtf file Mus_musculus.GRCm39.107.gtf from http://ftp.ensembl.org/pub/release-107/gtf/mus_musculus/Mus_musculus.GRCm39.107.gtf.gz
library(GenomicFeatures)
library(AnnotationDbi)
txdb = makeTxDbFromGFF('Mus_musculus.GRCm39.107.gtf')
library(org.Mm.eg.db)
library(GenomeInfoDb)
library(ChIPseeker)
peakAnno <- annotatePeak(SRR16936850Peaks_GR, tssRegion = c(-1000, 1000), TxDb=txdb, annoDb = "org.Mm.eg.db")
class(peakAnno)
peakAnno
peakAnno_GR <- as.GRanges(peakAnno)
peakAnno_DF <- as.data.frame(peakAnno)        
peakAnno_GR[1:2, ]
plotAnnoBar(peakAnno)
plotAnnoPie(peakAnno)
covplot(SRR16936850Peaks_GR, weightCol="fold_enrichment")
plotDistToTSS(peakAnno)
peakHeatmap(peakAnno_GR, weightCol="fold_enrichment", TxDb=txdb, upstream=10000, downstream=10000, xlab="genomic region (5' to 3')", color="blue")
SRR16936850_annot <- data.frame(peakAnno@anno)
SRR16936850_annot
peakAnno_GR[1, ]
annotatedPeaksGR_TSS <- peakAnno_GR[peakAnno_GR$annotation == "Promoter", ]
genesWithPeakInTSS <- unique(annotatedPeaksGR_TSS$geneId)
genesWithPeakInTSS[1:2]
allGeneGR <- genes(txdb)
allGeneGR[1:2, ]
allGeneIDs <- allGeneGR$gene_id
library(clusterProfiler)
GO_result <- enrichGO(gene = genesWithPeakInTSS, universe = allGeneIDs, keyType = "ENSEMBL", OrgDb = org.Mm.eg.db, ont = "BP")
GO_result_df <- data.frame(GO_result)
GO_result_df[1:5, ]
library(enrichplot)
GO_result_plot <- pairwise_termsim(GO_result)
dotplot(GO_result_plot, showCategory=20, font.size=10)



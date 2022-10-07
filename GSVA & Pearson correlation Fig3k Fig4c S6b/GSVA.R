library(DESeq2)
library(qusage)
library(GSVA)
library(org.Hs.eg.db)
library(magrittr)
library(msigdbr)
library(BaseSet)
library(AnnotationDbi)
hallmark_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") 
head(hallmark_gene_sets)
hallmarks_list <- split(hallmark_gene_sets$entrez_gene, hallmark_gene_sets$gs_name)
head(hallmarks_list, n = 2)
#download TCGA data
library(TCGAbiolinks)
query_seq <- GDCquery(project = "TCGA-LIHC", data.category = "Transcriptome Profiling", data.type="Gene Expression Quantification",  experimental.strategy = "RNA-Seq", sample.type= c("Primary Tumor"), access = "open")
GDCdownload(query_seq)
RnaseqSE <- GDCprepare(query_seq)
rnaseq <- assay(RnaseqSE)
rnaseq <- rnaseq[which(!duplicated(rownames(rnaseq))),]
rnaseqPatients <- cbind(colnames(rnaseq), gsub('\\.','-',substring(colnames(rnaseq),1,12)))
colnames(rnaseqPatients) <- c( "barcode","patient")
data_rnaseqPatients <- as.data.frame(rnaseqPatients)
data_rnaseqPatients$patient <- gsub('-','\\.', data_rnaseqPatients$patient)
colnames(rnaseq) <- gsub('-','\\.',colnames(rnaseq))
data_rnaseq <- as.data.frame(rnaseq)
data_rnaseq <- data_rnaseq[,colnames(data_rnaseq)[order(match(colnames(data_rnaseq), data_rnaseqPatients$patient))]]
write.table(data_rnaseq,"TCGA-LIHC_RNASeq_rawcounts.txt",col.name=TRUE,sep="\t",row.names=TRUE,quote=FALSE)
write.table(data_rnaseqPatients,"TCGA-LIHC_RNASeq_classdefinitions.txt",col.name=TRUE,sep="\t",row.names=TRUE,quote=FALSE)
data_rnaseq <- read.table('TCGA-LIHC_RNASeq_rawcounts.txt', row.names=1, header=TRUE, sep='\t')
data_rnaseqPatients<- read.table('TCGA-LIHC_RNASeq_classdefinitions.txt', header=TRUE, sep='\t')
data_rnaseqPatients$barcode <- gsub('-','\\.', data_rnaseqPatients$barcode)
library(textshape)
data_rnaseqPatients <- column_to_rownames(data_rnaseqPatients, "barcode")
data_rnaseqPatients$barcode <- gsub('-','\\.', data_rnaseqPatients$barcode)
library(magrittr)
data_rnaseqPatients <- data_rnaseqPatients %>%
  dplyr::select(-patient)
data_rnaseq <- data_rnaseq %>%
  dplyr::filter(rowMeans(.) >1) %>%
  round()
dds <- DESeqDataSetFromMatrix(
  countData = data_rnaseq, 
  colData = data_rnaseqPatients, 
  design = ~1)
dds_norm <- vst(dds)
vst_df <- assay(dds_norm) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_id") 
vst_df$ensembl_id <- gsub("\\..*","",vst_df$ensembl_id)
mapped_df <- data.frame(
  "entrez_id" = mapIds(org.Hs.eg.db, keys = vst_df$ensembl_id, keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first")) %>% 
  dplyr::filter(!is.na(entrez_id)) %>%
  tibble::rownames_to_column("Ensembl") %>%
  dplyr::inner_join(vst_df, by = c("Ensembl" = "ensembl_id"))
head(mapped_df)
sum(duplicated(mapped_df$entrez_id))
gene_means <- rowMeans(mapped_df %>% dplyr::select(-Ensembl, -entrez_id))
mapped_df <- mapped_df %>%
  dplyr::mutate(gene_means) %>%     
  dplyr::select(Ensembl, entrez_id, gene_means, dplyr::everything())
filtered_mapped_df <- mapped_df %>%
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  dplyr::distinct(entrez_id, .keep_all = TRUE)
sum(duplicated(filtered_mapped_df$entrez_id))
filtered_mapped_matrix <- filtered_mapped_df %>%
  dplyr::select(-Ensembl, -gene_means) %>%
  tibble::column_to_rownames("entrez_id") %>%  
  as.matrix()
gsva_results <- gsva(
  filtered_mapped_matrix,
  hallmarks_list,
  method = "gsva",
  kcdf = "Gaussian",
  min.sz = 10,
  max.sz = 9999,
  mx.diff = TRUE,
  verbose = FALSE)
head(gsva_results[, 1:10])
gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  readr::write_tsv("TCGA-LIHC_gsvaresults_C5_BP.tsv")
library(pheatmap)
db <- as.data.frame(filtered_mapped_matrix[,c("condition")])
rownames(db) <- colnames(gsva_results)
pheatmap::pheatmap(gsva_results, annotation_col = db, show_colnames = FALSE, fontsize_row = 2, height= 100,  filename="TCGA-LIHC_gsvaresults_C5_BP.pdf")

    
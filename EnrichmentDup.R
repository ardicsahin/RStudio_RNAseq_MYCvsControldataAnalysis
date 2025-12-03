library(readr)
library(dplyr)
library(DESeq2)
library(edgeR)
library(purrr)
library(readxl)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pheatmap)
library(ggplot2)
featureCountData <- read_csv("/Users/ardicsahin/Downloads/Galaxy42-[Column join on data 35, data 33, and others].csv")

featureCountData <- as.data.frame(featureCountData)
colnames(featureCountData) <- c("GeneID","MYC1","MYC2","CN1","CN2")
rownames(featureCountData) <- featureCountData$GeneID
featureCountData <- na.omit(featureCountData)

featureCountData <- featureCountData[ , -1]
entrez_rows <- featureCountData
entrez_rows <- entrez_rows[-28396,]

featureCountData <- data.frame(lapply(featureCountData, as.numeric))

featureCountData <- na.omit((featureCountData))
row.names(featureCountData) <- row.names(entrez_rows)
conditions <- c(rep("MYC",2),rep("CN",2))
sample.design <- data.frame(condition=conditions, row.names=colnames(featureCountData))

ddset <- DESeqDataSetFromMatrix(countData = featureCountData , colData = sample.design, design =~ condition)

dds <- estimateSizeFactors(ddset)
dds <- ddset[rowSums(counts(dds))>0 ,]

dds<- DESeq(dds)

MYCvsCN <- data.frame(results(dds, contrast= c("condition","MYC","CN")))
MYCvsCN_filter <- MYCvsCN %>% filter(MYCvsCN$pvalue< 0.05)
View(MYCvsCN_filter)

entrez_to_symbol <-bitr(row.names(MYCvsCN_filter), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
entrez_to_symbol <- entrez_to_symbol[!duplicated(entrez_to_symbol$ENTREZID),]

MYCvsCN_filter$ENTREZID <- row.names(MYCvsCN_filter)
MYCvsCN_filter <- merge(MYCvsCN_filter, entrez_to_symbol, by.x = "ENTREZID", by.y = "ENTREZID")

MYCvsCN_filter <- MYCvsCN_filter[!duplicated(MYCvsCN_filter$ENTREZID),]

rownames(MYCvsCN_filter) <- MYCvsCN_filter$SYMBOL

original_gene_list <- MYCvsCN_filter$log2FoldChange
names(original_gene_list) <- MYCvsCN_filter$SYMBOL

gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

res <- results(dds)
resultsNames(dds)

sig_genes <- rownames(res)[which(res$padj < 0.05)]
norm_counts <- counts(dds, normalized= T)

sig_gene_counts <- norm_counts[sig_genes,]


#Dot Plot
gse <- gseGO(geneList = gene_list, ont = "BP", keyType = "SYMBOL", nPerm = 1000, minGSSize = 1, maxGSSize = 1000,
             verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
p <- dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(.~.sign) #Sor

p + theme(axis.text.y = element_text(size=10))

#Enrichment
original_gene_list <- MYCvsCN$log2FoldChange
names(original_gene_list) <- row.names(MYCvsCN)
View(original_gene_list)

gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

sig_genes_df = subset(MYCvsCN, padj < 0.05)

genes <- sig_genes_df$log2FoldChange
names(genes) <- row.names(sig_genes_df)
genes <- na.omit(genes)

View(genes)

View(genes)
go_enrich <- enrichGO(gene = names(genes), universe = names(gene_list), OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
View(go_enrich)
upsetplot(go_enrich)

#Bar Plot
barplot(go_enrich, drop = TRUE, showCategory = 8, title = "GO Biological Pathways", font.size = 10)

#Dot Plot
dotplot(go_enrich)

ego_sim <- pairwise_termsim(go_enrich)
emapplot(ego_sim)

goplot(go_enrich, showCategory = 10)



cnetplot(go_enrich, categorySize = "geneNum", foldChange = gene_list, showCategory = 5)

cnetplot(go_enrich, categorySize = "pvalue", foldChange = gene_list, showCategory = 2)

enrich_df <- as.data.frame(go_enrich)

View(enrich_df)






library(readr)
library(DESeq2)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(pheatmap)
library(ggplot2)

featureCountData <- read_csv("/Users/ardicsahin/Downloads/Galaxy42-[Column join on data 35, data 33, and others].csv")

featureCountData <- as.data.frame(featureCountData)
colnames(featureCountData) <- c("GeneID","MYC1","MYC2","CN1","CN2")
rownames(featureCountData) <- featureCountData$GeneID
featureCountData <- featureCountData[ ,-1]
entrez_rows <- featureCountData
entrez_rows <- entrez_rows[-28396,]

featureCountData <- data.frame(lapply(featureCountData, as.numeric))

featureCountData <- na.omit(featureCountData)

conditions <- c(rep("MYC",2),rep("CN",2))
sample.design <- data.frame(condition=conditions, row.names=colnames(featureCountData))

ddset <- DESeqDataSetFromMatrix(countData = featureCountData , colData = sample.design, design =~ condition)

dds <- estimateSizeFactors(ddset)
dds <- ddset[rowSums(counts(dds))>0 ,]
dds<- DESeq(dds)

res <- results(dds)
resultsNames(dds)

sig_genes <- rownames(res)[which(res$padj < 0.05)]
norm_counts <- counts(dds, normalized= T)

sig_gene_counts <- norm_counts[sig_genes,]

pheatmap(sig_gene_counts, scale= "row", show_rownames= T, cluster_cols= T, main="Heatmap of Significant Genes")

res_t <- results(dds)
res_ordered <- res[order(res$padj),]
top_genes <- row.names(res_ordered)[1:20]

counts_t <- counts(dds, normalized= T)
counts_top <- counts_t[top_genes,]

log_counts_top <- log2(counts_top+1)

heatmap_20 <- pheatmap(log_counts_top, scale= "row", main="Top 20 Genes")




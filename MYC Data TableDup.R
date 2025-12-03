
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
featureCountData <- read_csv("/Users/ardicsahin/Downloads/Galaxy42-[Column join on data 35, data 33, and others].csv")

featureCountData <- as.data.frame(featureCountData)
colnames(featureCountData) <- c("GeneID","MYC1","MYC2","CN1","CN2")

View(featureCountData)
rownames(featureCountData) <- featureCountData$GeneID
featureCountData <- na.omit(featureCountData)
View(featureCountData)


featureCountData <- featureCountData[ , -1]
entrez_rows <- featureCountData
entrez_rows <- entrez_rows[-28396,]

featureCountData <- data.frame(lapply(featureCountData, as.numeric))

featureCountData <- na.omit(featureCountData)
row.names(featureCountData) <- row.names(entrez_rows)

conditions <- c(rep("MYC",2),rep("CN",2))
sample.design <- data.frame(condition=conditions, row.names=colnames(featureCountData))

ddset <- DESeqDataSetFromMatrix(countData = featureCountData , colData = sample.design, design =~ condition)

dds <- estimateSizeFactors(ddset)
dds <- ddset[rowSums(counts(dds))>0 ,]

dds<- DESeq(dds)

MYCvsCN <- data.frame(results(dds, contrast= c("condition","MYC","CN")))
MYCvsCN_filter <- MYCvsCN %>% filter(MYCvsCN$pvalue< 0.05)

entrez_to_symbol <-bitr(row.names(MYCvsCN_filter), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
entrez_to_symbol <- entrez_to_symbol[!duplicated(entrez_to_symbol$ENTREZID),]

View(entrez_to_symbol)
View(MYCvsCN_filter)

MYCvsCN_filter$ENTREZID <- row.names(MYCvsCN_filter)

View(MYCvsCN_filter)
MYCvsCN_filter <- merge(MYCvsCN_filter, entrez_to_symbol, by.x = "ENTREZID", by.y = "ENTREZID")

MYCvsCN_filter <- MYCvsCN_filter[!duplicated(MYCvsCN_filter$ENTREZID),]

rownames(MYCvsCN_filter) <- MYCvsCN_filter$SYMBOL


View(MYCvsCN)
View(MYCvsCN_filter)
MYCvsCN_filter <- MYCvsCN_filter[, -1]
MYCvsCN_filter <- MYCvsCN_filter[, -7]


# RStudio_RNAseq_MYCcsControldataAnalysis

This code performs a standard Differential Expression (DE) analysis pipeline, taking raw count data (likely from Galaxy/FeatureCounts), running it through DESeq2, and annotating the results.

![Flowchart of Data Pipeline](MYCFlowchart.jpg)

RNA-Seq Analysis & Heatmap Visualization PipelineOverviewThis R script acts as a downstream analysis and visualization tool for RNA-Seq data. It takes raw feature counts (typically from Galaxy/FeatureCounts), performs differential expression analysis using DESeq2, and generates high-quality heatmaps to visualize expression patterns between MYC (Treatment) and CN (Control) groups.Key FeaturesStatistical Analysis: Uses the Wald test (via DESeq2) to identify differentially expressed genes.Global Visualization: Generates a heatmap of all significantly differentially expressed genes ($P_{adj} < 0.05$) to show broad expression trends.Targeted Visualization: Generates a focused heatmap of the Top 20 most significant genes using Log2-transformed counts for precise comparison.System RequirementsR VersionR >= 4.0.0Required LibrariesEnsure the following packages are installed via CRAN or Bioconductor:R# Data Manipulation & plotting
install.packages(c("readr", "data.table", "ggplot2", "pheatmap"))

. Pre-processing: Reads the specific Galaxy output file.Cleaning:Sets Gene IDs as row names.Removes explicit outlier row 28396 (User-specific artifact).Forces all count data to numeric type and removes NA values.2. Differential Expression (DESeq2)Design: ~ condition (MYC vs CN).Normalization: Uses estimateSizeFactors to account for sequencing depth.

Filtering: Removes genes with zero counts across all samples.
Testing: Runs the negative binomial Wald test.3. 
Visualization 1: Global Significance 
HeatmapFilter: Selects genes with Adjusted P-value < 0.05.
Data: Uses normalized counts.
Scaling: Row-scaled (Z-score) to highlight relative up/down-regulation patterns regardless of absolute expression magnitude.4. Visualization 2: Top 20 Genes HeatmapSelection: Sorts all genes by adjusted p-value and selects the top 20.

Transformation: Applies $Log_2(Count + 1)$ transformation. This stabilizes variance, making high-expression and low-expression genes comparable on the color scale.UsageSet File Path:Modify the read_csv path in the script to point to your data:RfeatureCountData <- read_csv("/path/to/your/data.csv")
Verify Outlier Removal:Check line 17 (entrez_rows <- entrez_rows[-28396,]). If using a new dataset, ensure this row index corresponds to the row you intend to remove (e.g., a header or summary statistic), or remove this line if not needed.Run:Execute the script in RStudio. The heatmaps will appear in the "Plots" pane.OutputsDESeq2 

Object: Available in the R environment as dds.Plot 1: Heatmap of all significant genes (Row Z-score).Plot 2: Heatmap of top 20 significant genes (Log2 Expression).




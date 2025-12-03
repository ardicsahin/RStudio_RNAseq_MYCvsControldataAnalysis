# RStudio_RNAseq_MYCcsControldataAnalysis

 mermaid
graph TD
    %% Nodes
    Input[("Input CSV<br/>(Galaxy42.csv)")]
    
    subgraph Preprocessing [Data Cleaning]
        Import[read_csv]
        Clean1[Rename Columns]
        Clean2[Assign Rownames &<br/>Remove GeneID Col]
        Clean3[Convert to Numeric &<br/>Round to Integers]
    end

    subgraph Analysis [DESeq2 Pipeline]
        DDS_Obj[Create DESeqDataSet]
        Filter[Pre-filter low counts]
        Run[Run DESeq()]
    end

    subgraph Outputs [Results & Viz]
        Res(Results Extraction<br/>padj < 0.05)
        Norm(Normalization<br/>VST or rlog)
        Heatmap1[pheatmap:<br/>All Significant Genes]
        Heatmap2[pheatmap:<br/>Top 20 Genes]
    end

    %% Connections
    Input --> Import
    Import --> Clean1 --> Clean2 --> Clean3
    Clean3 --> DDS_Obj
    DDS_Obj --> Filter --> Run
    
    Run --> Res
    Run --> Norm
    
    Res -->|"Filter ID list"| Heatmap1
    Res -->|"Sort & Select Top 20"| Heatmap2
    Norm -->|"Data Matrix"| Heatmap1
    Norm -->|"Data Matrix"| Heatmap2
    
    %% Styling
    style Input fill:#f9f,stroke:#333,stroke-width:2px
    style Run fill:#bbf,stroke:#333,stroke-width:2px
    style Heatmap1 fill:#bfb,stroke:#333,stroke-width:2px
    style Heatmap2 fill:#bfb,stroke:#333,stroke-width:2px

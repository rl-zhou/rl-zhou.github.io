---
title:  "Differential Expression Analysis"
mathjax: true
layout: post
output: html_document
categories: project
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(apeglm)
library(gplots)
library(pheatmap)

#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install("apeglm")

#    if (!require("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
```


# Gene Quantification
Data loading
```{r}
# load data
rawCount <- read.delim("adar.rawcount.txt")
row.names(rawCount) <- rawCount[,1]
rawCount = rawCount[,-1]

# creating a meta table
name <- names(rawCount)
metaTable <- data.frame(row.names = name[1:length(name)], condition = as.factor(c(rep('adarsg1_ifnb', 3), rep('adarsg1_ns', 3), rep('ctrl_ifnb', 3), rep('ctrl_ns', 3))))
```

Count Normalization
```{r}
# DDS, Transforms count to log2FoldChange
dds <- DESeqDataSetFromMatrix(countData = rawCount, colData = metaTable, design = ~ condition)
keep <- rowSums(counts(dds)) >=10 ### keep genes that have at least 10 reads across all samples
dds <- dds[keep,] 
dds <- DESeq(dds) # the main DESeq function, including count  normalization

# extract the normalized counts
normCounts <- counts(dds, normalized=T)
write.csv(normCounts, 'normalizedCounts.csv')
```

Differential Gene Expression Analysis 
Creating Contrasts for comparisons
```{r}
# set the factor level
#dds$condition <- relevel(dds$condition, ref = 'ctrl') ###

adarsg1_ifnb_VS_adarsg1_ns <- c('condition', 'adarsg1_ifnb', 'adarsg1_ns')
adarsg1_ifnb_VS_ctrl_ifnb <- c('condition', 'adarsg1_ifnb', 'ctrl_ifnb')
adarsg1_ifnb_VS_ctrl_ns <- c('condition', 'adarsg1_ifnb', 'ctrl_ns')
adarsg1_ns_VS_ctrl_ifnb <- c('condition', 'adarsg1_ns', 'ctrl_ifnb')
adarsg1_ns_VS_ctrl_ns <- c('condition', 'adarsg1_ns', 'ctrl_ns')
ctrl_ifnb_VS_ctrl_ns <- c('condition', 'ctrl_ifnb', 'ctrl_ns')
```

The actual analysis data
```{r}
## unshrunken differential gene expression analysis results, alpha=0.05
### base mean, log2 fold change, p-value
res_adarsg1_ifnb_VS_adarsg1_ns <- results(dds, alpha=0.05, contrast=adarsg1_ifnb_VS_adarsg1_ns) 
res_adarsg1_ifnb_VS_ctrl_ifnb <- results(dds, alpha=0.05, contrast=adarsg1_ifnb_VS_ctrl_ifnb)
res_adarsg1_ifnb_VS_ctrl_ns <- results(dds, alpha=0.05, contrast=adarsg1_ifnb_VS_ctrl_ns) 
res_adarsg1_ns_VS_ctrl_ifnb <- results(dds, alpha=0.05, contrast=adarsg1_ns_VS_ctrl_ifnb) 
res_adarsg1_ns_VS_ctrl_ns <- results(dds, alpha=0.05, contrast=adarsg1_ns_VS_ctrl_ns) 
res_ctrl_ifnb_VS_ctrl_ns <- results(dds, alpha=0.05, contrast=ctrl_ifnb_VS_ctrl_ns) 
```

## PCA
```{r}
plotPCA(rlog(dds))
```





# Visualization 
## MA Plots
```{r}
#plotMA(res_adarsg1_ifnb_VS_adarsg1_ns)
plotMA(res_adarsg1_ifnb_VS_ctrl_ifnb)
#plotMA(res_adarsg1_ifnb_VS_ctrl_ns)
#plotMA(res_adarsg1_ns_VS_ctrl_ifnb)
#plotMA(res_adarsg1_ns_VS_ctrl_ns)
#plotMA(res_ctrl_ifnb_VS_ctrl_ns)
```
blue dots: significant

## Heatmap
y axis: log2fold change, avg expression, gene id
x axis: sample
- steps
  1. log2fold change and padj for res_adarsg1_ifnb_VS_ctrl_ifnb
  2. get the most differently expressed genes (negative log2 fold change?)
```{r}
# step1: extract the log2fold change and padj columns
heatmap_values <- res_adarsg1_ifnb_VS_ctrl_ifnb[, c(2,6)]

# step2: shrink the heatmap_values to contain the first 200 genes, alpha=0.01
topGenes <- subset(heatmap_values, padj<0.01 & (log2FoldChange>1 | log2FoldChange<(-1)))
head(topGenes[order(topGenes$log2FoldChange),], 100)
head(topGenes[order(topGenes$log2FoldChange, decreasing = T),])

topGenes.df <- as.data.frame(topGenes$log2FoldChange)
topGenes.df$genes <- rownames(normCounts)


```



# References
1. Youtube:
2. Evans, Ciaran et al. “Selecting between-sample RNA-Seq normalization methods from the perspective of their assumptions.” Briefings in bioinformatics vol. 19,5 (2018): 776-792. doi:10.1093/bib/bbx008


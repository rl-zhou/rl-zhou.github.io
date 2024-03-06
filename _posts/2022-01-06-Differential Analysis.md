---
title:  "Differential Expression Analysis"
mathjax: true
layout: post
output: html_document
categories: project
---


{% highlight r %}
knitr::opts_chunk$set(echo = TRUE)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(apeglm)
library(gplots)
library(pheatmap)
{% endhighlight %}


# Gene Quantification
Data loading
{% highlight r %}
# load data
rawCount <- read.delim("GSE110708_adar.rawcount.txt")
row.names(rawCount) <- rawCount[,1]
rawCount = rawCount[,-1]

# creating a meta table
name <- names(rawCount)
metaTable <- data.frame(row.names = name[1:length(name)], condition = as.factor(c(rep('adarsg1_ifnb', 3), rep('adarsg1_ns', 3), rep('ctrl_ifnb', 3), rep('ctrl_ns', 3))))
{% endhighlight %}

Count Normalization
{% highlight r %}
# DDS, Transforms count to log2FoldChange
dds <- DESeqDataSetFromMatrix(countData = rawCount, colData = metaTable, design = ~ condition)
keep <- rowSums(counts(dds)) >=10 ### keep genes that have at least 10 reads across all samples
dds <- dds[keep,] 
dds <- DESeq(dds) # the main DESeq function, including count  normalization

# extract the normalized counts
normCounts <- counts(dds, normalized=T)
write.csv(normCounts, 'normalizedCounts.csv')
{% endhighlight %}

Differential Gene Expression Analysis 
Creating Contrasts for comparisons
{% highlight r %}
# set the factor level
#dds$condition <- relevel(dds$condition, ref = 'ctrl') ###

adarsg1_ifnb_VS_adarsg1_ns <- c('condition', 'adarsg1_ifnb', 'adarsg1_ns')
adarsg1_ifnb_VS_ctrl_ifnb <- c('condition', 'adarsg1_ifnb', 'ctrl_ifnb')
adarsg1_ifnb_VS_ctrl_ns <- c('condition', 'adarsg1_ifnb', 'ctrl_ns')
adarsg1_ns_VS_ctrl_ifnb <- c('condition', 'adarsg1_ns', 'ctrl_ifnb')
adarsg1_ns_VS_ctrl_ns <- c('condition', 'adarsg1_ns', 'ctrl_ns')
ctrl_ifnb_VS_ctrl_ns <- c('condition', 'ctrl_ifnb', 'ctrl_ns')
{% endhighlight %}

The actual analysis data
{% highlight r %}
## unshrunken differential gene expression analysis results, alpha=0.05
### base mean, log2 fold change, p-value
res_adarsg1_ifnb_VS_adarsg1_ns <- results(dds, alpha=0.05, contrast=adarsg1_ifnb_VS_adarsg1_ns) 
res_adarsg1_ifnb_VS_ctrl_ifnb <- results(dds, alpha=0.05, contrast=adarsg1_ifnb_VS_ctrl_ifnb)
res_adarsg1_ifnb_VS_ctrl_ns <- results(dds, alpha=0.05, contrast=adarsg1_ifnb_VS_ctrl_ns) 
res_adarsg1_ns_VS_ctrl_ifnb <- results(dds, alpha=0.05, contrast=adarsg1_ns_VS_ctrl_ifnb) 
res_adarsg1_ns_VS_ctrl_ns <- results(dds, alpha=0.05, contrast=adarsg1_ns_VS_ctrl_ns) 
res_ctrl_ifnb_VS_ctrl_ns <- results(dds, alpha=0.05, contrast=ctrl_ifnb_VS_ctrl_ns) 
{% endhighlight %}

## PCA
{% highlight r %}
plotPCA(rlog(dds))
{% endhighlight %}
adar1_null and the control cells are significantly different under the stimulation. Without stimulation, they are not differentiable, but still different from any cells with stimulation.  


# Visualization 
## MA Plots
{% highlight r %}
plotMA(res_adarsg1_ifnb_VS_ctrl_ifnb)
{% endhighlight %}
All the blue dots represent differently expressed genes. We can see that there are quite amount of genes that are either upregulated or downregultaed in the adar1 null cells.    

## Heatmap
{% highlight r %}
# make a copy of res_adarsg1_ifnb_VS_ctrl_ifnb
res_adarsg1_ifnb_VS_ctrl_ifnb.df <- as.data.frame(res_adarsg1_ifnb_VS_ctrl_ifnb)

# only select genes that have high expression across all samples, 
# those with a high fold change, 
# and those  with significant p-values
sig.df <- res_adarsg1_ifnb_VS_ctrl_ifnb.df[(res_adarsg1_ifnb_VS_ctrl_ifnb.df$baseMean > 150) & (abs(res_adarsg1_ifnb_VS_ctrl_ifnb.df$log2FoldChange) > 3) & res_adarsg1_ifnb_VS_ctrl_ifnb.df$padj < 0.001, ]

# matrix for the heatmap
mat <- counts(dds, normalized = T)[rownames(sig.df),]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(metaTable)

heatmap.sig <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), name = "Z-score", row_labels = rownames(sig.df[rownames(mat.z), ]))

heatmap.sig

# print the figure
png('sigHeatmap.png', res = 250, width = 1500, height = 4100)
print(heatmap.sig)
dev.off()

{% endhighlight %}

To have a readable plot, I only select genes that have a high base mean, and a high log fold change. Those darker blocks represent higher fold of change. We may notice that, genes are likely to be differentially expressed under ifnb stimulation. And adar null cells have more differently expressed genes than control cells.  


# References
1. The raw count file is retrieved from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110708
2. Youtube: https://www.youtube.com/watch?v=OzNzO8qwwp0
3. Evans, Ciaran et al. “Selecting between-sample RNA-Seq normalization methods from the perspective of their assumptions.” Briefings in bioinformatics vol. 19,5 (2018): 776-792. doi:10.1093/bib/bbx008
4. http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow


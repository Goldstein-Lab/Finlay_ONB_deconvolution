#Load packages
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("gplots")
library('ggpubr')
library(viridis)
library(ggplot2)
library(ggrepel)
library(pheatmap)

#Set working directory
setwd("~/Library/CloudStorage/Box-Box/Goldstein Lab/Jack/ENB_Project/Malouf ENB Bulk RNA-Seq")

#Read bulk RNA-Seq samples
#This is .csv with genes in first column, each bulk sample in its own column
data <- read.csv("19_ENB_plus_3_control_filtered.csv", header = T, row.names = "Gene")

# identify conditions
#By grade
condition <- c(rep('Grade I/II', 9), rep("Grade III/IV", 10), rep("Control OE", 3))

# assign condition to each sample
my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(data)
my_colData

#Create DESeq2 Object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = my_colData,
                              design = ~condition)

#Run Deseq
#See here for code samples: https://erilu.github.io/bulk-rnaseq-analysis/
dds <- DESeq(dds)
#This normalizes entire dds object
dds <- estimateSizeFactors(dds)
dds

#Visualizations of data
vsd <- vst(dds, blind = TRUE)

#PCA plot
plot_PCA = function (vsd.obj) {
  pcaData <- plotPCA(vsd.obj,  intgroup = c("condition"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
         title = "PCA Plot colored by Tumor") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")} + theme_classic()+ theme(
      plot.title = element_text(hjust = 0.5, size=20)) 

plot_PCA(vsd)

#Boxplot for specific gene comparisons

#Set comparisons
my_comparisons <- list(c("Control OE", 'Grade I/II'), c('Control OE', 'Grade III/IV'), c('Grade I/II', 'Grade III/IV'))

Gene='EZH2'
plotCounts(dds, gene=Gene, intgroup = "condition", returnData = TRUE) %>%
  ggplot() + aes(condition, count) + geom_jitter() + geom_boxplot(aes(fill=condition), outlier.shape = NA, lwd=1) + theme_classic() + xlab(
    '') + ylab('Normalized Counts') + ggtitle(Gene) + theme(plot.title = element_text(size=20, hjust = 0.5),
   legend.position="none", text=element_text(size=25, color='black'), axis.text.x=element_text(angle=60, hjust=1)) + scale_fill_viridis(discrete=TRUE,
      alpha=0.3, option='C', direction=-1)  + stat_compare_means(comparisons = my_comparisons, 
bracket.size = 0.8, size=5, tip.length = 0.03, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1, 1.5), 
                                                                  symbols = c("****", "***", "**", "*", "ns", "ns")                                                                                                                                                                                                                                   )))
#  To add Kruskal Wallis on graph, use: '+  stat_compare_means(label.y = FALSE)' 

#To run with specific boxplot points labeled
plotCounts(dds, gene=Gene, intgroup = "condition", returnData = TRUE) %>%
  ggplot(aes(label=row.names(my_colData))) + geom_jitter(width = .1)+  aes(condition, count) + geom_boxplot(aes(fill=condition), outlier.shape = NA, lwd=1) + theme_classic() + xlab(
    'Tissue') + ylab('Normalized Counts') + ggtitle(Gene) + theme(plot.title = element_text(size=20, hjust = 0.5),
                                                                  legend.position="none", text=element_text(size=25, color='black')) + scale_fill_viridis(discrete=TRUE, alpha=0.3, option='E', 
                                                                                                                                                          direction=-1) + geom_label_repel()





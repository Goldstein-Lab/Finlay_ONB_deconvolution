#Load packages
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library('survMisc')
library(data.table)
library(dplyr)
library("DESeq2")
library("tidyverse")
library("gplots")
library('ggpubr')
library(ggplot2)
library(ggrepel)


#Setwd
setwd("~/Library/CloudStorage/Box-Box/Goldstein Lab/Jack/ENB_Project/Malouf ENB Bulk RNA-Seq")

#Read in data
#This is a .csv with overall survival and progression free survival for each patient
clin_df_ENB= read.csv('ENB_clin_df.csv', row.names=1)

head(clin_df_ENB)

#Can calculate quick kaplan Meier survival
Surv(clin_df_ENB$overall_survival, clin_df_ENB$deceased)

#Test with simple variable like Hyams_Grade_simple
Surv(clin_df_ENB$overall_survival, clin_df$deceased) ~ clin_df_ENB$Hyams_Grade_simple

#Fit survival model
fit = survfit(Surv(overall_survival, deceased) ~ Hyams_Grade_simple, data=clin_df_ENB)
print(fit)

#Plot Kaplan Meier
ggsurvplot(fit, data=clin_df_ENB, pval=T, conf.int=F)

#Optionally, can fit a Cox proportional hazards model
clin_df_ENB_dt=data.table(clin_df_ENB)

#Can specify which variables you want to regress
fit.coxph <- coxph(Surv(overall_survival, deceased) ~ Hyams_Grade_simple + Dulguerov_T_stage_simple, 
                   data = clin_df_ENB_dt)
#Show results
fit.coxph

#Plot results
ggforest(fit.coxph, data = clin_df_ENB_dt)

#####
#For survival based on gene expression data

#Load data
#This is unprocessed counts matrix
#No normal OE
data <- read.csv("19_samples_filtered.csv", header = T, row.names = "Gene")

head(data)

#Identify conditions, here by grade
condition <- c(rep('Grade I/II', 9), rep("Grade III/IV", 10))

# assign condition to each sample
my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(data)
my_colData

#Create DESeq2 Object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = my_colData,
                              design = ~condition)

#Run DeSeq2
dds <- DESeq(dds)
#This normalizes entire dds object
dds <- estimateSizeFactors(dds)
dds

#Normalize counts (this is internal normalization across each gene)
normalized_counts <- counts(dds, normalized = T)
head(normalized_counts)

#Check df
print(normalized_counts[1, ])
print(normalized_counts[ , 1])

#Specify gene
gene_id='TRPM5'

# get the expression values for the selected gene
clin_df_ENB$gene_value = normalized_counts[gene_id, rownames(clin_df_ENB)]

# find the median value of the gene and print it
median_value = median(clin_df_ENB$gene_value)
print(median_value)

#Divide patients into two groups based on gene expression
clin_df_ENB$gene = ifelse(clin_df_ENB$gene_value >= median_value, "UP", "DOWN")

# Fit survival model
fit = survfit(Surv(overall_survival, deceased) ~ clin_df_ENB$gene, data=clin_df_ENB)

#or if you want to use progession (instead of survival)
fit = survfit(Surv(overall_progression, progression) ~ clin_df_ENB$gene, data=clin_df_ENB)


#Extract the survival p-value
pval = surv_pvalue(fit, data=clin_df_ENB)$pval
print(pval)

#Plot with P-value
ggsurvplot(fit, data=clin_df_ENB, pval=F, risk.table=F, title=paste(gene_id), ylab='Survival Rate', xlab='Months',
           break.x.by= 25,
           legend.title='Expression', legend.labs=c('Low', 'High'), palette = 'lancet',
           conf.int = F, pval.coord=c(25,.2), pval.method = F, ggtheme = theme_classic2(base_size=20),
           font.tickslab = c(18))






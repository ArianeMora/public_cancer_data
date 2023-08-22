#install.packages('tidyverse')
#BiocManager::install("EnhancedVolcano")

#BiocManager::install("clusterProfiler")
#BiocManager::install("fgsea")
#BiocManager::install("enrichplot")
#BiocManager::install("GSEABase")
#BiocManager::install("missMethyl")
#BiocManager::install("minfi")
#BiocManager::install("DESeq2")

  
library("limma")
library("dplyr")
library(matrixStats)
library(DMRcate)
source("helper.R")


library(DESeq2)
library(edgeR)
require("impute")
require("remotes")
# You need to install via the below link
#install_github("WangLab-MSSM/DreamAI/Code")

library("DreamAI")


project_dir <- ''
cancer <- 'PanCan'
sample_file <- paste0( project_dir, cancer, '_Protein_samples.csv') 
imputed_file <- paste0(project_dir, cancer, '_Protein.csv')

output_file <- paste0(project_dir, cancer, '_filtered_DA_Protein.csv')


CondId <- 'CondID'
FullLabel <- 'Sample'
CaseId <- 'SafeCases'
paired=FALSE
gene = 'gene_name'

# Do protein analysis on the stage 4 sample

prot_data <- read.csv(imputed_file)
gene_names <- prot_data[['X']]
experimental_design <- read.csv(sample_file)
data_columns <- experimental_design[[FullLabel]] # get LFQ column numbers
prot_data_mat <- prot_data[, data_columns]

rownames(prot_data_mat) <- gene_names
rownames(prot_data) <- gene_names


# Get the column of interest from the experimental design
cond <- as.factor(experimental_design[[CondId]])
case_id <- as.factor(experimental_design[[CaseId]])
disease <- as.factor(experimental_design[['Disease.Type_protein']])

design <- model.matrix(~disease + cond) # We didn't have matching patient info
# Limma is good for detecting differentially abundent proteins
fit <- lmFit(prot_data_mat, design)
fit2 <- eBayes(fit, robust=TRUE) # Use robust:  https://www.biostars.org/p/496806/, https://support.bioconductor.org/p/118495/

# Keep in mind the condition of interest (tumour vs normal) is always the final column in the design
fit2_adj <- topTable(fit2, coef=length(colnames(design)), adjust="fdr", sort.by="none", number=1000000) 

# Create the dataframe to return
all_prot_df <- data.frame(prot_data[rownames(fit2_adj), colnames(prot_data)])
prot_data_mat <- prot_data_mat[rownames(fit2_adj), colnames(prot_data_mat)]
# Add in the statistics from the DA analysis
all_prot_df$gene_name <- rownames(fit2_adj)
all_prot_df$logFC_protein <- fit2_adj$logFC
all_prot_df$stat_protein <- fit2_adj$t
all_prot_df$pvalue_protein <- fit2_adj$P.Value
all_prot_df$padj_protein <- fit2_adj$adj.P.Val
all_prot_df$B_protein <- fit2_adj$B
all_prot_df$mean_protein <- fit2_adj$AveExpr

write.csv(all_prot_df, output_file, row.names = FALSE)




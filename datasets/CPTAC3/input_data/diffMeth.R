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


project_dir <- '../output_data/'
data_file <- paste0(project_dir, 'pancan_filtered_RNA.csv')
sample_file <- paste0( project_dir, 'pancan_filtered_samples_RNA.csv') 
output_file <- paste0(project_dir, 'pancan_filtered_DE_RNA.csv')
paired=FALSE


cat(paste("Differential Expression analysis for: \n", data_file, "\n"))

counts <- read.csv(data_file, header = TRUE, sep = ",")
rownames(counts) <- counts$gene_id

experimental_design <- read.csv(sample_file)
other_info <- counts[, c('gene_name', 'gene_id')] # Keep the metadata info
rownames(other_info) <- rownames(counts)
# Let's make sure our count data is in matrix format and is only the numeric columns i.e. everything but the genes
#nn <- counts[, !(names(counts) %in% experimental_design$Sample)]
experimental_design <- experimental_design[experimental_design$Sample %in% names(counts), ]

count_matrix <- as.matrix(counts[,  experimental_design$Sample])

count_matrix[is.na(count_matrix)] = 0
# We now set the row names to be the gene IDs
rownames(count_matrix) <- rownames(counts) 

# Separate out each of the sample names to be the different experiment conditions
condition_id <- as.factor(experimental_design$CondID)
case_id <- as.factor(experimental_design$SafeCases)
disease <- as.factor(experimental_design$Disease)
sample_df = data.frame(case_id = case_id, disease=disease, condition_id=condition_id)

# Before getting the means we want to normalise the counts (i.e. so our mean isn't stupidly big)
dge <- DGEList(counts=count_matrix)
dge <- calcNormFactors(dge, method="TMM")
tmm <- cpm(dge)


dds_mat <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = sample_df,
                                  design = ~disease + condition_id) # Just do on condition


dds <- estimateSizeFactors(dds_mat)

# Already did pre-filtering before
dds <- dds_mat
other_info_filtered <- other_info
counts_filtered <- counts
tmm_filtered <- tmm

# Log2 the TMM counts for better variance & mean estimation.
tmm_filtered <- log2(tmm_filtered + 1)
# Let's print the number of rows
cat(paste("Dataset dimensions: ", nrow(dds), ncol(dds), "\n"))
# Run DEseq2
dds <- DESeq(dds)
resultsNames(dds)
# Build results table
res <- results(dds, independentFiltering=FALSE)
other_info_filtered <- other_info_filtered[rownames(res), c('gene_name', 'gene_id')]

# Ensure ordering is correct 
tmm_filtered <- tmm_filtered[rownames(res), colnames(tmm_filtered)]
other_info_filtered$logFC_rna <- res$log2FoldChange
other_info_filtered$stat_rna <- res$stat
other_info_filtered$pvalue_rna <- res$pvalue
other_info_filtered$padj_rna <- res$padj
other_info_filtered$lfcSE_rna <- res$lfcSE
other_info_filtered$baseMean_rna <- res$baseMean
other_info_filtered$var_rna <- matrixStats::rowVars(as.matrix(tmm_filtered))

# Add in mean info
all_rna_df <- cbind(other_info_filtered, tmm_filtered)

write.csv(all_rna_df, file = output_file)


project_dir <- '../output_data/'
data_file <- paste0(project_dir, 'pancan_filtered_CpG.csv')
sample_file <- paste0( project_dir, 'pancan_filtered_samples_CpG.csv') 
output_file <- paste0(project_dir, 'pancan_filtered_DMC_CpG.csv')
paired=FALSE

cat(paste("Differential Methylation analysis for: \n", data_file, "\n"))

#### Import data ####
cpg_raw <- read.csv(data_file, header = TRUE, sep = ",")
sample_df <- read.csv(sample_file)

#### Change rownames ####
rowNames <- unlist(cpg_raw['id'])
sample_df <- sample_df[sample_df$Sample %in% names(cpg_raw), ]
cpg_data <- cpg_raw[, sample_df$Sample]
rownames(cpg_data) <- rowNames
rownames(cpg_raw) <- rowNames
#### QC/Filtering
# First remove all NA values
cpg_data[is.na(cpg_data)] <- 0
summary(cpg_data)

# Convert data to a matrix & plot the means for the rows
cpg_data_m <- as.matrix(cpg_data)

# The function model.matrix is used to generate the design matrix
cond_id <- as.factor(sample_df$CondID)
cases <- as.factor(sample_df$SafeCases)
disease <- as.factor(sample_df$Disease)
design = model.matrix(~disease + cond_id)  # Can't do paired in DNA methylation

# Before running limma we want to do two steps, (following the steps of Miss Methyl)
# 1) convert to M values 
# 2) perform limma
# References: 
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
# https://bioconductor.org/packages/release/bioc/vignettes/missMethyl/inst/doc/missMethyl.html 
# Add a very small amount so that if we have 0's we don't get an issue with NAs and log2
cpg_data_M_filtered <- log2((cpg_data_m) /(1-cpg_data_m))
# Normalise M values using https://github.com/regRCPqn/regRCPqn, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0229763
# Install: https://github.com/regRCPqn/regRCPqn
coef_id <- length(colnames(design))
fit <- lmFit(cpg_data_M_filtered, design)
fit2 <- eBayes(fit, robust=TRUE)

# Don't sort it otherwise we won't be able to match the data
fit2_adj <- as.data.frame(topTable(fit2, coef=coef_id, adjust="fdr", sort.by="none", number=1000000))
all_cpg_df <- cpg_raw[rownames(fit2_adj), c('id')]

# Add in the statistics from the CpG analysis
all_cpg_df$beta_logFC_meth <- fit2_adj$logFC
all_cpg_df$beta_stat_meth <- fit2_adj$t
all_cpg_df$beta_pvalue_meth <- fit2_adj$P.Value
all_cpg_df$beta_padj_meth <- fit2_adj$adj.P.Val
all_cpg_df$beta_B_meth <- fit2_adj$B
all_cpg_df$beta_mean_cpg <- fit2_adj$AveExpr

write.csv(fit2_adj, output_file)


# Do the ORA

GSEA_file_dir <- '../data/GSEA/'
#Import the pathways:
KEGG <- read.csv(file.path(GSEA_file_dir,  "c2.cp.kegg.v6.2.symbols.csv"))
Reactome <- read.csv(file.path(GSEA_file_dir,  "c2.cp.reactome.v6.2.symbols.csv"))
Biocarta <- read.csv(file.path(GSEA_file_dir, "c2.cp.biocarta.v6.2.symbols.csv"))
Hallmarks <- read.csv(file.path(GSEA_file_dir,  "h.all.v6.2.symbols.csv"))
GO_BP <- read.csv(file.path(GSEA_file_dir, "c5.go.bp.v7.2.symbols.csv"))
GO_CC <- read.csv(file.path(GSEA_file_dir, "c5.go.cc.v7.2.symbols.csv"))
GO_MF <- read.csv(file.path(GSEA_file_dir, "c5.go.mf.v7.2.symbols.csv"))
Metabolic <- read.csv(file.path(GSEA_file_dir, "41467_2016_BFncomms13041_MOESM340_ESM.csv"))
Correction <- read.csv(file.path(GSEA_file_dir, "41467_2016_BFncomms13041_MOESM340_ESM.csv"))

#Run the GSEA analysis
##1."KEGG", "Reactome", "Biocarta", "Hallmarks"
pathways <- rbind(KEGG, Reactome, Biocarta, Hallmarks)
pathway_list <- list()
for (pathway in unique(pathways$term)) {
  pathway_list[[pathway]] <- as.character(pathways[pathways$term == pathway, 1])
}

sircleORAHuman <- function(filename, entrezId, title, regLabels="RegulatoryLabels", emptyRegLabel="", fileType="pdf",
                           minGSSize=10, qvalueCutoff=0.2, pvalueCutoff=0.05, showCatagory=30, outputFolder='') {
  ## ------------ Setup and installs ----------- ##
  packages <- c("org.Hs.eg.db", "clusterProfiler", "svglite", "enrichplot")
  install.packages(setdiff(packages, rownames(installed.packages())))
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(svglite)
  library(enrichplot)
  ## ------------ Run ----------- ##
  # open the data
  df <- read.csv(filename)
  allGenes <- as.character(df[[entrezId]]) #
  clusterGenes <- subset(df, ! df[[regLabels]] == emptyRegLabel)
  grps_labels <- unlist(unique(clusterGenes[regLabels]))
  for(g in grps_labels) {
    grpGenes <- subset(df, df[[regLabels]] == g)
    print(g)
    print(dim(grpGenes))
    clusterGo <- enrichGO(gene = as.character(grpGenes[[entrezId]]),
                          universe = allGenes,
                          keyType = "ENTREZID",
                          OrgDb = org.Hs.eg.db,
                          ont = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1.0,
                          pvalueCutoff = 1.0,
                          minGSSize = minGSSize,
                          readable = TRUE)
    # We have a cutoff of all of them, and then only visualise the ones that the user wants...
    
    clusterGoSummary <- data.frame(clusterGo)
    write.csv(clusterGoSummary, paste(outputFolder, 'ClusterGoSummary_', g, title, '.csv', sep=""))#Export the ORA results as .csv
    
    if (!(dim(clusterGoSummary)[1] == 0)) {#exclude df's that have no observations
      Dotplot <- dotplot(clusterGo, showCategory=showCatagory) +
        ggtitle(paste("Dotplot ", g, sep=""))
      ggsave(file=paste(outputFolder, "SiRCle-ORA_Dotplot_Human_", g, title, ".", fileType, sep=""), plot=Dotplot, width=10, height=8)
      x2 <- pairwise_termsim(clusterGo)
      
      Emapplot <- emapplot(x2, pie_scale=1.5, layout = "nicely")+
        ggtitle(paste("Emapplot ", g, sep=""))
      ggsave(file=paste(outputFolder, "SiRCle-ORA_Emapplot_Human_", g, title, ".", fileType, sep="" ), plot=Emapplot, width=10, height=8)
      
      Heatplot <- heatplot(clusterGo,showCategory=showCatagory) +
        theme(axis.text.x =element_text(size=5), axis.text.y =element_text(size=8,face="bold"), axis.title=element_text(size=12,face="bold"))+
        ggtitle(paste("Heatplot ", g, sep=""))
      ggsave(file=paste(outputFolder, "SiRCle-ORA_Heatplot_Human_", g,title,  ".", fileType, sep="" ), plot=Heatplot, width=10, height=8)
      
    }
  }
}

input_data_dir <- '../manuscript_tables/'
output_data_dir <- '../ORA_figs/'
sircleFileName <- paste0(input_data_dir, 'Pred_genes_CPTAC.csv')
regLabel <- 'label'
test_title <- 'CPTAC'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, regLabel, emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder=output_data_dir)


# Do differential methylation analysis
project_dir <- '../filtered/'
output_dir <- '../DE/'


cptac_cancers = c(
  'Renal_cell_carcinoma__NOS',
  'Infiltrating_duct_carcinoma__NOS',
  'Adenocarcinoma__NOS',
  'Squamous_cell_carcinoma__NOS')

library("DreamAI")

# First impute the protein
for (c in cptac_cancers) {
  prot_data <- read.csv(paste0(project_dir, 'data_protein_', c, '.csv'), sep=',')
  prot_num_data <- prot_data[, 3:length(colnames(prot_data))]
  rownames(prot_num_data) <- prot_data$gene_name
  imputed_data <- DreamAI(prot_num_data, k = 10, maxiter_MF = 10, ntree = 100,
                          maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
                          gamma_ADMIN = NA, gamma = 50, CV = FALSE,
                          fillmethod = "row_mean", maxiter_RegImpute = 10,
                          conv_nrmse = 1e-06, iter_SpectroFM = 40, method = c("KNN",
                                                                              "MissForest", "Birnn", "SpectroFM", "RegImpute"),
                          out = c("Ensemble"))
  
  ens_data <- imputed_data$Ensemble
  # Save the imputed data to a csv for the remaining analysis
  write.csv(ens_data, paste0(project_dir, 'imputed_protein_', c, '.csv'))
}

for (c in cptac_cancers) {
  da <- pairedPatientDA(paste0(project_dir, 'imputed_protein_', c, '.csv'), 
                      paste0(project_dir, 'samples_Protein_', c, '.csv'), 
                      paste0(output_dir, 'DA_Protein_', c, '.csv'), 
                      gene='gene_name', FullLabel='FullLabel', CondId='CondID', CaseId='SafeCases', 
                      paired=FALSE)
}

for (c in cptac_cancers) {
  de <- pairedPatientDE(paste0(project_dir, 'data_RNA_', c, '.csv'), 
                        paste0(project_dir, 'samples_RNA_', c, '.csv'), 
                        paste0(output_dir, 'DE_RNA_', c, '.csv'),
                        paired=FALSE)
}

for (c in cptac_cancers) {
  dmc <- pairedPatientDMC(paste0(project_dir, 'data_CpG_', c, '.csv'), 
                          paste0(project_dir, 'samples_CpG_', c, '.csv'), 
                          paste0(output_dir, 'DCpG_CpG_', c, '.csv'), 
                          label=paste0('DMeth_', c),
                          project_dir = output_dir,
                          array_type='EPIC', # EPIC for CPTAC
                          paired=FALSE)
  
}


# Run ORA on each of the clusters
sircleFileName <- paste0('../data/sircle/PorMandR_Pancreas.csv')
regLabel <- 'RG2_Changes_filtered'
test_title <- '_Pancreas'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, regLabel, emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder='../data/sircle/Pancreas/')


# Run ORA on each of the clusters
sircleFileName <- paste0('../data/sircle/PorMandR_HeadNeck.csv')
regLabel <- 'RG2_Changes_filtered'
test_title <- '_HeadNeck'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, regLabel, emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder='../data/sircle/HeadNeck/')

# Run ORA on each of the clusters
sircleFileName <- paste0('../data/sircle/PorMandR_Kidney.csv')
regLabel <- 'RG2_Changes_filtered'
test_title <- '_Kidney'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, regLabel, emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder='../data/sircle/Kidney/')

# Run ORA on each of the clusters
sircleFileName <- paste0('../data/sircle/PorMandR_Lung.csv')
regLabel <- 'RG2_Changes_filtered'
test_title <- '_Lung'
sircleORAHuman(sircleFileName, "entrezgene_id", test_title, regLabel, emptyRegLabel="None", minGSSize=10, qvalueCutoff=0.1, outputFolder='../data/sircle/Lung/')



data_file <- '../data/CPTAC_Protein_Processed/Lung_data.csv'
sample_file <- '../data/CPTAC_Protein_Processed/Lung_samples.csv'
gene <- 'Gene'
FullLabel <- 'Full.Label'
CondId <- 'Sample.ID'
paired <- FALSE
cancers <- c( 'Liver')
project_dir <- '../data/CPTAC_Protein_Processed/'






project_dir <- '../data/Methylation/'
project_label <- 'CPTAC'
cancers <- c('Infiltrating_duct_carcinoma__NOS', 'Renal_cell_carcinoma__NOS', 'Squamous_cell_carcinoma__NOS', 'Adenocarcinoma__NOS')

for (cancer in cancers) {
  pairedPatientDMC(paste0(project_dir, 'CPTAC_matched_data_', cancer, '.csv'), 
                   paste0(project_dir, 'CPTAC_matched_samples_', cancer, '.csv'), 
                   paste0(project_dir, 'CPTAC_matched_dCpG_', cancer, '.csv'), 
                   label=paste0(project_label, '_Beta_', cancer),
                   project_dir = project_dir,
                   array_type='EPIC', # EPIC for CPTAC
                   paired=FALSE)
  
}

project_dir <- '../data/RNA/'
for (cancer in cancers) {
  pairedPatientDE(paste0(project_dir, 'CPTAC_matched_data_', cancer, '.csv'), 
                  paste0(project_dir, 'CPTAC_matched_samples_', cancer, '.csv'), 
                  paste0(project_dir, 'CPTAC_matched_DEseq2_', cancer, '.csv'), 
                  paired=FALSE)
}



project_dir <- '../data/Methylation/'
project_label <- 'TCGA'
cancers <- c('Clear_cell_adenocarcinoma__NOS', 'Papillary_adenocarcinoma__NOS', 'Hepatocellular_carcinoma__NOS',
             'Infiltrating_duct_carcinoma__NOS', 'Squamous_cell_carcinoma__NOS', 'Adenocarcinoma__NOS', 
             'Endometrioid_adenocarcinoma__NOS', 'Transitional_cell_carcinoma')

for (cancer in cancers) {
  pairedPatientDMC(paste0(project_dir, 'TCGA_matched_data_', cancer, '.csv'), 
                   paste0(project_dir, 'TCGA_matched_samples_', cancer, '.csv'), 
                   paste0(project_dir, 'TCGA_matched_dCpG_', cancer, '.csv'), 
                   label=paste0(project_label, '_', cancer),
                   project_dir = project_dir,
                   array_type='450K', # EPIC for CPTAC 450K for TCGA
                   paired=FALSE)
  
}

project_dir <- '../data/RNA/'
for (cancer in cancers) {
  pairedPatientDE(paste0(project_dir, 'TCGA_matched_data_', cancer, '.csv'), 
                  paste0(project_dir, 'TCGA_matched_samples_', cancer, '.csv'), 
                  paste0(project_dir, 'TCGA_matched_DEseq2_', cancer, '.csv'), 
                  paired=FALSE)
}


project_dir <- '../data/Methylation/'
project_label <- 'TCGA'
cancers <- c('TCGA-THCA', 'TCGA-ESCA', 'TCGA-LIHC', 'TCGA-BRCA', 'TCGA-LUAD',
             'TCGA-BLCA', 'TCGA-KIRP', 'TCGA-COAD', 'TCGA-UCEC', 'TCGA-HNSC',
             'TCGA-PRAD', 'TCGA-KIRC')

for (cancer in cancers) {
  pairedPatientDMC(paste0(project_dir, 'TCGA_matched_data_', cancer, '.csv'), 
                   paste0(project_dir, 'TCGA_matched_samples_', cancer, '.csv'), 
                   paste0(project_dir, 'TCGA_matched_dCpG_', cancer, '.csv'), 
                   label=paste0(project_label, '_', cancer),
                   project_dir = project_dir,
                   array_type='450K', # EPIC for CPTAC 450K for TCGA
                   paired=FALSE)
  
}

project_dir <- '../data/RNA/'
for (cancer in cancers) {
  pairedPatientDE(paste0(project_dir, 'TCGA_matched_data_', cancer, '.csv'), 
                  paste0(project_dir, 'TCGA_matched_samples_', cancer, '.csv'), 
                  paste0(project_dir, 'TCGA_matched_DEseq2_', cancer, '.csv'), 
                  paired=FALSE)
}








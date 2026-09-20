##############################################################################################################################
######################################################## FIGURE 2B ###########################################################
############################################ DEA HEATMAPS ACROSS DATASETS/LINEAGES ###########################################
##############################################################################################################################

library(Seurat)
library(limma)
library(tidyverse)
library(edgeR)

xpand <- readRDS(file = '/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/xpand_final_annot.rds')
ainci <- readRDS(file = '/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/ainci_final_annot.rds')
li <- readRDS(file = '/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/li_final_annot.rds')
hca <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_integrated_Marcal_analysis.rds')

# For HCA, keep only early progrenitor cells 
colnames(hca@meta.data)[colnames(hca@meta.data) == "Donor"] <- "sample"
hca <- subset(hca, subset = annotation_2 %in% c("HSC", "MPP", "CLP", "ProB", "PreB", "MEP", "ERP", "MDP", "GMP", "Eo/B/Mast"))
table(pseudobulk$annotation_2)
# Define H1 genes
features <- c('H1-0', 'H1-1', 'H1-2', 'H1-3', 'H1-4', 'H1-5', 'H1-10')

# Pseudobulking data
pseudobulk <- AggregateExpression(hca, return.seurat = TRUE, group.by = c("annotation_2", "sample"))

# Count number of cells per sample and annotation_2
cell_counts <- as.data.frame(table(hca$sample, hca$annotation_2))
cell_counts$pseudobulks <- paste(cell_counts$Var2, cell_counts$Var1, sep = '_')

# Filter out pseudobulks with less than 40 cells and cell populations with less than 3 valid pseudobulks
valid_pseudobulks <- subset(cell_counts, Freq >= 30) # Filter pseudobulks with at least 40 cells
valid_counts_per_population <- aggregate(Freq ~ Var2, data = valid_pseudobulks, FUN = length) # Count valid pseudobulks per cell population
valid_populations <- valid_counts_per_population$Var2[valid_counts_per_population$Freq >= 3] # Identify cell populations with at least 3 valid pseudobulks
final_pseudobulks <- subset(valid_pseudobulks, Var2 %in% valid_populations) # Filter valid cell populations for DEA

# filtering valid pseudobulks and extracting counts from the Seurat object 
pseudobulk_filt <- subset(pseudobulk, subset = orig.ident %in% final_pseudobulks$pseudobulks)
raw_counts <- GetAssayData(pseudobulk_filt, layer = "counts")

# Extract metadata
meta <- pseudobulk_filt@meta.data
meta <- meta[colnames(raw_counts), , drop = FALSE]  # Ensure metadata matches columns
meta$CellType <- factor(meta$annotation_2, levels = unique(meta$annotation_2))  # Convert to factor

# Create design matrix
design <- model.matrix(~ 0 + CellType, data = meta)  
colnames(design) <- levels(meta$CellType)  # Rename columns to cell type names
colnames(design) <- make.names(colnames(design))
cell_types <- colnames(design)

# limma preprocessing
dge <- DGEList(raw_counts, group = meta$CellType )
dge <- calcNormFactors(dge)
x <- new("EList")
x$E <- cpm(dge, log = TRUE, prior.count = 3)
# save DGE
# saveRDS(dge, file = "/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/DEA/DEA results/hca_dge_1v1.rds")
dge <- readRDS(file = "/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/DEA/DEA results/hca_dge_1v1.rds")

# Define the lineages
lineages <- list(
  "Lymphoid" = c("HSC", "MPP", "CLP", "ProB", "PreB"),
  "Erythroid" = c("HSC", "MPP", "MEP", "ERP"),
  "Monocyte-Dendritic" = c("HSC", "MPP", "MDP"),
  "Granulocyte" = c("HSC", "MPP", "GMP", "`Eo/B/Mast`")
)

# IMPORTANT
# lineages <- list( # RUN THIS ONLY FOR HCA DATASET, since Eo/B/Mast did not have enough pseudobulks to be included in the DEA analysis
#   "Lymphoid" = c("HSC", "MPP", "CLP", "ProB", "PreB"),
#   "Erythroid" = c("HSC", "MPP", "MEP", "ERP"),
#   "Monocyte-Dendritic" = c("HSC", "MPP", "MDP"),
#   "Granulocyte" = c("HSC", "MPP", "GMP")
# )

# Generate unique pairwise contrasts within each lineage
contrast_pairs <- unique(unlist(lapply(lineages, function(lineage) {
  combn(lineage, 2, FUN = function(pair) {
    paste(pair[2], pair[1], sep = "-")  # Ensure later cell type (more differentiated) is first
  })
})))

# Create contrast matrix
contrast_matrix <- makeContrasts(contrasts = contrast_pairs, levels = design)

# Running limma-trend
fit <- lmFit(x, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit, trend = TRUE, robust = TRUE)

# Extract results for all comparisons
dea_results <- lapply(colnames(contrast_matrix), function(contrast) {
  res <- topTable(fit2, coef = contrast, number = Inf, adjust.method = "fdr")  # All genes
  res$Gene <- rownames(res)
  res$Contrast <- contrast
  return(res)
}) %>% bind_rows()

# Filter for H1 genes
dea_h1_results <- dea_results %>%
  filter(Gene %in% features) %>%
  select(Gene, Contrast, logFC, adj.P.Val)


dea_results_filtered <- dea_h1_results %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.5)
# write.csv(dea_h1_results, file = "/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/DEA/DEA results/dea_limma_hca_AllvsAll.csv", row.names = FALSE)
# dea_h1_results <- read.table("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/DEA/DEA results/dea_limma_xpand_AllvsAll.csv", sep = ",", header = TRUE)
library(Seurat)

# 1. Merge datasets

integrated <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/integration_post_gene_fix.rds')
hca_dataset <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_integrated_Marcal_analysis.rds')
colnames(hca_dataset@meta.data)[colnames(hca_dataset@meta.data) == "percent_H1"] <- "percent.h1"
colnames(hca_dataset@meta.data)[colnames(hca_dataset@meta.data) == "Donor"] <- "sample"

# Fill NA values in annotation_2 with CellType values 
# (NAs  correspond to mature cell populations that had annotation only in the CellType column, that contains  original anotation by the authors)
hca_dataset$annotation_2 <- ifelse(
  is.na(hca_dataset$annotation_2),
  hca_dataset$CellType,
  hca_dataset$annotation_2
)

# Keep common columns and merge the two objects
common_cols <- intersect(colnames(integrated@meta.data), colnames(hca_dataset@meta.data))
integrated@meta.data <- integrated@meta.data[, common_cols]
hca_dataset@meta.data <- hca_dataset@meta.data[, common_cols]
all_datasets <- merge(integrated, hca_dataset)
all_datasets <- JoinLayers(all_datasets)

# 2. Compute percentage of individual H1

# Generating the percentage of H1 subtypes for each cell in the merged dataset
genes <- c("H1-0", "H1-1", "H1-2", "H1-3", "H1-4", "H1-5", "H1-10")
for (gene in genes) {
  percent_col <- paste0("percent.", gsub("-", ".", gene))
  all_datasets[[percent_col]] <- PercentageFeatureSet(all_datasets, features = gene)
}

# saveRDS(all_datasets, "/ibex/user/sotoorda/masterh1/public_data_analysis/merged_datasets/objects/all_datasets_merged.rds")
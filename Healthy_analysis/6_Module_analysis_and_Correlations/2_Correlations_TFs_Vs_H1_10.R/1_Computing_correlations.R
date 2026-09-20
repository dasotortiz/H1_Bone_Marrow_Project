# ================================================================================
# In this script we compute correlations between H1 genes and all other genes
# across different datasets and lineages, using pseudobulked expression data
# Then we identify highly correlated genes and potential TF regulators of H1 genes
# ================================================================================

library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(Seurat)

# Define H1 genes
h1_genes <- c('H1-0', 'H1-1', 'H1-2', 'H1-3', 'H1-4', 'H1-5', 'H1-10')
dataset_order <- c("Xpand", "Ainciburu", "Li", "HCA")
dataset_colors <- c(Xpand = "#1B9E77", Ainciburu = "#D95F02", Li = "#7570B3", HCA = "#E7298A")

# Define lineages
lineage_defs <- list(
  early_progenitors = c("HSC", "MPP"),
  Lymphocytes = c("CLP", "ProB", "PreB"),
  Erythrocytes = c('MEP', "ERP"),
  Monocytes_Dendritic_cells = c("MDP"),
  Eo_B_Mast = c("GMP", 'Eo/B/Mast')
)

compute_h1_correlations <- function(expr, h1_genes, all_genes) {
  # Keep only numeric columns
  expr <- expr[, sapply(expr, is.numeric), drop = FALSE]
  # All combinations of H1 and TFs 
  cors <- expand.grid(H1 = h1_genes, Gene = colnames(expr), stringsAsFactors = FALSE)
  # Compute correlation and p-value for each pair
  cors_stats <- mapply(function(h1, g) {
      if (h1 %in% colnames(expr) && g %in% colnames(expr)) { 
          test <- cor.test(expr[[h1]], expr[[g]], method = "spearman", exact = FALSE) 
          c(Correlation = as.numeric(test$estimate), PValue = test$p.value) 
        } else {
         c(Correlation = NA, PValue = NA)
        } 
    }, h1 = cors$H1, g = cors$Gene)
  # mapply returns a matrix; transpose to assign columns 
  cors$Correlation <- cors_stats["Correlation", ] 
  cors$PValue <- cors_stats["PValue", ] 
  # Adjust p-values for multiple testing (FDR) for each 'Dataset x Lineage x H1' combination
  cors <- cors %>%
    group_by(H1) %>%
    mutate(FDR = p.adjust(PValue, method = "BH")) %>%
    ungroup()
  return(cors)
}

# Load integrated Seurat object
all_datasets <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/merged_datasets/objects/all_datasets_merged.rds")
# List of all TFs
tf_list <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/TF_analysis/TF_List/JASPAR2024_TFs.rds")

# Define parameters for filtering during pseudobulk aggregation
datasets <- unique(all_datasets$Dataset)
min_cells_per_sample <- 30
min_samples <- 3

# RUNNING CORRELATIONS
lineage_cor_results <- lapply(names(lineage_defs), function(lin_name) { # For each lineage
  celltypes <- lineage_defs[[lin_name]] 
  dataset_results <- lapply(datasets, function(ds) { # For each dataset
     message("Analyzing ", ds, " / ", lin_name)
   
     # Subset Seurat object by dataset and lineage
     subset_obj <- subset(all_datasets, subset = Dataset == ds & annotation_2 %in% celltypes)
   
     #REMOVE SAMPLES WITH < 30 CELLS
     cell_counts <- subset_obj@meta.data %>% dplyr::count(sample, annotation_2) # Count cells per sample and cell type
     valid_counts <- cell_counts %>%dplyr::filter(n >= min_cells_per_sample) # Keep only sample-celltype pairs with >= 30 cells
     valid_samples <- valid_counts %>% dplyr::group_by(sample) %>% 
       dplyr::summarise(n_celltypes = n_distinct(annotation_2), .groups = "drop") %>%
       dplyr::filter(n_celltypes == length(celltypes)) %>%
       dplyr::pull(sample) # Identify cells to keep based on valid sample-celltype pairs
     subset_obj <- subset(subset_obj, subset = sample %in% valid_samples)
   
     if (length(valid_samples) < min_samples) {
         message("Skipping: not enough samples (", min_samples, ") with >=30 cells in ALL cell types in ", ds, " / ", lin_name)
         return(NULL)
        }
   
     # Pseudobulk
     if (length(celltypes) > 1) { # Some lineages (MDP) have only 1 cell type, so we can't group by annotation_2
         agg_obj <- AggregateExpression(subset_obj, group.by = c("annotation_2", "sample"), return.seurat = TRUE)
         #Expression matrix
         all_genes <- c(tf_list, h1_genes)
         expr_data <- FetchData(agg_obj,vars = c(all_genes, "annotation_2", "sample"))
         #Ensure ordering
         expr_data <- expr_data %>% mutate(annotation_2 = factor(annotation_2, levels = celltypes)) %>% arrange(annotation_2)
     } else {
         agg_obj <- AggregateExpression(subset_obj, group.by = "sample", return.seurat = TRUE)
         # Expression matrix
         all_genes <- c(tf_list, h1_genes)
         expr_data <- FetchData(agg_obj,vars = c(all_genes, "sample"))
     }
     # Correlations
     cors <- compute_h1_correlations(expr = expr_data, h1_genes = h1_genes, all_genes = all_genes)
     cors
  })

  names(dataset_results) <- datasets
  dataset_results
})

names(lineage_cor_results) <- names(lineage_defs)
# saveRDS(lineage_cor_results, "/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/Correlations/corr_early_progenitors_only_TFs.rds")



#============================================
# Filtering TFs from the correlation results
#============================================
str(lineage_cor_results)
# lineage_cor_results <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/Correlations/corr_early_progenitors_only_TFs.rds")

names(lineage_cor_results) <- names(lineage_defs)
lineage_cor_results <- lapply(lineage_cor_results, function(x) Filter(Negate(is.null), x)) # Remove NULLs from lineage results (In the earlier code, a dataset returns NULL when it has fewer than three valid samples)
# Convert nested list to data frame
cor_df <- imap_dfr(lineage_cor_results, function(lineage_list, lineage_name) {
  imap_dfr(lineage_list, function(df, dataset_name) {
    df %>%
      mutate(
        Lineage = lineage_name,
        Dataset = dataset_name
      )
  })
})


# Filter TFs in the correlation results
tf_list <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/TF_analysis/TF_List/JASPAR2024_TFs.rds")
tfs_cor <- cor_df[cor_df$Gene %in% tf_list, ]
# saveRDS(tfs_cor, "/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/Correlations/early_progenitors_out_TFs.rds")


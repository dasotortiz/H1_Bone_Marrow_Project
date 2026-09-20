library(dplyr)
library(purrr)
library(stringr)
library(ComplexHeatmap)
library(circlize)

# ====================================================
# CREATE HEATMAPS FOR SPECIFIC H1/LINEAGE COMBINATIONS
# ====================================================
rank_product_df <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/Correlations/rank_product_df.rds')
# Create separate data frames for each H1 × Lineage combination for visualization
# NOTE: The specific H1 × Lineage combinations chosen here were the ones containing most of the correlations
H1_10_df <- rank_product_df[rank_product_df$H1 == "H1-10" & rank_product_df$Lineage == "Lymphocytes", ]

# Get the original correlation data for filtering
tf_corr_orig <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/Correlations/early_progenitors_out_TFs.rds')
tf_corr_orig$Lineage <- ifelse(tf_corr_orig$Lineage == "Monocytes_Dendritic_cells", "Mono_DC",
                               ifelse(tf_corr_orig$Lineage == "early_progenitors", "EP", tf_corr_orig$Lineage))

# Define the H1-10 heatmap data and remove significant negative correlations
h1_target <- "H1-10"
lineage_target <- "Lymphocytes"
genes_to_remove <- tf_corr_orig %>%
  filter(
    H1 == h1_target,
    Lineage == lineage_target,
    !is.na(Correlation),
    !is.na(FDR),
    Correlation < 0,
    FDR < 0.05
  ) %>%
  distinct(Gene) %>%
  pull(Gene)

heatmap_df <- H1_10_df %>%
  filter(!Gene %in% genes_to_remove)

# Define annotation colors for datasets
ann_colors <- list(Dataset = dataset_colors)
font_family <- "Helvetica"
legend_font_size <- 12


# Get TFs in module H1-10 for Xpand to add them as additional annotation
tf_list <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/TF_analysis/TF_List/JASPAR2024_TFs.rds")
results <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/gene_modules_excluding_histones.rds')
xpand_lymphocytes <- results[['Xpand']]$gene_module_list$Lymphocytes
h1_location <- xpand_lymphocytes[xpand_lymphocytes$id %in% h1_genes,]
module_21 <- xpand_lymphocytes[xpand_lymphocytes$module == 21, ]
tfs_shared_in_modules <- intersect(module_21$id, tf_list)

# Loading TFs that have a Transcription Factor Binding Site (TFBS) in the promoter region of H1-10
TFBS_df <- read.csv("/ibex/user/sotoorda/masterh1/public_data_analysis/Accesibility_data/H1_10_exact_TF_motif_peak_positions.csv")
TFBS_list <- TFBS_df$TF

# Create PDF for all heatmaps
pdf('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_3/rank_product_heatmaps.pdf', 
    height = 90, width = 4.5)

genes_to_plot <- heatmap_df$Gene
  
  # Filter original correlation data for this H1 and lineage
  h1_lineage_data <- tf_corr_orig %>%
    filter(H1 == h1_target, Lineage == lineage_target, Gene %in% genes_to_plot)
  
  # Get sorted genes based on rank_product (best ranks first)
  genes_sorted <- heatmap_df %>%
    arrange(final_rank) %>%
    pull(Gene)
  
  # Create matrix
  combined_mat <- data.frame(Gene = genes_sorted)
  dataset_labels <- character()
  
  for (ds in dataset_order) {
    cor_vals <- h1_lineage_data %>%
      filter(Dataset == ds) %>%
      select(Gene, Correlation) %>%
      as.data.frame()
    
    col_data <- rep(NA, nrow(combined_mat))
    names(col_data) <- combined_mat$Gene
    
    if (nrow(cor_vals) > 0) {
      col_data[cor_vals$Gene] <- cor_vals$Correlation
    }
    
    combined_mat[[ds]] <- col_data
    dataset_labels <- c(dataset_labels, ds)
  }
  
  # Convert to matrix
  mat_combined <- as.matrix(combined_mat[, -1])
  rownames(mat_combined) <- genes_sorted
  
  # Build p-val matrix
  pval_mat <- matrix(NA, nrow = nrow(mat_combined), ncol = ncol(mat_combined), 
                     dimnames = dimnames(mat_combined))
  
  for (ds in dataset_order) {
    pvals <- h1_lineage_data %>%
      filter(Dataset == ds) %>%
      select(Gene, FDR)
    
    if (nrow(pvals) > 0) {
      pval_mat[pvals$Gene, ds] <- pvals$FDR
    }
  }
  
  sig_mat <- pval_mat < 0.05

  # Annotate genes shared with the module analysis and genes with an H1-10 TFBS.
  row_annotation <- rowAnnotation(
    Module_TF = factor(
      ifelse(genes_sorted %in% tfs_shared_in_modules, "Shared", "Not shared"),
      levels = c("Shared", "Not shared")
    ),
    H1_10_TFBS = factor(
      ifelse(genes_sorted %in% TFBS_df$TF, "TFBS", "No TFBS"),
      levels = c("TFBS", "No TFBS")
    ),
    col = list(
      Module_TF = c(Shared = "#414040", `Not shared` = "#e2e0e0"),
      H1_10_TFBS = c(TFBS = "#0077b6", `No TFBS` = "#e2e0e0")
    ),
    annotation_name_gp = gpar(
      fontsize = legend_font_size,
      fontface = "bold",
      fontfamily = font_family
    ),
    annotation_name_side = "top",
    annotation_name_rot = 90,
    annotation_legend_param = list(
      Module_TF = list(
        title = "Module TF",
        title_gp = gpar(fontsize = legend_font_size, fontface = "bold", fontfamily = font_family),
        labels_gp = gpar(fontsize = legend_font_size, fontfamily = font_family)
      ),
      H1_10_TFBS = list(
        title = "H1-10 TFBS",
        title_gp = gpar(fontsize = legend_font_size, fontface = "bold", fontfamily = font_family),
        labels_gp = gpar(fontsize = legend_font_size, fontfamily = font_family)
      )
    )
  )
  
  # Create heatmap
  hm <- Heatmap(
    mat_combined,
    name = "Correlation",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    show_row_names = TRUE,
    show_column_names = TRUE,
    show_row_dend = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    border = TRUE,
    left_annotation = row_annotation,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(sig_mat[i, j]) && sig_mat[i, j]) {
        grid.text(
          "*",
          x = x,
          y = y,
          gp = gpar(fontsize = 8, col = "black")
        )
      }
    },
    
    top_annotation = HeatmapAnnotation(
      Dataset = factor(dataset_labels, levels = dataset_order),
      col = ann_colors,
      annotation_legend_param = list(
        Dataset = list(
          at = dataset_order,
          labels = dataset_order,
          title_gp = gpar(fontsize = legend_font_size, fontface = "bold", fontfamily = font_family),
          labels_gp = gpar(fontsize = legend_font_size, fontfamily = font_family)
        )
      )
    ),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = legend_font_size, fontface = "bold", fontfamily = font_family),
      labels_gp = gpar(fontsize = legend_font_size, fontfamily = font_family)
    )
  )
  
  draw(hm, heatmap_legend_side = "right")

dev.off()

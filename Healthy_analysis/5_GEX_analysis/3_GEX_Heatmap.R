library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(tibble)

font_family <- "Helvetica"

# Genes to plot
features <- c("H1-0", "H1-10", "H1-1", "H1-2", "H1-3", "H1-4", "H1-5")

# Lineage definitions
lineages_by_dataset <- list(
  Xpand = list(
    Lymphoid = c("HSC", "MPP", "CLP", "ProB", "PreB"),
    Erythroid = c("HSC", "MPP", "MEP", "ERP"),
    `Monocyte-Dendritic` = c("HSC", "MPP", "MDP"),
    Granulocyte = c("HSC", "MPP", "GMP", "Eo/B/Mast")
  ),
  Ainciburu = list(
    Lymphoid = c("HSC", "MPP", "CLP", "ProB", "PreB"),
    Erythroid = c("HSC", "MPP", "MEP", "ERP"),
    `Monocyte-Dendritic` = c("HSC", "MPP", "MDP"),
    Granulocyte = c("HSC", "MPP", "GMP", "Eo/B/Mast")
  ),
  Li = list(
    Lymphoid = c("HSC", "MPP", "CLP", "ProB", "PreB"),
    Erythroid = c("HSC", "MPP", "MEP", "ERP"),
    `Monocyte-Dendritic` = c("HSC", "MPP", "MDP"),
    Granulocyte = c("HSC", "MPP", "GMP", "Eo/B/Mast")
  ),
  HCA = list(
    Lymphoid = c("HSC", "MPP", "CLP", "ProB", "PreB"),
    Erythroid = c("HSC", "MPP", "MEP", "ERP"),
    `Monocyte-Dendritic` = c("HSC", "MPP", "MDP"),
    Granulocyte = c("HSC", "MPP", "GMP")
  )
)

lineages <- lineages_by_dataset$Xpand

get_contrast_pairs <- function(lineage_set) {
  unique(unlist(lapply(lineage_set, function(lineage) {
    combn(lineage, 2, FUN = function(pair) paste(pair[2], pair[1], sep = "-"))
  })))
}

contrast_pairs_by_dataset <- lapply(lineages_by_dataset, get_contrast_pairs)
all_contrasts <- unique(unlist(contrast_pairs_by_dataset, use.names = FALSE))

# Dataset order and colors
dataset_order <- c("Xpand", "Ainciburu", "Li", "HCA")

# Load DEA results
h1_results <- list(
  Xpand = read.table("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/DEA/DEA results/dea_limma_xpand_AllvsAll.csv", sep = ",", header = TRUE),
  Ainciburu = read.table("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/DEA/DEA results/dea_limma_ainci_AllvsAll.csv", sep = ",", header = TRUE),
  Li = read.table("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/DEA/DEA results/dea_limma_li_AllvsAll.csv", sep = ",", header = TRUE),
  HCA = read.table("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/DEA/DEA results/dea_limma_hca_AllvsAll.csv", sep = ",", header = TRUE)
)

# Fix contrast label naming
h1_results <- lapply(h1_results, function(df) {
  df$Contrast <- gsub("Eo\\.B\\.Mast", "Eo/B/Mast", as.character(df$Contrast))
  df
})

# Helper to map contrasts to lineage
get_lineage <- function(contrast) {
  cells <- unlist(strsplit(contrast, "-"))
  if (contrast == "MPP-HSC") return("Pluripotent")
  for (lineage in names(lineages)) {
    if (all(cells %in% lineages[[lineage]])) return(lineage)
  }
  return(NA)
}

# Build long table with one row per Gene x Contrast x Dataset
sig_cutoff <- 0.05
lfc_cutoff <- 0.5

# Keep dataset-specific contrasts
dea_long <- bind_rows(lapply(dataset_order, function(dataset_name) {
  valid_contrasts <- contrast_pairs_by_dataset[[dataset_name]]
  h1_results[[dataset_name]] %>%
    filter(Gene %in% features, Contrast %in% valid_contrasts) %>%
    mutate(
      Dataset = dataset_name,
      Contrast = factor(Contrast, levels = valid_contrasts),
      is_significant = adj.P.Val < sig_cutoff & abs(logFC) > lfc_cutoff,
      column_id = paste(Contrast, Dataset, sep = "__")
    )
}))

# Use only columns that actually exist
column_metadata <- dea_long %>%
  distinct(Contrast, Dataset, column_id) %>%
  mutate(
    Lineage = vapply(as.character(Contrast), get_lineage, character(1)),
    Lineage = factor(Lineage, levels = names(lineage_colors)),
    Dataset = factor(Dataset, levels = dataset_order)
  ) %>%
  arrange(factor(Contrast, levels = unique(unlist(contrast_pairs_by_dataset, use.names = FALSE))), Dataset)

# Define column order and contrast split for heatmap
column_order <- column_metadata$column_id
contrast_split <- column_metadata$Contrast

# Create matrices for heatmap
lfc_wide <- dea_long %>%
  select(Gene, column_id, logFC) %>%
  pivot_wider(names_from = column_id, values_from = logFC, values_fill = NA_real_)

# Create a matrix of log fold changes with genes as rows and contrasts as columns
lfc_matrix <- as.matrix(column_to_rownames(lfc_wide, "Gene"))
lfc_matrix <- lfc_matrix[intersect(features, rownames(lfc_matrix)), column_order, drop = FALSE]

# Create a matrix of significance (TRUE/FALSE) for the asterics in the heatmap
sig_wide <- dea_long %>%
  select(Gene, column_id, is_significant) %>%
  pivot_wider(names_from = column_id, values_from = is_significant, values_fill = FALSE)
significance_matrix <- as.matrix(column_to_rownames(sig_wide, "Gene"))
significance_matrix <- significance_matrix[rownames(lfc_matrix), column_order, drop = FALSE]

# Add top annotations
slice_info <- column_metadata %>%
  distinct(Contrast) %>%
  mutate(Lineage = vapply(as.character(Contrast), get_lineage, character(1))) %>%
  arrange(factor(Contrast, levels = unique(column_metadata$Contrast)))

top_annotation <- HeatmapAnnotation(
  Lineage = column_metadata$Lineage,
  Dataset = column_metadata$Dataset,
  col = list(Lineage = lineage_colors, Dataset = dataset_colors),
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = font_family),
  annotation_legend_param = list(
    Lineage = list(
      title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = font_family),
      labels_gp = gpar(fontsize = 12, fontfamily = font_family)
    ),
    Dataset = list(
      title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = font_family),
      labels_gp = gpar(fontsize = 12, fontfamily = font_family)
    )
  )
)

pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_2/dea_limma_h1_three_datasets_complexheatmap.pdf", width = 15, height = 3)
Heatmap(
  lfc_matrix,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = contrast_split,
  column_order = column_order,
  column_gap = unit(1, "mm"),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  na_col = "grey95",
  row_names_gp = gpar(fontsize = 12, fontfamily = font_family),
  show_column_names = FALSE,
  column_title_gp = gpar(fontsize = 12, fontfamily = font_family),
  column_title_rot = 45,
  top_annotation = top_annotation,
  heatmap_legend_param = list(
    title = "LFC",
    title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = font_family),
    labels_gp = gpar(fontsize = 12, fontfamily = font_family)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (isTRUE(significance_matrix[i, j])) {
      asterisk_color <- ifelse(!is.na(lfc_matrix[i, j]) && lfc_matrix[i, j] < -1.5, "white", "black")
      grid.text("*", x = x, y = y, gp = gpar(fontsize = 10, fontface = "bold", fontfamily = font_family, col = asterisk_color))
    }
  }
)
dev.off()

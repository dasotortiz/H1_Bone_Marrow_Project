library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

colors_datasets <- c(
  Xpand = "#1B9E77",
  Ainciburu = "#D95F02",
  Li = "#7570B3",
  HCA = "#E7298A"
)

# Lineage colors
lineage_colors <- c(
  "Pluripotent" = "#6c757d",
  "Lymphoid" = "#1f77b4",
  "Erythroid" = "#d62728",
  "Monocyte-Dendritic" = "#2ca02c",
  "Granulocyte" = "#9467bd"
)

##############################################################################################################################
######################################################## FIGURE 2A ###########################################################
########################################### LINEPOTS OF H1 ACROSS DATASETS/LINEAGES ##########################################
##############################################################################################################################

lineages <- list(
  "Lymphoid" = c("HSC", "MPP", "CLP", "ProB", "PreB"),
  "Erythroid" = c("HSC", "MPP", "MEP", "ERP"),
  "Monocyte-Dendritic" = c("HSC", "MPP", "MDP"),
  "Granulocyte" = c("HSC", "MPP", "GMP", "Eo/B/Mast")
)

integrated <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/integration_post_gene_fix.rds')
hca_dataset <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_integrated_Marcal_analysis.rds')
colnames(hca_dataset@meta.data)[colnames(hca_dataset@meta.data) == "Donor"] <- "sample"

# Fill NA values in annotation_2 with CellType values (There were NAs previously in that column)
hca_dataset$annotation_2 <- ifelse(
  is.na(hca_dataset$annotation_2),
  hca_dataset$CellType,
  hca_dataset$annotation_2)

# Keep common columns and merge the two objects
common_cols <- intersect(colnames(integrated@meta.data), colnames(hca_dataset@meta.data))
integrated@meta.data <- integrated@meta.data[, common_cols]
hca_dataset@meta.data <- hca_dataset@meta.data[, common_cols]
all_datasets <- merge(integrated, hca_dataset)
all_datasets <- JoinLayers(all_datasets)
colnames(all_datasets@meta.data)

# Pseudobulk Aggregation
min_cells_per_sample <- 30
h1_genes <- c("H1-0", "H1-10", "H1-5")

# Filter out donor/celltype combinations with too few cells
cell_counts <- all_datasets@meta.data %>%
  count(sample, annotation_2) %>%
  filter(n >= min_cells_per_sample)

# Only keep cells that belong to a valid donor/celltype pair
all_datasets$keep <- paste0(all_datasets$sample, "_", all_datasets$annotation_2) %in%
                   paste0(cell_counts$sample, "_", cell_counts$annotation_2)
all_datasets <- subset(all_datasets, subset = keep == TRUE)

# Aggregate expression by annotation and sample
agg_obj <- AggregateExpression(all_datasets, group.by = c("annotation_2", "sample", "Dataset"), return.seurat = TRUE)

# Fetch the aggregated data
genes_to_fetch <- intersect(h1_genes, rownames(agg_obj))
pb_data <- FetchData(agg_obj, vars = c(genes_to_fetch, "annotation_2", "sample", "Dataset"))

# Convert to long format
df_long <- pb_data %>%
  pivot_longer(cols = any_of(h1_genes), names_to = "H1_variant", values_to = "expression") %>%
  rename(CellType = annotation_2)

# Create lookup table
lineage_df <- imap_dfr(lineages, ~ tibble(Lineage = .y,CellType = .x))

# Build ordered x-axis labels in the exact lineage order
plot_levels <- unlist(lapply(names(lineages), function(lineage_name) {
  paste(lineage_name, lineages[[lineage_name]], sep = "__")
}))

# Duplicate observations into every lineage they belong to
df_plot <- df_long %>%
  filter(H1_variant %in% h1_genes) %>%
  inner_join(lineage_df, by = "CellType") %>%
  mutate(
    H1_variant = factor(H1_variant, levels = c("H1-0", "H1-10", "H1-5")),
    Lineage = factor(Lineage, levels = names(lineages)),
    PlotCell = factor(paste(Lineage, CellType, sep = "__"), levels = plot_levels)
  )

summary_df <- df_plot %>%
  group_by(Dataset, Lineage, H1_variant, PlotCell, CellType) %>%
  summarize(
    mean_expr = mean(expression, na.rm = TRUE),
    sd_expr = sd(expression, na.rm = TRUE),
    n = n(),
    se = sd_expr / sqrt(n),
    ci = ifelse(n > 1, qt(0.975, df = n - 1) * se, NA_real_),
    .groups = "drop"
  )

p <- ggplot(summary_df,aes(x = PlotCell, y = mean_expr, color = Dataset, fill = Dataset, group = Dataset)) +
  geom_smooth(data = df_plot, aes(x = PlotCell, y = expression, color = Dataset, fill = Dataset, group = Dataset),
    method = "loess", se = TRUE, span = 0.7, linewidth = 1.4, alpha = 0.15) +
  geom_point(size = 3.5) +
  facet_grid(rows = vars(H1_variant), cols = vars(Lineage), scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_color_manual(values = colors_datasets) +
  scale_fill_manual(values = colors_datasets) +
  labs(x = NULL, y = "Normalized Expression") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 16),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))

pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_2/Line_plots_GEX_paper.pdf",
    width = 14, height = 10)
print(p)
dev.off()


##############################################################################################################################
######################################################## FIGURE 2B ###########################################################
############################################ DEA HEATMAPS ACROSS DATASETS/LINEAGES ###########################################
##############################################################################################################################

library(Seurat)
library(limma)
library(tidyverse)
library(ComplexHeatmap)
library(edgeR)
library(ggplotify)
library(circlize)
library(grid)

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

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tibble)
})

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



##############################################################################################################################
###################################################### FIGURE 2C-2D ##########################################################
############################################ DEA HEATMAPS ACROSS DATASETS/LINEAGES ###########################################
##############################################################################################################################


##### Line plots for B and T cell lineages in HCA dataset only, H1-10 variant

b_lineage <- c("HSC", "MPP", "CLP", "ProB", "PreB", "pre-PC", "Plasma Cell")
t_lineage <- c("HSC", "MPP", "CLP", "CD8 T-cell", "Naive T-cell")

# Run a function to build the lineage plots for HCA dataset only, H1-10 variant. Note: df_long is already defined in the previous code block and contains the pseudobulk expression data for all datasets.
build_hca_lineage_plot <- function(lineage_cells, output_path) {
  # Filter the data for HCA dataset, the specified lineage cells, and H1-10 variant
  df_plot_hca <- df_long %>%
    filter(
      Dataset == "HCA",
      CellType %in% lineage_cells,
      H1_variant == "H1-10"
    ) %>%
    mutate(Stage = as.numeric(factor(CellType, levels = lineage_cells)))

  # Create a summary dataframe for plotting by calculating mean expression and standard error for each stage
  summary_df_hca <- df_plot_hca %>%
    group_by(Dataset, Stage, CellType) %>%
    summarize(
      mean_expr = mean(expression, na.rm = TRUE),
      sd_expr = sd(expression, na.rm = TRUE),
      n = n(),
      se = sd_expr / sqrt(n),
      ci = ifelse(n > 1, qt(0.975, df = n - 1) * se, NA_real_),
      .groups = "drop"
    )
  # Build the ggplot for the HCA lineage plot
  p_hca <- ggplot(
    summary_df_hca,
    aes(x = Stage, y = mean_expr, color = Dataset, fill = Dataset, group = Dataset) +
    geom_smooth(data = df_plot_hca, aes(x = Stage, y = expression, color = Dataset, fill = Dataset), method = "loess", span = 0.5, degree = 1, se = TRUE, linewidth = 1.4, alpha = 0.15) +
    geom_point(size = 3.5) +
    scale_x_continuous(breaks = seq_along(lineage_cells), labels = lineage_cells) +
    scale_color_manual(values = colors_datasets) +
    scale_fill_manual(values = colors_datasets) +
    labs(x = "", y = "Normalized Expression") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20),
      axis.text.y = element_text(color = "black", size = 20),
      axis.title.y = element_text(size = 22),
      legend.title = element_blank(),
      legend.text = element_text(size = 20)
    )
  )
  pdf(output_path, width = 10, height = 6)
  print(p_hca)
  dev.off()
}

build_hca_lineage_plot(
  b_lineage,
  "/ibex/user/sotoorda/masterh1/public_data_analysis/Posters/EHA/draft_figures/Line_plots_GEX_HCA_B_lineage_H1-10.pdf"
)
build_hca_lineage_plot(
  t_lineage,
  "/ibex/user/sotoorda/masterh1/public_data_analysis/Posters/EHA/draft_figures/Line_plots_GEX_HCA_T_lineage_H1-10.pdf"
)
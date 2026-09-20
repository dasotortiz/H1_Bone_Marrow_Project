library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

integrated <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/integration_post_gene_fix.rds')
hca_dataset <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_integrated_Marcal_analysis.rds')
colnames(hca_dataset@meta.data)[colnames(hca_dataset@meta.data) == "Donor"] <- "sample"

# Fill NA values in annotation_2 with CellType values (There were NAs previously)
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

# Pseudobulk Aggregation
min_cells_per_sample <- 30
h1_genes <- c("H1-1", "H1-2", "H1-3", "H1-4")

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

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

lineages <- list(
  "Lymphoid" = c("HSC", "MPP", "CLP", "ProB", "PreB"),
  "Erythroid" = c("HSC", "MPP", "MEP", "ERP"),
  "Monocyte-Dendritic" = c("HSC", "MPP", "MDP"),
  "Granulocyte" = c("HSC", "MPP", "GMP", "Eo/B/Mast")
)

colors_datasets <- c(
  Xpand = "#1B9E77",
  Ainciburu = "#D95F02",
  Li = "#7570B3",
  HCA = "#E7298A"
)

# Create lookup table
lineage_df <- imap_dfr(lineages, ~
  tibble(
    Lineage = .y,
    CellType = .x
  )
)

# Build ordered x-axis labels in the exact lineage order
plot_levels <- unlist(lapply(names(lineages), function(lineage_name) {
  paste(lineage_name, lineages[[lineage_name]], sep = "__")
}))

# Duplicate observations into every lineage they belong to
df_plot <- df_long %>%
  filter(H1_variant %in% h1_genes) %>%
  inner_join(lineage_df, by = "CellType") %>%
  mutate(
    H1_variant = factor(H1_variant, levels = c("H1-1", "H1-2", "H1-3", "H1-4")),
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

p <- ggplot(
  summary_df,
  aes(
    x = PlotCell,
    y = mean_expr,
    color = Dataset,
    fill = Dataset,
    group = Dataset
  )
) +
  geom_smooth(
    data = df_plot,
    aes(
      x = PlotCell,
      y = expression,
      color = Dataset,
      fill = Dataset,
      group = Dataset
    ),
    method = "loess",
    se = TRUE,
    span = 0.7,
    linewidth = 1.4,
    alpha = 0.15
  ) +
  geom_point(size = 2) +
  facet_grid(
    rows = vars(H1_variant),
    cols = vars(Lineage),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_x_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_color_manual(values = colors_datasets) +
  scale_fill_manual(values = colors_datasets) +
  labs(x = NULL, y = "Normalized Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black",
      size = 16
    ),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )

pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_2/Line_plots_GEX_paper_supplementary.pdf",
    width = 14, height = 10)
print(p)
dev.off()
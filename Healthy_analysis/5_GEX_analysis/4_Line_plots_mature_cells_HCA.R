library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

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
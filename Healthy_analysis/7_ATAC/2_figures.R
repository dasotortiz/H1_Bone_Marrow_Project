# ==============================================================================
# 03. FIGURES - H1-10 coverage track and TF~H1-10 correlation panel
# ==============================================================================
# Input : xpand_atac_preprocessed_path (from 01_preprocessing.R),
#         common_tfs_path + xpand_pseudobulk_path (from 02_analysis.R)
# Output: coverage_H1-10_lymphocytes.pdf,
#         tf_h1_10_correlations_xpand_CLP_ProB_PreB.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
})

# ==============================================================================
# 0. CONFIG
# ==============================================================================
project_dir <- "/home/x_vazquem/workspace/h1_project/panel_3_C-E"
outdir      <- file.path(project_dir, "results")

xpand_atac_preprocessed_path <- file.path(outdir, "preprocessed_seurats/Xpand_atac_preprocessed.rds")
common_tfs_path              <- file.path(outdir, "common_tf_h1_10_correlations_min2_CLP_ProB_PreB.csv")
xpand_pseudobulk_path        <- file.path(outdir, "tf_h1_10_pseudobulk_xpand_CLP_ProB_PreB.csv")

# Must match tf_min_datasets in 02_analysis.R (only used in a message here)
tf_min_datasets <- 2

# Cell types shown in the coverage track (full lymphoid trajectory)
lymphocyte_lineage <- c("HSC", "MPP", "CLP", "ProB", "PreB")

# Cell types used for the TF~H1-10 expression correlations
correlation_celltypes <- c("CLP", "ProB", "PreB")

colors_celltypes <- c(
  HSC         = "#00441B", MPP         = "#00AF99",
  CLP         = "#98D9E9", ProB        = "#0081C9",
  PreB        = "#001588", MEP         = "#F6313E",
  ERP         = "#8F1336", MKP         = "#46A040",
  GMP         = "#FFC179", `Eo/B/Mast` = "#FFA300",
  MDP         = "#AFAFAF"
)

# ==============================================================================
# 1. COVERAGE PLOT - H1-10 (H1FX), full lymphoid trajectory
# ==============================================================================
xpand_atac <- readRDS(xpand_atac_preprocessed_path)
DefaultAssay(xpand_atac) <- "ATAC"

seurat_atac_sub <- subset(xpand_atac, subset = annotation_2 %in% lymphocyte_lineage)
seurat_atac_sub$annotation_2 <- factor(seurat_atac_sub$annotation_2, levels = lymphocyte_lineage)
Idents(seurat_atac_sub) <- "annotation_2"

p <- CoveragePlot(
  object = seurat_atac_sub,
  region = "H1FX",
  group.by = "annotation_2",
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000
) & scale_fill_manual(values = colors_celltypes)

ggsave(
  filename = file.path(outdir, "coverage_H1-10_lymphocytes.pdf"),
  plot = p,
  width = 8, height = 6, device = "pdf"
)

rm(xpand_atac, seurat_atac_sub); gc()

# ==============================================================================
# 2. CORRELATION PLOTS - PSEUDOBULKED H1-10 vs EACH FINAL TF (XPAND)
# ==============================================================================
message("--> Generating correlation plots for H1-10 vs final candidate TFs (Xpand, CLP/ProB/PreB)...")

final_results_min2  <- read_csv(common_tfs_path, show_col_types = FALSE)
xpand_tf_pseudobulk <- read_csv(xpand_pseudobulk_path, show_col_types = FALSE)

final_tf_list <- unique(final_results_min2$TF)

if (length(final_tf_list) == 0) {
  message("No TFs passed the |Rho| / FDR filter in >= ", tf_min_datasets,
          " datasets - skipping the TF correlation panel.")
} else {
  plot_data_final_tfs <- xpand_tf_pseudobulk %>%
    filter(TF %in% final_tf_list) %>%
    mutate(Cell_Type = factor(Cell_Type, levels = correlation_celltypes))

  stats_final_tfs <- plot_data_final_tfs %>%
    group_by(TF) %>%
    summarise(
      Rho = cor(pb_tf_expr, pb_h1_expr, method = "spearman", use = "complete.obs"),
      .groups = "drop"
    ) %>%
    mutate(stat_label = paste0("Rho = ", round(Rho, 2)))

  p_tf_h1_10_corr <- ggplot(plot_data_final_tfs, aes(x = pb_h1_expr, y = pb_tf_expr)) +
    geom_point(aes(color = Cell_Type), size = 3) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    geom_text(
      data = stats_final_tfs,
      aes(x = -Inf, y = Inf, label = stat_label),
      hjust = -0.1, vjust = 1.3, size = 5, inherit.aes = FALSE
    ) +
    facet_wrap(~ TF, scales = "free") +
    scale_color_manual(values = colors_celltypes, breaks = correlation_celltypes) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      strip.text = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold")
    ) +
    labs(x = "H1-10 expression", y = "TF expression", color = "")

  ggsave(
    filename = file.path(outdir, "tf_h1_10_correlations_xpand_CLP_ProB_PreB.pdf"),
    plot = p_tf_h1_10_corr,
    width = 10, height = 5, device = "pdf"
  )
}

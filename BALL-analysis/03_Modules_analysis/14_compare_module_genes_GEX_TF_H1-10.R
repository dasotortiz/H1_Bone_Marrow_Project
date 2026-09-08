# HEATMAP: H1-10 TFs, Healthy vs Hyperdiploid and Healthy vs TCF3::PBX1
# Row annotation is the sign of delta r (increased / decreased in disease).
# Input: ./output/delta_r_correlations.csv and the pseudobulk objects

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

delta_df     <- read.csv("./output/delta_r_correlations.csv")
expr_disease <- readRDS("../GEX-Plotting/output/pseudobulk_Ia.rds")
expr_healthy <- readRDS("../GEX-Plotting/output/pseudobulk_hc.rds")

lineage_cells <- c("HSC", "MPP", "CLP", "ProB", "PreB")
subtypes      <- c("Hyperdiploid", "TCF3::PBX1")

# JASPAR TF list, same file as the scatter script
jaspar_path <- "../data/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
lines      <- readLines(jaspar_path)
jaspar_tfs <- unique(gsub("^>\\S+\\s+", "", grep("^>", lines, value = TRUE)))

# H1-10 module TFs only, one row per Gene x Subtype
tf_anno <- delta_df %>%
  filter(H1_gene == "H1-10", Group %in% subtypes, Gene %in% jaspar_tfs) %>%
  mutate(
    change = case_when(
      delta_r > 0 ~ "Increased in disease",
      delta_r < 0 ~ "Decreased in disease",
      TRUE        ~ "No change"
    )
  ) %>%
  select(Gene, Subtype = Group, delta_r, change)

tf_genes <- unique(tf_anno$Gene)
message("TFs in H1-10 module: ", length(tf_genes))

# Healthy: average across donors, then across datasets
avail_hc <- intersect(tf_genes, rownames(expr_healthy))

expr_healthy_avg <- FetchData(expr_healthy,
                              vars  = c(avail_hc, "Celltype", "Dataset"),
                              layer = "data") %>%
  filter(Celltype %in% lineage_cells) %>%
  pivot_longer(cols = all_of(avail_hc), names_to = "Gene", values_to = "Expression") %>%
  group_by(Dataset, Celltype, Gene) %>%
  summarise(Mean = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  group_by(Gene, Celltype) %>%
  summarise(Healthy_mean = mean(Mean, na.rm = TRUE), .groups = "drop")

# Disease: average per subtype, per celltype
avail_ia <- intersect(tf_genes, rownames(expr_disease))

expr_disease_avg <- FetchData(expr_disease,
                              vars  = c(avail_ia, "Celltype", "Subtype"),
                              layer = "data") %>%
  filter(Celltype %in% lineage_cells, Subtype %in% subtypes) %>%
  pivot_longer(cols = all_of(avail_ia), names_to = "Gene", values_to = "Expression") %>%
  group_by(Subtype, Celltype, Gene) %>%
  summarise(Mean = mean(Expression, na.rm = TRUE), .groups = "drop")

ct_colors <- c(
  HSC  = "#00441B",
  MPP  = "#00AF99",
  CLP  = "#98D9E9",
  ProB = "#0081C9",
  PreB = "#001588"
)

change_colors <- c(
  "Increased in disease" = "#CC2222",
  "Decreased in disease" = "#1F78B4",
  "No change"            = "grey80"
)

# One gene x celltype matrix for healthy and one per subtype
mat_h <- expr_healthy_avg %>%
  pivot_wider(names_from = Celltype, values_from = Healthy_mean) %>%
  column_to_rownames("Gene") %>%
  as.matrix()
mat_h <- mat_h[, intersect(lineage_cells, colnames(mat_h)), drop = FALSE]
mat_h <- mat_h[complete.cases(mat_h), , drop = FALSE]

make_disease_mat <- function(subtype) {
  m <- expr_disease_avg %>%
    filter(Subtype == subtype) %>%
    select(Gene, Celltype, Mean) %>%
    pivot_wider(names_from = Celltype, values_from = Mean) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  m[, intersect(lineage_cells, colnames(m)), drop = FALSE]
}

mat_hyp <- make_disease_mat("Hyperdiploid")
mat_tcf <- make_disease_mat("TCF3::PBX1")

# Keep genes present and complete in all three
common_genes <- Reduce(intersect, list(rownames(mat_h), rownames(mat_hyp), rownames(mat_tcf)))
common_genes <- common_genes[
  complete.cases(mat_h[common_genes, ]) &
  complete.cases(mat_hyp[common_genes, ]) &
  complete.cases(mat_tcf[common_genes, ])
]

mat_h   <- mat_h[common_genes,   , drop = FALSE]
mat_hyp <- mat_hyp[common_genes, , drop = FALSE]
mat_tcf <- mat_tcf[common_genes, , drop = FALSE]
message("TFs retained after alignment: ", length(common_genes))

# Z-score per gene across all three matrices at once so they share one scale
mat_z_all <- t(scale(t(cbind(mat_h, mat_hyp, mat_tcf))))
idx_h   <- seq_len(ncol(mat_h))
idx_hyp <- ncol(mat_h) + seq_len(ncol(mat_hyp))
idx_tcf <- ncol(mat_h) + ncol(mat_hyp) + seq_len(ncol(mat_tcf))

mat_h_z   <- mat_z_all[, idx_h,   drop = FALSE]
mat_hyp_z <- mat_z_all[, idx_hyp, drop = FALSE]
mat_tcf_z <- mat_z_all[, idx_tcf, drop = FALSE]

colnames(mat_h_z)   <- lineage_cells[seq_len(ncol(mat_h_z))]
colnames(mat_hyp_z) <- lineage_cells[seq_len(ncol(mat_hyp_z))]
colnames(mat_tcf_z) <- lineage_cells[seq_len(ncol(mat_tcf_z))]

lim     <- quantile(abs(mat_z_all), 0.99, na.rm = TRUE)
col_fun <- colorRamp2(c(-lim, 0, lim), c("royalblue", "white", "firebrick"))

make_top_anno <- function(mat) {
  HeatmapAnnotation(
    Celltype = colnames(mat),
    col      = list(Celltype = ct_colors),
    annotation_name_gp   = gpar(fontsize = 8),
    annotation_name_side = "left",
    simple_anno_size     = unit(3.5, "mm")
  )
}

make_right_anno <- function(subtype) {
  anno_sub   <- tf_anno %>% filter(Gene %in% common_genes, Subtype == subtype)
  change_vec <- anno_sub$change[match(common_genes, anno_sub$Gene)]
  change_vec[is.na(change_vec)] <- "No change"
  rowAnnotation(
    Change = change_vec,
    col    = list(Change = change_colors),
    annotation_name_gp   = gpar(fontsize = 7),
    annotation_name_side = "top",
    simple_anno_size     = unit(4, "mm")
  )
}

# The healthy panel is always the one that gets clustered, the disease panels
# follow its row order. Both appear in two different figures, so they are wrapped.
healthy_ht <- function(row_names_side = "right") {
  Heatmap(
    mat_h_z,
    name             = "Z-score",
    col              = col_fun,
    top_annotation   = make_top_anno(mat_h_z),
    column_title     = "Healthy",
    column_title_gp  = gpar(fontsize = 10, fontface = "bold"),
    cluster_columns  = FALSE,
    column_names_gp  = gpar(fontsize = 7),
    column_names_rot = 45,
    cluster_rows       = TRUE,
    cluster_row_slices = FALSE,
    show_row_names   = TRUE,
    row_names_side   = row_names_side,
    row_names_gp     = gpar(fontsize = 7),
    rect_gp          = gpar(col = "white", lwd = 0.4),
    border           = TRUE,
    use_raster       = TRUE
  )
}

disease_ht <- function(mat_dis_z, subtype, show_row_names = FALSE) {
  Heatmap(
    mat_dis_z,
    name             = paste0("Z-score (", subtype, ")"),
    col              = col_fun,
    top_annotation   = make_top_anno(mat_dis_z),
    column_title     = subtype,
    column_title_gp  = gpar(fontsize = 10, fontface = "bold"),
    cluster_columns  = FALSE,
    column_names_gp  = gpar(fontsize = 7),
    column_names_rot = 45,
    cluster_rows     = FALSE,
    show_row_names   = show_row_names,
    row_names_side   = "right",
    row_names_gp     = gpar(fontsize = 7),
    show_heatmap_legend = FALSE,
    right_annotation = make_right_anno(subtype),
    rect_gp          = gpar(col = "white", lwd = 0.4),
    border           = TRUE,
    use_raster       = TRUE
  )
}

# All three side by side
pdf("./output/heatmap_TFGEX_H1-10_module.pdf",
    width = 16, height = max(6, length(common_genes) * 0.22))
draw(healthy_ht() + disease_ht(mat_hyp_z, "Hyperdiploid") + disease_ht(mat_tcf_z, "TCF3::PBX1"),
     column_title        = "H1-10 module TFs | row z-score",
     column_title_gp     = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "right",
     padding             = unit(c(5, 5, 5, 5), "mm"))
dev.off()
message("Saved: heatmap_TFGEX_H1-10_module.pdf")

# Pairwise version: heatmap then lineplots, one pair of pages per subtype
avail_common_hc <- intersect(common_genes, rownames(expr_healthy))

expr_healthy_common_avg <- FetchData(expr_healthy,
                                     vars  = c(avail_common_hc, "Celltype", "Dataset"),
                                     layer = "data") %>%
  filter(Celltype %in% lineage_cells) %>%
  pivot_longer(cols = all_of(avail_common_hc), names_to = "Gene", values_to = "Expression") %>%
  group_by(Dataset, Celltype, Gene) %>%
  summarise(Mean = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Celltype = factor(Celltype, levels = lineage_cells),
         Gene     = factor(Gene,     levels = common_genes))

expr_disease_common_avg <- expr_disease_avg %>%
  filter(Gene %in% common_genes) %>%
  mutate(Celltype = factor(Celltype, levels = lineage_cells),
         Gene     = factor(Gene,     levels = common_genes))

# 3 healthy datasets plus 1 disease subtype, one facet per TF
make_lineplot_by_dataset <- function(subtype) {
  healthy_sources <- unique(expr_healthy_common_avg$Dataset)
  plot_data <- bind_rows(
    expr_healthy_common_avg %>% mutate(Source = Dataset) %>% select(Source, Gene, Celltype, Mean),
    expr_disease_common_avg %>% filter(Subtype == subtype) %>%
      mutate(Source = Subtype) %>% select(Source, Gene, Celltype, Mean)
  ) %>%
    mutate(Source = factor(Source, levels = c(healthy_sources, subtype)))

  ggplot(plot_data, aes(x = Celltype, y = Mean, color = Source, group = Source)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_wrap(~ Gene, ncol = 5, scales = "free_y") +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = paste0("H1-10 module TFs | healthy datasets vs ", subtype),
      x     = "Cell type",
      y     = "Average expression",
      color = "Source"
    ) +
    theme_bw() +
    theme(
      strip.text       = element_text(size = 8),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.title       = element_text(size = 10),
      plot.title       = element_text(size = 12, face = "bold"),
      legend.position  = "bottom",
      panel.grid.major = element_line(color = "grey90")
    )
}

# Page tall enough for both the heatmap rows and the lineplot facets
page_height <- max(max(6, length(common_genes) * 0.22),
                   max(8, ceiling(length(common_genes) / 5) * 2.5))

pdf("./output/heatmap_TFGEX_H1-10_module-separate.pdf", width = 12, height = page_height)
for (st in subtypes) {
  mat_st <- if (st == "Hyperdiploid") mat_hyp_z else mat_tcf_z
  draw(healthy_ht(row_names_side = "left") + disease_ht(mat_st, st, show_row_names = TRUE),
       column_title        = paste0("H1-10 module TFs | Healthy vs ", st, " | row z-score"),
       column_title_gp     = gpar(fontsize = 12, fontface = "bold"),
       heatmap_legend_side = "right",
       padding             = unit(c(5, 5, 5, 5), "mm"))
  print(make_lineplot_by_dataset(st))
}
dev.off()
message("Saved: heatmap_TFGEX_H1-10_module-separate.pdf (4 pages)")

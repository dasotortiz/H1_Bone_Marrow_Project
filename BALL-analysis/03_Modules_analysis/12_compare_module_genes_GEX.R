# Pseudobulk expression of the module genes picked out by the delta r quadrants:
# healthy mean vs disease, as scatters and as heatmaps.
# Input: ./output/df_scatter_H1-0.rds, ./output/df_scatter_H1-10.rds (from 11)
#        pseudobulk objects from the GEX-Plotting folder

library(Seurat)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(grid)

ht_opt$message <- FALSE

lineage_cells <- c("HSC", "MPP", "CLP", "ProB", "PreB")

df_h10  <- readRDS("./output/df_scatter_H1-0.rds")
df_h110 <- readRDS("./output/df_scatter_H1-10.rds")

pb_ia_filt <- readRDS("../GEX-Plotting/output/pseudobulk_Ia.rds")
pb_hc_filt <- readRDS("../GEX-Plotting/output/pseudobulk_hc.rds")

# Gene sets, one per figure below
# Concordant genes: moved the same way in both subtypes
quad_h10 <- df_h10 %>%
  filter(quadrant %in% c("Both increased", "Both decreased")) %>%
  distinct(Gene, quadrant) %>%
  mutate(module = "H1-0")

quad_h110 <- df_h110 %>%
  filter(quadrant %in% c("Both increased", "Both decreased")) %>%
  distinct(Gene, quadrant) %>%
  mutate(module = "H1-10")

# Increased-correlation genes, one set per subtype
dux4_inc_anno <- df_h10 %>%
  filter(quadrant %in% c("Both increased", "DUX4-R only increased")) %>%
  distinct(Gene, quadrant)

znf_inc_anno <- df_h10 %>%
  filter(quadrant %in% c("Both increased", "ZNF384-R only increased")) %>%
  distinct(Gene, quadrant)

all_genes <- unique(c(quad_h10$Gene, quad_h110$Gene,
                      dux4_inc_anno$Gene, znf_inc_anno$Gene))

# Fetch + summarise expression from a pseudobulk Seurat object
fetch_module_expr <- function(pb_obj, genes, group_col, status_label) {
  available <- intersect(genes, rownames(pb_obj))
  missing   <- setdiff(genes, rownames(pb_obj))
  if (length(missing) > 0)
    message("Not found in ", status_label, " object: ", paste(missing, collapse = ", "))

  FetchData(pb_obj, vars = c(available, "Celltype", group_col), layer = "data") %>%
    filter(Celltype %in% lineage_cells) %>%
    pivot_longer(cols = all_of(available), names_to = "Gene", values_to = "Expression") %>%
    rename(Subtype = all_of(group_col)) %>%
    group_by(Subtype, Celltype, Gene) %>%
    summarise(Mean = mean(Expression, na.rm = TRUE),
              SD   = sd(Expression,   na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Status = status_label)
}

expr_disease <- fetch_module_expr(pb_ia_filt, all_genes, "Subtype", "Disease")

# Healthy is grouped by Dataset (Xpand / Ainciburu / Li), then averaged across the three
expr_healthy_raw <- fetch_module_expr(pb_hc_filt, all_genes, "Dataset", "Healthy")
expr_healthy <- expr_healthy_raw %>%
  group_by(Gene, Celltype) %>%
  summarise(Healthy_mean = mean(Mean, na.rm = TRUE), .groups = "drop")

saveRDS(expr_disease, "./output/expr_disease_module_genes.rds")
saveRDS(expr_healthy, "./output/expr_healthy_module_genes.rds")

message("Fetched: ", n_distinct(expr_healthy$Gene), " genes across ",
        n_distinct(expr_healthy$Celltype), " cell types")

# Scatter of healthy mean vs disease, concordant quadrant genes
quadrant_colors <- c("Both increased" = "#CC2222", "Both decreased" = "#1F78B4")

build_plot_df <- function(quad_df, subtypes) {
  expr_disease %>%
    filter(Gene %in% quad_df$Gene, Subtype %in% subtypes) %>%
    left_join(expr_healthy, by = c("Gene", "Celltype")) %>%
    left_join(quad_df %>% select(Gene, quadrant), by = "Gene") %>%
    filter(!is.na(quadrant), !is.na(Healthy_mean)) %>%
    mutate(
      Celltype   = factor(Celltype, levels = lineage_cells),
      delta_expr = Mean - Healthy_mean   # for labeling outliers
    )
}

plot_df_h10  <- build_plot_df(quad_h10,  c("DUX4-R", "ZNF384-R"))
plot_df_h110 <- build_plot_df(quad_h110, c("Hyperdiploid", "TCF3::PBX1"))

make_gex_scatter <- function(plot_df, subtype, h1_gene, label_thresh = 0.5) {
  sub_df <- plot_df %>% filter(Subtype == subtype)

  # Spearman rho per cell type, printed in the corner of each facet
  cor_df <- sub_df %>%
    group_by(Celltype) %>%
    summarise(r = cor(Healthy_mean, Mean, use = "pairwise.complete.obs", method = "spearman"),
              .groups = "drop") %>%
    mutate(label = paste0("rho = ", round(r, 2)))

  ggplot(sub_df, aes(x = Healthy_mean, y = Mean, color = quadrant)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey55", linewidth = 0.5) +
    geom_point(size = 2.2, alpha = 0.8) +
    geom_text_repel(
      data         = filter(sub_df, abs(delta_expr) > label_thresh),
      aes(label    = Gene),
      size         = 2.8,
      max.overlaps = 15,
      show.legend  = FALSE
    ) +
    geom_text(
      data        = cor_df,
      aes(label   = label),
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
      inherit.aes = FALSE, size = 3, color = "grey30"
    ) +
    scale_color_manual(values = quadrant_colors) +
    facet_wrap(~Celltype, nrow = 1, scales = "free") +
    labs(
      title    = paste0(h1_gene, " module | ", subtype, " vs Healthy"),
      subtitle = paste0(
        "Red = both increased, Blue = both decreased in delta-r scatter\n",
        "Dashed = y=x (no change); labeled if |disease - healthy mean| > ", label_thresh
      ),
      x     = "Healthy mean expression (avg across datasets)",
      y     = paste0(subtype, " mean expression"),
      color = "Quadrant"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold")
    )
}

p1 <- make_gex_scatter(plot_df_h10,  "DUX4-R",       "H1-0")
p2 <- make_gex_scatter(plot_df_h10,  "ZNF384-R",     "H1-0")
p3 <- make_gex_scatter(plot_df_h110, "Hyperdiploid", "H1-10")
p4 <- make_gex_scatter(plot_df_h110, "TCF3::PBX1",   "H1-10")

# One page per module, the two subtypes stacked
pdf("./output/GEX_scatter_healthy_vs_disease_module_genes.pdf", width = 14, height = 10)
print(plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 1)))
print(plot_grid(p3, p4, ncol = 1, rel_heights = c(1, 1)))
dev.off()
message("Saved: GEX_scatter_healthy_vs_disease_module_genes.pdf (2 pages)")

# Scatter of the increased-correlation H1-0 genes, one page per subtype
increased_quadrant_colors <- c(
  "Both increased"          = "#CC2222",
  "DUX4-R only increased"   = "#E69F00",
  "ZNF384-R only increased" = "#CC79A7"
)

make_increased_gex_scatter <- function(subtype, gene_anno) {
  plot_df <- expr_disease %>%
    filter(Subtype == subtype, Gene %in% gene_anno$Gene) %>%
    rename(Disease_mean = Mean) %>%
    left_join(expr_healthy, by = c("Gene", "Celltype")) %>%
    left_join(gene_anno, by = "Gene") %>%
    filter(!is.na(Healthy_mean), !is.na(Disease_mean)) %>%
    mutate(
      Celltype  = factor(Celltype, levels = lineage_cells),
      diff_expr = Disease_mean - Healthy_mean
    )

  ggplot(plot_df, aes(x = Disease_mean, y = Healthy_mean, color = quadrant)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey55", linewidth = 0.5) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text_repel(
      data = filter(plot_df, diff_expr > 0.3),
      aes(label = Gene),
      size = 2.6, max.overlaps = 20, show.legend = FALSE
    ) +
    facet_wrap(~Celltype, nrow = 1, scales = "free") +
    scale_color_manual(values = increased_quadrant_colors, drop = FALSE) +
    labs(
      title = paste0("H1-0 increased-correlation genes | ", subtype, " vs Healthy"),
      subtitle = "Points are Gene x Celltype pseudobulk means; labels shown per celltype if disease - healthy > 0.3",
      x = paste0(subtype, " pseudobulk GEX mean"),
      y = "Healthy pseudobulk GEX mean",
      color = "Quadrant"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold")
    )
}

pdf("./output/GEX_scatter_H1-0_increased_correlation_genes.pdf", width = 14, height = 5)
print(make_increased_gex_scatter("DUX4-R",   dux4_inc_anno))
print(make_increased_gex_scatter("ZNF384-R", znf_inc_anno))
dev.off()
message("Saved: GEX_scatter_H1-0_increased_correlation_genes.pdf (2 pages)")

# Same gene sets as heatmaps, columns fixed by cell type, rows clustered
ct_colors <- c(
  HSC  = "#00441B",
  MPP  = "#00AF99",
  CLP  = "#98D9E9",
  ProB = "#0081C9",
  PreB = "#001588"
)

group_colors <- c(
  "Healthy"  = "#228B22",
  "DUX4-R"   = "#E69F00",
  "ZNF384-R" = "#CC79A7"
)

make_increased_gex_heatmap <- function(subtype, gene_anno) {
  h_long <- expr_healthy %>%
    filter(Gene %in% gene_anno$Gene) %>%
    mutate(col = paste0(Celltype, "_Healthy")) %>%
    select(Gene, col, Mean = Healthy_mean)

  d_long <- expr_disease %>%
    filter(Subtype == subtype, Gene %in% gene_anno$Gene) %>%
    mutate(col = paste0(Celltype, "_", Subtype)) %>%
    select(Gene, col, Mean)

  col_order <- paste0(
    rep(lineage_cells, each = 2), "_",
    rep(c("Healthy", subtype), times = length(lineage_cells))
  )

  mat <- bind_rows(h_long, d_long) %>%
    group_by(Gene, col) %>%
    summarise(Mean = mean(Mean, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = col, values_from = Mean) %>%
    column_to_rownames("Gene") %>%
    as.matrix()

  # Pad any missing celltype/group column so all pages share the same layout
  missing_cols <- setdiff(col_order, colnames(mat))
  if (length(missing_cols) > 0) {
    mat <- cbind(mat, matrix(NA_real_, nrow = nrow(mat), ncol = length(missing_cols),
                             dimnames = list(rownames(mat), missing_cols)))
  }

  mat <- mat[, col_order, drop = FALSE]
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]

  message("H1-0 | ", subtype, " increased-correlation heatmap genes retained: ", nrow(mat))
  if (nrow(mat) == 0) {
    stop("No genes with pseudobulk GEX values available for ", subtype, " heatmap")
  }

  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0

  row_anno <- rowAnnotation(
    Quadrant = gene_anno$quadrant[match(rownames(mat_z), gene_anno$Gene)],
    col = list(Quadrant = increased_quadrant_colors),
    annotation_name_gp   = gpar(fontsize = 8),
    annotation_name_side = "top",
    simple_anno_size     = unit(4, "mm")
  )

  cols    <- colnames(mat_z)
  ct_vec  <- sub("_.*", "", cols)
  grp_vec <- sub("^[^_]+_", "", cols)

  top_anno <- HeatmapAnnotation(
    Celltype = ct_vec,
    Group    = grp_vec,
    col = list(Celltype = ct_colors, Group = group_colors),
    annotation_name_gp   = gpar(fontsize = 8),
    annotation_name_side = "left",
    simple_anno_size     = unit(3.5, "mm"),
    gap                  = unit(0.8, "mm")
  )

  lim <- quantile(abs(mat_z), 0.99, na.rm = TRUE)  # robust cap
  if (is.na(lim) || lim == 0) lim <- 1

  Heatmap(
    mat_z,
    name = "Row z-score",
    col  = colorRamp2(c(-lim, 0, lim), c("royalblue", "white", "firebrick")),
    top_annotation   = top_anno,
    column_split     = factor(ct_vec, levels = lineage_cells),
    column_gap       = unit(4, "mm"),
    column_title     = NULL,
    cluster_columns  = FALSE,
    column_labels    = grp_vec,
    column_names_gp  = gpar(fontsize = 7),
    column_names_rot = 45,
    left_annotation  = row_anno,
    cluster_rows     = TRUE,
    show_row_names   = TRUE,
    row_names_gp     = gpar(fontsize = 7),
    row_names_side   = "right",
    rect_gp    = gpar(col = "white", lwd = 0.4),
    border     = TRUE,
    use_raster = TRUE,
    na_col     = "grey90"
  )
}

pdf("./output/GEX_heatmap_H1-0_increased_correlation_genes_healthy_vs_subtype.pdf",
    width = 11, height = 8)
for (st in c("DUX4-R", "ZNF384-R")) {
  anno <- if (st == "DUX4-R") dux4_inc_anno else znf_inc_anno
  draw(make_increased_gex_heatmap(st, anno),
       column_title = paste0("H1-0 increased-correlation genes | Healthy vs ", st, " | row z-score"),
       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       padding = unit(c(5, 5, 5, 5), "mm"))
}
dev.off()
message("Saved: GEX_heatmap_H1-0_increased_correlation_genes_healthy_vs_subtype.pdf (2 pages)")

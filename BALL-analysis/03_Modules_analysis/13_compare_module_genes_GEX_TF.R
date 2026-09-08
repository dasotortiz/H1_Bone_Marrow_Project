# HEATMAP: H1-0 TFs, healthy vs DUX4-R vs ZNF384-R across lineage cell types
# Input: ./output/df_scatter_H1-0.rds (from 11) and the pseudobulk objects

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

df_h10       <- readRDS("./output/df_scatter_H1-0.rds")
expr_disease <- readRDS("../GEX-Plotting/output/pseudobulk_Ia.rds")
expr_healthy <- readRDS("../GEX-Plotting/output/pseudobulk_hc.rds")

lineage_cells <- c("HSC", "MPP", "CLP", "ProB", "PreB")
subtypes      <- c("DUX4-R", "ZNF384-R")

tf_df    <- df_h10 %>% filter(is_TF == TRUE) %>% select(Gene, quadrant)
tf_genes <- unique(tf_df$Gene)
message("TFs in H1-0 module: ", length(tf_genes)) # 7

# Healthy: average per dataset first, then across the three datasets
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

# Disease: average per subtype
avail_ia <- intersect(tf_genes, rownames(expr_disease))

expr_disease_avg <- FetchData(expr_disease,
                              vars  = c(avail_ia, "Celltype", "Subtype"),
                              layer = "data") %>%
  filter(Celltype %in% lineage_cells, Subtype %in% subtypes) %>%
  pivot_longer(cols = all_of(avail_ia), names_to = "Gene", values_to = "Expression") %>%
  group_by(Subtype, Celltype, Gene) %>%
  summarise(Mean = mean(Expression, na.rm = TRUE), .groups = "drop")

# Wide matrix: one column per celltype x group, ordered celltype-major
h_wide <- expr_healthy_avg %>%
  mutate(col = paste0(Celltype, "_Healthy")) %>%
  select(Gene, col, Mean = Healthy_mean)

d_wide <- expr_disease_avg %>%
  mutate(col = paste0(Celltype, "_", Subtype)) %>%
  select(Gene, col, Mean)

col_order <- paste0(
  rep(lineage_cells, each = length(c("Healthy", subtypes))), "_",
  rep(c("Healthy", subtypes), times = length(lineage_cells))
)

mat <- bind_rows(h_wide, d_wide) %>%
  group_by(Gene, col) %>%
  summarise(Mean = mean(Mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = col, values_from = Mean) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

mat <- mat[, intersect(col_order, colnames(mat)), drop = FALSE]
mat <- mat[complete.cases(mat), , drop = FALSE]
message("TFs retained after NA removal: ", nrow(mat))

mat_z <- t(scale(t(mat)))  # scale() works on columns, so transpose twice

# Row annotation: quadrant. Labels have to match what make_scatter() produces in 10.
quadrant_colors <- c(
  "Both increased"          = "#CC2222",
  "Both decreased"          = "#1F78B4",
  "DUX4-R only increased"   = "#B2DF8A",
  "ZNF384-R only increased" = "#e92584",
  "Near zero"               = "grey80"
)

row_quad  <- tf_df$quadrant[match(rownames(mat_z), tf_df$Gene)]
row_split <- factor(row_quad, levels = names(quadrant_colors))

row_anno <- rowAnnotation(
  Quadrant = row_quad,
  col = list(Quadrant = quadrant_colors),
  annotation_name_gp   = gpar(fontsize = 8),
  annotation_name_side = "top",
  simple_anno_size     = unit(4, "mm")
)

# Column annotation: cell type + group, split by cell type so
# Healthy | DUX4-R | ZNF384-R stay together
ct_colors <- c(
  HSC  = "#00441B",
  MPP  = "#00AF99",
  CLP  = "#98D9E9",
  ProB = "#0081C9",
  PreB = "#001588"
)

group_colors <- c(
  "Healthy"  = "#228B22",
  "DUX4-R"   = "#B2DF8A",
  "ZNF384-R" = "#e92584"
)

cols    <- colnames(mat_z)
ct_vec  <- sub("_.*", "", cols)        # HSC, MPP, ...
grp_vec <- sub("^[^_]+_", "", cols)    # Healthy, DUX4-R, ZNF384-R

top_anno <- HeatmapAnnotation(
  Celltype = ct_vec,
  Group    = grp_vec,
  col = list(Celltype = ct_colors, Group = group_colors),
  annotation_name_gp   = gpar(fontsize = 8),
  annotation_name_side = "left",
  simple_anno_size     = unit(3.5, "mm"),
  gap                  = unit(0.8, "mm")
)

lim     <- quantile(abs(mat_z), 0.99, na.rm = TRUE)  # robust cap
col_fun <- colorRamp2(c(-lim, 0, lim), c("royalblue", "white", "firebrick"))

ht <- Heatmap(
  mat_z,
  name = "Z-score",
  col  = col_fun,
  top_annotation   = top_anno,
  column_split     = factor(ct_vec, levels = lineage_cells),
  column_gap       = unit(4, "mm"),
  column_title     = NULL,
  cluster_columns  = FALSE,
  column_labels    = grp_vec,
  column_names_gp  = gpar(fontsize = 7),
  column_names_rot = 45,
  left_annotation    = row_anno,
  row_split          = row_split,
  cluster_rows       = TRUE,     # cluster within each quadrant block
  cluster_row_slices = FALSE,    # keep quadrant order fixed
  row_gap            = unit(2, "mm"),
  row_title_gp       = gpar(fontsize = 8, fontface = "bold"),
  row_title_rot      = 0,
  show_row_names     = TRUE,
  row_names_gp       = gpar(fontsize = 7),
  rect_gp    = gpar(col = "white", lwd = 0.4),
  border     = TRUE,
  use_raster = TRUE
)

pdf("./output/heatmap_TFGEX_H1-0_module.pdf",
    width = 12, height = max(6, nrow(mat_z) * 0.18))
draw(ht,
     column_title        = "H1-0 module TFs | row z-score expression",
     column_title_gp     = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "right",
     padding             = unit(c(5, 5, 5, 5), "mm"))
dev.off()
message("Saved: heatmap_TFGEX_H1-0_module.pdf")

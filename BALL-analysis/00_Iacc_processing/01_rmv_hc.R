# Find the healthy (non-malignant) clusters by marker expression and remove them.
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gtools)
library(tidyr)
set.seed(1234)

colors_celltypes <- c(
  # Stem
  'HSC'                  = "#00441B",
  'MPP'                  = "#00AF99",

  # B lineage (blue)
  'CLP'                  = "#98D9E9",
  'ProB'                 = "#0081C9",
  'PreB'                 = "#001588",
  'Follicular B cell'    = "#3B5BDB",
  'pre-PC'               = "#5E72E4",
  'Plasma Cell'          = "#7B8CFF",

  # T / NK lineage (purple)
  'pre-T'                = "#D4B9DA",
  'Naive T-cell'         = "#C994C7",
  'CD8 T-cell'           = "#9E4FB5",
  'NK cells'             = "#6A1B9A",

  # Erythroid (red)
  'MEP'                  = "#F6313E",
  'Early-Erythroblast'   = "#E95C68",
  'ERP'                  = "#8F1336",
  'Erythroblast'         = "#5C0B24",

  # Megakaryocyte (green)
  'MKP'                  = "#46A040",
  'Platelet'             = "#7BC96F",

  # Granulocyte (orange)
  'GMP'                  = "#FFF3B0",
  'Granulocytic-UNK'     = "#FFE066",
  'Immature-Neutrophil'  = "#FFB84D",
  'Neutrophil'           = "#F28E2B",
  'Eosinophil'           = "#FFD60A",
  'Eo/B/Mast'            = "#FFA300",

  # Monocyte / Dendritic (gray)
  'MDP'                  = "#AFAFAF",
  'Pre-Dendritic'        = "#C9C9C9",
  'Dendritic Cell'       = "#8C8C8C",
  'Monocyte'             = "#666666",

  # Other
  'Stromal'              = "#8B6BB8"
)

annotation_colors <- c(
  "B-Cells" = colors_celltypes[["Follicular B cell"]],
  "Erythroid" = colors_celltypes[["Erythroblast"]],
  "Monocytes" = colors_celltypes[["Monocyte"]],
  "Plasma Cells" = colors_celltypes[["Plasma Cell"]],
  "T+NK" = colors_celltypes[["NK cells"]]
)

subtype_colors <- c(
  "BCR::ABL1"        = "#8DD3C7",  # light turquoise
  "BCR::ABL1-like"   = "#56B4E9",  # sky blue
  "DUX4-R"           = "#0072B2",  # deep blue (focus)
  "ETV6::RUNX1-like" = "#CC79A7",  # muted magenta
  "Hyperdiploid"     = "#009E73",  # teal green (focus)
  "iAMP21"           = "#E69F00",  # amber
  "KMT2A-R"          = "#E31A1C",  # red
  "Low-hypodiploid"  = "#BDBDBD",  # gray
  "MEF2D-R"          = "#c6a3ec",  # purple
  "Near-haploid"     = "#4B5D73",  # slate
  "Other"            = "#8C564B",  # brown
  "PAX5alt"          = "#F4A3A8",  # soft pink
  "TCF3::PBX1"       = "#D55E00",  # burnt orange (focus)
  "ZNF384-R"         = "#7B3294",  # dark purple (focus)
  "ETV6::RUNX1"      = "#A65628",  # brown-orange

  "Healthy"          = "black",
  "Ainciburu"        = "black",
  "Li"               = "black",
  "Xpand"            = "black"
)

# I plot the same subtype composition barplot at every clustering step,
# so it lives here instead of being pasted four times.
composition_bar <- function(obj, group, title = "Composition of Clusters by Subtype") {
  obj[[]] %>%
    dplyr::count(.data[[group]], Subtype) %>%
    dplyr::group_by(.data[[group]]) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    ggplot(aes(x = .data[[group]], y = prop, fill = Subtype)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = subtype_colors, na.value = "grey80") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    labs(x = "Cluster", y = "Proportion", fill = "Subtype", title = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

combined <- readRDS("./output/combinedobj_raw_filtered_processed.rds")

umap <- DimPlot(combined, reduction = "umap.unintegrated", group.by = "orig.ident", label = TRUE, label.box = TRUE, repel = TRUE)
umapp <- DimPlot(combined, reduction = "umap.unintegrated", group.by = "Subtype") &
  scale_color_manual(values = subtype_colors, na.value = "grey80")
ggsave("./output/UMAP.png", umap | umapp , bg = "white", w = 17, h = 7)

combined <- FindClusters(combined, resolution = 0.1)
# Plot
cp11 <- DimPlot(combined, reduction = "umap.unintegrated",
              group.by = "Subtype") & NoAxes() &
  scale_color_manual(values = subtype_colors, na.value = "grey80")
cp1 <- DimPlot(combined, reduction = "umap.unintegrated",
              group.by = "RNA_snn_res.0.1", label = TRUE, label.box = TRUE, repel = TRUE) & NoAxes()
gp1 <- composition_bar(combined, "RNA_snn_res.0.1")

ggsave("./output/unintegrated_res0.1.png", cp11 | cp1 | gp1, bg = "white", w = 22, h = 8)

# Subcluster cluster 11
combined <- FindSubCluster(combined, cluster = "11", graph.name = "RNA_snn",
                            resolution = 0.1, subcluster.name = "res0.1_subc")

combined$res0.1_subc <- factor(combined$res0.1_subc, levels = mixedsort(unique(combined$res0.1_subc)))

cp2 <- DimPlot(combined, reduction = "umap.unintegrated",
              group.by = "res0.1_subc", label = TRUE, label.box = TRUE, repel = TRUE) & NoAxes()
gp2 <- composition_bar(combined, "res0.1_subc")

ggsave("./output/unintegrated_res0.1_subc.png", cp2 | gp2, bg = "white", w = 16, h = 8)

# combine 11_0, 11_1 --> 11: this will be 11 because it's all one subtype of BALL
# combine 11_2, 11_3, 11_4 and 11_5 --> will call this 11-H (healthy)
combined$res0.1_subc_final <- as.character(combined$res0.1_subc)
combined$res0.1_subc_final[combined$res0.1_subc %in% c("11_0", "11_1")] <- "11"
combined$res0.1_subc_final[combined$res0.1_subc %in% c("11_2", "11_3", "11_4", "11_5")] <- "11_H"
# Order
combined$res0.1_subc_final <- factor(combined$res0.1_subc_final, levels = mixedsort(unique(combined$res0.1_subc_final)))  # smart natural ordering
# Plot again
cp22 <- DimPlot(combined, reduction = "umap.unintegrated",
              group.by = "res0.1_subc_final", label = TRUE, label.box = TRUE, repel = TRUE) & NoAxes()
gp22 <- composition_bar(combined, "res0.1_subc_final")

ggsave("./output/unintegrated_res0.1_subc_relabeled.png", cp22 | gp22, bg = "white", w = 16, h = 8)

combined <- JoinLayers(combined)
Idents(combined) <- 'res0.1_subc_final'

# VlnPlots of markers, one set per lineage. Same plot every time so I wrote it once.
marker_violin <- function(obj, genes, set_name) {
  expr_data <- FetchData(obj, vars = genes)
  expr_data$cluster <- obj$res0.1_subc_final
  expr_data %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
    ggplot(aes(x = cluster, y = expression, fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE) + # scale = "width" means more comparable across groups
    facet_wrap(~gene, scales = "free_y", nrow = 1) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "Cluster", y = "Expression",
         title = paste(set_name, "Marker Expression per Cluster"))
}

marker_sets <- list(
  TNK       = c("CD3D", "CD3E", "CD3G", "GNLY"),
  BCell     = c("CD38", "IGLL1", "MME", "DNTT", "CD19", "MS4A1"),
  PC        = c("MZB1", "SDC1"),
  erythroid = c("HBG2", "HBM", "SLC4A1", "ACKR1", "CD34", "GATA1"),
  monocytes = c("CD14", "LYZ", "FCN1", "S100A8", "S100A9", "VCAN", "CSF1R")
)

for (set_name in names(marker_sets)) {
  genes <- marker_sets[[set_name]]
  pv <- marker_violin(combined, genes, set_name)
  pf <- FeaturePlot(combined, features = genes, order = TRUE) & scale_color_viridis_c()
  ggsave(paste0("./output/unintegrated_Vln_", set_name, ".png"), pv, bg = "white",
         w = max(10, 3.5 * length(genes)), h = 5)
  ggsave(paste0("./output/unintegrated_Feature_", set_name, ".png"), pf, bg = "white",
         w = 12, h = 4 * ceiling(length(genes) / 3))
}

# From those plots: cluster 7 is T/NK, 18 is B cells, 27 is plasma cells,
# 9 is erythroid and 11_H is monocytes :D

# Final Cluster annotations
healthy_map <- c("7"    = "T+NK",
                 "9"    = "Erythroid",
                 "18"   = "B-Cells",
                 "27"   = "Plasma Cells",
                 "11_H" = "Monocytes")

cl <- as.character(combined$res0.1_subc_final)
combined$annotation <- ifelse(cl %in% names(healthy_map), healthy_map[cl], cl)
combined$annotation <- factor(combined$annotation, levels = mixedsort(unique(combined$annotation)))  # smart natural ordering

# Keep my colours for the healthy classes, default hues for the rest
annotation_levels <- levels(combined$annotation)
annotation_plot_colors <- setNames(scales::hue_pal()(length(annotation_levels)), annotation_levels)
annotation_plot_colors[names(annotation_colors)] <- annotation_colors

# Median UMAP position per annotation, used to place the labels
annotation_umap <- as.data.frame(Embeddings(combined, "umap.unintegrated"))
annotation_umap$annotation <- combined$annotation
umap_cols <- colnames(annotation_umap)[1:2]
annotation_labels <- annotation_umap %>%
  dplyr::group_by(annotation) %>%
  dplyr::summarise(
    UMAP_1 = median(.data[[umap_cols[1]]]),
    UMAP_2 = median(.data[[umap_cols[2]]]),
    .groups = "drop"
  )

pf <- DimPlot(combined, reduction = "umap.unintegrated", group.by = "annotation", label = FALSE) +
      ggtitle("Final Healthy Annotations") +
      scale_color_manual(values = annotation_plot_colors) +
      ggrepel::geom_label_repel(
        data = annotation_labels,
        aes(x = UMAP_1, y = UMAP_2, label = annotation, color = annotation),
        fill = "white",
        fontface = "bold",
        label.size = 0.25,
        label.r = grid::unit(0.15, "lines"),
        size = 4,
        show.legend = FALSE
      ) +
      NoAxes()
ggsave("./output/final_annotation_healthy.png", pf, w = 10, h = 7)

saveRDS(combined, "./output/combinedobj_raw_filtered_processed_clustered.rds")

# Exclude healthy cells
combined_subset <- subset(combined,
                          subset = !annotation %in% c("B-Cells", "Erythroid", "Monocytes", "Plasma Cells", "T+NK"))
saveRDS(combined_subset, "./output/combinedobj_processed_noHC.rds")

# Plots
p1 <- DimPlot(combined, reduction = "umap.unintegrated",
              group.by = "Subtype") & NoAxes() &
  scale_color_manual(values = subtype_colors, na.value = "grey80")

gp1 <- composition_bar(combined, "annotation", "Subtype Composition of Clusters")

ggsave("./output/unintegrated_HC.png", p1 | pf | gp1, bg = "white", w = 26, h = 7)

m <- c("LYZ", "HBM", "MS4A1", "SDC1", "CD3E")
pcm <- FeaturePlot(combined, features = m, order = TRUE, ncol = 5) & scale_color_viridis_c() & NoAxes()
ggsave("./output/unintegrated_HC_markers.png", pcm, bg = "white", w = 26, h = 5)

p2 <- DimPlot(combined_subset, reduction = "umap.unintegrated",
              group.by = "Subtype") & NoAxes() & scale_color_manual(values = subtype_colors, na.value = "grey80")

ggsave("./output/unintegrated_post_HC_mv.png", p2, bg = "white", w = 7, h = 5)


# Plots for manuscript: reuses combined, combined_subset and the palettes above
p_all_clusters <- DimPlot(
  combined,
  reduction = "umap.unintegrated",
  group.by = "res0.1_subc_final",
  label = TRUE,
  label.box = TRUE,
  repel = TRUE
) +
  ggtitle("Unsupervised clustering of all cells") +
  NoAxes()

manuscript_marker_groups <- list(
  "B cells" = marker_sets$BCell,
  "Erythroid" = marker_sets$erythroid,
  "Monocytes" = marker_sets$monocytes,
  "Plasma cells" = marker_sets$PC,
  "T/NK" = marker_sets$TNK
)

p_marker_dotplot <- DotPlot(
  combined,
  features = manuscript_marker_groups,
  group.by = "annotation",
  dot.scale = 5
) +
  ggtitle("Canonical marker expression used for healthy-cell calls") +
  scale_color_gradientn(
    colors = c("#F0F0F0", "#6BAED6", "#08306B"),
    name = "Average expression"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing.x = grid::unit(0.15, "lines"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Marker genes", y = "Cluster / annotation")

p_post_removal <- DimPlot(combined_subset, reduction = "umap.unintegrated", group.by = "Subtype") &
  ggtitle("B-ALL cells after healthy-cell removal") &
  scale_color_manual(values = subtype_colors, na.value = "grey80") &
  NoAxes()

removed_counts <- combined[[]] %>%
  dplyr::filter(annotation %in% names(annotation_colors)) %>%
  dplyr::count(annotation) %>%
  dplyr::mutate(annotation = factor(annotation, levels = names(annotation_colors)))

p_removed_counts <- ggplot(removed_counts, aes(x = annotation, y = n, fill = annotation)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = annotation_colors) +
  scale_y_continuous(labels = scales::comma) +
  theme_bw() +
  labs(
    title = "Removed healthy-cell classes",
    x = NULL,
    y = "Cells removed"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

manuscript_design <- "
AABB
CCDE
"
manuscript_hc_panel <- p_all_clusters + p_marker_dotplot + pf + p_post_removal + p_removed_counts +
  plot_layout(design = manuscript_design, heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("./output/manuscript_healthy_cell_removal_panel.pdf", manuscript_hc_panel, bg = "white", w = 20, h = 12)

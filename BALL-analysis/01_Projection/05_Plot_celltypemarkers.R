# Marker dotplot across healthy annotations plus transferred labels, and the
# manuscript projection panel.
set.seed(1234)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Load data
ref <- readRDS("./data/integration_post_gene_fix.rds")
query <- readRDS("./output/query_predictions_final.rds")
queryOG <- readRDS("../Iacobucci/output/combinedobj_raw_filtered_processed_clustered.rds") # Iacabucci dataset with healthy

colors_celltypes <- c(
  "HSC"                 = "#00441B",
  "MPP"                 = "#00AF99",
  "CLP"                 = "#98D9E9",
  "ProB"                = "#0081C9",
  "PreB"                = "#001588",
  "Follicular B cell"   = "#3B5BDB",
  "pre-PC"              = "#5E72E4",
  "Plasma Cell"         = "#7B8CFF",
  "pre-T"               = "#D4B9DA",
  "Naive T-cell"        = "#C994C7",
  "CD8 T-cell"          = "#9E4FB5",
  "NK cells"            = "#6A1B9A",
  "MEP"                 = "#F6313E",
  "Early-Erythroblast"  = "#E95C68",
  "ERP"                 = "#8F1336",
  "Erythroblast"        = "#5C0B24",
  "MKP"                 = "#46A040",
  "Platelet"            = "#7BC96F",
  "GMP"                 = "#FFF3B0",
  "Granulocytic-UNK"    = "#FFE066",
  "Immature-Neutrophil" = "#FFB84D",
  "Neutrophil"          = "#F28E2B",
  "Eosinophil"          = "#FFD60A",
  "Eo/B/Mast"           = "#FFA300",
  "MDP"                 = "#AFAFAF",
  "Pre-Dendritic"       = "#C9C9C9",
  "Dendritic Cell"      = "#8C8C8C",
  "Monocyte"            = "#666666",
  "Stromal"             = "#8B6BB8"
)

subtype_colors <- c(
  "BCR::ABL1"        = "#8DD3C7",
  "BCR::ABL1-like"   = "#56B4E9",
  "DUX4-R"           = "#0072B2",
  "ETV6::RUNX1-like" = "#CC79A7",
  "Hyperdiploid"     = "#009E73",
  "iAMP21"           = "#E69F00",
  "KMT2A-R"          = "#E31A1C",
  "Low-hypodiploid"  = "#BDBDBD",
  "MEF2D-R"          = "#c6a3ec",
  "Near-haploid"     = "#4B5D73",
  "Other"            = "#8C564B",
  "PAX5alt"          = "#F4A3A8",
  "TCF3::PBX1"       = "#D55E00",
  "ZNF384-R"         = "#7B3294",
  "ETV6::RUNX1"      = "#A65628",
  "Healthy"          = "black",
  "Ainciburu"        = "black",
  "Li"               = "black",
  "Xpand"            = "black"
)

celltypeorder <- c("HSC","MPP","ERP","MEP","MKP","GMP","MDP","Eo/B/Mast","CLP","ProB","PreB")
celltype_colors <- colors_celltypes[celltypeorder]
ref$annotation_2 <- factor(ref$annotation_2, levels = celltypeorder)
query$predicted.celltype <- factor(query$predicted.celltype, levels = celltypeorder)

# queryOG formatting
# Label HC categories + MALIGNANT
cell_types <- c("B-Cells", "Erythroid", "Monocytes", "Plasma Cells", "T+NK")
queryOG$annotation2 <- ifelse(as.character(queryOG$annotation) %in% cell_types,
                              as.character(queryOG$annotation),
                              "MALIGNANT")
table(queryOG$annotation2)
# For malignant cells, copy the predicted celltypes from the projections
# annotation3 has to stay a character vector while I fill it in
queryOG$annotation3 <- as.character(queryOG$annotation2)
malignant_cells <- colnames(queryOG)[queryOG$annotation3 == "MALIGNANT"]
predicted_labels <- as.character(query$predicted.celltype[malignant_cells])
queryOG$annotation3[malignant_cells] <- predicted_labels
table(queryOG$annotation3)
# Make sure we can distinguish healthy cells from disease
queryOG$annotation3 <- ifelse(
  as.character(queryOG$annotation3) %in% cell_types,
  paste0(as.character(queryOG$annotation3), "(H)"),
  as.character(queryOG$annotation3)
)
table(queryOG$annotation3)
# order: (H) groups first
current_levels <- sort(unique(queryOG$annotation3))
h_groups <- current_levels[grepl("\\(H\\)$", current_levels)]
other_groups <- current_levels[!grepl("\\(H\\)$", current_levels)]
queryOG$annotation3 <- factor(queryOG$annotation3, levels = c(h_groups, other_groups))

saveRDS(queryOG, "./output/query_HCann_projectionann.rds")

# Marker genes used to annotate the healthy reference, extended with the
# healthy populations I removed in 01_rmv_hc.R
markers_list <- list(
  HSC = c("CD34", "CD164", "BEX1", "BEX2", "AVP", "CRHBP", "HLF"),
  MPP = c("CD33", "MPO", "FLT3"),
  GMP = c("AZU1", "PRTN3", "ELANE"),
  MDP = c("IRF8", "LY86", "RUNX2", "LILRB4"),
  MEP = c("GATA2", "FCER1A", "ITGA2B", "CSF2RB", "HBG2", "HBM", "SLC4A1", "ACKR1"),
  ERP = c("GATA1", "EPOR", "CA1", "CA2", "EPCAM", "KLF1", "BLVRB", "APOC1", "APOE"),
  MKP = c("PLEK", "PPBP", "PF4", "GP9"),
  `Eo/B/Mast` = c("MS4A2", "MS4A3", "TPSAB1", "TPSB2", "HDC", "CLC", "PRG2"),
  CLP = c("TRBC2", "LTB", "JCHAIN", "ADA", "BCL2"),
  ProB = c("DNTT", "MME", "PAX5", "RAG1", "RAG2"),
  PreB = c("VPREB1", "IGLL1", "TCL1A"),
  Plasma = c("MZB1","SDC1"),
  B_cells = c("CD38", "CD19", "MS4A1"),
  T_NK = c("CD3D", "CD3E", "CD3G", "GNLY"),
  Monocytes = c("CD14", "FCN1", "S100A8", "S100A9", "VCAN", "CSF1R")
)

# Removed duplicates:
# CD34 (was in HSC and MPP - kept in HSC)
# MPO (was in MPP, GMP, MDP - kept in MPP)
# FLT3 (was in MPP and CLP - kept in MPP)
# MZB1 (was in MPP and CLP - kept in MPP)
# LYZ (was in GMP and MDP - kept in GMP)
# GATA2 (was in MEP and MKP - kept in MEP)
# FCER1A (was in MEP and MKP - kept in MEP)
# ITGA2B (was in MEP and MKP - kept in MEP)
# JCHAIN (was in CLP and PreB - kept in CLP)

# Extract the markers in order
markers <- unlist(markers_list, use.names = FALSE)

# Map each marker back to its cell type, keeping the order defined above
marker_celltype_map <- data.frame(
  gene = markers,
  celltype = rep(names(markers_list), times = lengths(markers_list))
)
marker_celltype_map$celltype <- factor(marker_celltype_map$celltype, levels = names(markers_list))

# Build the DotPlot, then rebuild it by hand so the markers can be faceted by group
Idents(queryOG) <- "annotation3"
p <- DotPlot(queryOG, features = markers) +
  RotatedAxis() +
  scale_colour_viridis_c(option = "turbo")

plot_data <- p$data %>%
  left_join(marker_celltype_map, by = c("features.plot" = "gene"))
plot_data$features.plot <- factor(plot_data$features.plot, levels = markers)

p_grouped <- ggplot(plot_data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_colour_viridis_c(option = "turbo") +
  scale_size_continuous(range = c(0, 6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Markers", y = "Cell Type", color = "Average\nExpression", size = "Percent\nExpressed") +
  facet_grid(. ~ celltype, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("./output/query_marker_dotplot_extensivelist_wHCcells.png", plot = p_grouped, width = 16, height = 8, bg = "white")

### Manuscript projection figure
# This panel supports the text describing projection of malignant B-ALL cells
# onto the healthy hematopoietic atlas and their developmental label transfer.
present_subtypes <- unique(as.character(query$Subtype))
subtype_order <- c(
  intersect(names(subtype_colors), present_subtypes),
  setdiff(sort(present_subtypes), names(subtype_colors))
)
query$Subtype <- factor(query$Subtype, levels = subtype_order)

p_ref_atlas <- DimPlot(
  ref,
  reduction = "umap_integrated_rpca",
  group.by = "annotation_2",
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  cols = celltype_colors
) +
  ggtitle("Healthy hematopoietic reference atlas") +
  NoAxes() +
  NoLegend()

p_query_transfer <- DimPlot(
  query,
  reduction = "ref.umap",
  group.by = "predicted.celltype",
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  cols = celltype_colors
) +
  ggtitle("Projected B-ALL cells by transferred healthy label") +
  NoAxes() +
  NoLegend()

transfer_summary <- query@meta.data %>%
  count(orig.ident, Subtype, predicted.celltype) %>%
  group_by(orig.ident, Subtype) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

transfer_summary$orig.ident <- factor(
  transfer_summary$orig.ident,
  levels = transfer_summary %>%
    distinct(Subtype, orig.ident) %>%
    arrange(Subtype, orig.ident) %>%
    pull(orig.ident)
)

p_transfer_proportions <- ggplot(
  transfer_summary,
  aes(x = orig.ident, y = proportion, fill = predicted.celltype)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~Subtype, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(
    title = "Transferred healthy labels by sample, grouped by B-ALL subtype",
    x = "Sample",
    y = "Proportion of malignant cells",
    fill = "Transferred label"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = grid::unit(0.3, "lines"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

projection_marker_panel <- p_grouped +
  ggtitle("Canonical marker expression supports transferred healthy labels")

projection_main_panel <- (p_ref_atlas | p_query_transfer) / p_transfer_proportions / projection_marker_panel +
  plot_layout(heights = c(1, 1.1, 1.35)) +
  plot_annotation(tag_levels = "A")

ggsave("./output/manuscript_projection_main_panel.pdf", projection_main_panel, bg = "white", w = 18, h = 22)

# Project the Iacobucci B-ALL cells onto the healthy hematopoietic reference and
# transfer cell type labels.
# https://satijalab.org/seurat/articles/integration_mapping
# nohup Rscript 02_Projection_toHealthy.R > projection.log 2>&1
set.seed(1234)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# Query: Iacobucci dataset without the healthy cells (from 01_rmv_hc.R)
query <- readRDS("../Iacobucci/output/combinedobj_processed_noHC.rds")

# Healthy reference for the projection
ref <- readRDS("./data/integration_post_gene_fix.rds")

# Process reference
# return.model = TRUE matters here, MapQuery needs the stored UMAP model
ref <- FindNeighbors(ref, reduction = "integrated.rpca", dims = 1:30)
ref <- FindClusters(ref, resolution = 0.1)
ref <- RunUMAP(ref, reduction = "integrated.rpca", dims = 1:30,
               reduction.name = "umap_integrated_rpca", seed.use = 123,
               return.model = TRUE)
saveRDS(ref, "./data/integration_post_gene_fix.rds")

# confirm the model is stored
length(SeuratObject::Misc(ref[["umap_integrated_rpca"]], slot = "model"))

# Plot reference
p <- DimPlot(ref, reduction = "umap_integrated_rpca", group.by = c("annotation_2", "Dataset"),
             label = TRUE, label.box = TRUE, ncol = 2)
ggsave("./output/healthyRef_int_rpca.png", p, w = 13, h = 6)

# Transfer data:
anchors <- FindTransferAnchors(reference = ref,
                               query = query, dims = 1:30,
                               reference.reduction = "pca") # From Seurat tutorial
saveRDS(anchors, "./output/anchors.rds")
predictions <- TransferData(anchorset = anchors, refdata = ref$annotation_2, dims = 1:30)
query <- AddMetaData(query, metadata = predictions)
saveRDS(query, "./output/query_predictions.rds")

table(query$predicted.id)
    #   CLP Eo/B/Mast       ERP       GMP       HSC       MDP       MEP       MKP
    # 47768        31      2860      7041     14767     34736      1067         3
    #   MPP      PreB      ProB
    # 23767     78917     91063

# Map query
query <- MapQuery(anchorset = anchors, reference = ref, query = query,
    refdata = list(celltype = "annotation_2"), reference.reduction = "pca",
    reduction.model = "umap_integrated_rpca")

saveRDS(query, "./output/query_predictions_final.rds")

# Plotting:
ref <- readRDS("./data/integration_post_gene_fix.rds")
query <- readRDS("./output/query_predictions_final.rds")

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

celltypeorder <- c( "HSC",        # hematopoietic stem cell
                    "MPP",        # multipotent progenitor
                    "ERP",        # early erythroid progenitor
                    "MEP",        # megakaryocyte-erythroid progenitor
                    "MKP",        # megakaryocyte progenitor
                    "GMP",        # granulocyte-monocyte progenitor
                    "MDP",        # monocyte-dendritic progenitor
                    "Eo/B/Mast",  # eosinophil/B/mast mixed lineage
                    "CLP",        # common lymphoid progenitor
                    "ProB",       # pro-B cell
                    "PreB"        # pre-B cell
                                )
celltype_colors <- colors_celltypes[celltypeorder]

ref$annotation_2 <- factor(ref$annotation_2, levels = celltypeorder)
query$predicted.celltype <- factor(query$predicted.celltype, levels = celltypeorder)

p1 <- DimPlot(ref, reduction = "umap_integrated_rpca", group.by = "annotation_2", label = TRUE, label.box = TRUE, cols = celltype_colors,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")

p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.box = TRUE,
    repel = TRUE, cols = celltype_colors) + NoLegend() + ggtitle("Query transferred labels")

ggsave("./output/projected_labels_Ia.png", p1 | p2, bg = "white", w = 14, h = 7)

# Proportions per Sample grouped by Subtype
celltype_summary <- query@meta.data %>%
  count(orig.ident, predicted.celltype, Subtype) %>%  # include Subtype
  group_by(orig.ident, Subtype) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(Subtype, orig.ident)
# Preserve ordering by Subtype
celltype_summary$orig.ident <- factor(
  celltype_summary$orig.ident,
  levels = celltype_summary %>%
    distinct(Subtype, orig.ident) %>%
    arrange(Subtype, orig.ident) %>%
    pull(orig.ident)
)

p3 <- ggplot(celltype_summary, aes(x = orig.ident, y = proportion, fill = predicted.celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = celltype_colors) +
  theme_bw() +
  labs(
    y = "Proportion",
    x = "Sample",
    fill = "Predicted Celltype"
  ) +
  facet_grid(~Subtype, scales = "free_x", space = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0.3, "lines"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))  # keep the legend on one row

ggsave("./output/predicted_celltype_proportions_per_sample_GroupedbySubtype.png", p3, bg = "white", w = 12, h = 7)

# Proportions per SUBTYPE:
celltype_summary <- query@meta.data %>%
  count(Subtype, predicted.celltype) %>%
  group_by(Subtype) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(Subtype)

p4 <- ggplot(celltype_summary, aes(x = Subtype, y = proportion, fill = predicted.celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = celltype_colors) +
  theme_bw() +
  labs(
    y = "Proportion",
    x = "Subtype",
    fill = "Predicted Celltype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0.3, "lines"),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))  # consistent legend layout

ggsave("./output/predicted_celltype_proportions_per_subtype.png",
  p4,
  bg = "white",
  w = 10,
  h = 6
)

# Same projection, but onto the BoneMarrowMap reference instead of the healthy
# reference used in 02, so the two label sets can be compared.
# https://satijalab.org/seurat/articles/integration_mapping
# nohup Rscript 04_Projection_toBoneMarrowMap_SeuratPipeline.R > projection.log 2>&1
set.seed(1234)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

out_dir <- "./output/BoneMarrowMapResults"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Query
query <- readRDS("../Iacobucci/output/combinedobj_processed_noHC.rds") # Iacabucci dataset w/o healthy

# Reference for projection
ref <- readRDS("./data/BoneMarrowMap_Annotated_Dataset_expandedFeatures.rds")
ref

# Plot reference
p <- DimPlot(ref, reduction = "umap", group.by = "CellType_Broad", label = TRUE, label.box = TRUE) &
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "BoneMarrowMap_umap.png"), p, w = 13, h = 13)

# Transfer data:
anchors <- FindTransferAnchors(reference = ref,
                               query = query, dims = 1:30,
                               reference.reduction = "pca") # From Seurat tutorial
saveRDS(anchors, file.path(out_dir, "anchors.rds"))
predictions <- TransferData(anchorset = anchors, refdata = ref$CellType_Broad, dims = 1:30)
query <- AddMetaData(query, metadata = predictions)
saveRDS(query, file.path(out_dir, "query_predictions.rds"))

table(query$predicted.id)

# The reference already carries a UMAP model, so MapQuery can use it directly
length(SeuratObject::Misc(ref[["umap"]], slot = "model"))

# Map query
query <- MapQuery(anchorset = anchors, reference = ref, query = query,
    refdata = list(celltype = "CellType_Broad"), reference.reduction = "pca",
    reduction.model = "umap")

saveRDS(query, file.path(out_dir, "query_predictions_final.rds"))

# Plotting:
ref <- readRDS("./data/BoneMarrowMap_Annotated_Dataset_expandedFeatures.rds")
query <- readRDS(file.path(out_dir, "query_predictions_final.rds"))

celltype_colors <- c(
  "HSC MPP" = "#005AB5FF",
  "LMPP" = "#579D1CFF",
  "Early Lymphoid" = "#7AD151FF",
  "Pro-B" = "#4B1F6FFF",
  "Pre-B" = "#FF950EFF",
  "B" = "#FFB347FF",
  "Plasma_Cell" = "#FF7F0EFF",
  "Early GMP" = "#A63D00FF",
  "Late GMP" = "#7E0021FF",
  "Monocyte" = "#B55D00FF",
  "Pro-Monocyte" = "#D86E00FF",
  "cDC" = "#83CAFFFF",
  "pDC" = "#5FB6D9FF",
  "Cycling Progenitor" = "#808080FF",
  "Early Erythroid" = "#E6194BFF",
  "Late Erythroid" = "#C5000BFF",
  "MEP" = "#314004FF",
  "Megakaryocyte Precursor" = "#AECF00FF",
  "EoBasoMast Precursor" = "#0084D1FF",
  "Naive T" = "#58A94EFF",
  "CD4 Memory T" = "#2E7A38FF",
  "CD8 Memory T" = "#1E5E2BFF",
  "NK" = "#3C9F6EFF",
  "Stromal" = "#C4A484FF"
)

celltypeorder <- c(
  # Myeloid lineage
  "HSC MPP",                   # hematopoietic stem / multipotent progenitor
  "LMPP",                      # lymphoid-primed multipotent progenitor
  "Early GMP",                 # early granulocyte-monocyte progenitor
  "Late GMP",                  # late granulocyte-monocyte progenitor
  "Pro-Monocyte",              # pro-monocyte
  "Monocyte",                  # monocyte
  "cDC",                       # conventional dendritic cell
  "pDC",                       # plasmacytoid dendritic cell
  "EoBasoMast Precursor",      # eosinophil / basophil / mast precursor

  # Erythroid & megakaryocytic lineage
  "Early Erythroid",           # early erythroid progenitor
  "Late Erythroid",            # late erythroid progenitor
  "MEP",                       # megakaryocyte-erythroid progenitor
  "Megakaryocyte Precursor",   # megakaryocyte precursor

  # Lymphoid lineage
  "Early Lymphoid",            # early lymphoid progenitor
  "Pro-B",                     # pro-B cell
  "Pre-B",                     # pre-B cell
  "B",                         # mature B cell
  "Plasma_Cell",               # plasma cell
  "Naive T",                   # naive T cell
  "CD4 Memory T",              # CD4+ memory T cell
  "CD8 Memory T",              # CD8+ memory T cell
  "NK",                        # natural killer cell

  # Other
  "Cycling Progenitor",        # cycling progenitor
  "Stromal"                    # stromal cell
)

ref$CellType_Broad <- factor(ref$CellType_Broad, levels = celltypeorder)
query$predicted.celltype <- factor(query$predicted.celltype, levels = celltypeorder)

p1 <- DimPlot(ref, reduction = "umap", group.by = "CellType_Broad", label = TRUE, label.box = TRUE, cols = celltype_colors,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")

p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.box = TRUE,
    repel = TRUE, cols = celltype_colors) + NoLegend() + ggtitle("Query transferred labels")

ggsave(file.path(out_dir, "projected_labels.png"), p1 | p2, bg = "white", w = 14, h = 7)

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
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))  # 24 labels, so 4 rows

ggsave(file.path(out_dir, "predicted_celltype_proportions_per_sample_GroupedbySubtype.png"), p3, bg = "white", w = 12, h = 7)

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
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))  # consistent legend layout

ggsave(file.path(out_dir, "predicted_celltype_proportions_per_subtype.png"),
  p4,
  bg = "white",
  w = 10,
  h = 6
)

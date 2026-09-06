library(Seurat)
library(scCustomize)
library(ggplot2)

integrated <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/integrated_datasets_filtered_final_annotation.rds")
hca_dataset <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_merging_ready.rds')

colors_celltypes <- c(
  HSC        = "#00441B",
  MPP        = "#00AF99",
  CLP        = "#98D9E9",
  ProB       = "#0081C9",
  PreB       = "#001588",
  MEP        = "#F6313E",
  ERP        = "#8F1336",
  MKP        = "#46A040",
  GMP        = "#FFC179",
  `Eo/B/Mast`= "#FFA300",
  MDP        = "#AFAFAF"
)

markers <- c(
  #MKP
  'ITGA2B', 'HACD1', 'SELP',
  #MEP
 'GATA2','FCER1A', 'PBX1',
  # MDP
  'IRF8','IRF7','SPIB',
  # GMP
  'ELANE', 'CTSG', 'PRTN3',
  # Eo/B/Mast
 'CSF2RB', 'MS4A2', 'MS4A3', 'CLC',
  # ERP
  'GATA1', 'CA1', 'APOC1', 'HBB',
  # PreB
  'CD24', 'SMC4', 'VPREB1', 'TOP2A', 'IGHM',
  # ProB
  'RAG1', 'AKAP12', 'ARPP21', 'TOP2B',
  #CLP
  'ADA', 'JCHAIN', 'DNTT',
  #MPP
  'SPINK2', 'SMIM24',
  # HSC
  'AVP', 'CD164', 'CRHBP'
)
# Dotplot markers
library(viridis)
pal <- viridis(n = 10, option = "D")
Idents(integrated) <- "annotation_2"
pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Supplementary_1/dotplot_markers.pdf",width = 8, height = 8)
DotPlot_scCustom(seurat_object = integrated, features = markers, x_lab_rotate = TRUE,colors = pal) & coord_flip() & theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

# Dotplot HCA

cells <- c("HSC", "MPP", "CLP", "ProB", "PreB", "ERP", "Eo/B/Mast", "GMP", "MDP", "MEP", "MKP")
hca <- subset(hca_dataset, subset = annotation_2 %in% cells)
hca$annotation_2 <- factor(hca$annotation_2, levels = cells)
Idents(hca) <- "annotation_2"
pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Supplementary_1/dotplot_markers_hca.pdf",width = 8, height = 8)
DotPlot_scCustom(seurat_object = hca, features = markers, x_lab_rotate = TRUE,colors = pal) & coord_flip() & theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
DotPlot_scCustom(seurat_object = hca, features = markers, x_lab_rotate = TRUE,colors = pal) & coord_flip() & theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

#UMAP HCA
pdf('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Supplementary_1/umap_hca.pdf', width = 6, height = 6)
DimPlot(hca, reduction = "umap.integrated", group.by = 'annotation_2', cols = colors_celltypes, raster = TRUE, label = TRUE, repel = TRUE, pt.size = 1.2)
DimPlot(hca, reduction = "umap.integrated", group.by = 'annotation_2', cols = colors_celltypes, raster = TRUE, pt.size = 1.2, label = TRUE, repel = TRUE) + NoLegend()
dev.off()

# UMAP Integrated
pdf('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Supplementary_1/umap_integrated.pdf', width = 6, height = 6)
DimPlot(integrated, reduction = "umap_integrated_rpca", group.by = 'annotation_2', cols = colors_celltypes, raster = TRUE, label = TRUE, repel = TRUE)
DimPlot(integrated, reduction = "umap_integrated_rpca", group.by = 'annotation_2', cols = colors_celltypes, raster = TRUE, pt.size = 1.2, label = TRUE, repel = TRUE) + NoLegend()
dev.off()

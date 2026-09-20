library(ggplot2)
library(dplyr)
library(Seurat)


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

dataset_colors <- c(
	"Xpand" 	= "#1B9E77", 
	"Ainciburu" = "#D95F02", 
	"Li" 		= "#7570B3", 
	"HCA" 		= "#E7298A")

colors_h1 <- c(
  `H1-1` = "#eb5e28", `H1-2` = "#e09f3e",
  `H1-3` = "#bf0603", `H1-4` = "#f4acb7",
  `H1-5` = "#40916c", `H1-0` = "#5a189a",
  `H1-10` = "#0077b6"
)

########################################################################################################
############################################# FIGURE 1A-1B #############################################
########################################################################################################
table(all_datasets$Dataset, all_datasets$annotation_2)
# UMAPS
library(Seurat)
integrated <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/integration_post_gene_fix.rds')
hca_dataset <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_integrated_Marcal_analysis.rds')
colnames(hca_dataset@meta.data)[colnames(hca_dataset@meta.data) == "percent_H1"] <- "percent.h1"
colnames(hca_dataset@meta.data)[colnames(hca_dataset@meta.data) == "Donor"] <- "sample"

# Fill NA values in annotation_2 with CellType values 
# (NAs  correspond to mature cell populations that had annotation only in the CellType column, that contains  original anotation by the authors)
hca_dataset$annotation_2 <- ifelse(
  is.na(hca_dataset$annotation_2),
  hca_dataset$CellType,
  hca_dataset$annotation_2
)

# Keep common columns and merge the two objects
common_cols <- intersect(colnames(integrated@meta.data), colnames(hca_dataset@meta.data))
integrated@meta.data <- integrated@meta.data[, common_cols]
hca_dataset@meta.data <- hca_dataset@meta.data[, common_cols]
all_datasets <- merge(integrated, hca_dataset)
all_datasets <- JoinLayers(all_datasets)
# saveRDS(all_datasets, "/ibex/user/sotoorda/masterh1/public_data_analysis/merged_datasets/objects/all_datasets_merged.rds")

# Ordering the labels
v1 <- as.vector(names(colors_celltypes))
v2 <- c('HSC', 'MPP', 'CLP', 'ProB', 'PreB', 'MEP', 'ERP', 'MKP', 'GMP', 'Eo/B/Mast', 'MDP')
v1 <- c(v2[v2 %in% v1], setdiff(v1, v2))


hca_dataset <- RunUMAP(hca_dataset, dims = 1:30, reduction = "integrated.rpca", reduction.name = "umap.integrated")
hca_dataset$annotation_2 <- factor(hca_dataset$annotation_2, levels = v1)

# saveRDS(hca_dataset, '/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_merging_ready.rds')

# UMAP HCA
pdf('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_1/umap_hca.pdf', width = 6, height = 6)
DimPlot(hca_dataset, reduction = "umap.integrated", group.by = 'annotation_2', cols = colors_celltypes, raster = TRUE, label = TRUE, repel = TRUE, pt.size = 1.2)
DimPlot(hca_dataset, reduction = "umap.integrated", group.by = 'annotation_2', cols = colors_celltypes, raster = TRUE, pt.size = 1.2, label = TRUE, repel = TRUE) + NoLegend()
dev.off()

# UMAP Integrated
pdf('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_1/umap_integrated.pdf', width = 6, height = 6)
DimPlot(integrated, reduction = "umap_integrated_rpca", group.by = 'annotation_2', cols = colors_celltypes, raster = TRUE, label = TRUE, repel = TRUE)
DimPlot(integrated, reduction = "umap_integrated_rpca", group.by = 'annotation_2', cols = colors_celltypes, raster = TRUE, pt.size = 1.2, label = TRUE, repel = TRUE) + NoLegend()
dev.off()
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
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


########################################################################################################
############################################# FIGURE 1A-1B #############################################
########################################################################################################

# UMAPS
library(Seurat)
integrated <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/integration_post_gene_fix.rds')
hca_dataset <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_integrated_Marcal_analysis.rds')
colnames(hca_dataset@meta.data)
colnames(hca_dataset@meta.data)[colnames(hca_dataset@meta.data) == "Donor"] <- "sample"

# Fill NA values in annotation_2 with CellType values (There were NAs previously)
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

####################################################################################
################################### FIGURE 1C ######################################
############################# CELL POPULATION BARPLOTS #############################
####################################################################################

int <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/integration_post_gene_fix.rds')
hca <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_merging_ready.rds')
all_datasets <- merge(int, hca)
all_datasets <- JoinLayers(all_datasets)
all_datasets

# Ordering the labels
v1 <- as.vector(names(colors_celltypes))
v2 <- c('HSC', 'MPP', 'CLP', 'ProB', 'PreB', 'MEP', 'ERP', 'MKP', 'GMP', 'Eo/B/Mast', 'MDP')
v1 <- c(v2[v2 %in% v1], setdiff(v1, v2))

celltype_counts <- all_datasets@meta.data %>%
  group_by(Dataset, annotation_2) %>%
  summarise(cell_count = n()) %>%  # Count the number of cells for each CellType in each sample
  ungroup()
celltype_counts$annotation_2 <- factor(celltype_counts$annotation_2, levels = v1)
celltype_counts$Dataset <- factor(celltype_counts$Dataset, levels = c("Xpand", "Ainciburu", "Li", "HCA"))
pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_1/stackedbarplot_cell_composition_vertical.pdf", width = 5, height = 10)
ggplot(celltype_counts, aes(x = annotation_2, y = cell_count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = dataset_colors) +
  labs(title = "Cell Type Composition Across Datasets",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_1/stackedbarplot_cell_composition_horizontal.pdf", width = 12, height = 4)
ggplot(celltype_counts, aes(x = annotation_2, y = cell_count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "stack") +  # Use dodge to split the bars
  scale_fill_manual(values = dataset_colors) +
  labs(title = "Cell Type Composition Across Datasets",
       x = "Cell Type",
       y = "Number of Cells") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

####################################################################################
################################### FIGURE 1D ######################################
############################# CELL POPULATION BARPLOTS #############################
####################################################################################

##############################################################################################################################
######################################################## FIGURE 1E ###########################################################
########################################### HEATMAP OF FEATURE IMPORTANCE RANKINGS ###########################################
##############################################################################################################################

input_dir <- "/ibex/user/sotoorda/masterh1/public_data_analysis/ML_analysis/Feature_importance_order_csv"
dataset_files <- c(
  Ainciburu = file.path(input_dir, "Ainciburu_feature_importance_order.csv"),
  Li        = file.path(input_dir, "Li_feature_importance_order.csv"),
  Xpand     = file.path(input_dir, "Xpand_feature_importance_order.csv"), 
  HCA       = file.path(input_dir, "HCA_feature_importance_order.csv")
)
output_dir <- "/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_1"

h1_genes <- c("H1-0", "H1-1", "H1-2", "H1-3", "H1-4", "H1-5", "H1-10")

# Merge tables
read_rank_table <- function(file_path, dataset_name) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  df <- df[, c("gene", "importance_order", "mean_importance")]
  names(df) <- c("gene", dataset_name, paste0(dataset_name, "_mean_importance"))
  df
}
rank_tables <- Map(read_rank_table, dataset_files, names(dataset_files))
rank_table <- Reduce(function(l, r) merge(l, r, by = "gene", all = TRUE), rank_tables)

setdiff(rank_tables$HCA$gene, rank_tables$Li$gene)

dataset_cols <- names(dataset_files)
# compute consensus rank and sort
rank_table$consensus_rank <- rowMeans(rank_table[, dataset_cols], na.rm = TRUE)
rank_table <- rank_table[order(rank_table$consensus_rank, rank_table$gene), ]

# build matrix
mat <- as.matrix(rank_table[, dataset_cols])
rownames(mat) <- rank_table$gene
# numeric matrix (should already be numeric from CSV); ensure numeric
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- rank_table$gene

vals <- mat[!is.na(mat)]
rng <- range(vals, na.rm = TRUE)
mid <- median(vals, na.rm = TRUE)
col_fun <- colorRamp2(c(1, 40, 80), c("#ee3e32", "white", "steelblue"))

# cell label function
cell_text <- function(j, i, x, y, width, height, fill) {
  v <- mat[rownames(mat)[i], colnames(mat)[j]]
  if (!is.na(v)) grid::grid.text(as.character(round(v)), x = x, y = y, gp = grid::gpar(fontsize = 8))
}

# Short heatmap (top non-H1 + H1 genes like original script)
non_h1_table <- rank_table[!rank_table$gene %in% h1_genes, ]
top15_non_h1 <- head(non_h1_table$gene, 15)
selected_genes <- unique(c(top15_non_h1, intersect(h1_genes, rank_table$gene)))
short_mat <- mat[selected_genes, , drop = FALSE]
short_file <- file.path(output_dir, "feature_ranking_heatmap_short_complexheatmap.pdf")
n_genes_short <- nrow(short_mat)

pdf(short_file, width = 3, height = 6)
ht_short <- Heatmap(
  short_mat,
  name = "Ranking",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 10),
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- short_mat[rownames(short_mat)[i], colnames(short_mat)[j]]
    if (!is.na(v)) grid::grid.text(as.character(round(v)), x = x, y = y, gp = grid::gpar(fontsize = 10))
  },
  heatmap_legend_param = list(
    title = "Ranking",
    at = c(rng[1], round(mid), rng[2]),
    direction = "horizontal",
    labels_gp = grid::gpar(fontsize = 10),
    title_gp = grid::gpar(fontsize = 10),
    legend_height = grid::unit(6, "mm")
  )
)
draw(ht_short, heatmap_legend_side = "bottom")
dev.off()


####################################################################################
################################### FIGURE 1F ######################################
############################# PERCENTAGE OF H1 BARPLOTS ############################
####################################################################################

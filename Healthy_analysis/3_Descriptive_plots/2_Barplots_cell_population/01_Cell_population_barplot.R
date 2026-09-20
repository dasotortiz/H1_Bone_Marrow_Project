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

colors_h1 <- c(
  `H1-1` = "#eb5e28", `H1-2` = "#e09f3e",
  `H1-3` = "#bf0603", `H1-4` = "#f4acb7",
  `H1-5` = "#40916c", `H1-0` = "#5a189a",
  `H1-10` = "#0077b6"
)

####################################################################################
################################### FIGURE 1C ######################################
############################# CELL POPULATION BARPLOTS #############################
####################################################################################

all_datasets <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/merged_datasets/objects/all_datasets_merged.rds")

# Ordering the labels
v1 <- as.vector(names(colors_celltypes))
v2 <- c('HSC', 'MPP', 'CLP', 'ProB', 'PreB', 'MEP', 'ERP', 'MKP', 'GMP', 'Eo/B/Mast', 'MDP')
mature <- setdiff(v1, v2)

all_datasets_early <- subset(all_datasets, subset = annotation_2 %in% v2)
all_datasets_mature <- subset(all_datasets, subset = annotation_2 %in% mature)

# Barplot for early cell populations

celltype_counts <- all_datasets_early@meta.data %>%
  group_by(Dataset, annotation_2) %>%
  summarise(cell_count = n()) %>%  # Count the number of cells for each CellType in each sample
  ungroup()
celltype_counts$annotation_2 <- factor(celltype_counts$annotation_2, levels = v2)
celltype_counts$Dataset <- factor(celltype_counts$Dataset, levels = c("Xpand", "Ainciburu", "Li", "HCA"))


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

# Mature cell populations

unique(all_datasets_mature$annotation_2)
celltype_counts <- all_datasets_mature@meta.data %>%
  group_by(annotation_2) %>%
  summarise(cell_count = n()) %>%  # Count the number of cells for each CellType in each sample
  ungroup()
celltype_counts$annotation_2 <- factor(celltype_counts$annotation_2, levels = mature)


pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_1/stackedbarplot_cell_composition_vertical.pdf", width = 5, height = 10)
ggplot(celltype_counts, aes(x = annotation_2, y = cell_count)) +
  geom_col(fill = dataset_colors[["HCA"]]) +
  labs(title = "Cell Type Composition Across Datasets",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none")
dev.off()
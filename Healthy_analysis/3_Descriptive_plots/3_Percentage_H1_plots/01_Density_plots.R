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
################################### FIGURE 1D ######################################
############################### H1 PERCENTAGE DENSITIES ############################
####################################################################################

library(ggridges)
library(tidyr)

h1_genes <- c("H1-0", "H1-1", "H1-2", "H1-3",
              "H1-4", "H1-5", "H1-10")


all_datasets <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/merged_datasets/objects/all_datasets_merged.rds")

# One page will be generated for each lineage below
density_lineages <- list(
  Lymphoid = c("HSC", "MPP", "CLP", "ProB", "PreB"),
  Erythroid = c("HSC", "MPP", "MEP", "ERP"),
  `Monocyte-Dendritic` = c("HSC", "MPP", "MDP"),
  Granulocyte = c("HSC", "MPP", "GMP", "Eo/B/Mast")
)

percent_cols <- paste0("percent.", gsub("-", ".", h1_genes))
percent_cols_present <- intersect(percent_cols, colnames(all_datasets@meta.data))

# Use the per-cell H1 percentages already stored in the Seurat metadata
density_df <- all_datasets@meta.data %>%
  tibble::rownames_to_column("Cell") %>%
  dplyr::select(Cell, Dataset, annotation_2, dplyr::all_of(percent_cols_present)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(percent_cols_present),
    names_to = "percent_col",
    values_to = "percentage"
  ) %>%
  dplyr::mutate(
    H1_gene = sub("^percent\\.", "", percent_col),
    H1_gene = gsub("\\.", "-", H1_gene)
  ) %>%
  dplyr::filter(
    !is.na(Dataset),
    !is.na(annotation_2),
    is.finite(percentage)
  )

density_df$H1_gene <- factor(density_df$H1_gene, levels = h1_genes)
density_df$Dataset <- factor(
  density_df$Dataset,
  levels = c("Xpand", "Ainciburu", "Li", "HCA")
)

# One PDF page per lineage
pdf(
  "/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_1/H1_percentage_density_by_lineage.pdf",
  width = 14,
  height = 7
)

for (lineage_name in names(density_lineages)) {

  lineage_data <- density_df %>%
    dplyr::filter(annotation_2 %in% density_lineages[[lineage_name]]) %>%
    dplyr::mutate(
      annotation_2 = factor(
        annotation_2,
        levels = rev(density_lineages[[lineage_name]])
      )
    )

  # Use the 95th percentile within each dataset to limit displayed outliers.
  x_max <- lineage_data %>%
    dplyr::summarise(
      x_max = quantile(percentage, probs = 0.95, na.rm = TRUE),
      .by = Dataset
    ) %>%
    dplyr::summarise(x_max = max(x_max, na.rm = TRUE)) %>%
    dplyr::pull(x_max)

  p <- ggplot(
    lineage_data,
    aes(
      x = percentage,
      y = annotation_2,
      fill = H1_gene,
      color = H1_gene,
      group = interaction(annotation_2, H1_gene)
    )
  ) +
    ggridges::geom_density_ridges(
      alpha = 0.45,
      scale = 1,
      rel_min_height = 0.01,
      na.rm = TRUE,
      # color = "white",
      size = 0.2
    ) +
    facet_wrap(~Dataset, nrow = 1, drop = FALSE) +
    scale_fill_manual(values = colors_h1, drop = FALSE) +
    scale_color_manual(values = colors_h1, drop = FALSE) +
    coord_cartesian(xlim = c(0, x_max)) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    labs(
      title = paste(lineage_name, "lineage"),
      x = "H1 counts as percentage of total RNA counts",
      y = "Cell population",
      fill = "H1 gene",
      color = "H1 gene"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  print(p)
}

dev.off()
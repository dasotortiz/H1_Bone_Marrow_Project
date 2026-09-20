library(dplyr)
library(ComplexHeatmap)
library(circlize)


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
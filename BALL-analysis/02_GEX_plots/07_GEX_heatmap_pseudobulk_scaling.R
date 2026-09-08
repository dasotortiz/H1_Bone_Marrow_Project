# Per-gene heatmaps of pseudobulk mean expression along the lymphoid lineage,
# z-scored against the celltypes that healthy and disease groups have in common.
# Input: ./output/i_df.rds and ./output/h_df.rds (from 06_build_pseudobulk_summaries.R)

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

lineage_cells <- c("HSC", "MPP", "CLP", "ProB", "PreB")
common_celltypes <- c("CLP", "ProB", "PreB")  # the only ones present in both healthy and disease

# h1_df.rds lives outside this repo, set the path to wherever you keep it
h1_path <- "../h1_df.rds"
h1_df <- readRDS(h1_path)
gene_order <- h1_df$Symbol

# Starting from pbs + mean + SD per subtype dfs
i_df <- readRDS("./output/i_df.rds")
h_df <- readRDS("./output/h_df.rds")
df <- bind_rows(h_df, i_df)
write.csv(as.data.frame(df), "./output/df_healthy_ia.csv", row.names = FALSE)  # the lineplot script reads this

df$Gene     <- factor(df$Gene,     levels = gene_order)
df$Celltype <- factor(df$Celltype, levels = lineage_cells)

genes_present <- gene_order[gene_order %in% as.character(df$Gene)]

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

hm_colors     <- subtype_colors
status_colors <- c("Disease" = "red", "Healthy" = "#156415")
col_fun       <- colorRamp2(c(-2, 0, 2), c("royalblue", "white", "firebrick"))

dev_stage_mapping <- c(
  "ZNF384-R" = "stem", "DUX4-R" = "stem",
  "BCR::ABL1" = "early_lymph", "BCR::ABL1-like" = "early_lymph",
  "KMT2A-R" = "early_lymph", "PAX5alt" = "early_lymph",
  "ETV6::RUNX1-like" = "early_lymph", "ETV6::RUNX1" = "early_lymph",
  "Hyperdiploid" = "early_lymph", "Low-hypodiploid" = "early_lymph",
  "Near-haploid" = "early_lymph", "iAMP21" = "early_lymph",
  "TCF3::PBX1" = "diff", "MEF2D-R" = "diff",
  # Healthy datasets mapped individually
  "Xpand" = "healthy", "Ainciburu" = "healthy", "Li" = "healthy"
)
dev_stage_colors_full <- c(
  "stem" = "#1b9e77", "early_lymph" = "#d95f02",
  "diff" = "#7570b3", "healthy" = "grey60"
)

# Custom scaling: mean and SD come from the 3 common celltypes only,
# but the correction is applied to every celltype in the row.
scale_by_common_cols <- function(mat, common_cols) {
  cols_for_stats <- intersect(common_cols, colnames(mat))

  t(apply(mat, 1, function(row) {
    ref_values <- row[cols_for_stats]
    ref_mean   <- mean(ref_values, na.rm = TRUE)
    ref_sd     <- sd(ref_values,   na.rm = TRUE)
    # Avoid division by zero (e.g. if all ref values are identical)
    if (is.na(ref_sd) || ref_sd == 0) ref_sd <- 1
    (row - ref_mean) / ref_sd
  }))
}

# Scaled matrix plus its row annotations for one gene
gene_heatmap_data <- function(g) {
  gene_df <- df %>% filter(Gene == g)

  gene_mat <- gene_df %>%
    select(Subtype, Celltype, Mean) %>%
    pivot_wider(names_from = Celltype, values_from = Mean) %>%
    column_to_rownames("Subtype") %>%
    as.matrix()
  # Enforce lineage_cells column order
  gene_mat <- gene_mat[, intersect(lineage_cells, colnames(gene_mat))]

  scaled_mat <- scale_by_common_cols(gene_mat, common_celltypes)
  dimnames(scaled_mat) <- dimnames(gene_mat)

  anno_df <- gene_df %>%
    distinct(Subtype, Status) %>%
    column_to_rownames("Subtype")
  anno_df <- anno_df[rownames(scaled_mat), , drop = FALSE]

  list(mat = scaled_mat, anno = anno_df)
}

# show_devstage adds the developmental stage annotation column,
# show_legends off is for the grid PDFs further down where legends are shared
gene_heatmap <- function(g, show_devstage = TRUE, show_legends = TRUE) {
  hm <- gene_heatmap_data(g)
  anno_args <- list(Subtype = rownames(hm$mat), Status = hm$anno$Status)
  anno_cols <- list(Subtype = hm_colors, Status = status_colors)

  if (show_devstage) {
    stages <- dev_stage_mapping[rownames(hm$anno)]
    present_stages <- unique(na.omit(stages))
    anno_args$DevStage <- factor(stages, levels = present_stages)
    anno_cols$DevStage <- dev_stage_colors_full[present_stages]
  }

  row_ha <- do.call(rowAnnotation, c(
    anno_args,
    list(col = anno_cols, show_annotation_name = TRUE, show_legend = show_legends)
  ))

  Heatmap(
    hm$mat,
    name             = "Z-score",
    column_title     = paste("Gene:", g),
    col              = col_fun,
    cluster_rows     = TRUE,
    cluster_columns  = FALSE,
    left_annotation  = row_ha,
    rect_gp          = gpar(col = "white", lwd = 0.5),
    row_names_side   = "right",
    column_names_rot = 45,
    show_heatmap_legend = show_legends
  )
}

# One heatmap per page, with and without the DevStage annotation
write_gene_heatmap_pdf <- function(file, show_devstage) {
  pdf(file, width = 8, height = 6)
  for (g in genes_present) {
    draw(gene_heatmap(g, show_devstage = show_devstage), merge_legend = TRUE)
  }
  dev.off()
  message("Saved: ", file)
}

write_gene_heatmap_pdf("./output/pseudobulk_gene_heatmaps_CustomScale.pdf", show_devstage = TRUE)
write_gene_heatmap_pdf("./output/pseudobulk_gene_heatmaps_CustomScale_noDevStage.pdf", show_devstage = FALSE)

# Manuscript PDFs: several no-DevStage heatmaps on one page with shared legends
make_shared_legends <- function(genes) {
  hm_list <- lapply(genes, gene_heatmap_data)
  present_subtypes <- unique(unlist(lapply(hm_list, function(x) rownames(x$mat))))
  present_status <- unique(na.omit(unlist(lapply(hm_list, function(x) as.character(x$anno$Status)))))

  lgd_z <- Legend(title = "Z-score", col_fun = col_fun,
                  direction = "horizontal", legend_width = unit(35, "mm"))
  lgd_status <- Legend(title = "Status", labels = present_status,
                       legend_gp = gpar(fill = status_colors[present_status]))
  lgd_subtype <- Legend(title = "Subtype", labels = present_subtypes,
                        legend_gp = gpar(fill = hm_colors[present_subtypes]), nrow = 2)

  packLegend(lgd_z, lgd_status, lgd_subtype, direction = "horizontal", gap = unit(6, "mm"))
}

draw_gene_heatmap_grid_pdf <- function(genes, file, ncol, nrow, width, height) {
  genes <- genes[genes %in% genes_present]
  if (length(genes) == 0) {
    warning("No requested genes were found in df$Gene. Skipping: ", file)
    return(invisible(NULL))
  }

  heatmap_grobs <- lapply(genes, function(g) {
    grid.grabExpr(draw(gene_heatmap(g, show_devstage = FALSE, show_legends = FALSE),
                       merge_legend = FALSE))
  })
  shared_legends <- make_shared_legends(genes)

  pdf(file, width = width, height = height)
  grid.newpage()
  # extra row at the bottom holds the shared legends
  pushViewport(viewport(
    layout = grid.layout(
      nrow = nrow + 1,
      ncol = ncol,
      heights = unit.c(unit(rep(1, nrow), "null"), unit(3, "cm"))
    )
  ))

  for (i in seq_along(heatmap_grobs)) {
    pushViewport(viewport(layout.pos.row = ceiling(i / ncol),
                          layout.pos.col = ((i - 1) %% ncol) + 1))
    grid.draw(heatmap_grobs[[i]])
    popViewport()
  }

  pushViewport(viewport(layout.pos.row = nrow + 1, layout.pos.col = 1:ncol))
  draw(shared_legends, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center")
  popViewport(2)
  dev.off()

  message("Saved: ", file)
}

h1_pair_genes <- c("H1-0", "H1-10")
h1_rest_genes <- setdiff(as.character(gene_order), h1_pair_genes)

draw_gene_heatmap_grid_pdf(
  genes = h1_pair_genes,
  file = "./output/pseudobulk_H1_0_H1_10_side_by_side_noDevStage.pdf",
  ncol = 2, nrow = 1, width = 14, height = 7
)

draw_gene_heatmap_grid_pdf(
  genes = h1_rest_genes,
  file = "./output/pseudobulk_remaining_H1_genes_grid_noDevStage.pdf",
  ncol = 2, nrow = 3, width = 14, height = 14
)

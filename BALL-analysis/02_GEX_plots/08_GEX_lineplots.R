# Lineplots of pseudobulk mean expression along the lymphoid lineage,
# healthy references next to selected B-ALL subtypes.
# Input: ./output/df_healthy_ia.csv, written by 07_GEX_heatmap_pseudobulk_scaling.R

library(dplyr)
library(ggplot2)

df <- read.csv("./output/df_healthy_ia.csv")
lineage_cells <- c("HSC", "MPP", "CLP", "ProB", "PreB")
df$Celltype <- factor(df$Celltype, levels = lineage_cells)

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

colors_datasets <- subtype_colors[c("Xpand", "Ainciburu", "Li")]
dataset_shapes  <- c("Xpand" = 15, "Ainciburu" = 17, "Li" = 18)

# Healthy references on the left panel, the chosen disease subtypes on the right
gene_lineplot <- function(gene, target_disease) {
  plot_df <- df %>%
    filter(Gene == gene) %>%
    filter(Status == "Healthy" | (Status == "Disease" & Subtype %in% target_disease)) %>%
    mutate(Panel = factor(ifelse(Status == "Healthy", "Healthy Reference", "Selected Disease Subtypes"),
                          levels = c("Healthy Reference", "Selected Disease Subtypes")))

  current_colors <- c(colors_datasets, subtype_colors[target_disease])
  current_shapes <- c(dataset_shapes, setNames(rep(16, length(target_disease)), target_disease))

  ggplot(plot_df, aes(x = Celltype, y = Mean, group = Subtype,
                      color = Subtype, fill = Subtype, shape = Subtype)) +
    facet_wrap(~ Panel, scales = "free_y") +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.15, color = NA) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Cell Type",
      y = "Mean Expression",
      title = paste("Gene Expression Profile:", gene)
    ) +
    scale_color_manual(values = current_colors) +
    scale_fill_manual(values = current_colors) +
    scale_shape_manual(values = current_shapes)
}

# H1-0: the two subtypes projecting to the stem end
# target_disease was c("DUX4-R", "ZNF384-R", "Hyperdiploid") before dropping Hyperdiploid
p_h1_0 <- gene_lineplot("H1-0", c("DUX4-R", "ZNF384-R"))
ggsave("./output/pseudobulk_H1-0_lineplot.pdf", p_h1_0, width = 16, height = 4)

# H1-10: the more differentiated subtypes
p_h1_10 <- gene_lineplot("H1-10", c("TCF3::PBX1", "Hyperdiploid"))
ggsave("./output/pseudobulk_H1-10_lineplot.pdf", p_h1_10, width = 16, height = 4)

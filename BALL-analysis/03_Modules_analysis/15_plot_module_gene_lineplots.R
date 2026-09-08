# Lineplots of individual module genes along the lymphoid lineage,
# the 3 healthy datasets next to a chosen disease subtype.

library(Seurat)
library(tidyverse)
library(patchwork)

pb_hc_filt <- readRDS("../GEX-Plotting/output/pseudobulk_hc.rds")
pb_ia_filt <- readRDS("../GEX-Plotting/output/pseudobulk_Ia.rds")

# Genes to plot, edit these two lines per figure
module_genes   <- c("SOX4")
target_disease <- "TCF3::PBX1"

genes_to_use <- Reduce(intersect, list(module_genes, rownames(pb_hc_filt), rownames(pb_ia_filt)))
missing <- setdiff(module_genes, genes_to_use)
if (length(missing) > 0) message("Not found in both objects, skipping: ", paste(missing, collapse = ", "))

lineage_cells <- c("HSC", "MPP", "CLP", "ProB", "PreB")

# Healthy, grouped by Dataset so each of the three shows as its own line
h_df <- FetchData(pb_hc_filt, vars = c(genes_to_use, "orig.ident", "Celltype", "Dataset"),
                  layer = "data") %>%
  pivot_longer(cols = all_of(genes_to_use), names_to = "Gene", values_to = "Expression") %>%
  rename(Subtype = Dataset) %>%
  group_by(Subtype, Celltype, Gene) %>%
  summarise(Mean = mean(Expression, na.rm = TRUE),
            SD   = sd(Expression,   na.rm = TRUE), .groups = "drop") %>%
  mutate(Status = "Healthy")

# Disease, grouped by Subtype
i_df <- FetchData(pb_ia_filt, vars = c(genes_to_use, "orig.ident", "Celltype", "Subtype"),
                  layer = "data") %>%
  pivot_longer(cols = all_of(genes_to_use), names_to = "Gene", values_to = "Expression") %>%
  group_by(Subtype, Celltype, Gene) %>%
  summarise(Mean = mean(Expression, na.rm = TRUE),
            SD   = sd(Expression,   na.rm = TRUE), .groups = "drop") %>%
  mutate(Status = "Disease")

expr_data <- bind_rows(h_df, i_df) %>%
  filter(Celltype %in% lineage_cells) %>%
  mutate(Celltype = factor(Celltype, levels = lineage_cells))

subtype_colors <- c(
  "Hyperdiploid"     = "#56B4E9",
  "DUX4-R"           = "#E69F00",
  "ZNF384-R"         = "#CC79A7",
  "TCF3::PBX1"       = "#F0E442",
  "BCR::ABL1"        = "#0072B2",
  "BCR::ABL1-like"   = "#6B4C9A",
  "ETV6::RUNX1"      = "#C4A882",
  "ETV6::RUNX1-like" = "#44AA99",
  "iAMP21"           = "#B2DF8A",
  "KMT2A-R"          = "#000000",
  "Low-hypodiploid"  = "#4B5D73",
  "Near-haploid"     = "#C1703A",
  "PAX5alt"          = "#F4A3B0",
  "MEF2D-R"          = "#FDB462",
  "Other"            = "#BEBEBE"
)

colors_datasets <- c(Xpand = "#1B9E77", Ainciburu = "#D95F02", Li = "#7570B3")
current_colors  <- c(colors_datasets, subtype_colors[target_disease])

# One panel per gene, healthy on the left facet and the subtype on the right
plots <- lapply(genes_to_use, function(gene) {
  plot_df <- expr_data %>%
    filter(Gene == gene) %>%
    filter(Status == "Healthy" | (Status == "Disease" & Subtype %in% target_disease)) %>%
    mutate(Panel = factor(ifelse(Status == "Healthy", "Healthy Reference", "Selected Disease Subtypes"),
                          levels = c("Healthy Reference", "Selected Disease Subtypes")))

  ggplot(plot_df, aes(x = Celltype, y = Mean, group = Subtype,
                      color = Subtype, fill = Subtype)) +
    facet_wrap(~ Panel, scales = "fixed") +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = current_colors) +
    scale_fill_manual(values = current_colors) +
    labs(x = "Cell Type", y = "Mean Expression", title = gene) +
    theme_bw() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    )
})

ggsave("./output/GEX_lineplot_TCF-SOX.pdf", wrap_plots(plots, ncol = 1),
       width = 16, height = 4 * length(plots))

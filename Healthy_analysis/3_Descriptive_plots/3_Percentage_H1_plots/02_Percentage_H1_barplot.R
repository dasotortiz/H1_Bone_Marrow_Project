library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(Seurat)



####################################################################################
################################### FIGURE 1F ######################################
############################# PERCENTAGE OF H1 BARPLOTS ############################
####################################################################################

genes <- c("H1-0", "H1-1", "H1-2", "H1-3", "H1-4", "H1-5", "H1-10")
lineages <- list(
  "Lymphoid" = c("HSC", "MPP", "CLP", "ProB", "PreB"),
  "Erythroid" = c("HSC", "MPP", "MEP", "ERP"),
  "Monocyte-Dendritic" = c("HSC", "MPP", "MDP"),
  "Granulocyte" = c("HSC", "MPP", "GMP", "Eo/B/Mast")
)

all_datasets <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/merged_datasets/objects/all_datasets_merged.rds")

# Generating the percentage of H1 subtypes for each cell in the merged dataset
for (gene in genes) {
  percent_col <- paste0("percent.", gsub("-", ".", gene))
  all_datasets[[percent_col]] <- PercentageFeatureSet(all_datasets, features = gene)
}



lineage_df <- purrr::imap_dfr(
  lineages,
  ~ tibble::tibble(Lineage = .y, annotation_2 = .x)
)

percent_cols <- paste0("percent.", gsub("-", ".", genes))

# Create a data frame with the mean percentage of each H1 subtype for each Dataset and annotation_2 combination
plot_data <- all_datasets@meta.data %>%
  dplyr::select(Dataset, annotation_2, dplyr::all_of(percent_cols)) %>%
  dplyr::group_by(Dataset, annotation_2) %>%
  # Compute the mean percentage for each H1 subtype within cells of the same Dataset and annotation_2 group
  dplyr::summarise( 
    dplyr::across(dplyr::all_of(percent_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(percent_cols),
    names_to = "percent_col",
    values_to = "percentage"
  ) %>%
  dplyr::mutate(
    H1_subtype = sub("^percent\\.", "", percent_col),
    H1_subtype = gsub("\\.", "-", H1_subtype)
  ) %>%
  dplyr::select(-percent_col) %>%
  # Join together with the lineage information to get the Lineage column
  dplyr::inner_join(lineage_df, by = "annotation_2") %>%
  dplyr::mutate(
    H1_subtype = factor(H1_subtype, levels = genes),
    Lineage = factor(Lineage, levels = names(lineages)),
    annotation_2 = factor(annotation_2, levels = unique(unlist(lineages)))
  )

pdf(
  "/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/plots/Manuscript_figure/Panel_1/H1_subtype_percentages_by_lineage.pdf",
  width = 9,
  height = 3
)

for (lineage_name in names(lineages)) {
  lineage_plot_data <- plot_data %>%
    dplyr::filter(Lineage == lineage_name) %>%
    dplyr::mutate(annotation_2 = factor(annotation_2, levels = lineages[[lineage_name]]))

  h1_percentage_plot <- ggplot(
    lineage_plot_data,
    aes(x = annotation_2, y = percentage, fill = H1_subtype)
  ) +
    geom_col(width = 0.8) +
    facet_wrap(~Dataset, nrow = 1) +
    scale_fill_manual(values = colors_h1, drop = FALSE) +
    labs(
      title = lineage_name,
      x = "Cell population",
      y = "H1 subtype percentage",
      fill = "H1 subtype"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )

  print(h1_percentage_plot)
}

dev.off()
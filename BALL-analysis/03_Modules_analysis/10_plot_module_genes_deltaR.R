# Delta r per module gene: disease correlation minus the healthy mean correlation.
# Input:  ./output/module_gene_correlations.csv (from 09)
# Output: ./output/delta_r_correlations.csv, read by 11 and 14.

library(tidyverse)
library(patchwork)

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

cor_df <- read.csv("./output/module_gene_correlations.csv")

# Mean healthy r per gene (across the 3 datasets), then subtract from each subtype
delta_df <- cor_df %>%
  group_by(H1_gene, Gene) %>%
  summarise(r_healthy = mean(Correlation[Status == "Healthy"], na.rm = TRUE), .groups = "drop") %>%
  left_join(
    cor_df %>%
      filter(Status == "Disease") %>%
      select(H1_gene, Gene, Group, Correlation) %>%
      rename(r_disease = Correlation),
    by = c("H1_gene", "Gene")
  ) %>%
  mutate(delta_r = r_disease - r_healthy)

write.csv(delta_df, "./output/delta_r_correlations.csv", row.names = FALSE)

# The raw r dot plot and density are drawn twice below (all subtypes, then a
# selected pair), so they live in these two helpers.
raw_r_dotplot <- function(sub_raw, h1_name, subtitle, plot_colors) {
  ggplot(sub_raw, aes(x = Gene, y = Correlation, color = color_group, shape = Status)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_jitter(data = filter(sub_raw, Status == "Healthy"),
                width = 0.25, size = 5, alpha = 0.15, color = subtype_colors[["Healthy"]]) +
    geom_jitter(width = 0.25, size = 2.5, alpha = 0.85) +
    scale_color_manual(values = plot_colors) +
    scale_shape_manual(values = c("Healthy" = 16, "Disease" = 17)) +
    guides(color = guide_legend(title = "Group",  nrow = 2, byrow = TRUE, override.aes = list(size = 3)),
           shape = guide_legend(title = "Status", override.aes = list(size = 3))) +
    labs(
      title    = paste("Raw Spearman r - module genes vs", h1_name),
      subtitle = subtitle,
      x        = NULL,
      y        = "Spearman r"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.x = element_line(color = "grey92"),
      legend.position    = "bottom",
      axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

raw_r_density <- function(sub_raw, h1_name, subtitle, plot_colors) {
  ggplot(sub_raw, aes(x = Correlation, color = color_group, fill = color_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_density(alpha = 0.15, linewidth = 0.8) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values  = plot_colors) +
    guides(color = guide_legend(title = "Group", nrow = 2, byrow = TRUE),
           fill  = guide_legend(title = "Group", ncol = 2)) +
    labs(
      title    = paste("Density of raw r - module genes vs", h1_name),
      subtitle = subtitle,
      x        = "Spearman r",
      y        = "Density"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

# PDF 1: one H1 gene per 4 pages, all subtypes overlaid
pdf("./output/delta_r_overall.pdf", width = 16, height = 8)

for (h1_name in c("H1-10", "H1-0")) {

  sub_delta <- delta_df %>% filter(H1_gene == h1_name)
  sub_raw   <- cor_df   %>% filter(H1_gene == h1_name)

  # Rank genes by mean delta r across subtypes, keeps the x axis consistent
  gene_order <- sub_delta %>%
    group_by(Gene) %>%
    summarise(mean_delta = mean(delta_r, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_delta) %>%
    pull(Gene)

  sub_delta$Gene <- factor(sub_delta$Gene, levels = gene_order)
  sub_raw$Gene   <- factor(sub_raw$Gene,   levels = gene_order)

  present_subtypes <- unique(sub_delta$Group)
  disease_colors   <- subtype_colors[names(subtype_colors) %in% present_subtypes]
  plot_colors      <- c("Healthy" = subtype_colors[["Healthy"]], disease_colors)

  sub_raw <- sub_raw %>% mutate(color_group = if_else(Status == "Healthy", "Healthy", Group))

  subtitle_all <- "Healthy (black) + all disease subtypes"

  print(raw_r_dotplot(sub_raw, h1_name, subtitle_all, plot_colors))

  # delta r dot plot, disease only
  print(
    ggplot(sub_delta, aes(x = Gene, y = delta_r, color = Group)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_jitter(width = 0.25, size = 2.5, alpha = 0.85) +
      scale_color_manual(values = disease_colors) +
      guides(color = guide_legend(title = "Subtype", nrow = 2, byrow = TRUE, override.aes = list(size = 3))) +
      labs(
        title    = paste("Delta Spearman r - module genes vs", h1_name),
        subtitle = "delta r = r(disease) - r(healthy mean)",
        x        = NULL,
        y        = "Delta Spearman r"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.major.x = element_line(color = "grey92"),
        legend.position    = "bottom",
        axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  )

  print(raw_r_density(sub_raw, h1_name, subtitle_all, plot_colors))

  # delta r density, disease only
  print(
    ggplot(sub_delta, aes(x = delta_r, color = Group, fill = Group)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      geom_density(alpha = 0.15, linewidth = 0.8) +
      scale_color_manual(values = disease_colors) +
      scale_fill_manual(values  = disease_colors) +
      guides(color = guide_legend(title = "Subtype", nrow = 2, byrow = TRUE),
             fill  = guide_legend(title = "Subtype", ncol = 2)) +
      labs(
        title    = paste("Density of delta r - module genes vs", h1_name),
        subtitle = "Distribution of r(disease) - r(healthy mean) across module genes",
        x        = "Delta Spearman r",
        y        = "Density"
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom")
  )
}

dev.off()
message("Saved: delta_r_overall.pdf")

# PDF 2: one page per H1 gene x subtype, ranked dot plot next to its density
disease_subtypes <- unique(delta_df$Group)

for (h1_name in c("H1-10", "H1-0")) {

  fname <- paste0("./output/delta_r_per_subtype_", gsub("::", "-", h1_name), ".pdf")
  pdf(fname, width = 22, height = 7)

  for (st in disease_subtypes) {

    sub_df <- delta_df %>% filter(H1_gene == h1_name, Group == st)
    if (nrow(sub_df) == 0) next

    st_color <- subtype_colors[st]
    sub_df <- sub_df %>% arrange(delta_r) %>% mutate(Gene = factor(Gene, levels = Gene))

    p_dot <- ggplot(sub_df, aes(x = Gene, y = delta_r)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_segment(aes(xend = Gene, y = 0, yend = delta_r), color = "grey70", linewidth = 0.4) +
      geom_point(color = st_color, size = 2.5, alpha = 0.9) +
      labs(
        title    = paste(h1_name, "-", st),
        subtitle = "Genes ranked by delta r",
        x        = NULL,
        y        = "Delta Spearman r"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_line(color = "grey92")
      )

    # delta r = 0 is the null here, so the rug shows how far genes sit from healthy
    p_dens <- ggplot(sub_df, aes(x = delta_r)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      geom_density(fill = st_color, color = st_color, alpha = 0.3, linewidth = 0.9) +
      geom_rug(color = st_color, alpha = 0.6) +
      labs(
        title    = paste(h1_name, "-", st),
        subtitle = "Density of delta r across module genes",
        x        = "Delta Spearman r",
        y        = "Density"
      ) +
      theme_bw(base_size = 11)

    print(p_dot + p_dens + plot_layout(widths = c(4, 1)))
  }

  dev.off()
  message("Saved: ", fname)
}

# PDF 3: raw r only, healthy plus the two subtypes of interest per H1 gene
subtype_selection <- list(
  "H1-0"  = c("DUX4-R", "ZNF384-R"),
  "H1-10" = c("Hyperdiploid", "TCF3::PBX1")
)

pdf("./output/raw_r_H1-0_H1-10_SpecificSubtypes.pdf", width = 16, height = 8)

for (h1_name in c("H1-10", "H1-0")) {

  selected_subtypes <- subtype_selection[[h1_name]]

  sub_raw <- cor_df %>%
    filter(H1_gene == h1_name, Status == "Healthy" | Group %in% selected_subtypes)

  # Here genes are ranked by mean raw r across all included groups, not by delta r
  gene_order <- sub_raw %>%
    group_by(Gene) %>%
    summarise(mean_cor = mean(Correlation, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_cor) %>%
    pull(Gene)

  sub_raw$Gene <- factor(sub_raw$Gene, levels = gene_order)

  present_subtypes <- unique(sub_raw$Group[sub_raw$Status == "Disease"])
  plot_colors <- c("Healthy" = subtype_colors[["Healthy"]],
                   subtype_colors[names(subtype_colors) %in% present_subtypes])

  sub_raw <- sub_raw %>% mutate(color_group = if_else(Status == "Healthy", "Healthy", Group))

  subtitle_sel <- paste("Healthy (black) +", paste(selected_subtypes, collapse = ", "))

  print(raw_r_dotplot(sub_raw, h1_name, subtitle_sel, plot_colors))
  print(raw_r_density(sub_raw, h1_name, subtitle_sel, plot_colors))
}

dev.off()
message("Saved: raw_r_H1-0_H1-10_SpecificSubtypes.pdf")

# Delta r scatters for the H1-0 and H1-10 modules, with JASPAR TFs highlighted,
# plus the rewiring summary statistics that go into the text.
# Input: ./output/delta_r_correlations.csv (from 10_plot_module_genes_deltaR.R)

library(tidyverse)
library(ggrepel)

# JASPAR2024 TF list, used to flag which module genes are TFs
# Sitting outside the repo, set this to wherever you keep the file
jaspar_path <- "../data/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
lines <- readLines(jaspar_path)
jaspar_tfs <- unique(gsub("^>\\S+\\s+", "", grep("^>", lines, value = TRUE)))
# "IRF1" %in% jaspar_tfs  # TRUE, IRF1 is in the JASPAR2024 TF list

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

quadrant_colors_h10 <- c(
  "Both increased"          = "#6B3F1D",
  "Both decreased"          = "#3F3F3F",
  "DUX4-R only increased"   = subtype_colors[["DUX4-R"]],
  "ZNF384-R only increased" = subtype_colors[["ZNF384-R"]],
  "Near zero"               = "grey70"
)

delta_df <- read.csv("./output/delta_r_correlations.csv")

# Wide df of two subtypes' delta r, one row per gene, with quadrant and TF flag
make_scatter <- function(delta_df, h1_gene, subtype_a, subtype_b) {
  delta_df %>%
    filter(H1_gene == h1_gene, Group %in% c(subtype_a, subtype_b)) %>%
    select(Gene, Group, delta_r) %>%
    pivot_wider(names_from = Group, values_from = delta_r) %>%
    rename(A = all_of(subtype_a), B = all_of(subtype_b)) %>%
    mutate(
      quadrant = case_when(
        A > 0 & B > 0 ~ "Both increased",
        A < 0 & B < 0 ~ "Both decreased",
        A > 0 & B < 0 ~ paste0(subtype_a, " only increased"),
        A < 0 & B > 0 ~ paste0(subtype_b, " only increased"),
        TRUE           ~ "Near zero"
      ),
      is_TF = Gene %in% jaspar_tfs
    )
}

# All genes labelled, TFs drawn bigger and boxed so they stand out
make_plot <- function(df, subtype_a, subtype_b, h1_gene, quadrant_colors) {
  r_val <- cor(df$A, df$B, use = "pairwise.complete.obs", method = "spearman")

  ggplot(df, aes(x = B, y = A, color = quadrant, shape = is_TF)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
    geom_point(data = filter(df, !is_TF), size = 2,   alpha = 0.75) +
    geom_point(data = filter(df,  is_TF), size = 3.5, alpha = 0.95) +
    geom_text_repel(
      data          = filter(df, !is_TF),
      aes(label     = Gene, color = quadrant),
      fontface      = "plain",
      size          = 2.4,
      max.overlaps  = Inf,
      segment.color = "grey75",
      segment.size  = 0.25,
      show.legend   = FALSE
    ) +
    geom_label_repel(
      data          = filter(df, is_TF),
      aes(label     = Gene, color = quadrant, fill = quadrant),
      fontface      = "bold.italic",
      size          = 3,
      label.size    = 0.25,
      label.r       = grid::unit(0.12, "lines"),
      box.padding   = 0.35,
      point.padding = 0.25,
      max.overlaps  = Inf,
      segment.color = "grey40",
      segment.size  = 0.35,
      show.legend   = FALSE
    ) +
    scale_color_manual(values = quadrant_colors) +
    scale_fill_manual(values = scales::alpha(quadrant_colors, 0.14), guide = "none") +
    scale_shape_manual(
      values = c("FALSE" = 16, "TRUE" = 18),
      labels = c("FALSE" = "Other gene", "TRUE" = "JASPAR2024 TF")
    ) +
    guides(
      color = guide_legend(title = "Quadrant",  nrow = 2, byrow = TRUE,
                           override.aes = list(size = 3)),
      shape = guide_legend(title = "Gene type", override.aes = list(size = c(2.5, 4)))
    ) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("Spearman r = ", round(r_val, 2)),
             size = 4, color = "grey30") +
    labs(
      title    = paste0("delta r: ", subtype_a, " vs ", subtype_b),
      subtitle = paste0(
        h1_gene, " module genes | ",
        "delta r = r(subtype) - mean(r(healthy))\n",
        "All genes labeled; TFs (JASPAR2024) shown as boxed diamond labels"
      ),
      x     = paste0(subtype_b, " delta r"),
      y     = paste0(subtype_a, " delta r"),
      color = "Quadrant",
      shape = "Gene type"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.box = "vertical")
}

# Scripts 11, 12 and ORA.R read these two RDS files
df_h10  <- make_scatter(delta_df, "H1-0",  "DUX4-R",       "ZNF384-R")
df_h110 <- make_scatter(delta_df, "H1-10", "Hyperdiploid", "TCF3::PBX1")
saveRDS(df_h10,  "./output/df_scatter_H1-0.rds")
saveRDS(df_h110, "./output/df_scatter_H1-10.rds")

# Only the H1-0 scatter is drawn this way, H1-10 is done per subtype further down
p1 <- make_plot(df_h10, "DUX4-R", "ZNF384-R", "H1-0", quadrant_colors_h10)

pdf("./output/scatter_H1modules_TF_highlighted.pdf", width = 8, height = 9)
print(p1)
dev.off()

# Rewiring summary statistics
# Set to 0 for any directional change, or 0.3 for stronger rewiring only.
REWIRING_SUMMARY_THRESH <- 0

summarise_module_rewiring <- function(scatter_df, h1_gene, subtype_a, subtype_b,
                                      threshold = 0) {
  long_df <- bind_rows(
    scatter_df %>% transmute(Gene, Subtype = subtype_a, delta_r = A),
    scatter_df %>% transmute(Gene, Subtype = subtype_b, delta_r = B)
  )

  subtype_summary <- long_df %>%
    group_by(Subtype) %>%
    summarise(
      H1_gene = h1_gene,
      n_genes = sum(!is.na(delta_r)),
      n_reduced = sum(delta_r < -threshold, na.rm = TRUE),
      pct_reduced = round(100 * n_reduced / n_genes, 1),
      n_gained = sum(delta_r > threshold, na.rm = TRUE),
      pct_gained = round(100 * n_gained / n_genes, 1),
      median_delta_r = round(median(delta_r, na.rm = TRUE), 3),
      mean_delta_r = round(mean(delta_r, na.rm = TRUE), 3),
      .groups = "drop"
    )

  n_pairwise_genes <- sum(!is.na(scatter_df$A) & !is.na(scatter_df$B))
  both_lost <- sum(scatter_df$A < -threshold & scatter_df$B < -threshold, na.rm = TRUE)
  both_gained <- sum(scatter_df$A > threshold & scatter_df$B > threshold, na.rm = TRUE)
  subtype_a_only_gained <- sum(scatter_df$A > threshold & scatter_df$B <= threshold, na.rm = TRUE)
  subtype_b_only_gained <- sum(scatter_df$B > threshold & scatter_df$A <= threshold, na.rm = TRUE)
  subtype_a_only_lost <- sum(scatter_df$A < -threshold & scatter_df$B >= -threshold, na.rm = TRUE)
  subtype_b_only_lost <- sum(scatter_df$B < -threshold & scatter_df$A >= -threshold, na.rm = TRUE)
  pct_pairwise <- function(n) round(100 * n / n_pairwise_genes, 1)

  quadrant_summary <- tibble(
    H1_gene = h1_gene,
    comparison = paste(subtype_a, "vs", subtype_b),
    n_pairwise_genes = n_pairwise_genes,
    both_lost = both_lost,
    both_lost_pct = pct_pairwise(both_lost),
    both_gained = both_gained,
    both_gained_pct = pct_pairwise(both_gained),
    subtype_a_only_gained = subtype_a_only_gained,
    subtype_a_only_gained_pct = pct_pairwise(subtype_a_only_gained),
    subtype_b_only_gained = subtype_b_only_gained,
    subtype_b_only_gained_pct = pct_pairwise(subtype_b_only_gained),
    subtype_a_only_lost = subtype_a_only_lost,
    subtype_a_only_lost_pct = pct_pairwise(subtype_a_only_lost),
    subtype_b_only_lost = subtype_b_only_lost,
    subtype_b_only_lost_pct = pct_pairwise(subtype_b_only_lost)
  )

  subtype_a_stats <- subtype_summary %>% filter(Subtype == subtype_a)
  subtype_b_stats <- subtype_summary %>% filter(Subtype == subtype_b)
  overall_median <- round(median(c(scatter_df$A, scatter_df$B), na.rm = TRUE), 3)

  sentences <- c(
    sprintf(
      "%s%% of %s module genes exhibited reduced coordination in %s; the median delta r was %s.",
      subtype_a_stats$pct_reduced, h1_gene, subtype_a, subtype_a_stats$median_delta_r
    ),
    sprintf(
      "%s%% of %s module genes exhibited reduced coordination in %s; the median delta r was %s.",
      subtype_b_stats$pct_reduced, h1_gene, subtype_b, subtype_b_stats$median_delta_r
    ),
    sprintf(
      "%s%% of %s module genes lost coordination in both %s and %s (%s/%s genes), whereas %s%% and %s%% gained coordination uniquely in %s and %s, respectively (%s and %s genes).",
      quadrant_summary$both_lost_pct, h1_gene, subtype_a, subtype_b,
      quadrant_summary$both_lost, quadrant_summary$n_pairwise_genes,
      quadrant_summary$subtype_a_only_gained_pct, quadrant_summary$subtype_b_only_gained_pct,
      subtype_a, subtype_b,
      quadrant_summary$subtype_a_only_gained, quadrant_summary$subtype_b_only_gained
    ),
    sprintf(
      "Across the %s module and both subtypes, the median delta r was %s.",
      h1_gene, overall_median
    )
  )

  list(
    subtype_summary = subtype_summary,
    quadrant_summary = quadrant_summary,
    sentences = sentences
  )
}

# Same three outputs per module, so the writing is wrapped
write_rewiring <- function(rewiring, tag) {
  print(rewiring$subtype_summary)
  print(rewiring$quadrant_summary)
  cat(paste(rewiring$sentences, collapse = "\n"), "\n")

  write.csv(rewiring$subtype_summary,
            paste0("./output/rewiring_summary_", tag, "_by_subtype.csv"), row.names = FALSE)
  write.csv(rewiring$quadrant_summary,
            paste0("./output/rewiring_summary_", tag, "_quadrants.csv"), row.names = FALSE)
  writeLines(rewiring$sentences,
             paste0("./output/rewiring_summary_", tag, "_sentences.txt"))
}

rewiring_h10 <- summarise_module_rewiring(
  df_h10, h1_gene = "H1-0", subtype_a = "DUX4-R", subtype_b = "ZNF384-R",
  threshold = REWIRING_SUMMARY_THRESH
)
write_rewiring(rewiring_h10, "H1-0_DUX4R_ZNF384R")

rewiring_h110 <- summarise_module_rewiring(
  df_h110, h1_gene = "H1-10", subtype_a = "Hyperdiploid", subtype_b = "TCF3::PBX1",
  threshold = REWIRING_SUMMARY_THRESH
)
write_rewiring(rewiring_h110, "H1-10_Hyperdiploid_TCF3PBX1")

# Replotting scatter for H1-10, one panel per subtype instead of one shared panel
# x = subtype r, y = healthy mean r, colored by whether delta r crosses the threshold
DELTA_THRESH <- 0.3

change_colors <- c(
  "Increased cor in disease" = "#6B3F1D",
  "Decreased cor in disease" = "#3F3F3F",
  "No change"                = "grey70"
)

h110_df <- delta_df %>%
  filter(H1_gene == "H1-10", Group %in% c("Hyperdiploid", "TCF3::PBX1")) %>%
  mutate(
    is_TF  = Gene %in% jaspar_tfs,
    change = case_when(
      delta_r >  DELTA_THRESH ~ "Increased cor in disease",
      delta_r < -DELTA_THRESH ~ "Decreased cor in disease",
      TRUE                    ~ "No change"
    )
  )

make_subtype_df <- function(df, subtype) {
  df %>%
    filter(Group == subtype) %>%
    select(Gene, r_disease, r_healthy, delta_r, is_TF, change)
}

make_subtype_plot <- function(df, subtype, h1_gene = "H1-10") {
  r_val <- cor(df$r_disease, df$r_healthy,
               use = "pairwise.complete.obs", method = "spearman")

  ggplot(df, aes(x = r_disease, y = r_healthy, color = change, shape = is_TF)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    geom_point(data = filter(df, !is_TF), size = 2,   alpha = 0.75) +
    geom_point(data = filter(df,  is_TF), size = 3.5, alpha = 0.95) +
    geom_text_repel(
      data          = filter(df, !is_TF),
      aes(label     = Gene, color = change),
      fontface      = "plain",
      size          = 2.4,
      max.overlaps  = Inf,
      segment.color = "grey75",
      segment.size  = 0.25,
      show.legend   = FALSE
    ) +
    geom_label_repel(
      data          = filter(df, is_TF),
      aes(label     = Gene, color = change, fill = change),
      fontface      = "bold.italic",
      size          = 3,
      label.size    = 0.25,
      label.r       = grid::unit(0.12, "lines"),
      box.padding   = 0.35,
      point.padding = 0.25,
      max.overlaps  = Inf,
      segment.color = "grey40",
      segment.size  = 0.35,
      show.legend   = FALSE
    ) +
    scale_color_manual(values = change_colors) +
    scale_fill_manual(values = scales::alpha(change_colors, 0.14), guide = "none") +
    scale_shape_manual(
      values = c("FALSE" = 16, "TRUE" = 18),
      labels = c("FALSE" = "Other gene", "TRUE" = "JASPAR2024 TF")
    ) +
    guides(
      color = guide_legend(
        title        = paste0("Change (|delta r| > ", DELTA_THRESH, ")"),
        override.aes = list(shape = 16, size = 3)
      ),
      shape = guide_legend(title = "Gene type", override.aes = list(size = c(2.5, 4)))
    ) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("Spearman r = ", round(r_val, 2)),
             size = 4, color = "grey30") +
    labs(
      title    = paste0(subtype, " - ", h1_gene, " module genes"),
      subtitle = paste0(
        "x = Spearman r (", subtype, " vs ", h1_gene, ") | y = mean Spearman r (healthy vs ", h1_gene, ")\n",
        "Dashed diagonal = identity (y = x) | ",
        "Color = change (delta r threshold +/-", DELTA_THRESH, ") | ",
        "TFs (JASPAR2024) shown as boxed diamond labels"
      ),
      x     = paste0(subtype, " r"),
      y     = "Healthy mean r",
      color = "Change",
      shape = "Gene type"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.box = "vertical")
}

p_hyper <- make_subtype_plot(make_subtype_df(h110_df, "Hyperdiploid"), "Hyperdiploid")
p_tcf3  <- make_subtype_plot(make_subtype_df(h110_df, "TCF3::PBX1"),   "TCF3::PBX1")

pdf("./output/scatter_H1-10_subtype_vs_healthy_r.pdf", width = 16, height = 9)
print(
  patchwork::wrap_plots(p_hyper, p_tcf3, ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
)
dev.off()

message("Saved: scatter_H1-10_subtype_vs_healthy_r.pdf (1 page, side by side)")

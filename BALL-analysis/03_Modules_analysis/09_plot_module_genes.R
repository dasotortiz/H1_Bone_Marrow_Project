# Spearman correlation of each module gene with its H1 gene, computed per healthy
# dataset and per B-ALL subtype across pseudobulks.
# Output: ./output/module_gene_correlations.csv, read by 10.

library(Seurat)
library(tidyverse)

target_genes <- list(
  "H1-10" = list(healthy = "H1-10",  disease = "H1FX"),
  "H1-0"  = list(healthy = "H1-0",   disease = "H1F0")
)

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

pbh  <- readRDS("./output/pseudobulk_hc_GSVA_modules.rds")
pbia <- readRDS("./output/pseudobulk_Ia_GSVA_modules.rds")

modules <- read.delim("./data/modules.tsv", header = TRUE, stringsAsFactors = FALSE)

module_genes <- list(
  "H1-10" = modules$id[modules$module == "21"],
  "H1-0"  = modules$id[modules$module == "28"]
)

# The disease object still uses the old HGNC symbols for these genes
std_to_old <- c(
  "H1-10"      = "H1FX",
  "SPACDR"     = "C7orf61",
  "CFAP251"    = "WDR66",
  "H3-3B"      = "H3F3B",
  "ZNF407-AS1" = "LINC00909",
  "H1-0"       = "H1F0",
  "H2AC25"     = "HIST3H2A",
  "STING1"     = "TMEM173",
  "LDAF1"      = "TMEM159"
)
old_to_std <- setNames(names(std_to_old), std_to_old)

# Fetch expression plus the grouping column, renaming old symbols back to standard
get_expr_with_group <- function(obj, genes_standard, group_col, remap = NULL) {
  fetch_names <- ifelse(genes_standard %in% names(remap),
                        remap[genes_standard],
                        genes_standard)

  present <- fetch_names[fetch_names %in% rownames(obj)]
  missing <- fetch_names[!fetch_names %in% rownames(obj)]
  if (length(missing) > 0) message("Missing: ", paste(missing, collapse = ", "))

  mat <- FetchData(obj, vars = c(present, group_col), layer = "data")

  if (!is.null(remap)) {
    colnames(mat) <- ifelse(colnames(mat) %in% names(old_to_std),
                            old_to_std[colnames(mat)],
                            colnames(mat))
  }
  mat
}

# Correlation of every module gene against the H1 gene, within one group
cor_within_groups <- function(mat, group_col, h1_col, mod_genes) {
  mat %>%
    group_by(.data[[group_col]]) %>%
    group_modify(~ tibble(
      Gene = mod_genes,
      Correlation = sapply(mod_genes, function(g) {
        if (!g %in% colnames(.x)) return(NA_real_)
        cor(.x[[h1_col]], .x[[g]], use = "pairwise.complete.obs", method = "spearman")
      })
    )) %>%
    ungroup() %>%
    rename(Group = all_of(group_col))
}

compute_correlations <- function(h1_name) {
  mod_genes  <- module_genes[[h1_name]]
  h1_healthy <- target_genes[[h1_name]]$healthy
  h1_disease <- target_genes[[h1_name]]$disease

  # Healthy: one correlation per Dataset
  h_mat <- get_expr_with_group(pbh, c(h1_healthy, mod_genes), "Dataset")
  cor_healthy <- cor_within_groups(h_mat, "Dataset", h1_healthy, mod_genes) %>%
    mutate(Status = "Healthy")

  # Disease: one per Subtype. The H1 gene is fetched under its old symbol
  # but comes back renamed, so index it by the standard name.
  d_mat <- get_expr_with_group(pbia, c(h1_disease, mod_genes), "Subtype", remap = std_to_old)
  cor_disease <- cor_within_groups(d_mat, "Subtype", h1_name, mod_genes) %>%
    mutate(Status = "Disease")

  bind_rows(cor_healthy, cor_disease) %>%
    mutate(H1_gene = h1_name) %>%
    filter(!is.na(Correlation))
}

cor_df <- bind_rows(
  compute_correlations("H1-10"),
  compute_correlations("H1-0")
)

write.csv(cor_df, "./output/module_gene_correlations.csv", row.names = FALSE)
message("Saved: module_gene_correlations.csv")

# Ranked dot plot, genes ordered by mean healthy correlation.
# subtype = NULL draws all subtypes at once, otherwise healthy vs that one subtype.
plot_correlation_dotplot <- function(h1_name, df, subtype = NULL) {
  sub_df <- df %>% filter(H1_gene == h1_name)
  if (!is.null(subtype)) {
    sub_df <- sub_df %>% filter(Status == "Healthy" | Group == subtype)
  }

  gene_order <- sub_df %>%
    filter(Status == "Healthy") %>%
    group_by(Gene) %>%
    summarise(mean_cor = mean(Correlation, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_cor) %>%
    pull(Gene)

  sub_df$Gene <- factor(sub_df$Gene, levels = gene_order)

  # Healthy datasets share one color, disease keeps its subtype color
  present_subtypes <- unique(sub_df$Group[sub_df$Status == "Disease"])
  plot_colors <- c("Healthy" = "#156415",
                   subtype_colors[names(subtype_colors) %in% present_subtypes])

  sub_df <- sub_df %>% mutate(color_group = if_else(Status == "Healthy", "Healthy", Group))

  halo_size  <- if (is.null(subtype)) 10 else 5
  legend_col <- if (is.null(subtype)) 6 else 2
  base_size  <- if (is.null(subtype)) 13 else 11
  subtitle   <- if (is.null(subtype)) "Spearman r across pseudobulks, per subtype / dataset"
                else paste("Healthy vs", subtype)

  ggplot(sub_df, aes(y = Correlation, x = Gene, color = color_group, shape = Status)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    # halo behind the healthy points so they read as a background band
    geom_jitter(data = filter(sub_df, Status == "Healthy"),
                width = 0.25, size = halo_size, alpha = 0.15, color = "#156415") +
    geom_jitter(width = 0.25, size = 2.5, alpha = 0.85) +
    scale_color_manual(values = plot_colors) +
    scale_shape_manual(values = c("Healthy" = 16, "Disease" = 17)) +
    guides(
      color = guide_legend(title = "Group",  ncol = legend_col, override.aes = list(size = 3)),
      shape = guide_legend(title = "Status", override.aes = list(size = 3))
    ) +
    labs(
      title    = paste("Module gene correlations with", h1_name),
      subtitle = subtitle,
      x        = NULL,
      y        = "Spearman r"
    ) +
    theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.x = element_line(color = "grey92"),
      legend.position    = "bottom",
      axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

pdf("./output/module_gene_correlations.pdf", width = 21, height = 9)
print(plot_correlation_dotplot("H1-10", cor_df))
print(plot_correlation_dotplot("H1-0",  cor_df))
dev.off()
message("Saved: module_gene_correlations.pdf")

# One PDF per H1 gene, one page per disease subtype
disease_subtypes <- unique(cor_df$Group[cor_df$Status == "Disease"])

for (h1_name in c("H1-10", "H1-0")) {
  fname <- paste0("./output/module_gene_correlations_", gsub("::", "-", h1_name), ".pdf")
  pdf(fname, width = 21, height = 7)
  for (st in disease_subtypes) {
    if (!any(cor_df$Group == st & cor_df$H1_gene == h1_name)) next
    print(plot_correlation_dotplot(h1_name, cor_df, subtype = st))
  }
  dev.off()
  message("Saved: ", fname)
}

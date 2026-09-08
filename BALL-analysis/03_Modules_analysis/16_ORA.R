# Over-representation analysis of the module gene sets defined by the delta r
# quadrants, against all genes detected in the healthy pseudobulk as background.
# Input: ./output/df_scatter_H1-0.rds and ./output/delta_r_correlations.csv

# Background was built once, outside this script:
# pb_hc_filt <- readRDS("../GEX-Plotting/output/pseudobulk_hc.rds")
# saveRDS(rownames(pb_hc_filt), "./data/hc_genes.rds")

# Run in Docker:
# docker start ORA_H1
# docker exec -it ORA_H1 /bin/bash

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(tidytable)
library(tidyverse)

run_ora_go <- function(genes, background, label, ont = "BP") {
  message("Running GO:", ont, " ORA for: ", label, " (n=", length(genes), ")")
  enrichGO(
    gene          = genes,
    universe      = background,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,   # keep everything, filtering happens at plotting
    qvalueCutoff  = 1,
    readable      = TRUE
  )
}

run_ora_hallmarks <- function(genes, background, label) {
  message("Running MSigDB Hallmarks ORA for: ", label, " (n=", length(genes), ")")
  h_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
    dplyr::select(gs_name, gene_symbol)
  enricher(
    gene          = genes,
    universe      = background,
    TERM2GENE     = h_sets,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    minGSSize     = 1,
    maxGSSize     = 500
  )
}

empty_plot <- function(label, text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = text, size = 5, color = "grey50") +
    theme_void() +
    labs(title = label)
}

plot_ora <- function(ego, label, fill_color, top_n = 15) {
  if (is.null(ego)) return(empty_plot(label, paste0("No terms found\n", label)))

  df <- ego@result %>%
    arrange(pvalue) %>%
    slice_head(n = top_n) %>%
    mutate(
      Count       = as.numeric(Count),
      Description = fct_reorder(Description, Count),
      gene_label  = str_replace_all(geneID, "/", ", ")
    )

  if (nrow(df) == 0) return(empty_plot(label, paste0("No terms found\n", label)))

  ggplot(df, aes(x = Count, y = Description, fill = pvalue)) +
    geom_col() +
    geom_text(aes(x = 0.1, label = gene_label),
              hjust = 0, size = 2.5, color = "white", fontface = "italic") +
    geom_text(aes(label = ifelse(p.adjust < 0.05, "*", "")),
              hjust = -0.3, size = 4, color = "black") +
    scale_fill_gradient(low = fill_color, high = "grey85", name = "raw p-value") +
    labs(
      title    = label,
      subtitle = paste0("Top ", top_n, " terms (ORA, ranked by raw p-value)\n* = adj. p < 0.05"),
      x        = "Gene count",
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")
}

# Run one analysis over a named list of gene sets, write the CSVs and one PDF.
# Each set carries its plotting color and the short tag used in filenames.
run_ora_panel <- function(sets, module, ont, suffix, title_tag) {
  plots <- lapply(names(sets), function(nm) {
    s <- sets[[nm]]
    label <- paste(module, "|", nm)
    ego <- if (ont == "Hallmarks") {
      run_ora_hallmarks(s$genes, background, label)
    } else {
      run_ora_go(s$genes, background, label, ont = ont)
    }
    write.csv(as.data.frame(ego),
              paste0("./output/ORA_", module, "_", s$tag, "_", suffix, "-bgHCgenes.csv"),
              row.names = FALSE)
    plot_ora(ego, paste0(label, " | ", title_tag), fill_color = s$color)
  })

  pdf(paste0("./output/ORA_", suffix, "_", module, "-bgHCgenes.pdf"), width = 10, height = 7)
  for (p in plots) print(p)
  dev.off()
}

background <- readRDS("./data/hc_genes.rds")
message("Background: n=", length(background)) # n=14170

# H1-0: gene sets come from the delta r quadrants
df_h10 <- readRDS("./output/df_scatter_H1-0.rds")

h10_sets <- list(
  "DUX4-R" = list(
    genes = filter(df_h10, quadrant %in% c("Both increased", "DUX4-R only increased"))$Gene,
    color = "#CC2222", tag = "DUX4-R"),
  "ZNF384-R" = list(
    genes = filter(df_h10, quadrant %in% c("Both increased", "ZNF384-R only increased"))$Gene,
    color = "#7B2D8B", tag = "ZNF384-R"),
  "Both decreased" = list(
    genes = filter(df_h10, quadrant == "Both decreased")$Gene,
    color = "#1F78B4", tag = "BLUE")
)

for (nm in names(h10_sets)) message(nm, " gene set: n=", length(h10_sets[[nm]]$genes))

run_ora_panel(h10_sets, "H1-0", ont = "BP",        suffix = "GOBP",      title_tag = "GO:BP")
run_ora_panel(h10_sets, "H1-0", ont = "MF",        suffix = "GOMF",      title_tag = "GO:MF")
run_ora_panel(h10_sets, "H1-0", ont = "Hallmarks", suffix = "Hallmarks", title_tag = "Hallmarks")
message("Saved: H1-0 ORA results to ./output/")

# H1-10: gene sets come from the sign of delta r, per subtype
DELTA_THRESH <- 0.3

h110_df <- read.csv("./output/delta_r_correlations.csv") %>%
  filter(H1_gene == "H1-10", Group %in% c("Hyperdiploid", "TCF3::PBX1")) %>%
  mutate(
    change = case_when(
      delta_r >  DELTA_THRESH ~ "Increased cor in disease",
      delta_r < -DELTA_THRESH ~ "Decreased cor in disease",
      TRUE                    ~ "No change"
    )
  )

pick <- function(group, change) filter(h110_df, Group == group, change == !!change)$Gene

h110_sets <- list(
  "Hyperdiploid - Decreased" = list(genes = pick("Hyperdiploid", "Decreased cor in disease"),
                                    color = "#1F78B4", tag = "hypDOWN"),
  "Hyperdiploid - Increased" = list(genes = pick("Hyperdiploid", "Increased cor in disease"),
                                    color = "#CC2222", tag = "hypUP"),
  "TCF3::PBX1 - Decreased"   = list(genes = pick("TCF3::PBX1", "Decreased cor in disease"),
                                    color = "#1F78B4", tag = "tcfDOWN"),
  "TCF3::PBX1 - Increased"   = list(genes = pick("TCF3::PBX1", "Increased cor in disease"),
                                    color = "#CC2222", tag = "tcfUP")
)

for (nm in names(h110_sets)) message(nm, ": n=", length(h110_sets[[nm]]$genes))
# Hyperdiploid DOWN: n=127, UP: n=7 | TCF3::PBX1 DOWN: n=127, UP: n=9

run_ora_panel(h110_sets, "H1-10", ont = "BP",        suffix = "GOBP",      title_tag = "GO:BP")
run_ora_panel(h110_sets, "H1-10", ont = "Hallmarks", suffix = "Hallmarks", title_tag = "Hallmark")
message("Saved: H1-10 ORA results to ./output/")

# Replot the saved H1-10 GO:BP tables without rerunning ORA, keeping only
# significant terms and coloring by adjusted p-value.
DECREASED_COR_COLOR <- "#3F3F3F"

plot_sig_ora_from_csv <- function(file, label, fill_color = DECREASED_COR_COLOR,
                                  top_n = 15, padj_cutoff = 0.05, show_genes = TRUE) {
  if (!file.exists(file)) return(empty_plot(label, paste0("Missing ORA table\n", file)))

  df <- read.csv(file, check.names = FALSE)
  if (nrow(df) == 0 || !"p.adjust" %in% colnames(df)) {
    return(empty_plot(label, paste0("No ORA results found\n", label)))
  }

  df <- df %>%
    filter(!is.na(p.adjust), p.adjust < padj_cutoff) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      Count = as.numeric(Count),
      Description = factor(Description, levels = rev(Description)),
      gene_label = if ("geneID" %in% colnames(.)) str_replace_all(geneID, "/", ", ") else ""
    )

  if (nrow(df) == 0) {
    return(empty_plot(label, paste0("No significant terms at FDR < ", padj_cutoff, "\n", label)))
  }

  p <- ggplot(df, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_col() +
    scale_fill_gradient(low = fill_color, high = "grey85",
                        limits = c(0, padj_cutoff), name = "adj. p-value") +
    labs(
      title    = label,
      subtitle = "Significant GO:BP terms only (ranked by adjusted p-value)",
      x        = "Gene count",
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")

  if (show_genes) {
    p <- p + geom_text(aes(x = 0.1, label = gene_label),
                       hjust = 0, size = 2.5, color = "#FFD60A", fontface = "italic")
  }
  p
}

sig_files <- c(
  "H1-10 | Hyperdiploid - Decreased | GO:BP" = "./output/ORA_H1-10_hypDOWN_GOBP-bgHCgenes.csv",
  "H1-10 | TCF3::PBX1 - Decreased | GO:BP"   = "./output/ORA_H1-10_tcfDOWN_GOBP-bgHCgenes.csv"
)

# Same figure twice, with and without the gene names printed inside the bars
for (show_genes in c(TRUE, FALSE)) {
  fname <- if (show_genes) {
    "./output/ORA_GOBP_H1-10-bgHCgenes_padj_sig_decreased.pdf"
  } else {
    "./output/ORA_GOBP_H1-10-bgHCgenes_padj_sig_decreased_noGeneLabels.pdf"
  }

  panels <- lapply(names(sig_files), function(lbl) {
    plot_sig_ora_from_csv(sig_files[[lbl]], lbl, show_genes = show_genes)
  })

  pdf(fname, width = 20, height = 7)
  print(
    patchwork::wrap_plots(panels, ncol = 2, guides = "collect") &
      theme(legend.position = "bottom")
  )
  dev.off()
  message("Saved: ", fname)
}

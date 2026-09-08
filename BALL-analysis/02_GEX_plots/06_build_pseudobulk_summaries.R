# Summarise H1 gene expression from the pseudobulk objects into the mean + SD
# tables that the heatmap and lineplot scripts start from.
# Output: ./output/h_df.rds (healthy, grouped by dataset)
#         ./output/i_df.rds (Iacobucci, grouped by subtype)
#
# pseudobulk_hc.rds and pseudobulk_Ia.rds are inputs, not built here.

library(Seurat)
library(tidyverse)

pb_hc <- readRDS("./output/pseudobulk_hc.rds")
pb_ia <- readRDS("./output/pseudobulk_Ia.rds")

# h1_df.rds lives outside this repo, set the path to wherever you keep it
h1_path <- "../h1_df.rds"
h1_df <- readRDS(h1_path)

lineage_cells <- c("HSC", "MPP", "CLP", "ProB", "PreB")

# The healthy object carries the standard symbols, the Iacobucci object still
# uses the alternative ones, so the disease side gets renamed back after fetching.
h1_std <- h1_df$Symbol
h1_alt <- h1_df$Alternative_symbol

summarise_h1 <- function(obj, genes, group_col, status_label) {
  available <- intersect(genes, rownames(obj))
  missing   <- setdiff(genes, rownames(obj))
  if (length(missing) > 0)
    message("Not found in ", status_label, " object: ", paste(missing, collapse = ", "))

  FetchData(obj, vars = c(available, "Celltype", group_col), layer = "data") %>%
    pivot_longer(cols = all_of(available), names_to = "Gene", values_to = "Expression") %>%
    rename(Subtype = all_of(group_col)) %>%
    group_by(Subtype, Celltype, Gene) %>%
    summarise(Mean = mean(Expression, na.rm = TRUE),
              SD   = sd(Expression,   na.rm = TRUE), .groups = "drop") %>%
    mutate(Status = status_label)
}

h_df <- summarise_h1(pb_hc, h1_std, "Dataset", "Healthy")

i_df <- summarise_h1(pb_ia, h1_alt, "Subtype", "Disease") %>%
  mutate(Gene = h1_std[match(Gene, h1_alt)])  # back to the standard symbols

# Keep the lymphoid lineage only, in developmental order
h_df <- h_df %>% filter(Celltype %in% lineage_cells) %>%
  mutate(Celltype = factor(Celltype, levels = lineage_cells))
i_df <- i_df %>% filter(Celltype %in% lineage_cells) %>%
  mutate(Celltype = factor(Celltype, levels = lineage_cells))

saveRDS(h_df, "./output/h_df.rds")
saveRDS(i_df, "./output/i_df.rds")

message("Saved: h_df.rds (", nrow(h_df), " rows) and i_df.rds (", nrow(i_df), " rows)")

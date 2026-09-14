# ==============================================================================
# 0. SETUP & LIBRARIES
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(JASPAR2024)
  library(TFBSTools)
  library(BSgenome.Hsapiens.NCBI.GRCh38)
  library(motifmatchr)
  library(ensembldb)
  library(readr)
})

# ==============================================================================
# 1. PATHS
# ==============================================================================
seurat_path     <- "/home/x_vazquem/workspace/h1_project/panel_3_C-E/data/Xpand_atac.rds"
base_dir        <- "/home/x_vazquem/workspace/h1_project/panel_3_C-E/data/atac/xpand_20260223"
samples         <- c("MO291", "MO292", "MO294", "MO296", "MO298")
outdir          <- "/home/x_vazquem/workspace/h1_project/panel_3_C-E/results/"

# Other preprocessed datasets + TF analysis paths
seurat_ainci_path <- "/home/x_vazquem/workspace/h1_project/panel_3_C-E/results/preprocessed_seurats/seurat_ainci_final.rds"
seurat_li_path    <- "/home/x_vazquem/workspace/h1_project/panel_3_C-E/results/preprocessed_seurats/seurat_li_final.rds"
seurat_hca_path   <- "/home/x_vazquem/workspace/h1_project/panel_3_C-E/results/preprocessed_seurats/seurat_hca_final.rds"
tfs_rds_path      <- "/home/x_vazquem/workspace/h1_project/panel_3_C-E/data/JASPAR2024_TFs.rds"

# TF~H1-10 correlation thresholds
tf_corr_threshold  <- 0.50
tf_fdr_threshold   <- 0.05   # significance is called on the BH-adjusted p-value
tf_min_datasets    <- 2

# ==============================================================================
# 2. GLOBAL VARIABLES & PALETTES
# ==============================================================================
colors_celltypes <- c(
  HSC         = "#00441B", MPP         = "#00AF99",
  CLP         = "#98D9E9", ProB        = "#0081C9",
  PreB        = "#001588", MEP         = "#F6313E",
  ERP         = "#8F1336", MKP         = "#46A040",
  GMP         = "#FFC179", `Eo/B/Mast` = "#FFA300",
  MDP         = "#AFAFAF"
)

gene_name_map <- c(
  "H1F0"     = "H1-0", "HIST1H1A" = "H1-1",
  "HIST1H1C" = "H1-2", "HIST1H1D" = "H1-3",
  "HIST1H1E" = "H1-4", "HIST1H1B" = "H1-5",
  "H1FX"     = "H1-10"
)

target_gene <- "H1-10"

# Cell types shown in the coverage track (full lymphoid trajectory)
lymphocyte_lineage <- c("HSC", "MPP", "CLP", "ProB", "PreB")

# Cell types used for the TF~H1-10 expression correlations
correlation_celltypes <- c("CLP", "ProB", "PreB")

# Pseudobulk QC settings for the TF~H1-10 correlations
min_cells_per_sample     <- 30
min_samples_per_celltype <- 3

# H1 gene-name aliases across datasets -> harmonized name
gene_lookup <- data.frame(
  Original = c(
    "H1F0",     "HIST1H1A", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H1B", "H1FX",
    "H1-0",     "H1-1",     "H1-2",     "H1-3",     "H1-4",     "H1-5",     "H1-10",
    "H1.0",     "H1.1",     "H1.2",     "H1.3",     "H1.4",     "H1.5",     "H1.10"
  ),
  Final    = c(
    "H1-0",     "H1-1",     "H1-2",     "H1-3",     "H1-4",     "H1-5",     "H1-10",
    "H1-0",     "H1-1",     "H1-2",     "H1-3",     "H1-4",     "H1-5",     "H1-10",
    "H1-0",     "H1-1",     "H1-2",     "H1-3",     "H1-4",     "H1-5",     "H1-10"
  )
)

# ==============================================================================
# 3. LOAD SEURAT OBJECT
# ==============================================================================
xpand_atac <- readRDS(seurat_path)

# ==============================================================================
# 4. RNA PRE-PROCESSING (only if RNA assay present and PCA not yet computed)
# ==============================================================================
if ("RNA" %in% Assays(xpand_atac)) {
  DefaultAssay(xpand_atac) <- "RNA"

  if (!inherits(xpand_atac[["RNA"]], "Assay5")) {
    xpand_atac[["RNA"]] <- as(xpand_atac[["RNA"]], Class = "Assay5")
  }

  if (!"pca" %in% Reductions(xpand_atac)) {
    message("Running RNA preprocessing...")
    xpand_atac <- NormalizeData(xpand_atac)
    xpand_atac <- FindVariableFeatures(xpand_atac)
    xpand_atac <- ScaleData(xpand_atac)
    xpand_atac <- RunPCA(xpand_atac)
    xpand_atac <- RunUMAP(xpand_atac, dims = 1:30,
                          reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
  } else {
    message("PCA already found — skipping RNA preprocessing.")
  }
}

# ==============================================================================
# 5. ATAC PRE-PROCESSING
# ==============================================================================
DefaultAssay(xpand_atac) <- "ATAC"
Idents(xpand_atac) <- "annotation_2"

# Exclude first LSI dimension (correlated with sequencing depth)
xpand_atac <- RunTFIDF(xpand_atac)
xpand_atac <- FindTopFeatures(xpand_atac, min.cutoff = "q0")
xpand_atac <- RunSVD(xpand_atac)
xpand_atac <- RunUMAP(xpand_atac, reduction = "lsi", dims = 2:30,
                      reduction.name = "umap.atac", reduction.key = "atacUMAP_")
xpand_atac <- FindNeighbors(xpand_atac, reduction = "lsi", dims = 2:30)
xpand_atac <- FindClusters(xpand_atac, verbose = FALSE, algorithm = 3)

DefaultAssay(xpand_atac) <- "ATAC"
Idents(xpand_atac) <- "annotation_2"

# ==============================================================================
# 6. FRAGMENT RE-BINDING
# ==============================================================================
Fragments(xpand_atac) <- NULL

# Re-bind fragments by matching sample prefixes directly from cell names
frag_list <- lapply(samples, function(s) {
  frag_path <- file.path(base_dir, s, "atac_fragments.tsv.gz")

  # Pattern to match leading prefix (e.g., "_MO291_")
  prefix_pattern <- paste0("^_", s, "_")

  # Extract cell names belonging to sample 's'
  seurat_cells <- grep(prefix_pattern, colnames(xpand_atac), value = TRUE)

  # Strip "_MO291_" prefix to match raw fragment barcodes (e.g., "AAACAGCCACCACAAC-1")
  raw_barcodes <- sub(prefix_pattern, "", seurat_cells)

  # Map Seurat cell names (vector names) to fragment file barcodes (vector values)
  cells_mapped <- setNames(raw_barcodes, seurat_cells)

  CreateFragmentObject(
    path = frag_path,
    cells = cells_mapped
  )
})

# Assign updated fragment list to object
Fragments(xpand_atac) <- frag_list

# ==============================================================================
# 7. H1 PROMOTER PEAKS (ATAC)
# ==============================================================================
# Defines the peak set used downstream for JASPAR motif matching. Computed on
# the full object: peak calling is cell-type agnostic and motif matching is
# purely sequence-based.
DefaultAssay(xpand_atac) <- "ATAC"
gene_coords <- Annotation(xpand_atac)
h1_coords <- gene_coords[gene_coords$gene_name %in% names(gene_name_map)]
h1_promoters <- promoters(h1_coords, upstream = 2000, downstream = 500)

peak_coords <- StringToGRanges(rownames(xpand_atac), sep = c("-", "-"))
overlaps <- findOverlaps(query = h1_promoters, subject = peak_coords)

clean_map <- data.frame(
  Gene_Original = h1_coords$gene_name[queryHits(overlaps)],
  Peak_ID = rownames(xpand_atac)[subjectHits(overlaps)]
) %>% distinct() %>% mutate(Gene_Final = gene_name_map[Gene_Original])

# ==============================================================================
# 8. COVERAGE PLOT — H1-10 (H1FX), full lymphoid trajectory
# ==============================================================================
# Shows the whole lymphoid trajectory (HSC -> PreB) for context; only the
# TF~H1-10 correlations below are restricted to CLP/ProB/PreB.

seurat_atac_sub <- subset(xpand_atac, subset = annotation_2 %in% lymphocyte_lineage)
seurat_atac_sub$annotation_2 <- factor(seurat_atac_sub$annotation_2, levels = lymphocyte_lineage)
Idents(seurat_atac_sub) <- "annotation_2"

p <- CoveragePlot(
  object = seurat_atac_sub,
  region = "H1FX",
  group.by = "annotation_2",
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000
) & scale_fill_manual(values = colors_celltypes)

ggsave(
  filename = paste0(outdir, "coverage_H1-10_lymphocytes.pdf"),
  plot = p,
  width = 8, height = 6, device = "pdf"
)

# ==============================================================================
# 9. CANDIDATE TFs BOUND TO THE H1-10 PROMOTER (XPAND, JASPAR MOTIF MATCH)
# ==============================================================================
message("--> Loading target TF list and matching JASPAR motifs on Xpand ATAC...")

tfs_rds  <- readRDS(tfs_rds_path)
all_tfs  <- if (is.data.frame(tfs_rds)) {
  gene_col <- intersect(c("Gene", "TF", "gene", "tf"), colnames(tfs_rds))[1]
  unique(tfs_rds[[gene_col]])
} else {
  unique(as.character(tfs_rds))
}

DefaultAssay(xpand_atac) <- "ATAC"

db_connection <- JASPAR2024()@db
pfm <- getMatrixSet(
  x = db_connection,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)

if (is.null(Motifs(xpand_atac[["ATAC"]]))) {
  xpand_atac <- AddMotifs(object = xpand_atac, genome = BSgenome.Hsapiens.NCBI.GRCh38, pfm = pfm)
}

motif_names_map <- xpand_atac[["ATAC"]]@motifs@motif.names
tf_to_id_map    <- setNames(names(motif_names_map), unlist(motif_names_map))

valid_tfs        <- intersect(all_tfs, names(tf_to_id_map))
target_motif_ids <- tf_to_id_map[valid_tfs]

h1_peaks <- unique(clean_map$Peak_ID)

motif_matrix    <- xpand_atac[["ATAC"]]@motifs@data
h1_motif_subset <- motif_matrix[h1_peaks, target_motif_ids, drop = FALSE]

h1_motif_df <- as.data.frame(as.matrix(h1_motif_subset))
colnames(h1_motif_df) <- names(target_motif_ids)
h1_motif_df <- h1_motif_df %>% rownames_to_column(var = "Peak_ID")

tf_h1_bindings <- h1_motif_df %>%
  pivot_longer(cols = -Peak_ID, names_to = "TF", values_to = "Motif_Present") %>%
  filter(Motif_Present == 1) %>%
  left_join(clean_map, by = "Peak_ID") %>%
  select(Gene_Final, Peak_ID, TF) %>%
  distinct()

tf_h1_10_bindings <- tf_h1_bindings %>%
  filter(Gene_Final == target_gene)

candidate_tfs_h1_10 <- unique(tf_h1_10_bindings$TF)
message(sprintf("Identified %d unique TFs from your list bound to open H1-10 promoter peaks in Xpand.", length(candidate_tfs_h1_10)))

write.table(
  candidate_tfs_h1_10,
  file = file.path(outdir, "candidate_tfs_H1_10_open_regions.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

# ==============================================================================
# 10. EXACT TF MOTIF COORDINATES WITHIN H1-10 PEAKS
# ==============================================================================
message("--> Finding exact genomic coordinates for TFs within H1-10 peaks...")

h1_10_peak_ids <- unique(tf_h1_10_bindings$Peak_ID)
h1_10_peaks_gr <- StringToGRanges(h1_10_peak_ids, sep = c("-", "-"))

candidate_motif_ids <- target_motif_ids[candidate_tfs_h1_10]
h1_10_pfms          <- pfm[candidate_motif_ids]

motif_positions <- matchMotifs(
  pwms    = h1_10_pfms,
  subject = h1_10_peaks_gr,
  genome  = BSgenome.Hsapiens.NCBI.GRCh38,
  out     = "positions"
)

tf_exact_locations <- as.data.frame(motif_positions) %>%
  rename(Motif_ID = group_name) %>%
  mutate(TF = names(target_motif_ids)[match(Motif_ID, target_motif_ids)]) %>%
  select(TF, Motif_ID, seqnames, start, end, strand, score) %>%
  arrange(seqnames, start)

message(sprintf("Identified %d exact motif occurrences within the H1-10 ATAC peak boundaries.", nrow(tf_exact_locations)))

write_csv(tf_exact_locations, file.path(outdir, "H1_10_exact_TF_motif_peak_positions.csv"))

# ==============================================================================
# 11. RNA PSEUDOBULK TF~H1-10 CORRELATION FUNCTION (per dataset)
# ==============================================================================
# Pseudobulk uses AggregateExpression (sum counts, then normalize), grouped by
# annotation_2 x sample, restricted to CLP/ProB/PreB and QC-filtered
# (min_cells_per_sample / min_samples_per_celltype).
#
# Significance is called on the BH-adjusted p-value (FDR < tf_fdr_threshold),
# adjusted across all candidate TFs tested within that dataset, rather than on
# the raw p-value.

compute_dataset_correlations <- function(seurat_obj, dataset_name, sample_col) {
  message(sprintf("--> Processing RNA pseudobulk TF~H1-10 correlations for %s (CLP/ProB/PreB)...", dataset_name))

  if ("RNA" %in% Assays(seurat_obj)) {
    DefaultAssay(seurat_obj) <- "RNA"
  } else if ("SCT" %in% Assays(seurat_obj)) {
    DefaultAssay(seurat_obj) <- "SCT"
  }
  rna_assay_ds <- DefaultAssay(seurat_obj)

  available_rna_genes <- rownames(seurat_obj)
  valid_rna_tfs        <- intersect(candidate_tfs_h1_10, available_rna_genes)

  h1_10_original <- gene_lookup %>%
    filter(Final == target_gene, Original %in% available_rna_genes) %>%
    pull(Original) %>%
    unique()

  if (length(h1_10_original) == 0 || length(valid_rna_tfs) == 0) {
    warning(sprintf("Skipping %s: H1-10 gene or candidate TFs not found in its RNA assay.", dataset_name))
    return(list(stats = tibble(), pseudobulk = tibble()))
  }
  h1_10_original <- h1_10_original[1]

  # QC filter: min cells/sample, min samples/celltype, CLP/ProB/PreB only
  raw_meta <- FetchData(seurat_obj, vars = c(sample_col, "annotation_2"))

  valid_pairs <- raw_meta %>%
    filter(annotation_2 %in% correlation_celltypes) %>%
    count(!!sym(sample_col), annotation_2, name = "n_cells") %>%
    filter(n_cells >= min_cells_per_sample)

  valid_celltypes <- valid_pairs %>%
    count(annotation_2, name = "n_samples") %>%
    filter(n_samples >= min_samples_per_celltype) %>%
    pull(annotation_2)

  valid_pairs <- valid_pairs %>% filter(annotation_2 %in% valid_celltypes)

  seurat_obj$keep <- paste(seurat_obj@meta.data[[sample_col]], seurat_obj$annotation_2, sep = "__") %in%
    paste(valid_pairs[[sample_col]], valid_pairs$annotation_2, sep = "__")
  seurat_sub <- subset(seurat_obj, subset = keep == TRUE)

  if (ncol(seurat_sub) == 0) {
    warning(sprintf("Skipping %s: no cells passed the QC filter.", dataset_name))
    return(list(stats = tibble(), pseudobulk = tibble()))
  }

  agg_obj_ds <- AggregateExpression(
    seurat_sub,
    assays   = rna_assay_ds,
    group.by = c("annotation_2", sample_col),
    return.seurat = TRUE
  )

  pb_expr <- FetchData(agg_obj_ds, vars = unique(c(h1_10_original, valid_rna_tfs))) %>%
    rownames_to_column("pb_sample")

  pb_meta_ds <- agg_obj_ds@meta.data %>%
    rownames_to_column("pb_sample") %>%
    select(pb_sample, Cell_Type = annotation_2, Sample_ID = !!sym(sample_col))

  h1_pb <- pb_expr %>%
    select(pb_sample, pb_h1_expr = all_of(h1_10_original)) %>%
    mutate(Gene_Final = target_gene)

  tf_pb <- pb_expr %>%
    select(pb_sample, all_of(valid_rna_tfs)) %>%
    pivot_longer(cols = -pb_sample, names_to = "TF", values_to = "pb_tf_expr")

  corr_input <- tf_pb %>%
    inner_join(h1_pb, by = "pb_sample") %>%
    inner_join(pb_meta_ds, by = "pb_sample")

  tf_h1_stats <- corr_input %>%
    group_by(TF, Gene_Final) %>%
    summarise(
      n_points = n(),
      Rho      = cor(pb_tf_expr, pb_h1_expr, method = "spearman", use = "complete.obs"),
      P_val    = tryCatch(cor.test(pb_tf_expr, pb_h1_expr, method = "spearman", exact = FALSE)$p.value, error = function(e) NA_real_),
      .groups  = "drop"
    ) %>%
    filter(!is.na(Rho) & !is.na(P_val)) %>%
    mutate(FDR = p.adjust(P_val, method = "BH")) %>%
    filter(abs(Rho) > tf_corr_threshold & FDR < tf_fdr_threshold) %>%
    arrange(desc(Rho))

  message(sprintf("    %s: %d TF~H1-10 pairs pass |Rho| > %.2f & FDR < %.2f.",
                  dataset_name, nrow(tf_h1_stats), tf_corr_threshold, tf_fdr_threshold))

  list(stats = tf_h1_stats, pseudobulk = corr_input)
}

# ==============================================================================
# 12. RUN TF~H1-10 CORRELATIONS ACROSS ALL 4 DATASETS
# ==============================================================================
# Xpand reuses the in-memory preprocessed object (xpand_atac) instead of
# re-reading it from disk.
xpand_result        <- compute_dataset_correlations(xpand_atac, "Xpand", "sample")
res_xpand           <- xpand_result$stats
xpand_tf_pseudobulk <- xpand_result$pseudobulk

seurat_ainci <- readRDS(seurat_ainci_path)
ainci_result <- compute_dataset_correlations(seurat_ainci, "Ainci", "sample")
res_ainci    <- ainci_result$stats
rm(seurat_ainci); gc()

seurat_li <- readRDS(seurat_li_path)
li_result <- compute_dataset_correlations(seurat_li, "Li", "sample")
res_li    <- li_result$stats
rm(seurat_li); gc()

seurat_hca <- readRDS(seurat_hca_path)
hca_result <- compute_dataset_correlations(seurat_hca, "HCA", "orig.ident")
res_hca    <- hca_result$stats
rm(seurat_hca); gc()

# ==============================================================================
# 13. TFs SIGNIFICANT (|Rho| > tf_corr_threshold, FDR < tf_fdr_threshold) IN >= tf_min_datasets DATASETS
# ==============================================================================
message(sprintf("--> Finding TF~H1-10 pairs significant (|Rho| > %.2f, FDR < %.2f) in at least %d of 4 datasets...",
                tf_corr_threshold, tf_fdr_threshold, tf_min_datasets))

pairs_xpand <- res_xpand %>% distinct(TF) %>% mutate(in_xpand = TRUE)
pairs_ainci <- res_ainci %>% distinct(TF) %>% mutate(in_ainci = TRUE)
pairs_li    <- res_li    %>% distinct(TF) %>% mutate(in_li = TRUE)
pairs_hca   <- res_hca   %>% distinct(TF) %>% mutate(in_hca = TRUE)

all_pairs <- pairs_xpand %>%
  full_join(pairs_ainci, by = "TF") %>%
  full_join(pairs_li, by = "TF") %>%
  full_join(pairs_hca, by = "TF") %>%
  replace_na(list(in_xpand = FALSE, in_ainci = FALSE, in_li = FALSE, in_hca = FALSE)) %>%
  mutate(n_datasets = in_xpand + in_ainci + in_li + in_hca)

final_results_min2 <- all_pairs %>%
  filter(n_datasets >= tf_min_datasets) %>%
  left_join(res_xpand %>% rename(Rho_Xpand = Rho, P_val_Xpand = P_val, FDR_Xpand = FDR) %>% select(TF, Rho_Xpand, P_val_Xpand, FDR_Xpand), by = "TF") %>%
  left_join(res_ainci %>% rename(Rho_Ainci = Rho, P_val_Ainci = P_val, FDR_Ainci = FDR) %>% select(TF, Rho_Ainci, P_val_Ainci, FDR_Ainci), by = "TF") %>%
  left_join(res_li    %>% rename(Rho_Li = Rho, P_val_Li = P_val, FDR_Li = FDR)          %>% select(TF, Rho_Li, P_val_Li, FDR_Li), by = "TF") %>%
  left_join(res_hca   %>% rename(Rho_HCA = Rho, P_val_HCA = P_val, FDR_HCA = FDR)       %>% select(TF, Rho_HCA, P_val_HCA, FDR_HCA), by = "TF") %>%
  mutate(Gene_Final = target_gene) %>%
  arrange(desc(n_datasets), TF)

message("--- Summary: TFs significant (FDR < ", tf_fdr_threshold, ") in >= ", tf_min_datasets, " datasets ---")
print(as_tibble(final_results_min2))
write_csv(final_results_min2, file.path(outdir, "common_tf_h1_10_correlations_min2_CLP_ProB_PreB.csv"))

# ==============================================================================
# 14. INDIVIDUAL DATASET CSV EXPORTS
# ==============================================================================
message("Xpand TFs: ", length(unique(res_xpand$TF)))
message("Ainci TFs: ", length(unique(res_ainci$TF)))
message("Li TFs:    ", length(unique(res_li$TF)))
message("HCA TFs:   ", length(unique(res_hca$TF)))

write_csv(res_xpand, file.path(outdir, "tf_h1_10_correlations_xpand_CLP_ProB_PreB.csv"))
write_csv(res_ainci, file.path(outdir, "tf_h1_10_correlations_ainci_CLP_ProB_PreB.csv"))
write_csv(res_li, file.path(outdir, "tf_h1_10_correlations_li_CLP_ProB_PreB.csv"))
write_csv(res_hca, file.path(outdir, "tf_h1_10_correlations_hca_CLP_ProB_PreB.csv"))

# ==============================================================================
# 15. CORRELATION PLOTS — PSEUDOBULKED H1-10 vs EACH FINAL TF (XPAND, CLP/ProB/PreB)
# ==============================================================================
message("--> Generating correlation plots for H1-10 vs final candidate TFs (Xpand, CLP/ProB/PreB)...")

final_tf_list <- unique(final_results_min2$TF)

if (length(final_tf_list) == 0) {
  message("No TFs passed the |Rho| / FDR filter in >= ", tf_min_datasets,
          " datasets — skipping the TF correlation panel.")
} else {
  plot_data_final_tfs <- xpand_tf_pseudobulk %>%
    filter(TF %in% final_tf_list) %>%
    mutate(Cell_Type = factor(Cell_Type, levels = correlation_celltypes))

  stats_final_tfs <- plot_data_final_tfs %>%
    group_by(TF) %>%
    summarise(
      Rho   = cor(pb_tf_expr, pb_h1_expr, method = "spearman", use = "complete.obs"),
      .groups = "drop"
    ) %>%
    mutate(stat_label = paste0("Rho = ", round(Rho, 2)))

  p_tf_h1_10_corr <- ggplot(plot_data_final_tfs, aes(x = pb_h1_expr, y = pb_tf_expr)) +
    geom_point(aes(color = Cell_Type), size = 3) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    geom_text(
      data = stats_final_tfs,
      aes(x = -Inf, y = Inf, label = stat_label),
      hjust = -0.1, vjust = 1.3, size = 5, inherit.aes = FALSE
    ) +
    facet_wrap(~ TF, scales = "free") +
    scale_color_manual(values = colors_celltypes, breaks = correlation_celltypes) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      strip.text = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold")
    ) +
    labs(
      x = "H1-10 expression",
      y = "TF expression",
      color = ""
    )

  ggsave(
    filename = file.path(outdir, "tf_h1_10_correlations_xpand_CLP_ProB_PreB.pdf"),
    plot = p_tf_h1_10_corr,
    width = 10, height = 5, device = "pdf"
  )
}

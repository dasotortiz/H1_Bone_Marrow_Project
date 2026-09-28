# ==============================================================================
# 02. ANALYSIS - TFs binding the H1-10 promoter and TF~H1-10 correlations
# ==============================================================================
# Input : xpand_atac_preprocessed_path (from 01_preprocessing.R),
#         Ainci / Li / HCA preprocessed Seurat objects, TF lists
# Output: candidate TF list, exact motif positions, per-dataset correlation
#         tables, TF summary table, Xpand pseudobulk table (for 03_figures.R)

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(JASPAR2024)
  library(TFBSTools)
  library(BSgenome.Hsapiens.NCBI.GRCh38)
  library(motifmatchr)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(readr)
})

# ==============================================================================
# 0. CONFIG
# ==============================================================================
project_dir <- "/home/x_vazquem/workspace/h1_project/panel_3_C-E"
outdir      <- file.path(project_dir, "results")

# ---- Inputs ------------------------------------------------------------------
xpand_atac_preprocessed_path <- file.path(outdir, "preprocessed_seurats/Xpand_atac_preprocessed.rds")
seurat_ainci_path <- file.path(outdir, "preprocessed_seurats/seurat_ainci_final.rds")
seurat_li_path    <- file.path(outdir, "preprocessed_seurats/seurat_li_final.rds")
seurat_hca_path   <- file.path(outdir, "preprocessed_seurats/seurat_hca_final.rds")
tfs_rds_path      <- file.path(project_dir, "data/JASPAR2024_TFs.rds")
rank_tf_list_path <- file.path(project_dir, "data/list_rank_product_tfs.txt")

# ---- Outputs read by 03_figures.R --------------------------------------------
common_tfs_path       <- file.path(outdir, "common_tf_h1_10_correlations_min2_CLP_ProB_PreB.csv")
xpand_pseudobulk_path <- file.path(outdir, "tf_h1_10_pseudobulk_xpand_CLP_ProB_PreB.csv")

# ---- Parameters --------------------------------------------------------------
# TF~H1-10 correlation thresholds
tf_corr_threshold  <- 0.50
tf_fdr_threshold   <- 0.05   # significance is called on the BH-adjusted p-value
tf_min_datasets    <- 2

target_gene <- "H1-10"

# TFs belonging to the H1-10 co-expression module
tfs_h1_10 <- c(
  "MEF2D", "RARB", "IRF2", "EBF1", "IRF4", "SOX4",
  "ZSCAN16", "BACH2", "MYB", "KLF10", "TCF3", "ETV2"
)

# Cell types used for the TF~H1-10 expression correlations
correlation_celltypes <- c("CLP", "ProB", "PreB")

# Pseudobulk QC settings for the TF~H1-10 correlations
min_cells_per_sample     <- 30
min_samples_per_celltype <- 3

# ---- Gene names --------------------------------------------------------------
# Old symbols (used in the ATAC gene annotation) -> harmonized name
gene_name_map <- c(
  "H1F0"     = "H1-0", "HIST1H1A" = "H1-1",
  "HIST1H1C" = "H1-2", "HIST1H1D" = "H1-3",
  "HIST1H1E" = "H1-4", "HIST1H1B" = "H1-5",
  "H1FX"     = "H1-10"
)

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

xpand_atac <- readRDS(xpand_atac_preprocessed_path)
DefaultAssay(xpand_atac) <- "ATAC"
Idents(xpand_atac) <- "annotation_2"

# ==============================================================================
# 1. H1 PROMOTER PEAKS (ATAC)
# ==============================================================================
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
# 2. CANDIDATE TFs BOUND TO THE H1-10 PROMOTER (XPAND, JASPAR MOTIF MATCH)
# ==============================================================================
message("--> Loading target TF list and matching JASPAR motifs on Xpand ATAC...")

tfs_rds  <- readRDS(tfs_rds_path)
all_tfs  <- if (is.data.frame(tfs_rds)) {
  gene_col <- intersect(c("Gene", "TF", "gene", "tf"), colnames(tfs_rds))[1]
  unique(tfs_rds[[gene_col]])
} else {
  unique(as.character(tfs_rds))
}

pfm <- getMatrixSet(
  x = JASPAR2024()@db,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)

motif_names_map <- xpand_atac[["ATAC"]]@motifs@motif.names

# ---- Motif -> TF mapping rule ----------------------------------------------
# 1. Dimer matrices (e.g. MA0091.2 "TAL1::TCF3") are DROPPED. A composite
#    heterodimer site is not evidence that either constituent TF binds on its
#    own, and the partner may not even be expressed in CLP/ProB/PreB.
# 2. The TF universe is the UNION of JASPAR2024_TFs.rds and the rank-product
#    list, so a TF present in either source is evaluated.
# 3. A TF counts as bound if ANY of its monomer matrices hits. The original
#    `tf_to_id_map[TF]` returned only the FIRST matrix, silently missing hits
#    from a TF's other matrices - that was a bug, not a definition.

rank_tfs <- unique(read_tsv(rank_tf_list_path, show_col_types = FALSE)$Gene)
tf_universe <- sort(union(all_tfs, rank_tfs))

message(sprintf("TF universe: %d TFs = %d from JASPAR2024_TFs.rds U %d from the rank-product list (%d shared).",
                length(tf_universe), length(all_tfs), length(rank_tfs),
                length(intersect(all_tfs, rank_tfs))))

motif_tf_tbl <- tibble(
  Motif_ID   = rep(names(motif_names_map), lengths(motif_names_map)),
  Motif_Name = unlist(motif_names_map, use.names = FALSE)
) %>%
  filter(!grepl("::", Motif_Name, fixed = TRUE)) %>%   # drop dimer matrices
  mutate(TF = trimws(Motif_Name)) %>%
  filter(TF %in% tf_universe) %>%
  distinct(Motif_ID, Motif_Name, TF)

message(sprintf("Motif map: %d monomer matrices covering %d TFs of the universe.",
                n_distinct(motif_tf_tbl$Motif_ID), n_distinct(motif_tf_tbl$TF)))

# TFs in the universe that JASPAR has no monomer matrix for: they can never be
# called TFBS-positive, so their has_TFBS = FALSE means "not testable".
tfs_without_motif <- setdiff(tf_universe, motif_tf_tbl$TF)
message(sprintf("  %d TF(s) in the universe have no monomer JASPAR matrix.",
                length(tfs_without_motif)))

# Module TFs outside the universe would be silently called FALSE - flag them.
module_outside_universe <- setdiff(tfs_h1_10, tf_universe)
if (length(module_outside_universe) > 0) {
  warning("H1-10 module TFs absent from the TF universe (has_TFBS forced FALSE): ",
          paste(module_outside_universe, collapse = ", "))
}
module_without_motif <- setdiff(intersect(tfs_h1_10, tf_universe), motif_tf_tbl$TF)
if (length(module_without_motif) > 0) {
  message("H1-10 module TFs with no monomer JASPAR matrix: ",
          paste(module_without_motif, collapse = ", "))
}

h1_peaks     <- unique(clean_map$Peak_ID)
motif_matrix <- xpand_atac[["ATAC"]]@motifs@data

# ---- 2a. Motif presence in H1 promoter peaks --------------------------------
# Restricted to the monomer matrices of the TF universe defined above.
h1_motif_all <- as.data.frame(as.matrix(motif_matrix[h1_peaks, , drop = FALSE])) %>%
  rownames_to_column("Peak_ID") %>%
  pivot_longer(cols = -Peak_ID, names_to = "Motif_ID", values_to = "Motif_Present") %>%
  filter(Motif_Present > 0) %>%
  inner_join(motif_tf_tbl, by = "Motif_ID") %>%
  left_join(clean_map, by = "Peak_ID") %>%
  select(Gene_Final, Peak_ID, Motif_ID, Motif_Name, TF) %>%
  distinct()

# TFBS table restricted to the H1-10 promoter peaks
tfbs_h1_10 <- h1_motif_all %>% filter(Gene_Final == target_gene)

tfbs_summary <- tfbs_h1_10 %>%
  group_by(TF) %>%
  summarise(
    has_TFBS     = TRUE,
    n_TFBS_peaks = n_distinct(Peak_ID),
    TFBS_peaks   = paste(sort(unique(Peak_ID)), collapse = ";"),
    TFBS_motifs  = paste(sort(unique(Motif_ID)), collapse = ";"),
    .groups = "drop"
  )

# ---- 2b. Candidate list (TFs from the active TF universe with a motif hit) --
tf_h1_bindings <- h1_motif_all %>%
  filter(TF %in% tf_universe) %>%
  select(Gene_Final, Peak_ID, TF) %>%
  distinct()

tf_h1_10_bindings   <- tf_h1_bindings %>% filter(Gene_Final == target_gene)
candidate_tfs_h1_10 <- unique(tf_h1_10_bindings$TF)

message(sprintf("Identified %d unique TFs from your list bound to open H1-10 promoter peaks in Xpand.",
                length(candidate_tfs_h1_10)))

write.table(
  candidate_tfs_h1_10,
  file = file.path(outdir, "candidate_tfs_H1_10_open_regions2.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

# ---- 2c. Full TF set to be correlated --------------------------------------
# Union of motif-matched candidates and the H1-10 module, so the output table
# also covers module TFs without a TFBS and TFBS TFs outside the module.
tfs_to_test <- sort(union(candidate_tfs_h1_10, tfs_h1_10))
message(sprintf("Testing %d TFs in total (%d motif candidates U %d module TFs).",
                length(tfs_to_test), length(candidate_tfs_h1_10), length(tfs_h1_10)))

module_tfs_no_tfbs <- setdiff(tfs_h1_10, tfbs_summary$TF)
if (length(module_tfs_no_tfbs) > 0) {
  message("Module TFs with no JASPAR motif hit in the H1-10 promoter peaks: ",
          paste(module_tfs_no_tfbs, collapse = ", "))
}

# ==============================================================================
# 3. EXACT TF MOTIF COORDINATES WITHIN H1-10 PEAKS
# ==============================================================================
message("--> Finding exact genomic coordinates for TFs within H1-10 peaks...")

h1_10_peak_ids <- unique(tfbs_h1_10$Peak_ID)
h1_10_peaks_gr <- StringToGRanges(h1_10_peak_ids, sep = c("-", "-"))

# All motifs whose TF is in the tested set
candidate_motif_ids <- tfbs_h1_10 %>%
  filter(TF %in% tfs_to_test) %>%
  pull(Motif_ID) %>% unique()
h1_10_pfms <- pfm[candidate_motif_ids]

motif_positions <- matchMotifs(
  pwms    = h1_10_pfms,
  subject = h1_10_peaks_gr,
  genome  = BSgenome.Hsapiens.NCBI.GRCh38,
  out     = "positions"
)

motif_id_to_tf <- motif_tf_tbl %>%
  group_by(Motif_ID) %>%
  summarise(TF = paste(sort(unique(TF)), collapse = "::"), .groups = "drop")

tf_exact_locations <- as.data.frame(motif_positions) %>%
  rename(Motif_ID = group_name) %>%
  left_join(motif_id_to_tf, by = "Motif_ID") %>%
  select(TF, Motif_ID, seqnames, start, end, strand, score) %>%
  arrange(seqnames, start)

message(sprintf("Identified %d exact motif occurrences within the H1-10 ATAC peak boundaries.",
                nrow(tf_exact_locations)))

write_csv(tf_exact_locations, file.path(outdir, "H1_10_exact_TF_motif_peak_positions.csv"))

# ==============================================================================
# 4. RNA PSEUDOBULK TF~H1-10 CORRELATION FUNCTION (per dataset)
# ==============================================================================
# Returns the UNFILTERED statistics for every tested TF, plus a Sig flag.
# BH adjustment is done across all TFs successfully tested in that dataset.

compute_dataset_correlations <- function(seurat_obj, dataset_name, sample_col,
                                         tf_set = tfs_to_test) {
  message(sprintf("--> Processing RNA pseudobulk TF~H1-10 correlations for %s (CLP/ProB/PreB)...", dataset_name))

  if ("RNA" %in% Assays(seurat_obj)) {
    DefaultAssay(seurat_obj) <- "RNA"
  } else if ("SCT" %in% Assays(seurat_obj)) {
    DefaultAssay(seurat_obj) <- "SCT"
  }
  rna_assay_ds <- DefaultAssay(seurat_obj)

  available_rna_genes <- rownames(seurat_obj)
  valid_rna_tfs       <- intersect(tf_set, available_rna_genes)

  missing_tfs <- setdiff(tf_set, available_rna_genes)
  if (length(missing_tfs) > 0) {
    message(sprintf("    %s: %d TF(s) absent from the RNA assay -> NA: %s",
                    dataset_name, length(missing_tfs), paste(missing_tfs, collapse = ", ")))
  }

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
      Rho      = tryCatch(cor(pb_tf_expr, pb_h1_expr, method = "spearman", use = "complete.obs"),
                          error = function(e) NA_real_),
      P_val    = tryCatch(cor.test(pb_tf_expr, pb_h1_expr, method = "spearman", exact = FALSE)$p.value,
                          error = function(e) NA_real_),
      .groups  = "drop"
    ) %>%
    mutate(FDR = p.adjust(P_val, method = "BH")) %>%
    mutate(Sig = !is.na(Rho) & !is.na(FDR) & abs(Rho) > tf_corr_threshold & FDR < tf_fdr_threshold) %>%
    arrange(desc(Rho))

  message(sprintf("    %s: %d/%d TF~H1-10 pairs pass |Rho| > %.2f & FDR < %.2f.",
                  dataset_name, sum(tf_h1_stats$Sig), nrow(tf_h1_stats),
                  tf_corr_threshold, tf_fdr_threshold))

  list(stats = tf_h1_stats, pseudobulk = corr_input)
}

# ==============================================================================
# 5. RUN TF~H1-10 CORRELATIONS ACROSS ALL 4 DATASETS
# ==============================================================================
xpand_result        <- compute_dataset_correlations(xpand_atac, "Xpand", "sample")
res_xpand           <- xpand_result$stats
xpand_tf_pseudobulk <- xpand_result$pseudobulk
rm(xpand_atac); gc()

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

# Xpand pseudobulk values, used by 03_figures.R for the correlation panel
write_csv(xpand_tf_pseudobulk, xpand_pseudobulk_path)

# ==============================================================================
# 6. SUMMARY TABLE
# ==============================================================================
# One row per TF: Rho + FDR + significance in each of the 4 datasets,
# H1-10 module membership, and TFBS presence.
message("--> Assembling the TF summary table...")

rename_ds <- function(df, suffix) {
  if (nrow(df) == 0) {
    return(tibble(TF = character()))
  }
  df %>%
    select(TF, Rho, FDR, Sig) %>%
    rename_with(~ paste0(.x, "_", suffix), -TF)
}

tf_summary <- tibble(TF = tfs_to_test) %>%
  mutate(
    in_H1_10_module = TF %in% tfs_h1_10,
    has_TFBS        = TF %in% tfbs_summary$TF
  ) %>%
  left_join(rename_ds(res_xpand, "Xpand"), by = "TF") %>%
  left_join(rename_ds(res_ainci, "Ainci"), by = "TF") %>%
  left_join(rename_ds(res_li,    "Li"),    by = "TF") %>%
  left_join(rename_ds(res_hca,   "HCA"),   by = "TF")

# Make sure every expected column exists even if a dataset returned nothing
for (sfx in c("Xpand", "Ainci", "Li", "HCA")) {
  for (col in c("Rho", "FDR", "Sig")) {
    nm <- paste0(col, "_", sfx)
    if (!nm %in% names(tf_summary)) {
      tf_summary[[nm]] <- if (col == "Sig") NA else NA_real_
    }
  }
}

sig_cols <- paste0("Sig_", c("Xpand", "Ainci", "Li", "HCA"))

tf_summary <- tf_summary %>%
  mutate(.n_sig = rowSums(across(all_of(sig_cols), ~ !is.na(.x) & .x))) %>%
  arrange(desc(.n_sig), desc(in_H1_10_module), desc(has_TFBS), TF) %>%
  select(
    TF, in_H1_10_module, has_TFBS,
    Rho_Xpand, FDR_Xpand, Sig_Xpand,
    Rho_Ainci, FDR_Ainci, Sig_Ainci,
    Rho_Li,    FDR_Li,    Sig_Li,
    Rho_HCA,   FDR_HCA,   Sig_HCA
  )

print(as_tibble(tf_summary), n = nrow(tf_summary), width = Inf)

tf_summary_path <- file.path(outdir, "H1_10_TF_summary_module_TFBS_correlations2.csv")
write_csv(tf_summary, tf_summary_path)
message("Wrote: ", tf_summary_path)

# Cross-tab overview
message("--- Module membership x TFBS status ---")
print(with(tf_summary, table(in_H1_10_module = in_H1_10_module, has_TFBS = has_TFBS)))

# ==============================================================================
# 7. LEGACY OUTPUTS (significant-only tables, unchanged semantics)
# ==============================================================================
final_results_min2 <- tf_summary %>%
  filter(rowSums(across(all_of(sig_cols), ~ !is.na(.x) & .x)) >= tf_min_datasets)

message("--- Summary: TFs significant (FDR < ", tf_fdr_threshold, ") in >= ", tf_min_datasets, " datasets ---")
print(as_tibble(final_results_min2))
write_csv(final_results_min2, common_tfs_path)

message("Xpand TFs: ", sum(res_xpand$Sig, na.rm = TRUE))
message("Ainci TFs: ", sum(res_ainci$Sig, na.rm = TRUE))
message("Li TFs:    ", sum(res_li$Sig,    na.rm = TRUE))
message("HCA TFs:   ", sum(res_hca$Sig,   na.rm = TRUE))

write_csv(res_xpand, file.path(outdir, "tf_h1_10_correlations_xpand_CLP_ProB_PreB.csv"))
write_csv(res_ainci, file.path(outdir, "tf_h1_10_correlations_ainci_CLP_ProB_PreB.csv"))
write_csv(res_li,    file.path(outdir, "tf_h1_10_correlations_li_CLP_ProB_PreB.csv"))
write_csv(res_hca,   file.path(outdir, "tf_h1_10_correlations_hca_CLP_ProB_PreB.csv"))

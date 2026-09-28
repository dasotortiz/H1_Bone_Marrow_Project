# ==============================================================================
# 01. PREPROCESSING - Xpand multiome (RNA + ATAC)
# ==============================================================================
# Input : seurat_path (raw Xpand multiome object) + per-sample fragment files
# Output: xpand_atac_preprocessed_path (RNA/ATAC reductions, fragments, motifs)

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(JASPAR2024)
  library(TFBSTools)
  library(BSgenome.Hsapiens.NCBI.GRCh38)
  library(motifmatchr)
})

# ==============================================================================
# 0. CONFIG
# ==============================================================================
project_dir <- "/home/x_vazquem/workspace/h1_project/panel_3_C-E"

seurat_path <- file.path(project_dir, "data/Xpand_atac.rds")
base_dir    <- file.path(project_dir, "data/atac/xpand_20260223")
samples     <- c("MO291", "MO292", "MO294", "MO296", "MO298")

xpand_atac_preprocessed_path <- file.path(project_dir, "results/preprocessed_seurats/Xpand_atac_preprocessed.rds")

# ==============================================================================
# 1. LOAD SEURAT OBJECT
# ==============================================================================
xpand_atac <- readRDS(seurat_path)

# ==============================================================================
# 2. RNA PRE-PROCESSING (only if RNA assay present and PCA not yet computed)
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
    message("PCA already found - skipping RNA preprocessing.")
  }
}

# ==============================================================================
# 3. ATAC PRE-PROCESSING
# ==============================================================================
DefaultAssay(xpand_atac) <- "ATAC"

xpand_atac <- RunTFIDF(xpand_atac)
xpand_atac <- FindTopFeatures(xpand_atac, min.cutoff = "q0")
xpand_atac <- RunSVD(xpand_atac)
xpand_atac <- RunUMAP(xpand_atac, reduction = "lsi", dims = 2:30,
                      reduction.name = "umap.atac", reduction.key = "atacUMAP_")
xpand_atac <- FindNeighbors(xpand_atac, reduction = "lsi", dims = 2:30)
xpand_atac <- FindClusters(xpand_atac, verbose = FALSE, algorithm = 3)

# FindClusters resets the identities
Idents(xpand_atac) <- "annotation_2"

# ==============================================================================
# 4. FRAGMENT RE-BINDING
# ==============================================================================
Fragments(xpand_atac) <- NULL

frag_list <- lapply(samples, function(s) {
  frag_path <- file.path(base_dir, s, "atac_fragments.tsv.gz")
  prefix_pattern <- paste0("^_", s, "_")
  seurat_cells <- grep(prefix_pattern, colnames(xpand_atac), value = TRUE)
  raw_barcodes <- sub(prefix_pattern, "", seurat_cells)
  cells_mapped <- setNames(raw_barcodes, seurat_cells)
  CreateFragmentObject(path = frag_path, cells = cells_mapped)
})

Fragments(xpand_atac) <- frag_list

# ==============================================================================
# 5. JASPAR MOTIF ANNOTATION OF ATAC PEAKS
# ==============================================================================
if (is.null(Motifs(xpand_atac[["ATAC"]]))) {
  message("--> Annotating ATAC peaks with JASPAR2024 motifs...")
  pfm <- getMatrixSet(
    x = JASPAR2024()@db,
    opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
  )
  xpand_atac <- AddMotifs(object = xpand_atac, genome = BSgenome.Hsapiens.NCBI.GRCh38, pfm = pfm)
}

# ==============================================================================
# 6. SAVE
# ==============================================================================
saveRDS(xpand_atac, xpand_atac_preprocessed_path)
message("Wrote: ", xpand_atac_preprocessed_path)

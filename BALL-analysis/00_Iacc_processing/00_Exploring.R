# Build one Seurat object from the Iacobucci et al. samples, QC it, and run the
# standard pipeline. Dataset: Nature Cancer 2025, doi:10.1038/s43018-025-00987-2
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readxl)
library(scDblFinder)
library(BiocParallel)
set.seed(1234)

# 1. Get sample prefixes (remove just the ".matrix.mtx.gz")
files <- list.files("./data/samples", pattern = "matrix.mtx.gz", full.names = FALSE)
samples <- gsub(".matrix.mtx.gz", "", files, fixed = TRUE)
samples

# 2. Read each sample
seu_list <- lapply(samples, function(s) {
  counts   <- readMM(file.path("./data/samples", paste0(s, ".matrix.mtx.gz")))
  features <- read.delim(file.path("./data/samples", paste0(s, ".features.tsv.gz")),
                         header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(file.path("./data/samples", paste0(s, ".barcodes.tsv.gz")),
                         header = FALSE, stringsAsFactors = FALSE)

  # Use gene symbols for rownames
  rownames(counts) <- make.unique(features[, 2])
  colnames(counts) <- barcodes[, 1]

  seu <- CreateSeuratObject(counts, project = s, min.cells = 0, min.features = 0)
  # Prefix cell IDs with sample name
  RenameCells(seu, add.cell.id = s)
})

# 3. Merge all samples
combined <- merge(seu_list[[1]], y = seu_list[-1])
combined # 890387 cells across 89 samples!
saveRDS(combined, "./output/combinedobj_raw.rds")

# 4. Load metadata and filter samples
# "We analyzed children, adolescents and young adults (AYA) with newly diagnosed
# B-ALL (n = 84) and matched relapsed B-ALL (n = 5)"
metadata <- read_excel("./data/md.xlsx")
dim(metadata) # 89 x 21
table(metadata$'Disease phase') # only keep Initial diagnosis samples = 84 samples
metadata <- metadata[metadata$'Disease phase' == "Initial diagnosis",]
table(metadata$'Tissue') # only keep Bone marrow samples = 81 samples
metadata <- metadata[metadata$'Tissue' == "Bone marrow",]
table(metadata$"RNA") # all Yes
table(metadata$"Risk group") # only keeping childhood related entries (to be sure)
metadata <- metadata[metadata$'Risk group' %in% c("Childhood HR", "Childhood SR"),]
dim(metadata) # 51 x 21

# Subset the object to those samples:
# Remove everything before first underscore
new_names <- sub("^[^_]+_", "", colnames(combined))
# Extract sample IDs with the underscore preserved (SJALL040053_D1)
sample_ids <- sub("^([^_]+_[^_]+).*", "\\1", new_names)
# Keep cells with matching sample IDs
keep_cells <- colnames(combined)[sample_ids %in% metadata$'Sample ID']
combinedf <- subset(combined, cells = keep_cells)

# 5. Attach metadata to object:
colnames(combinedf) <- sub("^[^_]+_", "", colnames(combinedf)) # clean colnames()
combinedf$orig.ident <- sub("^[^_]+_", "", combinedf$orig.ident) # clean orig.ident
# QC:
length(unique(combinedf$orig.ident)) # 51! Perfect :D
table(combinedf$orig.ident) # number of cells per sample (or patient)
# Attach
combinedf$Subtype  <- metadata$Subtype[match(combinedf$orig.ident, metadata$'Sample ID')]
combinedf$Subgroup <- metadata$Subgroup[match(combinedf$orig.ident, metadata$'Sample ID')]
combinedf$Age      <- metadata$'Age (years)'[match(combinedf$orig.ident, metadata$'Sample ID')]

table(combinedf$orig.ident %in% metadata$'Sample ID')
saveRDS(combinedf, "./output/combinedfiltered.rds")

combinedfinal <- readRDS("./output/combinedfiltered.rds")
# QC and plot:
# 1) % mitochondrial
combinedfinal[["percent.mt"]] <- PercentageFeatureSet(combinedfinal, pattern = "^MT-")
# Extract metadata
df <- combinedfinal@meta.data
pmt <- ggplot(df, aes(x = percent.mt, fill = orig.ident)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  labs(
    title = "Mitochondrial Content Distribution",
    x = "Percent Mitochondrial Genes",
    y = "Density"
  ) +
  theme_minimal()
pmt
ggsave("./output/QC_percent_mt.png", pmt, bg = "white")

# 2) nFeature and nCount
p2 <- FeatureScatter(combinedfinal, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
ggsave("./output/QC_ncount_nfeature_before.png", p2, bg = "white")

# 3) Doublets
x <- JoinLayers(combinedfinal)
x <- as.SingleCellExperiment(x)
# default is clusters(the mode)=FALSE, run per sample based on orig.ident, multithread
x <- scDblFinder(x, samples = "orig.ident", BPPARAM = MulticoreParam(12))
table(x$scDblFinder.class)
# singlet doublet
# 449,523   59,578
combinedfinal$dbl.class <- x$scDblFinder.class
combinedfinal$dbl.score <- x$scDblFinder.score
saveRDS(combinedfinal, "./output/combinedobj_qc.rds")

# Filtering: initial 509,101
# Note: these are my own thresholds, the paper's methods used
# percent.mt <= 8, nFeature >= 500 and nCount >= 2500.
f1 <- subset(combinedfinal,
              subset = dbl.class == "singlet" &
              percent.mt < 10 &
              nFeature_RNA > 300 &
              nFeature_RNA < 7000) # keep singlets and less than 10% mito and cells with more than 300 genes
f1 # 349,769, around 69% of total kept
saveRDS(f1, "./output/combinedobj_raw_filtered.rds")

pb <- VlnPlot(combinedfinal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", pt.size = 0)
pa <- VlnPlot(f1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", pt.size = 0)
ggsave("./output/QC_before_after.png", pb / pa , bg = "white", w = 24, h = 10)

# Standard pipeline
# Normalize data
combined <- f1
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
# Cell cycle: score (Tirosh et al. 2015 markers, shipped with Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
combined <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes,
                             search = TRUE, set.ident = FALSE)
combined$CC.Difference <- combined$S.Score - combined$G2M.Score
combined <- RunPCA(combined)
pca <- DimPlot(combined, reduction = "pca", group.by = "orig.ident")
pcap <- DimPlot(combined, reduction = "pca", group.by = "Phase")
ggsave("./output/PCA.png", pca | pcap , bg = "white")

combined <- FindNeighbors(combined, dims = 1:30, reduction = "pca")
combined <- FindClusters(combined, resolution = 0.2)
combined <- RunUMAP(combined, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
saveRDS(combined, "./output/combinedobj_raw_filtered_processed.rds")

umap <- DimPlot(combined, reduction = "umap.unintegrated", group.by = "orig.ident", label = TRUE, label.box = TRUE, repel = TRUE)
umapp <- DimPlot(combined, reduction = "umap.unintegrated", group.by = "Phase")
ggsave("./output/UMAP.png", umap | umapp , bg = "white", w = 17, h = 7)

# H1 genes across subtypes (quick look, before healthy cells are removed)
# h1_df.rds lives outside this repo, set the path to wherever you keep it
h1_path <- "../h1_df.rds"
h1_df <- readRDS(h1_path)
h1_sym <- h1_df$Symbol
h1_sym2 <- h1_df$Alternative_symbol

# object uses the alternative symbols, so subset on those then rename back
h1obj <- subset(combined, features = h1_sym2)
rownames(h1obj) <- h1_df$Symbol[match(rownames(h1obj), h1_df$Alternative_symbol)]

pv <- VlnPlot(h1obj, features = h1_sym, group.by = "Subtype", stack = TRUE)
pv
ggsave("./output/VlnPlot_ballSC.png", pv , width = 20, height = 20, units = "cm", bg = "white")

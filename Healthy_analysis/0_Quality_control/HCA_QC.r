library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)

index <- readxl::read_xlsx("/ibex/user/sotoorda/masterh1/project_raw_data/hca_raw/index.xlsx")

seurat_objects <- list()
for (i in 1:nrow(index)){
    ProjectName <- index[i,]$sample_id
    counts <- Read10X(index[i,]$path)
    x <- CreateSeuratObject(counts, min.cells = 1, min.features = 300, project = ProjectName)
    seurat_objects <- c(seurat_objects, x)
}

for(i in 1:length(seurat_objects)){
  seurat_objects[[i]][["log10GenesPerUMI"]] <- log10(seurat_objects[[i]]$nFeature_RNA) / log10(seurat_objects[[i]]$nCount_RNA) 
  # %mt 
  seurat_objects[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_objects[[i]], pattern = c("^MT-"), assay ='RNA')
  # %rb 
  seurat_objects[[i]][["percent.rb"]] <- PercentageFeatureSet(seurat_objects[[i]], pattern = c("^RPS"), assay ='RNA')
  # identify doublets using scDblFinder
  sce <- as.SingleCellExperiment(seurat_objects[[i]])
  sce <- scDblFinder(sce)
  seurat_objects[[i]]$scDblFinder.class <- sce$scDblFinder.class
}

saveRDS(seurat_objects, file = "/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/raw_hca.rds")
seurat_objects_merged <- readRDS(file = "/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/raw_hca.rds")
seurat_objects <- SplitObject(seurat_objects_merged, split.by = 'orig.ident')

# Merge for QC plots
merged_seurat <- merge(seurat_objects[[1]], y = seurat_objects[2:8], add.cell.ids = c("manton1", "manton2", "manton3", "manton4", "manton5", "manton6", "manton7", 'manton8'))
merged_seurat <- JoinLayers(merged_seurat)
saveRDS(merged_seurat, file = "/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/raw_hca.rds")

# Joint QC plots
pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/hca/QC/mergedQC.pdf")
VlnPlot(merged_seurat, features = "nFeature_RNA", ncol = 1, split.by = "orig.ident", alpha = 0, group.by = "orig.ident") + geom_boxplot(width=0.1,fill="white") +
  theme(legend.position = 'none') 
VlnPlot(merged_seurat, features = "nCount_RNA", ncol = 1, split.by = "orig.ident", alpha = 0, group.by = "orig.ident") + geom_boxplot(width=0.1,fill="white") +
  theme(legend.position = 'none') 
VlnPlot(merged_seurat, features = "percent.mt", ncol = 1, split.by = "orig.ident", alpha = 0, group.by = "orig.ident") + geom_boxplot(width=0.1,fill="white") +
  theme(legend.position = 'none') 
dev.off()

# Scatter plot to check what to discard
filter_ncount <- 500
filter_nfeature_quantile <- 0.9
filter_ncount_quantile <- 0.9
filter_mt <- 10

# Assessing the number of cells to discard
filter_nfeature_q <- quantile(merged_seurat$nFeature_RNA, filter_nfeature_quantile)
filter_ncount_q <- quantile(merged_seurat$nCount_RNA, filter_ncount_quantile)
merged_seurat$mt_discard10 <- ifelse(merged_seurat$percent.mt > 10 | merged_seurat$nFeature_RNA < 500 | merged_seurat$nFeature_RNA > filter_ncount_q |  
  merged_seurat$nCount_RNA > filter_ncount_q , "discard", "keep")
table(merged_seurat$mt_discard10)
table(merged_seurat$scDblFinder.class)

# What cells will be discarded? Scatter plot showing this
merged_seurat_df <- as.data.frame(merged_seurat@meta.data)
pdf('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/QC/unfiltered_scatter_plot_qc.pdf')
ggplot(merged_seurat_df, aes(x = nCount_RNA, y = nFeature_RNA, color = mt_discard10)) + 
  facet_wrap(~ orig.ident, ncol = 4) +
  geom_point(size = 0.005) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic()
  dev.off()

# Definitive filtering parameters
filter_nfeature <- 500
filter_nfeature_quantile <- 0.9
filter_ncount_quantile <- 0.9
filter_mt <- 10

# Filtering
seurat_filtered <- list()
for (i in 1:length(seurat_objects)){
filter_nfeature_q <- quantile(seurat_objects[[i]]$nFeature_RNA, filter_nfeature_quantile)
filter_ncount_q <- quantile(seurat_objects[[i]]$nCount_RNA, filter_ncount_quantile)
seurat_filtered[[i]] <- subset(seurat_objects[[i]], subset = nCount_RNA < filter_ncount_q & nFeature_RNA < filter_nfeature_q)
seurat_filtered[[i]] <- subset(seurat_filtered[[i]], subset = percent.mt < filter_mt)
seurat_filtered[[i]] <- subset(seurat_filtered[[i]], subset = nFeature_RNA > filter_nfeature)
seurat_filtered[[i]] <- subset(seurat_filtered[[i]], subset = scDblFinder.class == "singlet")
}
summary(seurat_filtered[[1]]$percent.mt)

# Plotting filtered data
merged_seurat_filt <- merge(seurat_filtered[[1]], y = seurat_filtered[2:8], add.cell.ids = c("manton1", "manton2", "manton3", "manton4", "manton5", "manton6", "manton7", 'manton8'))
pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/hca/QC/mergedQC_filtered.pdf")
VlnPlot(merged_seurat_filt, features = "nFeature_RNA", ncol = 1, split.by = "orig.ident", alpha = 0, group.by = "orig.ident") + geom_boxplot(width=0.1,fill="white") +
  theme(legend.position = 'none') 
VlnPlot(merged_seurat_filt, features = "nCount_RNA", ncol = 1, split.by = "orig.ident", alpha = 0, group.by = "orig.ident") + geom_boxplot(width=0.1,fill="white") +
  theme(legend.position = 'none') 
VlnPlot(merged_seurat_filt, features = "percent.mt", ncol = 1, split.by = "orig.ident", alpha = 0, group.by = "orig.ident") + geom_boxplot(width=0.1,fill="white") +
  theme(legend.position = 'none') 
dev.off()

# Scatter plot of filtered data
merged_seurat_df <- as.data.frame(merged_seurat_filt@meta.data)
pdf('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/QC/scatter_plot_qc.pdf')
ggplot(merged_seurat_df, aes(x = nCount_RNA, y = nFeature_RNA)) + 
  facet_wrap(~ orig.ident, ncol = 4) +
  geom_point(size = 0.005) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic()
  dev.off()

# Saving filtered data
saveRDS(seurat_filtered, file = "/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_filtered.rds")
seurat_filtered <- readRDS(file = "/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_filtered.rds")







######## PLOTTING THE HCA ORIGINAL DATA ########

merged_seurat_filt <- readRDS(file = "/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/MantonBM.rds")
saveRDS(merged_seurat_filt, file = "/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/MantonBM.rds")

  merged_seurat_filt[["log10GenesPerUMI"]] <- log10(merged_seurat_filt$nFeature_RNA) / log10(merged_seurat_filt$nCount_RNA) 
  # %mt 
  merged_seurat_filt[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_filt, pattern = c("^MT-"), assay ='RNA')
  # %rb 
  merged_seurat_filt[["percent.rb"]] <- PercentageFeatureSet(merged_seurat_filt, pattern = c("^RPS"), assay ='RNA')
colnames(merged_seurat_filt@meta.data) 

pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/hca/QC/mergedQC_OriginalData.pdf")
VlnPlot(merged_seurat_filt, features = "nFeature_RNA", ncol = 1, split.by = "orig.ident",  group.by = "orig.ident", pt.size = 0) + geom_boxplot(width=0.1,fill="white") +
  theme(legend.position = 'none') 
VlnPlot(merged_seurat_filt, features = "nCount_RNA", ncol = 1, split.by = "orig.ident",  group.by = "orig.ident", pt.size = 0) + geom_boxplot(width=0.1,fill="white") +
  theme(legend.position = 'none') 
VlnPlot(merged_seurat_filt, features = "percent.mt", ncol = 1, split.by = "orig.ident",  group.by = "orig.ident", pt.size = 0) + geom_boxplot(width=0.1,fill="white") +
  theme(legend.position = 'none') 
dev.off()


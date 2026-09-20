library(dplyr)
library(readxl)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)

#### Reading Metadata ####
metadata <- read.csv("/ibex/user/sotoorda/masterh1/GSE189161/GSE189161_metadata.csv")

#### STEP 1: Creating object list ####
counts <- Read10X('/ibex/user/sotoorda/masterh1/GSE189161')
se <- CreateSeuratObject(counts, min.cells = 3, min.features = 300)
se <- AddMetaData(object = se, metadata = metadata)

#### STEP 2: QC and filtering ####
sce <- as.SingleCellExperiment(se)
samples <- as.factor(sce$orig.ident)
sce <- scDblFinder(sce, samples = samples)
se$scDblFinder.class <- sce$scDblFinder.class

# Removing samples that are not of interest
excluded_samples <- c('GR97', 'GR72', 'GR73', 'GR74', 'CB1', 'CB2', 'Samp22', 'Samp23', 'Samp41', 'Samp42', 'Samp21')
for (i in 1:length(excluded_samples)){
  se <- subset(se, subset = orig.ident == excluded_samples[i], invert = TRUE)
}

# Mithocondrial and ribosomal percentage
se[["log10GenesPerUMI"]] <- log10(se$nFeature_RNA) / log10(se$nCount_RNA) 
se[["percent.mt"]] <- PercentageFeatureSet(se, pattern = c("^MT\\."), assay ='RNA')
se[["percent.rb"]] <- PercentageFeatureSet(se, pattern = c("^RPL|^RPS"), assay ='RNA')

# Replacing sample names for plotting
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "GR95", "2yr")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "Samp44", "4yr")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "GR96", "10y")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "Samp43", "12yr")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "GR90", "17yr")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "BM1", "25yr_a")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "BM3", "25yr_b")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "BM2", "32yr")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "Samp24", "35yr")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "Samp25", "45yr")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "Samp26", "53yr")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "GR91", "62yr_a")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "GR94", "62yr_b")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "GR92", "76yr")
se@meta.data$orig.ident <- replace(se@meta.data$orig.ident, se@meta.data$orig.ident == "GR93", "77yr")

# Defining the order of my samples in the plots
se@meta.data$orig.ident <- factor(se@meta.data$orig.ident, levels = c("2yr", "4yr", "10y", "12yr", "17yr", "25yr_a", "25yr_b", 
"32yr", "35yr", "45yr", "53yr", "62yr_a", "62yr_b", "76yr", "77yr"))

# Splitting object by sample
se_split <- SplitObject(se, split.by = "orig.ident")

#Filtering parameters
filter_ncount <- 500
filter_nfeature_quantile <- 0.1
filter_mt <- 10

# Filtering
for (i in 1:length(se_split)){
filter_nfeature <- quantile(se_split[[i]]@meta.data$nFeature_RNA, filter_nfeature_quantile)
se_split[[i]] <- subset(se_split[[i]], subset = scDblFinder.class == "singlet")
se_split[[i]] <- subset(se_split[[i]], subset = nCount_RNA > filter_ncount & percent.mt < filter_mt & nFeature_RNA > filter_nfeature)
}

# Merging filtered objects
se_merged <- merge(se_split[[1]], y = se_split[2:length(se_split)])
se_merged <- JoinLayers(se_merged)

saveRDS(se_merged, file = "/ibex/user/sotoorda/masterh1/public_data_analysis/li_2021/objects/li_filtered.rds")
library(Seurat)
library(Signac)
library(Matrix)
library(batchelor)
library(SeuratWrappers)

ainci <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/ainci_final_annot.rds')
xpand <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/xpand_final_annot.rds')
li <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/li_final_annot.rds')

# Merging datasets and splitting the layers by dataset
merged_datasets <- merge(x = xpand, y = c(ainci, li))
merged_datasets <- JoinLayers(merged_datasets)
merged_datasets <- split(merged_datasets, f = merged_datasets$Dataset)

# Preprocessing functions
npcs<-30
merged_datasets <- NormalizeData(merged_datasets)
merged_datasets <- FindVariableFeatures(merged_datasets)
merged_datasets <- ScaleData(merged_datasets)
merged_datasets <- RunPCA(merged_datasets, npcs = npcs)
# merged_datasets <-RunUMAP(object = merged_datasets, dims = 1:npcs, reduction = "pca", reduction.name = "umap_unintegrated", seed.use = 123)

# Integration
options(future.globals.maxSize = 64 * 1024^3) # Increase memory allocation for the integration
# integrated <- IntegrateLayers(object = merged_datasets, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE, assay = "RNA")
# integrated <- IntegrateLayers(object = integrated, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = TRUE)
integrated <- IntegrateLayers(object = integrated, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", verbose = FALSE, assay = "RNA")
# integrated <- IntegrateLayers(object = merged_datasets, method = FastMNNIntegration, new.reduction = "integrated.mnn", verbose = FALSE, assay = "RNA")
# integrated <- FastMNNIntegration(object = integrated, reduction = "integrated.mnn", dims = 1:npcs, verbose = FALSE)
# re-join layers after integration
# integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

saveRDS(integrated, file = '/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/integration_post_gene_fix.rds')
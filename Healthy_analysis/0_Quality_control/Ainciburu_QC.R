library(dplyr)
library(readxl)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)

create_se <- source('/ibex/user/sotoorda/masterh1/scripts/R/Functions/create_se_h5.R')$value
QC_matrices <- source('/ibex/user/sotoorda/masterh1/scripts/R/Functions/QC_matrices.R')$value
QC_filter_cells <- source('/ibex/user/sotoorda/masterh1/scripts/R/Functions/QC_filter_cells.R')$value
preprocess_seurat <- source('/ibex/user/sotoorda/masterh1/scripts/R/Functions/normalization.R')$value
cluster_seurat <- source('/ibex/user/sotoorda/masterh1/scripts/R/Functions/clustering.R')$value


#### Reading Metadata ####
metadata <- readxl::read_excel('/ibex/user/sotoorda/masterh1/GSE180298/raw_data/metadata.xlsx')
#### STEP 1: Creating seurat object ####
object_list <- list()
for (i in 1:nrow(metadata)){
  if (metadata[i,]$Processed == "NO"){
    x <- create_se(dir = metadata[i,]$filtered_feature_bc_matrix_PATH, min.cells = 3, min.features = 100, ProjectName = metadata[i,]$SampleID)
    object_list <- c(object_list, x)}
}
head(object_list[[1]]@meta.data)
object <- merge(object_list[[1]], y = object_list[2:8])

#-------------------------------------------------------------------------------------------------------------------------------------------------------

#### STEP 2: compute QC ####
object_list2 <- list()
for (i in 1:length(object_list)){
  x <- QC_matrices(object_list[[i]])
  object_list2 <- c(object_list2, x)
}

#------------------------------------------------------------------------------------------------------------------------------------------------------- # nolint

#### STEP 3: remove cells ####

object_list_filt <- list()
object_list_filt[[1]] <- subset(object_list2[[1]], subset = percent.mt < 5 & nFeature_RNA > 500 & nCount_RNA < 30000)
object_list_filt[[2]] <- subset(object_list2[[2]], subset = percent.mt < 5 & nFeature_RNA > 500 & nCount_RNA < 21000)
object_list_filt[[3]] <- subset(object_list2[[3]], subset = percent.mt < 5 & nFeature_RNA > 500 & nCount_RNA < 40000)
object_list_filt[[4]] <- subset(object_list2[[4]], subset = percent.mt > 0.3 & percent.mt < 5 & nFeature_RNA > 500 & nFeature_RNA < 4500)
object_list_filt[[5]] <- subset(object_list2[[5]], subset = percent.mt < 10 & nFeature_RNA > 700 & nCount_RNA < 40000)
object_list_filt[[6]] <- subset(object_list2[[6]], subset = percent.mt > 0.3 & percent.mt < 5 & nFeature_RNA > 500 & nCount_RNA < 40000)
object_list_filt[[7]] <- subset(object_list2[[7]], subset = percent.mt > 0.5 & percent.mt < 5 & nFeature_RNA > 500 & nCount_RNA < 35000)
object_list_filt[[8]] <- subset(object_list2[[8]], subset = percent.mt > 0.5 & percent.mt < 10 & nFeature_RNA > 500 & nCount_RNA < 35000)

# change file name/path
saveRDS(object_list_filt, file = "/ibex/user/sotoorda/masterh1/public_data_analysis/ainciburu_2023/objects/replicate_ainciburu_filtered.rds")
#object_list_filt <- readRDS(file = "/ibex/user/sotoorda/masterh1/public_data_analysis/ainciburu_2023/objects/replicate_ainciburu_filtered.rds")
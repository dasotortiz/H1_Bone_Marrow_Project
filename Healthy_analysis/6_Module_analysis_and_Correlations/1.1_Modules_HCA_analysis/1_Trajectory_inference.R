library(Seurat)
library(SeuratWrappers)
library(monocle3)
# library(units) requieres module 'units' to be loaded 
library(Signac)


get_earliest_principal_node <- function(cds, timeStart){
  cell_ids <- which(colData(cds)[, "annotation_2"] == timeStart)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}

cells <- c("HSC", "MPP", "CLP", "ProB", "PreB", "ERP", "Eo/B/Mast", "GMP", "MDP", "MEP", "MKP")
hca_dataset <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_merging_ready.rds')
hca_subset <- subset(hca_dataset, subset = annotation_2 %in% cells)

embeddings <- Embeddings(hca_subset,reduction = "umap.integrated")
hca_subset@reductions <- list() # Clean the embeddings slot to avoid conflicts
hca_subset@reductions <- list(UMAP = CreateDimReducObject(embeddings = embeddings, key = "UMAP_", assay = DefaultAssay(hca_subset)))

# From seurat to cds
cds <- as.cell_data_set(hca_subset,group.by = "annotation_2", reduction_key = "UMAP_")


# cds <- estimate_size_factors(cds)
cds <- cluster_cells(cds)
colData(cds)$cluster_seurat <- hca_subset$annotation_2
rowData(cds)$gene_short_name <- rownames(cds)
colData(cds)$partition <- 1  # If you want all cells in one partition
colData(cds)$cluster <- colData(cds)$cluster_seurat

cds_list <- list(HCA = cds[, colData(cds)$Dataset == 'HCA'])

# Learn the principal graph and order cells for the HCA dataset
cds_list$HCA <- learn_graph(cds_list$HCA,use_partition = TRUE, learn_graph_control = list(ncenter = 1200))
cds_list$HCA <- order_cells(cds_list$HCA, root_pr_nodes=get_earliest_principal_node(cds, 'HSC'))

lineage_defs <- list(
  Lymphocytes = c("HSC", "MPP", "ProB", "PreB", "CLP"),
  Erythrocytes = c("HSC", "MPP", 'MEP', "ERP"),
  Monocytes_Dendritic_cells = c("HSC", "MPP", "MDP"),
  Eo_B_Mast = c("HSC", "MPP", "GMP", 'Eo/B/Mast'))

outdir <- "/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/pr_results"

# Run graph_test for each dataset and lineage, and save the results
for (ds_name in names(cds_list)) {
  cds_ds <- cds_list[[ds_name]]
  
  for (lin_name in names(lineage_defs)) {
    lineage_cells <- colnames(cds_ds)[
      colData(cds_ds)$annotation_2 %in% lineage_defs[[lin_name]]
    ]
    
    if (length(lineage_cells) > 0) {
      cds_subset <- cds_ds[, lineage_cells]
      
      message("Running graph_test for ", ds_name, " - ", lin_name)
      res <- graph_test(cds_subset, neighbor_graph = "principal_graph", cores = 20)
      
      saveRDS(
        res,
        file = file.path(
          outdir,
          paste0("ds_pr_test_res_", ds_name, "_", lin_name, ".rds")
        )
      )
    } else {
      message("Skipping ", ds_name, " - ", lin_name, " (no cells)")
    }
  }
}

# saveRDS(cds_list$HCA, '/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/cds_HCA.rds')
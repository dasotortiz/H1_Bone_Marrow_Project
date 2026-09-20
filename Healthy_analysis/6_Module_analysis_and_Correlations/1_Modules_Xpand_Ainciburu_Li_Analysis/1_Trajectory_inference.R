library(Seurat)
library(Signac)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(tidyr)
library(dplyr)
library(purrr)
library(clusterProfiler)
library(org.Hs.eg.db)

get_earliest_principal_node <- function(cds, timeStart){
  cell_ids <- which(colData(cds)[, "annotation_2"] == timeStart)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}

all_cells <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/integration_post_gene_fix.rds')

all_cells <- JoinLayers(all_cells)
# Trasnferring UMAP coordinates to monocle3 cds object
all_cells[["umap"]] <- all_cells[["umap_integrated_rpca"]]
cds <- as.cell_data_set(all_cells, group.by='annotation_2', reduction_key = "umap_",)

cds <- estimate_size_factors(cds)
colData(cds)$cluster_seurat <- all_cells$annotation_2
rowData(cds)$gene_short_name <- rownames(cds)

# copy back the consensus pseudotime to the cds object
cds <- cluster_cells(cds)

colData(cds)$partition <- 1  # If you want all cells in one partition
colData(cds)$cluster <- colData(cds)$cluster_seurat

# saveRDS(cds, '/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/all_cells.rds')
cds <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/all_cells.rds')
cds <- learn_graph(cds, learn_graph_control=list(ncenter=1200))
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, 'HSC'))

# Run the graph test for each dataset and lineage, and save the results
cds_list <- list(
  Xpand = cds[, colData(cds)$Dataset == 'Xpand'],
  Ainciburu = cds[, colData(cds)$Dataset == 'Ainciburu'],
  Li = cds[, colData(cds)$Dataset == 'Li']
)

lineage_defs <- list(
  Lymphocytes = c("HSC", "MPP", "ProB", "PreB", "CLP"),
  Erythrocytes = c("HSC", "MPP", 'MEP', "ERP"),
  Monocytes_Dendritic_cells = c("HSC", "MPP", "MDP"),
  Eo_B_Mast = c("HSC", "MPP", "GMP", 'Eo/B/Mast'))

outdir <- "/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/pr_results"

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

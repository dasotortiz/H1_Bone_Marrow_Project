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


# Load the cds data and subset it into datasets
cds <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/all_cells.rds')
cds_xpand <- cds[, colData(cds)$Dataset == 'Xpand']
cds_ainci <- cds[, colData(cds)$Dataset == 'Ainciburu']
cds_li <- cds[, colData(cds)$Dataset == 'Li']

# Run the preprocessing step neccesary for downstream analysis
cds_xpand <- preprocess_cds(cds_xpand, method = "PCA", num_dim = 50)
cds_ainci <- preprocess_cds(cds_ainci, method = "PCA", num_dim = 50)
cds_li <- preprocess_cds(cds_li, method = "PCA", num_dim = 50)

# Define  datasets and lineages
datasets <- c(Xpand = cds_xpand, Ainciburu = cds_ainci, Li = cds_li)
lineage_defs <- list(
  Lymphocytes = c("HSC", "MPP", "ProB", "PreB", "CLP"),
  Erythrocytes = c("HSC", "MPP", 'MEP', "ERP"),
  Monocytes_Dendritic_cells = c("HSC", "MPP", "MDP"),
  Eo_B_Mast = c("HSC", "MPP", "GMP", "Eo/B/Mast")
)
lineage_names <- names(lineage_defs)
# gene_module_df <- find_gene_modules(cds_lineage[pr_deg_ids,], resolution=1e-2)

# Load previous graph test results into a list
outdir <- "/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/pr_results"
pr_results <- list()

for (ds_name in names(datasets)) {
  for (lin_name in lineage_names) {
    file_path <- file.path(outdir, paste0("ds_pr_test_res_", ds_name, "_", lin_name, ".rds"))
    pr_results[[ds_name]][[lin_name]] <- readRDS(file_path)
  }
}
str(pr_results)
# Define lineages
lineage_defs <- list(
  Lymphocytes = c("HSC", "MPP", "ProB", "PreB", "CLP"),
  Erythrocytes = c("HSC", "MPP", "MEP", "ERP"),
  Monocytes_Dendritic_cells = c("HSC", "MPP", "MDP"),
  Eo_B_Mast = c("HSC", "MPP", "GMP", 'Eo/B/Mast'))

# Define H1 genes
h1_genes <- c('H1-0', 'H1-1', 'H1-2', 'H1-3', 'H1-4', 'H1-5', 'H1-10')

# Running the gene module analysis for each dataset and lineage
results <- list()
for (dataset_name in names(datasets)) {
  cds_dataset <- datasets[[dataset_name]]
  gene_module_list <- list()
  agg_mat_list <- list()
  lineage_names <- names(lineage_defs)
  for (lineage in lineage_names) {
    lineage_cells <- colData(cds_dataset)$annotation_2 %in% lineage_defs[[lineage]]
    cds_lineage <- cds_dataset[, lineage_cells]
    cds_lineage <- preprocess_cds(cds_lineage, method = "PCA", num_dim = 50)
    
    # Use the appropriate graph test results for each dataset and lineage
    pr_deg_ids <- row.names(subset(pr_results[[dataset_name]][[lineage]], q_value < 0.05))

    # Group genes into modules
    gene_module_df <- find_gene_modules(cds_lineage[pr_deg_ids,], resolution=1e-2)
    cell_group_df <- tibble::tibble(cell=row.names(colData(cds_lineage)), 
                                    cell_group=colData(cds_lineage)$annotation_2)

    # Aggregate the expression of the gene modules to plot them into heatmaps
    agg_mat <- aggregate_gene_expression(cds_lineage, gene_module_df, cell_group_df)
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat)) 
    
    gene_module_list[[lineage]] <- gene_module_df
    agg_mat_list[[lineage]] <- agg_mat

    # Check what modules H1 genes were assigned to
    print(paste("Dataset:", dataset_name, "Lineage:", lineage))
    print(gene_module_df[gene_module_df$id %in% h1_genes,], n=Inf)
  results[[dataset_name]] <- list(
    gene_module_list = gene_module_list,
    agg_mat_list = agg_mat_list
  )
}}
# saveRDS(results, '/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/gene_modules.rds')



#################################
# Check the location of H1 genes#
#################################

h1_genes <- c('H1-0', 'H1-1', 'H1-2', 'H1-3', 'H1-4', 'H1-5', 'H1-10')
# Function to extract H1 gene info for a given dataset and lineage
get_h1_info <- function(dataset, lineage) {
  tab <- results[[dataset]]$gene_module_list[[lineage]]
  tab_h1 <- tab[tab$id %in% h1_genes, ]
  if (nrow(tab_h1) > 0) {
    tab_h1$Dataset <- dataset
    tab_h1$Lineage <- lineage
    return(tab_h1[, c("id", "module", "Dataset", "Lineage")])
  } else {
    return(NULL)
  }
}

# Use lapply over all combinations
h1_module_info_list <- lapply(datasets, function(ds) {
  lapply(lineages, function(ln) get_h1_info(ds, ln))
})

# Flatten and combine into a dataframe
h1_module_info_df <- do.call(rbind, unlist(h1_module_info_list, recursive = FALSE))
colnames(h1_module_info_df)[1] <- "Gene"
print(h1_module_info_df, n=Inf)


#######################################################
# Compute ORA for all the modules containing H1 genes #
#######################################################

ora_results <- list()
str(ora_results)
# Get unique combinations of Dataset, Lineage, and module
  unique_combos <- unique(h1_module_info_df[, c("Dataset", "Lineage", "module")])
unique_combos$names <- paste(unique_combos$Dataset, unique_combos$Lineage, unique_combos$module, sep = "_")

for (i in seq_len(nrow(unique_combos))) {
  ds <- unique_combos$Dataset[i]
  ln <- unique_combos$Lineage[i]
  mod <- unique_combos$module[i]
  # Get all genes in this module
  tab <- results[[ds]]$gene_module_list[[ln]]
  genes_in_module <- tab$id[as.character(tab$module) == as.character(mod)]
  # Run ORA (example: enrichGO, you can change to your preferred function)
  ora <- tryCatch({
    enrichGO(gene = genes_in_module,
             OrgDb = org.Hs.eg.db,
             keyType = "SYMBOL",
             ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2)
  }, error = function(e) NULL)
  ora_results[[paste(ds, ln, mod, sep = "_")]] <- ora
}
names(ora_results) <- unique_combos$names

# saveRDS(ora_results, "/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/H1_related_genes_modules_ORA_results.rds")
ora_results <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/H1_related_genes_modules_ORA_results.rds")


###############################################
########### Plotting the ORA results ##########
###############################################

# Selecting only ORA results relevant for my research question

str(ora_results)
pdf("public_data_analysis/pseudotime/Monocle3/Co-regulation/enrichment_analysis/ORA_plots_CCI_H1_related_modules.pdf")
for (nm in names(ora_results)) {
  res <- ora_results[[nm]]
  if (!is.null(res) && "result" %in% slotNames(res) && nrow(res@result) > 0) {
    tryCatch({
      print(dotplot(res, showCategory = 15) + ggtitle(nm))
    }, error = function(e) {
      message("Skipping ", nm, ": ", e$message)
    })
  }
}
dev.off()
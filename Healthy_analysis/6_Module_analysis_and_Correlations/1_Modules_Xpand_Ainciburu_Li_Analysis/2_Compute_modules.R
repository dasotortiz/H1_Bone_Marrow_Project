library(monocle3)

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
# saveRDS(results, '/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/gene_modules_excluding_histones.rds')

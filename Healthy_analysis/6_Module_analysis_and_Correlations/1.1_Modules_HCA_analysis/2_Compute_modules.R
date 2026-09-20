library(monocle3)

##############################################################
### Excluding other Histone genes from the module analysis ###
##############################################################

cds <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/cds_HCA.rds')

genes_2_exclude <- c(
  "H2BC18", "H2BC21", "H2AC20", "H2AC21", "H3C1", "H4C1", "H4C2",
  "H3C2", "H2AC4", "H2BC3", "H3C3", "H4C3", "H2AC6",
  "H4C4", "H4C5", "H2AC8", "H4C6", "H3C7", "H2BC9", "H3C8",
  "H4C8", "H2BC11", "H2AC11", "H4C9", "H2BC12", "H2AC12", "H2BC13", "H2AC13",
  "H3C10", "H2AC14", "H2BC14", "H2AC15", "H2BC15", "H2AC16", "H3C11",
  "H4C13", "H3C12", "H2AC17", "H2BC17"
)
lineage_names <- names(lineage_defs)

cds <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/cds_HCA.rds')

# Define  datasets and lineages
datasets <- c(HCA = cds)
lineages <- names(lineage_defs)

# Load previous results into a list
outdir <- "/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/pr_results"
pr_results <- list()
for (ds_name in names(datasets)) {
  for (lin_name in lineage_names) {
    file_path <- file.path(outdir, paste0("ds_pr_test_res_", ds_name, "_", lin_name, ".rds"))
    pr_results[[ds_name]][[lin_name]] <- readRDS(file_path)
  }
}

gene_module_list <- list()
agg_mat_list <- list()
lineage_names <- names(lineage_defs)

# Define H1 genes depending on the dataset (h1_genes for Xpand, h1_genes_2 for Ainciburu and Li)
h1_genes <- c('H1-0', 'H1-1', 'H1-2', 'H1-3', 'H1-4', 'H1-5', 'H1-10')

# Find modules of co-regulared genes for each dataset and lineage
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
    pr_deg_ids <- pr_deg_ids[!pr_deg_ids %in% genes_2_exclude]

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
  }
  results[[dataset_name]] <- list(
    gene_module_list = gene_module_list,
    agg_mat_list = agg_mat_list
  )
}
str(results)
# saveRDS(results, '/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/HCA_gene_modules_excluding_histones.rds')

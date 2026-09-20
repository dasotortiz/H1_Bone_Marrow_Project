library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

#################################
# Check the location of H1 genes#
#################################
results <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/gene_modules_excluding_histones.rds')

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
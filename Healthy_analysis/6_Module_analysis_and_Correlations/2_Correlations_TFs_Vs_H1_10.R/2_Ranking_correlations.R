library(dplyr)
library(purrr)
library(stringr)

# ======================================================================================================================================================
# # Compute the rank product across datasets for each H1 × Lineage combination, prioritizing significant correlations and penalizing inconsistent signs
# ======================================================================================================================================================
tf_correlated_df <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/Correlations/early_progenitors_out_TFs.rds')
# Rename lineages for clarity
tf_correlated_df$Lineage <- ifelse(tf_correlated_df$Lineage == "Monocytes_Dendritic_cells", "Mono_DC",
                            ifelse(tf_correlated_df$Lineage == "early_progenitors", "EP", tf_correlated_df$Lineage))

# Define a small epsilon to avoid division by zero in later calculations
epsilon <- 1e-6

ranked_df <- tf_correlated_df %>%
  # 1. Create ranking score prioritizing significant correlations
  mutate(
    sig_priority = ifelse(FDR < 0.05, 0, 1),
    rank_metric = sig_priority * 100 + (1 - abs(Correlation))
  ) %>%
  
  # 2. Rank within H1 × Dataset × Lineage
  group_by(H1, Dataset, Lineage) %>%
  mutate(rank = rank(rank_metric, ties.method = "average")) %>%
  ungroup()

# 3. Compute rank product across datasets
rank_product_df <- ranked_df %>%
  group_by(Gene, H1, Lineage) %>%
  summarise(
    RP = exp(mean(log(rank))),  # geometric mean (rank product)
    n_datasets = n(),
    
    # 4. compute sign consistency
    sign_sum = abs(sum(sign(Correlation))),
    sign_consistency = sign_sum / n_datasets,
    .groups = "drop"
  ) %>%
  
  # 5. apply penalty
  mutate(
    # Add a small value to the sign consistency to avoid division by zero
    final_rank = RP / (sign_consistency + epsilon)
  )

# 6. final ranking within each H1 × Lineage
final_ranking <- rank_product_df %>%
  group_by(H1, Lineage) %>%
  arrange(final_rank) %>%
  mutate(final_rank = rank(final_rank, ties.method = "average")) %>%
  ungroup()

saveRDS(final_ranking, '/ibex/user/sotoorda/masterh1/public_data_analysis/pseudotime/Monocle3/Correlations/rank_product_df.rds')

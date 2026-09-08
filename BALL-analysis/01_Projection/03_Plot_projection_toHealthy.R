# H1 gene expression along the healthy developmental trajectory, reference vs query.
set.seed(1234)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
ref <- readRDS("./data/integration_post_gene_fix.rds")
query <- readRDS("./output/query_predictions_final.rds")

subtype_colors <- c(
  "BCR::ABL1"        = "#8DD3C7",
  "BCR::ABL1-like"   = "#56B4E9",
  "DUX4-R"           = "#0072B2",
  "ETV6::RUNX1-like" = "#CC79A7",
  "Hyperdiploid"     = "#009E73",
  "iAMP21"           = "#E69F00",
  "KMT2A-R"          = "#E31A1C",
  "Low-hypodiploid"  = "#BDBDBD",
  "MEF2D-R"          = "#c6a3ec",
  "Near-haploid"     = "#4B5D73",
  "Other"            = "#8C564B",
  "PAX5alt"          = "#F4A3A8",
  "TCF3::PBX1"       = "#D55E00",
  "ZNF384-R"         = "#7B3294",
  "ETV6::RUNX1"      = "#A65628",
  "Healthy"          = "black",
  "Ainciburu"        = "black",
  "Li"               = "black",
  "Xpand"            = "black"
)

celltypeorder <- c("HSC","MPP","ERP","MEP","MKP","GMP","MDP","Eo/B/Mast","CLP","ProB","PreB")
lymphoid_cells <- c("HSC", "MPP", "CLP", "ProB", "PreB")
ref$annotation_2 <- factor(ref$annotation_2, levels = celltypeorder)
query$predicted.celltype <- factor(query$predicted.celltype, levels = celltypeorder)

# Load H1 genes
# h1_df.rds lives outside this repo, set the path to wherever you keep it
h1_path <- "../h1_df.rds"
h1_df <- readRDS(h1_path)
h1_genes <- h1_df$Alternative_symbol   # naming used in the query
h1_genes2 <- h1_df$Symbol              # naming used in the reference

# Extract expression and metadata
ref <- JoinLayers(ref)
ref_expr <- as.data.frame(t(as.matrix(GetAssayData(ref, layer = "data")[h1_genes2, ])))
ref_expr$CellType <- ref$annotation_2
ref_expr$Dataset <- ref$Dataset
ref_expr$Source <- "Reference"

query_expr <- as.data.frame(t(as.matrix(GetAssayData(query, layer = "data")[h1_genes, ])))
query_expr$CellType <- query$predicted.celltype
query_expr$Dataset <- "Iacobucci"
query_expr$Source <- "Query"
colnames(query_expr) <- colnames(ref_expr)  # query uses the alternative symbols, force the ref naming

# Add LineGroup column for coloring
ref_expr$LineGroup <- ref_expr$Dataset
query_expr$LineGroup <- query$Subtype

# Pivot to long format
ref_long <- ref_expr %>%
  pivot_longer(cols = starts_with("H1-"), names_to = "Gene", values_to = "Expression")
query_long <- query_expr %>%
  pivot_longer(cols = starts_with("H1-"), names_to = "Gene", values_to = "Expression")

combined_long <- rbind(ref_long, query_long)

# Factor levels
combined_long$CellType <- factor(combined_long$CellType, levels = celltypeorder)
combined_long$Gene <- factor(combined_long$Gene, levels = c("H1-0","H1-10","H1-1","H1-2","H1-3","H1-4","H1-5"))
first_levels <- c("Ainciburu","Li","Xpand")
other_levels <- setdiff(unique(combined_long$LineGroup), first_levels)
linegroup_levels <- c(first_levels, sort(other_levels))
combined_long$LineGroup <- factor(combined_long$LineGroup, levels = linegroup_levels)
combined_long$LegendGroup <- combined_long$LineGroup

# Aesthetics: colors and shapes
combined_long <- combined_long %>%
  mutate(
    PlotColor = ifelse(Source=="Reference","black",subtype_colors[as.character(LineGroup)]),
    PointShape = factor(case_when(
      Source=="Reference" & LineGroup=="Ainciburu" ~ "Ainciburu",
      Source=="Reference" & LineGroup=="Li" ~ "Li",
      Source=="Reference" & LineGroup=="Xpand" ~ "Xpand",
      TRUE ~ "Query"
    ), levels=c("Ainciburu","Li","Xpand","Query"))
  )

shape_values <- c("Ainciburu"=15, "Li"=17, "Xpand"=18, "Query"=16)
color_values <- c(subtype_colors,"Ainciburu"="black","Li"="black","Xpand"="black")

# Same line plot for every version below, so it is written once here
h1_lineplot <- function(dat, title) {
  ggplot(dat, aes(x=CellType, y=Expression, group=LineGroup, color=LegendGroup, shape=PointShape)) +
    stat_summary(geom="line", fun=mean, size=0.6) +
    stat_summary(geom="point", fun=mean, size=2) +
    facet_wrap(~Gene, scales="free_y", ncol=2) +
    scale_color_manual(values=color_values, name="Dataset / Subtype") +
    scale_shape_manual(values=shape_values, name="Reference / Query") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          strip.background=element_rect(fill="lightgrey",color=NA)) +
    labs(x="Cell Type", y="Expression", title=title)
}

p1 <- h1_lineplot(combined_long, "H1 Genes Expression Across Cell Types")
ggsave("./output/H1_ref_query.png", p1)

###### only lymphoid lineage
combined_lymphoid <- combined_long %>%
  filter(CellType %in% lymphoid_cells)

p_lymphoid <- h1_lineplot(combined_lymphoid, "H1 Genes Expression Across Lymphoid Lineage")
ggsave("./output/H1_ref_query_lymphoid.png", p_lymphoid)

##### PDF: each page is a BALL subtype, references kept on every page
query_subtypes <- unique(query$Subtype)

pdf("./output/H1_ref_query_by_subtype.pdf")
for (sub in query_subtypes) {
  plot_data <- combined_long %>%
    filter(Source == "Reference" | (Source == "Query" & LineGroup == sub))
  plot_data$CellType <- factor(plot_data$CellType, levels = celltypeorder)
  print(h1_lineplot(plot_data, paste("H1 Genes Expression - Subtype:", sub)))
}
dev.off()

## PDF: only lymphoid lineage thoooooo
pdf("./output/H1_lymphoid_ref_query_by_subtype.pdf")
for (sub in query_subtypes) {
  plot_data <- combined_long %>%
    filter(CellType %in% lymphoid_cells &
           (Source == "Reference" | (Source == "Query" & LineGroup == sub)))
  plot_data$CellType <- factor(plot_data$CellType, levels = lymphoid_cells)
  print(h1_lineplot(plot_data, paste("H1 Genes Expression - Subtype:", sub, "(Lymphoid Lineage)")))
}
dev.off()

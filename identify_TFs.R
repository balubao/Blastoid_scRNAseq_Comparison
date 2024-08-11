## identify TF in DEG lists
# acquire TF list from https://www.tfcheckpoint.org/index.php/download-data
TFDB.df = read_excel('data/TFC2_16102023b.xlsx')
# TFDB = read.table('data/TFC2_16102023b.tsv', sep = '\t', header = TRUE, fill = TRUE)
TFDB.df

# load DEG list
aggregated_ranks_celltype_pe.ls = readRDS('06_deg_analysis/results/aggregated_rank_gene_celltype_pe_list.rds')
aggregated_ranks_celltype_pe.ls = lapply(aggregated_ranks_celltype_pe.ls, function(x)x[x$Score<0.01,])

## Assign TF label
aggregated_ranks_celltype_pe_tf.ls = lapply(aggregated_ranks_celltype_pe.ls, function(aggregated_ranks_celltype_pe){
  aggregated_ranks_celltype_pe$TF = aggregated_ranks_celltype_pe$Name %in% TFDB.df$Associated.Gene.Name
  return(aggregated_ranks_celltype_pe)
})

lapply(aggregated_ranks_celltype_pe_tf.ls, function(x) table(x$TF[x$Score < 0.01]))
lapply(aggregated_ranks_celltype_pe_tf.ls, function(x) head(x))

# Write each dataframe to a separate sheet in an Excel file
# Load necessary library
if (!require(writexl)) install.packages("writexl")
library(writexl)
write_xlsx(aggregated_ranks_celltype_pe_tf.ls, "06_deg_analysis/results/aggregated_rank_gene_celltype_pe_TFs.xlsx")

################################################################################
# load DEG list
aggregated_rank_gene_samplelc.ls = readRDS('06_deg_analysis/results/aggregated_rank_gene_samplelc_list.rds')
aggregated_rank_gene_samplelc.ls = lapply(aggregated_rank_gene_samplelc.ls, function(x)x[x$Score<0.01,])
names(aggregated_rank_gene_samplelc.ls) = paste0("cluster_", names(aggregated_rank_gene_samplelc.ls))

## Assign TF label
aggregated_rank_gene_samplelc_tf.ls = lapply(aggregated_rank_gene_samplelc.ls, function(aggregated_rank_gene_samplelc){
  aggregated_rank_gene_samplelc$TF = aggregated_rank_gene_samplelc$Name %in% TFDB.df$Associated.Gene.Name
  return(aggregated_rank_gene_samplelc)
})

lapply(aggregated_rank_gene_samplelc_tf.ls, function(x) table(x$TF[x$Score < 0.01]))
lapply(aggregated_rank_gene_samplelc_tf.ls, function(x) head(x))
lapply(aggregated_rank_gene_samplelc_tf.ls, function(x) head(x[x$TF,]))

# Write each dataframe to a separate sheet in an Excel file
# Load necessary library
if (!require(writexl)) install.packages("writexl")
library(writexl)
write_xlsx(aggregated_rank_gene_samplelc_tf.ls, "06_deg_analysis/results/aggregated_rank_gene_samplelc_TFs.xlsx")
################################################################################
# load DEG list
aggregated_ranks_cluster.ls = readRDS('06_deg_analysis/results/aggregated_rank_gene_rpca_clust_list.rds')
aggregated_ranks_cluster.ls = lapply(aggregated_ranks_cluster.ls, function(x)x[x$Score<0.01,])
names(aggregated_ranks_cluster.ls) = paste0("cluster_", names(aggregated_ranks_cluster.ls))

## Assign TF label
aggregated_ranks_cluster_tf.ls = lapply(aggregated_ranks_cluster.ls, function(aggregated_rank_gene_df){
  aggregated_rank_gene_df$TF = aggregated_rank_gene_df$Name %in% TFDB.df$Associated.Gene.Name
  return(aggregated_rank_gene_df)
})

lapply(aggregated_ranks_cluster_tf.ls, function(x) table(x$TF[x$Score < 0.01]))
lapply(aggregated_ranks_cluster_tf.ls, function(x) head(x))
lapply(aggregated_ranks_cluster_tf.ls, function(x) head(x[x$TF,]))

aggregated_ranks_cluster_tf.ls = aggregated_ranks_cluster_tf.ls[sapply(aggregated_ranks_cluster_tf.ls, nrow) > 0]

# Write each dataframe to a separate sheet in an Excel file
# Load necessary library
if (!require(writexl)) install.packages("writexl")
library(writexl)
write_xlsx(aggregated_ranks_cluster_tf.ls, "06_deg_analysis/results/aggregated_rank_gene_rpca_clust_TFs.xlsx")

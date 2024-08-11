# similarity between local clusters


seurat.list$Petropoulos2016_ref$seurat_clusters

seurat.list = lapply(seurat.list, function(x) {
  Idents(x) = x$seurat_clusters
  return(x)})
pseudobulk.list = lapply(seurat.list, AggregateExpression)

# use shared genes across datasets
keep_genes = lapply(pseudobulk.list, function(x) rownames(x$RNA))
length(Reduce(union, keep_genes))
length(Reduce(intersect, keep_genes))
keep_genes = Reduce(intersect, keep_genes)
pseudobulk_comparable.list = lapply(pseudobulk.list, function(x) x$RNA = x$RNA[rownames(x$RNA) %in% keep_genes,])
pseudobulk.mat = do.call('cbind', pseudobulk_comparable.list)
colnames(pseudobulk.mat) = make.unique(rep(names(pseudobulk.list), times = sapply(pseudobulk_comparable.list, ncol)))

geneexp_cor_mat = cor(as.matrix(pseudobulk.mat), method = 'spearman')
cor_plot(geneexp_cor_mat, tl.cex = 0.4)

################################################################################
## Try again with integrated embeddings
for(i in seq_along(seurat.list)){
  select_c = plotting.data$sample == names(seurat.list)[i]
  plotting.data$seurat_clusters[select_c] = seurat.list[[i]]$seurat_clusters
}
#create rpca_embedding dataframe
integrated_rpca.df = split_obj@reductions$integrated.rpca@cell.embeddings %>% as.data.frame()
dim(integrated_rpca.df)
# create celltype dataframe
integrated_rpca.df$sample_clusters = paste0(split_obj$sample, '_', split_obj$seurat_clusters)
integrated_rpca_celltypes.df = data.frame(sample_clusters = integrated_rpca.df$sample_clusters,
                                          celltype_pe = plotting.data$celltype_pe)
# generate pseudobulk
integrated_rpca_pseudobulk.df = integrated_rpca.df %>%
  group_by(sample_clusters) %>%
  summarise(across(colnames(integrated_rpca.df)[!colnames(integrated_rpca.df) %in% c('sample_clusters')], sum))
  # column_to_rownames(var = "sample_clusters")

# Assign a celltype label to each sample_cluster based on the majority celltype
majority_celltype_df <- integrated_rpca_celltypes.df %>%
  group_by(sample_clusters) %>%
  summarize(celltype_pe = names(sort(table(celltype_pe), decreasing = TRUE))[1])

# Merge the celltype labels with the pseudobulk dataframe
pseudobulk_df <- integrated_rpca_pseudobulk.df %>%
  left_join(majority_celltype_df, by = "sample_clusters")

pseudobulk_df

# Compute similarity matrix (e.g., correlation matrix) between pseudobulks
pseudobulk_matrix <- pseudobulk_df %>%
  column_to_rownames('sample_clusters') %>%
  dplyr::select(-celltype_pe) %>%
  as.matrix()

similarity_matrix <- cor(t(pseudobulk_matrix), method = 'spearman')

similarity_matrix

# Load necessary libraries for visualization
if (!require(pheatmap)) install.packages("pheatmap")
library(pheatmap)

# Step 4: Visualize clustered heatmap and annotate with celltype and sample labels
# Prepare annotation for the heatmap
annotation <- data.frame(celltype = pseudobulk_df$celltype_pe, sample = sub("_[^_]+$", "", pseudobulk_df$sample_clusters))
rownames(annotation) <- pseudobulk_df$sample_clusters


annotation_colors <- list(celltype = celltype_pe_color)

# Visualize clustered heatmap
pheatmap(similarity_matrix,
         annotation_col = annotation,
         annotation_row = annotation,
         annotation_colors = annotation_colors,
         main = "Clustered Heatmap of Pseudobulk Similarity")

################################################################################
sample_clusters_dist = as.dist(1 - similarity_matrix)
sample_clusters_hc = hclust(sample_clusters_dist)
sim_clusters = cutree(sample_clusters_hc, 10) #number of datasets

# Assign the new cluster labels to each sample_cluster
pseudobulk_df <- pseudobulk_df %>%
  mutate(new_cluster = sim_clusters[as.character(sample_clusters)])

# Step 5: Propagate the new cluster labels back to the celltype dataframe
integrated_rpca_celltypes.df <- integrated_rpca_celltypes.df %>%
  left_join(pseudobulk_df %>% dplyr::select(sample_clusters, new_cluster), by = "sample_clusters")

# View the updated celltype_df with new cluster labels
plotting.data$samplelc_clusters = integrated_rpca_celltypes.df$new_cluster
table(plotting.data$samplelc_clusters)
plotting.data$samplelc_clusters = factor(plotting.data$samplelc_clusters)

ggplot()+
  geom_point(aes(x=UMAP_1, y=UMAP_2, color = samplelc_clusters), size = 0.5, plotting.data)+
  geom_label_repel(aes(x=UMAP_1, y=UMAP_2, label=samplelc_clusters), plotting.data %>% group_by(samplelc_clusters) %>% summarise(across(c('UMAP_1','UMAP_2'), mean)))

heatmap(table(plotting.data$rpca_clusters, plotting.data$samplelc_clusters))
heatmap(table(plotting.data$sample, plotting.data$samplelc_clusters))
heatmap(table(plotting.data$celltype_pe, plotting.data$samplelc_clusters))

################################################################################

################################################################################
# transfer cluster labels to each object (from Sample Similarity clusters)
# i=1
seurat.list = lapply(seurat.list, function(seurat.obj) {
  # seurat.obj = seurat.list[[i]]
  plotting.data_sub = plotting.data[plotting.data$sample == unique(seurat.obj$sample),]
  seurat.obj$clusters = plotting.data_sub$samplelc_clusters
  return(seurat.obj)
})


split_obj$celltype_pe = plotting.data$celltype_pe
split_obj$samplelc_clusters = plotting.data$samplelc_clusters
DimPlot(split_obj, reduction = 'umap.rpca', group.by = 'samplelc_clusters')#QC

# get DEGs for each dataset using integrated.clusters  (RPCA)
topmarkers.all.ls = lapply(seurat.list, function(seurat.obj) {
  Idents(seurat.obj) = seurat.obj$clusters
  top_markers = FindAllMarkers(seurat.obj)
  return(top_markers)
})

# saveRDS(topmarkers.all.ls, "./06_deg_analysis/results/wilcox_deg_list_samplelcclust_unfiltered_all.rds")
topmarkers.all.ls = readRDS("./06_deg_analysis/results/wilcox_deg_list_samplelcclust_unfiltered_all.rds")
################################################################################

# label each df by dataset
topmarkers.all.ls = lapply(seq_along(topmarkers.all.ls), function(i) {
  topmarkers.all.ls[[i]]$sample = names(topmarkers.all.ls)[i]
  return(topmarkers.all.ls[[i]])
})
names(topmarkers.all.ls) = sapply(topmarkers.all.ls, function(x)
  unique(x$sample))

# 1. keep genes appearing in at least 2 datasets
cluster_levels = sort(as.numeric(unique(do.call('c',lapply(topmarkers.all.ls, function(x) levels(x$cluster))))))
keep_genes.ls = lapply(cluster_levels, function(i)
{
  cluster_genes = lapply(topmarkers.all.ls, function(topmarkers) {
    gene_names = topmarkers$gene[topmarkers$cluster == i]
    return(gene_names)
  })
  cluster_genes = do.call("c", cluster_genes)
  keep_genes = cluster_genes[duplicated(cluster_genes)]
  return(keep_genes)
})
names(keep_genes.ls) = cluster_levels
sapply(keep_genes.ls, length)

n_dup_genes = sapply(keep_genes.ls, length)
plotting.data$w_dup_genes = plotting.data$samplelc_clusters %in% names(n_dup_genes)
ggplot()+
  geom_point(aes(x=UMAP_1, y=UMAP_2, color=w_dup_genes), plotting.data)+
  geom_label_repel(aes(x=UMAP_1, y=UMAP_2, label=samplelc_clusters), plotting.data %>% group_by(samplelc_clusters) %>% summarise(across(c('UMAP_1','UMAP_2'), mean)))

# 2. filter DEG tables to only retain duplicated genes
topmarkers.all_dupfilt.ls = lapply(topmarkers.all.ls, function(topmarkers.all)
{
  topmarkers.all_sub.ls = lapply(seq_along(cluster_levels), function(i) {
    topmarkers = topmarkers.all[topmarkers.all$cluster == cluster_levels[i],]
    topmarkers = topmarkers[topmarkers$gene %in% keep_genes.ls[[i]],]
    return(topmarkers)
  })

  topmarkers.all_sub = do.call("rbind", topmarkers.all_sub.ls)
  return(topmarkers.all_sub)
})

sapply(topmarkers.all_dupfilt.ls, dim)
sapply(topmarkers.all_dupfilt.ls, function(x)
  table(x$cluster))

# 3. filter duplicated DEG by LFC>.25 and pvalue<.01
topmarkers.all_filt_dup.ls = lapply(topmarkers.all_dupfilt.ls, function(topmarkers.all)
{
  topmarkers.all = topmarkers.all[(topmarkers.all$avg_log2FC > 0.25) &
                                    (topmarkers.all$p_val_adj <= 0.01), ]
  return(topmarkers.all)
})

sapply(topmarkers.all_filt_dup.ls, dim)
sapply(topmarkers.all_filt_dup.ls, function(x)
  table(x$cluster))

# filter DEG by LFC>.25 and pvalue<.01
topmarkers.all_filt.ls = lapply(topmarkers.all.ls, function(topmarkers.all)
{
  topmarkers.all = topmarkers.all[(topmarkers.all$avg_log2FC>0.25) & (topmarkers.all$p_val_adj<=0.01),]
  return(topmarkers.all)
})

sapply(topmarkers.all_filt.ls, dim)
sapply(topmarkers.all_filt.ls, function(x)
  table(x$cluster))

# # filter universal DEG by LFC>.25 and pvalue<.01
# lapply(topmarkers.all.ls, function(x)
#   table(x$cluster))
# tmp_df = do.call('rbind',topmarkers.all.ls)
# tmp_df = tmp_df %>%
#   group_by(cluster, gene) %>%
#   count() %>%
#   filter(n==10)
# table(tmp_df$cluster)
# universal_genes = lapply(topmarkers.all.ls, function(topmarkers.all) topmarkers.all$gene)
# universal_genes = Reduce(intersect, universal_genes)
# topmarkers.all_filt.ls = lapply(topmarkers.all.ls, function(topmarkers.all)
# {
#   topmarkers.all = topmarkers.all[(topmarkers.all$avg_log2FC>0.25) & (topmarkers.all$p_val_adj<=0.01),]
#   topmarkers.all = topmarkers.all[topmarkers.all$gene %in% universal_genes,]
#   return(topmarkers.all)
# })
# 
# sapply(topmarkers.all_filt.ls, dim)
# sapply(topmarkers.all_filt.ls, function(x)
#   table(x$cluster))

################################################################################
# Function to merge data frames on Var1 and fill missing values with 0
merge_dfs <- function(df_list) {
  merged_df <- reduce(df_list, function(x, y) {
    full_join(x, y, by = "Var1")
  })
  merged_df[is.na(merged_df)] <- 0
  return(merged_df)
}
GetNGenesPerDatasetPerClust = function(topmarkers.all.ls) {
  ndegs_per_clust.ls = lapply(topmarkers.all.ls, function(x)
    as.data.frame(table(x$cluster)))
  names(ndegs_per_clust.ls) = sapply(topmarkers.all.ls, function(x)
    unique(x$sample))
  
  # Find all unique Var1 values across all data frames
  all_vars <-
    unique(unlist(lapply(ndegs_per_clust.ls, function(df)
      df$Var1)))
  
  # Add all unique Var1 values to each data frame with missing Freq values set to 0
  filled_dfs <- lapply(ndegs_per_clust.ls, function(df) {
    df <- full_join(data.frame(Var1 = all_vars), df, by = "Var1")
    df[is.na(df)] <- 0
    return(df)
  })
  
  # Merge all data frames
  combined_df = do.call("cbind", lapply(filled_dfs, function(x)
    x[, 2]))
  rownames(combined_df) = filled_dfs$Petropoulos2016_ref[, 1]
  combined_df = as.data.frame(combined_df)
  # combined_df <- filled_dfs %>%
  #   reduce(function(x, y)
  #     full_join(x, y, by = "Var1", suffix = c(".x", ".y"))) %>%
  #   replace(is.na(.), 0)
  
  # Rename columns for clarity (optional)
  # colnames(combined_df) <- c("Var1", names(ndegs_per_clust.ls))
  
  # View the combined data frame
  return(combined_df)
  
}
# Function to create the merged data frame with "value1/value2" format
merge_data_frames <- function(df1, df2) {
  merged_df <- df1 %>%
    mutate(across(everything(), ~ paste0(.x, "/", df2[[cur_column()]])))
  return(merged_df)
}

# only duplicates versus all detected DEGs
merge_data_frames(
  GetNGenesPerDatasetPerClust(topmarkers.all_dupfilt.ls),
  GetNGenesPerDatasetPerClust(topmarkers.all.ls)
)

# only DEGs with LFC>0.25 and pval<0.01 versus all detected DEGs
merge_data_frames(
  GetNGenesPerDatasetPerClust(topmarkers.all_filt.ls),
  GetNGenesPerDatasetPerClust(topmarkers.all.ls)
)

# only duplicated DEGs with LFC>0.25 and pval<0.01 versus all detected DEGs
merge_data_frames(
  GetNGenesPerDatasetPerClust(topmarkers.all_filt_dup.ls),
  GetNGenesPerDatasetPerClust(topmarkers.all.ls)
)


################################################################################
# rank genes by LFC for each cluster
topmarkers.all_ranked.ls = lapply(topmarkers.all_filt_dup.ls, function(topmarkers) {
  topmarkers = topmarkers %>%
    group_by(cluster) %>%
    mutate(rank = rank(-avg_log2FC, ties.method = "average")) #rank from largest to smallest
  
  return(topmarkers)
})

CheckGeneDist = function(topmarkers.clust.df) {
  require(dplyr)
  pre_m = topmarkers.clust.df %>%
    dplyr::select(sample, gene) %>%
    split(.$sample) %>%
    map( ~ .x %>% dplyr::select(-sample))
  pre_m = sapply(pre_m, function(x)
    x$gene)
  m = make_comb_mat(pre_m)
  UpSet(m)
}

# define variables
topmarkers.all_ranked.df = do.call("rbind", topmarkers.all_ranked.ls)
clusters = levels(topmarkers.all_ranked.df$cluster)

CheckGeneDist(topmarkers.all_ranked.df)

require(RobustRankAggreg)
aggregated_ranks.ls = lapply(seq_along(clusters), function(i) {
  print(i)
  topmarkers.clust.df = topmarkers.all_ranked.df[topmarkers.all_ranked.df$cluster == clusters[i], ]
  pre_ranked_list = topmarkers.clust.df %>%
    dplyr::select(sample, gene, rank, p_val_adj) %>%
    split(.$sample) %>%
    map( ~ .x %>% dplyr::select(-sample))
  
  ranked_list = lapply(pre_ranked_list, function(x)
    x$gene[order(x$rank)])
  
  # Aggregate ranks
  aggregated_ranks <- aggregateRanks(ranked_list)
  
  return(aggregated_ranks)
})
names(aggregated_ranks.ls) = clusters
# CheckGeneDist(topmarkers.clust.df)

# sanity checks
sapply(aggregated_ranks.ls, nrow) #number of ranked genes
sapply(aggregated_ranks.ls, function(x)
  sum(x$Score < 0.01)) #number of ranked genes w/ p-value>0.01

n_sig_aggr_genes = sapply(aggregated_ranks.ls, function(x)
  sum(x$Score < 0.01)) #number of ranked genes w/ p-value>0.01
plotting.data$w_sig_aggr_genes = plotting.data$samplelc_clusters %in% names(n_sig_aggr_genes)[n_sig_aggr_genes>0]
ggplot()+
  geom_point(aes(x=UMAP_1, y=UMAP_2, color=w_sig_aggr_genes), plotting.data)+
  geom_label_repel(aes(x=UMAP_1, y=UMAP_2, label=samplelc_clusters), plotting.data %>% group_by(samplelc_clusters) %>% summarise(across(c('UMAP_1','UMAP_2'), mean)))

saveRDS(aggregated_ranks.ls, "./06_deg_analysis/results/aggregated_rank_gene_samplelc_list.rds")
################################################################################

aggregated_ranks.ls = lapply(aggregated_ranks.ls, function(genes_df){
  genes_df = genes_df[genes_df$Name %in% rownames(split_obj),]
  return(genes_df)
})

# sort by score (decreasing)
aggregated_ranks.ls = lapply(aggregated_ranks.ls, function(x) x[order(x[,2], decreasing = FALSE),])

TopN = 5

#select top N genes
aggregated_topN.ls = lapply(aggregated_ranks.ls, function(x) x[seq(TopN),])

aggregated_topN = do.call("rbind",aggregated_topN.ls)
aggregated_topN$cluster = rep(seq(length(aggregated_topN.ls)), each = TopN)

# visualize by heatmap (mean over cluster vs dataset)
split_obj = ScaleData(split_obj, features = aggregated_topN$Name, vars.to.regress = "sample")
# DoHeatmap(split_obj, features = aggregated_topN$Name)
# p1 = FeaturePlot(split_obj, features = aggregated_topN$Name,
#             reduction = "umap.rpca", order = TRUE) &
#   theme_minimal() +
#   theme(
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     legend.position = "none"
#   )


# split_obj = ScaleData(split_obj, features = aggregated_celltype_pe_topN$Name)
# DoHeatmap(split_obj, features = aggregated_celltype_pe_topN$Name)

PlotHeatmap_samplelc_cluster = function(aggregated_cluster_topN = aggregated_cluster_topN) {
  # library(groupdata2)
  set.seed(2024)
  split_obj = ScaleData(split_obj, features = aggregated_cluster_topN$Name)
  heatmap_data = split_obj@assays$RNA$scale.data[aggregated_cluster_topN$Name, ]
  sampling_df = split_obj@meta.data[, c("sample", "celltype_pe", "samplelc_clusters", "rpca_clusters")]
  sampling_df$sample = factor(sampling_df$sample)
  sampling_df$sample = fct_inorder(sampling_df$sample)
  
  sampling_df$idx = seq(nrow(sampling_df))
  sampling_df = sampling_df %>%
    group_by(samplelc_clusters, sample, rpca_clusters) %>%
    sample_n(20, replace = TRUE)
  subsample_idx = sampling_df$idx
  # subsample_idx = sample(seq(ncol(heatmap_data)), size = 5000)
  subsample_idx = subsample_idx[order(split_obj$samplelc_clusters[subsample_idx])]
  subsample_idx = subsample_idx[order(split_obj$sample[subsample_idx])]
  sampling_df = sampling_df[match(subsample_idx, sampling_df$idx), ]
  
  samplelc_cluster_color = pal_brewer(type = 'qual', palette = 'Set3')(nlevels(sampling_df$samplelc_clusters))
  names(samplelc_cluster_color) = levels(sampling_df$samplelc_clusters)
  
  heatmap_data_sub = heatmap_data[, subsample_idx]
  
  row_ha = rowAnnotation(clusters = aggregated_cluster_topN$cluster,
                         col = list(clusters = samplelc_cluster_color))
  col_ha = HeatmapAnnotation(
    # dataset = sampling_df$sample,
    clusters = sampling_df$rpca_clusters,
    celltype = sampling_df$celltype_pe,
    # clusters = sampling_df$samplelc_clusters,
    col = list(
      # clusters = samplelc_cluster_color,
      celltype = celltype_pe_color
      ))
  
  h1 = Heatmap(
    heatmap_data_sub,
    name = "Expression",
    top_annotation = col_ha,
    right_annotation = row_ha,
    column_split = split_obj$rpca_clusters[subsample_idx],
    show_column_names = FALSE,
    cluster_columns = TRUE,
    show_column_dend = FALSE,
    cluster_rows = FALSE,
    use_raster = FALSE
  ) %>%
    draw()
}
PlotHeatmap_samplelc_cluster(aggregated_topN)


PlotHeatmap_cluster = function(aggregated_cluster_topN = aggregated_cluster_topN) {
  # library(groupdata2)
  set.seed(2024)
  split_obj = ScaleData(split_obj, features = aggregated_cluster_topN$Name)
  heatmap_data = split_obj@assays$RNA$scale.data[aggregated_cluster_topN$Name, ]
  sampling_df = split_obj@meta.data[, c("sample", "rpca_clusters")]
  sampling_df$sample = factor(sampling_df$sample)
  sampling_df$sample = fct_inorder(sampling_df$sample)

  sampling_df$idx = seq(nrow(sampling_df))
  sampling_df = sampling_df %>%
    group_by(rpca_clusters, sample) %>%
    sample_n(20, replace = TRUE)
  subsample_idx = sampling_df$idx
  # subsample_idx = sample(seq(ncol(heatmap_data)), size = 5000)
  subsample_idx = subsample_idx[order(split_obj$rpca_clusters[subsample_idx])]
  subsample_idx = subsample_idx[order(split_obj$sample[subsample_idx])]
  sampling_df = sampling_df[match(subsample_idx, sampling_df$idx), ]


  heatmap_data_sub = heatmap_data[, subsample_idx]

  row_ha = rowAnnotation(clusters = aggregated_cluster_topN$cluster)
  col_ha = HeatmapAnnotation(
    dataset = sampling_df$sample,
    clusters = sampling_df$rpca_clusters)
  # row_ha@anno_list$celltype@color_mapping = col_ha@anno_list$celltype@color_mapping

  h1 = Heatmap(
    heatmap_data_sub,
    name = "Expression",
    top_annotation = col_ha,
    right_annotation = row_ha,
    column_split = split_obj$rpca_clusters[subsample_idx],
    show_column_names = FALSE,
    cluster_columns = TRUE,
    show_column_dend = FALSE,
    cluster_rows = FALSE,
    use_raster = FALSE
  ) %>%
    draw()
}
PlotHeatmap_cluster(aggregated_topN)

# show celltype specifc gene expression across clusters with >.5 unknown cells
PlotHeatmap_cluster = function(aggregated_cluster_topN = aggregated_cluster_topN) {
  set.seed(2024)
  
  # get expression scaled data
  split_obj = ScaleData(split_obj, features = aggregated_celltype_pe_topN$Name)
  heatmap_data = split_obj@assays$RNA$scale.data[aggregated_celltype_pe_topN$Name, ]
  
  # retrieve cluster and celltype information
  sampling_df = split_obj@meta.data[, c("sample", "rpca_clusters", "celltype_pe")]
  sampling_df$sample = factor(sampling_df$sample)
  sampling_df$sample = fct_inorder(sampling_df$sample)
  
  #filter out clusters with <.5 unknown cells
  sampling_df_coarse = sampling_df %>%
    group_by(rpca_clusters, celltype_pe) %>%
    count() %>% #how many cells of type in cluster
    ungroup(celltype_pe) %>%
    mutate(fraction = n/sum(n)) %>% #fraction of celltype in cluster
    filter((celltype_pe == "unknown") & (fraction >=0.5)) #filter clusters with unknown>.5
  
  sampling_df$idx = seq(nrow(sampling_df)) #assign index to full df
  sampling_df = sampling_df[sampling_df$rpca_clusters %in% sampling_df_coarse$rpca_clusters,] #filter for clusters of interest
  
  sampling_df = sampling_df %>% #to reduce load on plotting, subsample n cells from each dataset for each cluster
    group_by(rpca_clusters, sample) %>%
    sample_n(30, replace = TRUE)
  
  subsample_idx = sampling_df$idx #retrieve sampling index
  # subsample_idx = sample(seq(ncol(heatmap_data)), size = 5000)
  subsample_idx = subsample_idx[order(split_obj$rpca_clusters[subsample_idx])]
  subsample_idx = subsample_idx[order(split_obj$sample[subsample_idx])]
  sampling_df = sampling_df[match(subsample_idx, sampling_df$idx), ]
  
  heatmap_data_sub = heatmap_data[, subsample_idx]
  
  row_ha = rowAnnotation(celltype = aggregated_celltype_pe_topN$cluster,
                         col = list(celltype = celltype_pe_color))
  col_ha = HeatmapAnnotation(
    dataset = sampling_df$sample,
    celltype = sampling_df$celltype_pe,
    col = list(celltype = celltype_pe_color)
  )
  
  h1 = Heatmap(
    heatmap_data_sub,
    name = "Expression",
    top_annotation = col_ha,
    right_annotation = row_ha,
    column_split = split_obj$rpca_clusters[subsample_idx],
    show_column_names = FALSE,
    cluster_columns = TRUE,
    show_column_dend = FALSE,
    cluster_rows = FALSE,
    use_raster = FALSE
  ) %>%
    draw()
}
PlotHeatmap_cluster(aggregated_cluster_topN)

p2 = FeaturePlot(split_obj, features = aggregated_topN$Name,
                 reduction = "umap.rpca", order = TRUE) &
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )
p2

# DimPlot(split_obj, reduction = "umap.rpca", group.by = "celltype_pe") &
#   theme_minimal() +
#   theme(
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     legend.position = "none"
#   )

















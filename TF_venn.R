excel_sheet_reader <- function(filename) {
  sheets <- excel_sheets(filename)
  x <- lapply(sheets, function(X) read_excel(filename, sheet = X))
  names(x) <- sheets
  x
}
aggregated_ranks_celltype_pe_tf.ls = excel_sheet_reader("06_deg_analysis/results/aggregated_rank_gene_celltype_pe_TFs.xlsx")
aggregated_rank_gene_samplelc_tf.ls = excel_sheet_reader("06_deg_analysis/results/aggregated_rank_gene_samplelc_TFs.xlsx")

library("ggVennDiagram")
# Default plot
ggVennDiagram(lapply(aggregated_ranks_celltype_pe_tf.ls, function(x) x$Name))
ggVennDiagram(lapply(aggregated_rank_gene_samplelc_tf.ls, function(x) x$Name))



plts_tf_venn = lapply(seq_along(aggregated_rank_gene_samplelc_tf.ls), function(i) {
  inp_tf_list = c(
    lapply(aggregated_ranks_celltype_pe_tf.ls[c(2, 4, 5, 1)], function(x)
      x$Name),
    list(aggregated_rank_gene_samplelc_tf.ls[order(as.numeric(names(aggregated_rank_gene_samplelc_tf.ls)))][[i]]$Name)
  )
  names(inp_tf_list)[length(inp_tf_list)] = names(aggregated_rank_gene_samplelc_tf.ls)[order(as.numeric(names(aggregated_rank_gene_samplelc_tf.ls)))][i]
  p_cell_v_clust = ggVennDiagram(inp_tf_list,
                                 set_color = c("black", "black", "black", "black", "red"), label = "count", label_alpha = 0) + 
    scale_x_continuous(expand = expansion(mult = .2)) + 
    scale_fill_distiller(palette = "Reds", direction = 1)
  p_cell_v_clust
})
plts_tf_venn[[1]]
wrap_plots(plts_tf_venn, ncol = 3)
ggsave('06_deg_analysis/TF_celltype_pe_vs_samplelc_venn.jpeg', height = 16, width = 20)

aggregated_rank_gene_samplelc_tf_only.ls = lapply(aggregated_rank_gene_samplelc_tf.ls, function(x) x[x$TF,])
aggregated_ranks_celltype_pe_tf_only.ls = lapply(aggregated_ranks_celltype_pe_tf.ls, function(x) x[x$TF,])

plts_tf = lapply(seq_along(aggregated_rank_gene_samplelc_tf_only.ls), function(j){
  plts_tf_single = lapply(c(2,4,5,1), function(i){ 
    inp_tf_list = list(Cluster = aggregated_rank_gene_samplelc_tf_only.ls[[j]]$Name, Celltype = aggregated_ranks_celltype_pe_tf_only.ls[[i]]$Name)
    names(inp_tf_list) = c(names(aggregated_rank_gene_samplelc_tf_only.ls)[j], names(aggregated_ranks_celltype_pe_tf_only.ls)[i])
    plt_tf_single = ggVennDiagram(inp_tf_list, label_alpha = 0) + scale_x_continuous(expand = expansion(mult = .2))
    return(plt_tf_single)
  })
})
length(plts_tf)
wrap_plots(plts_tf[[1]])
plts_tf_wrapped = lapply(plts_tf, wrap_plots) 

# test
plts_tf[[1]][[1]]+ 
  theme(
    plot.title = element_blank()
  )
################################################################################
# Function to calculate percentage overlap
# Function to calculate percentage overlap, including non-overlapping TFs
calculate_overlap <- function(cluster_tfs, celltype_tfs) {
  # Extract the TF names from the cluster and cell type dataframes
  cluster_tf_names <- lapply(cluster_tfs, function(df) df$Name[df$TF])
  celltype_tf_names <- lapply(celltype_tfs, function(df) df$Name[df$TF])
  
  # Initialize a list to store the results
  overlap_results <- list()
  
  # Calculate the percentage overlap for each cluster with each cell type
  for (cluster_name in names(cluster_tf_names)) {
    overlaps <- sapply(celltype_tf_names, function(celltype_tfs) {
      cluster_tfs <- cluster_tf_names[[cluster_name]]
      length(intersect(cluster_tfs, celltype_tfs)) / length(cluster_tfs) * 100
    })
    
    # Calculate non-overlapping TFs
    total_tfs <- length(cluster_tf_names[[cluster_name]])
    overlapping_tfs <- length(unique(unlist(lapply(celltype_tf_names, function(celltype_tfs) {
      intersect(cluster_tf_names[[cluster_name]], celltype_tfs)
    }))))
    non_overlapping_percentage <- (total_tfs - overlapping_tfs) / total_tfs * 100
    
    # Store results
    overlap_results[[cluster_name]] <- c(overlaps, None = non_overlapping_percentage)
  }
  
  # Convert the results to a data frame for easier interpretation
  overlap_df <- do.call(rbind, lapply(names(overlap_results), function(cluster_name) {
    do.call(rbind, lapply(names(overlap_results[[cluster_name]]), function(celltype_name) {
      data.frame(Cluster = cluster_name, CellType = celltype_name, OverlapPercentage = overlap_results[[cluster_name]][[celltype_name]])
    }))
  }))
  
  return(overlap_df)
}

# Calculate the overlap percentages
overlap_df <- calculate_overlap(aggregated_rank_gene_samplelc_tf.ls[order(as.numeric(names(aggregated_rank_gene_samplelc_tf.ls)))], 
                                aggregated_ranks_celltype_pe_tf.ls[c(2,4,5,1)])

# Print the results
print(overlap_df)
overlap_df$CellType = factor(overlap_df$CellType)
overlap_df$CellType = fct_inorder(overlap_df$CellType)
overlap_df$Cluster = factor(overlap_df$Cluster)
overlap_df$Cluster = fct_inorder(overlap_df$Cluster)

# Visualize the results with a heatmap
heatmap_plot <- ggplot(overlap_df, aes(x = Cluster, y = CellType, fill = OverlapPercentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", OverlapPercentage)), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  labs(title = "Percentage Overlap of TFs Between Similarity Clusters and Cell Types",
       x = "Cluster",
       y = "Cell Type",
       fill = "Overlap (%)") +
  theme_minimal() +
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the heatmap
print(heatmap_plot)

################################################################################
################################################################################
aggregated_ranks_cluster_tf.ls = excel_sheet_reader("06_deg_analysis/results/aggregated_rank_gene_rpca_clust_TFs.xlsx")
names(aggregated_ranks_cluster_tf.ls) = as.character(seq_along(aggregated_ranks_cluster_tf.ls)-1)

aggregated_ranks_cluster_tf.ls = aggregated_ranks_cluster_tf.ls[sapply(aggregated_ranks_cluster_tf.ls, function(x) sum(x$TF)>0)]

library("ggVennDiagram")
aggregated_ranks_cluster_tf.ls = aggregated_ranks_cluster_tf.ls[sapply(aggregated_ranks_cluster_tf.ls, nrow)>0]
plts_tf_venn = lapply(seq_along(aggregated_ranks_cluster_tf.ls), function(i) {
  inp_tf_list = c(
    lapply(aggregated_ranks_celltype_pe_tf.ls[c(2, 4, 5, 1)], function(x)
      x$Name),
    list(aggregated_ranks_cluster_tf.ls[order(as.numeric(names(aggregated_ranks_cluster_tf.ls)))][[i]]$Name)
  )
  names(inp_tf_list)[length(inp_tf_list)] = names(aggregated_ranks_cluster_tf.ls)[order(as.numeric(names(aggregated_ranks_cluster_tf.ls)))][i]
  p_cell_v_clust = ggVennDiagram(inp_tf_list,
                                 set_color = c("black", "black", "black", "black", "red"), label = "count", label_alpha = 0) + 
    scale_x_continuous(expand = expansion(mult = .2)) + 
    scale_fill_distiller(palette = "Reds", direction = 1)
  p_cell_v_clust
})
plts_tf_venn[[1]]
wrap_plots(plts_tf_venn, ncol = 3)
ggsave('06_deg_analysis/TF_celltype_pe_vs_rpca_venn.jpeg', height = 16, width = 20)


# Calculate the overlap percentages
overlap_df <- calculate_overlap(aggregated_ranks_cluster_tf.ls[order(as.numeric(names(aggregated_ranks_cluster_tf.ls)))], 
                                aggregated_ranks_celltype_pe_tf.ls[c(2,4,5,1)])

# Print the results
print(overlap_df)
overlap_df$CellType = factor(overlap_df$CellType)
overlap_df$CellType = fct_inorder(overlap_df$CellType)
overlap_df$Cluster = factor(overlap_df$Cluster)
overlap_df$Cluster = fct_inorder(overlap_df$Cluster)

# Visualize the results with a heatmap
heatmap_plot <- ggplot(overlap_df, aes(x = Cluster, y = CellType, fill = OverlapPercentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", OverlapPercentage)), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  labs(
    # title = "Percentage Overlap of TFs Between Integrated Clusters and Cell Types",
    x = "Cluster",
    y = "Cell Type",
    fill = "Overlap (%)"
  ) +
  theme_minimal() +
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the heatmap
print(heatmap_plot)

p2 = plotting.data %>% 
  filter(rpca_clusters %in% unique(overlap_df$Cluster)) %>%
  ggplot(aes(x=rpca_clusters, fill=celltype_pe))+
  geom_bar(position = 'fill')+
  scale_fill_manual(values = celltype_pe_color)

(p2+coord_fixed())/heatmap_plot
wrap_plots(p2, heatmap_plot, design = 'A\nB\nB\nB')

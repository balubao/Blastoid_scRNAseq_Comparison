## This script is to evaluate the abundance of cell types between datasets
##

# load packages and functions
source('./script/Initiation.R', chdir = TRUE)

# load data
split_obj = readRDS('04_integration/results/annotated_integrated_blastoids_revision2')
plotting.data = readRDS('results/plotting_data_revision_2.RDS')

print(split_obj)

##
### -------------------------------------------------------------------------------------------------------------
plotting.data$celltype = plotting.data$annot_major
celltype_prop.df = plotting.data[,c("sample", "celltype")]
celltype_prop.df$celltype = factor(celltype_prop.df$celltype)

celltype_prop.df$source = factor(celltype_prop.df$sample)
celltype_prop.df$source = fct_inorder(celltype_prop.df$source)
levels(celltype_prop.df$source) = c("Blastocyst",
                                    "Blastocyst",
                                    "Blastocyst",
                                    "nPSC",
                                    "nPSC",
                                    "Fibr.",
                                    "EPSC",
                                    "EPSC",
                                    "nPSC",
                                    "nPSC")
# levels(celltype_prop.df$source) = c("EPSC",
#                                     "nPSC",
#                                     "Fibroblast",
#                                     "Blastocyst",
#                                     "EPSC",
#                                     "Blastocyst",
#                                     "nPSC",
#                                     "Blastocyst",
#                                     "nPSC",
#                                     "nPSC")

## run test for all datasets
chisq_test(celltype_prop.df$celltype, celltype_prop.df$source)
# # A tibble: 1 × 6
# n statistic     p    df method          p.signif
# * <int>     <dbl> <dbl> <int> <chr>           <chr>   
#   1 35051     3504.     0     9 Chi-square test ****    

### - significance of annotation with respect to samples
# Create a contingency table
contingency_table <- table(celltype_prop.df$sample, celltype_prop.df$celltype)

# Perform the chi-square test
chi_square_result <- chisq.test(contingency_table)

# View the result
print(chi_square_result)

### - significance of annotation with respect to source cell line
# Create a contingency table
contingency_table <- table(celltype_prop.df$source, celltype_prop.df$celltype)

# Perform the chi-square test
chi_square_result <- chisq.test(contingency_table)

# View the result
print(chi_square_result)

### ---------------------------------------------------------------------------
celltype_prop_wounk.df = celltype_prop.df[celltype_prop.df$celltype!="unknown",]
celltype_prop_wounk.df$celltype = factor(celltype_prop_wounk.df$celltype)

### - significance of annotation with respect to samples
# Create a contingency table
contingency_table <- table(celltype_prop_wounk.df$sample, celltype_prop_wounk.df$celltype)

# Perform the chi-square test
chi_square_result <- chisq.test(contingency_table)

# View the result
print(chi_square_result)

### - significance of annotation with respect to source cell line
# Create a contingency table
contingency_table <- table(celltype_prop_wounk.df$source, celltype_prop_wounk.df$celltype)

# Perform the chi-square test
chi_square_result <- chisq.test(contingency_table)

# View the result
print(chi_square_result)

require(scales)
hex3 <- hue_pal()(3)
celltype_color <- c(hex3, "lightgrey")
names(celltype_color) <- c("EPI_ICM", "PE", "TE", "unknown")

celltype_prop_nounk.plotting_data = celltype_prop_wounk.df %>%
  group_by(sample, source, celltype) %>%
  count(celltype) %>%
  ungroup(celltype) %>%
  mutate(value = round(n/sum(n)*1700))

require(ggplot2)
p1 = ggplot(celltype_prop_nounk.plotting_data) +
  geom_col(aes(x = sample, y=value,fill=celltype), position = position_fill(reverse = TRUE)) +
  coord_flip() + 
  scale_fill_manual(values = celltype_color[seq(3)], name = "Cell Type")+
  facet_grid(source ~ ., scales = "free_y", space = "free_y", switch = "y")+
  theme(strip.placement = "outside",
        # strip.switch.pad.grid = unit(1, "cm"),
        axis.title.y = element_blank()) #element_text(vjust = -15))
p1
ggsave("./03_annotation/figures/celltype_nounk_prop_by_source_sample_bar.jpeg", width = 6, height = 5)

celltype_prop.plotting_data = celltype_prop.df %>%
  group_by(sample, source, celltype) %>%
  count(celltype)
# ungroup(celltype) %>%
# mutate(value = round(n/sum(n)*1700))
# celltype_prop.plotting_data$celltype = as.character(celltype_prop.plotting_data$celltype)

require(ggplot2)
p2 = ggplot(celltype_prop.plotting_data) +
  geom_col(aes(x = sample, y=n, fill=celltype), position = position_fill(reverse = TRUE)) +
  coord_flip() + 
  scale_fill_manual(values = celltype_color, name = "Cell Type") +
  facet_grid(source ~ ., scales = "free_y", space = "free_y", switch = "y")+
  theme(strip.placement = "outside",
        # strip.switch.pad.grid = unit(1, "cm"),
        axis.title.y = element_blank()) #element_text(vjust = -15))
p2
ggsave("./03_annotation/figures/celltype_prop_by_source_sample_bar.jpeg", width = 6, height = 5)

p1/p2
ggsave(plot = p1/p2, "./03_annotation/figures/celltype_prop_by_source_sample_bar_both.jpeg", width = 8, height = 8)
### - significance of annotation with respect to source cell line + dataset
# Create a contingency table
# contingency_table <- table(celltype_prop.df$source, celltype_prop.df$sample)
# contingency_table = CountNormalize(contingency_table)
# 
# # Perform the chi-square test
# chi_square_result <- chisq.test(contingency_table)
# 
# # View the result
# print(chi_square_result)

### -------------------------------------------------------------------------------------------------------------
## Color PE subtypes
pe_list = lapply(seurat.list, function(seurat.obj){
  seurat.obj = subset(seurat.obj, annot_major == "PE")
})
sapply(pe_list, ncol)
pe_list = lapply(pe_list, function(seurat.obj){
  seurat.obj = seurat.obj %>%
    # RunPCA(seurat.obj) %>%
    # FindNeighbors() %>%
    # FindClusters() %>%
    RunUMAP(dims = 1:10, n.neighbors = 8)
  return(seurat.obj)
})
plts = lapply(pe_list, function(seurat.obj){
  plt = DimPlot(seurat.obj, group.by = "annot_major") + ggtitle(unique(seurat.obj$sample))
  return(plt)
})
wrap_plots(plts)

# split_obj$celltype_pe = split_obj$celltype
# select_pe = split_obj$celltype_pe == "PE"
# split_obj$celltype_pe[select_pe] = paste0(split_obj$celltype_pe[select_pe], ifelse(split_obj@reductions$umap.rpca@cell.embeddings[select_pe,2]>0, "I", "II"))
# DimPlot(split_obj, group.by = "celltype_pe", reduction = "umap.rpca")

umap_vec = grep("umap", names(split_obj@reductions), value = TRUE)

require(scales)
hex3 <- hue_pal()(3)
hex3_sub <- hue_pal()(7)
celltype_pe_color <- c(hex3[1], hex3_sub[c(3,4)], hex3[3], "lightgrey")
names(celltype_pe_color) <- c("EPI_ICM", "PEI", "PEII", "TE", "unknown")

DimPlot(split_obj, group.by = "celltype_pe", reduction = "umap.rpca", label = TRUE)+
  scale_color_manual(values = celltype_pe_color)+
  theme(legend.position = "none")+
  theme_minimal()

plts = lapply(seq_along(umap_vec), function(i){
  DimPlot(split_obj, group.by = "celltype_pe", reduction = umap_vec[i])+
    scale_color_manual(values = celltype_pe_color)
})
wrap_plots(plts)

plotting.data$celltype_pe = split_obj$celltype_pe

celltype_prop.df = plotting.data[,c("sample", "celltype_pe")]
celltype_prop.df$celltype = factor(celltype_prop.df$celltype)

celltype_prop.df$source = factor(celltype_prop.df$sample)
celltype_prop.df$source = fct_inorder(celltype_prop.df$source)
levels(celltype_prop.df$source) = c("Blastocyst",
                                    "Blastocyst",
                                    "Blastocyst",
                                    "nPSC",
                                    "nPSC",
                                    "Fibr.",
                                    "EPSC",
                                    "EPSC",
                                    "nPSC",
                                    "nPSC")
celltype_prop_wounk.df = celltype_prop.df[plotting.data$celltype_pe!="unknown",]
celltype_prop_wounk.df$celltype = factor(celltype_prop_wounk.df$celltype_pe)

celltype_prop_nounk.plotting_data = celltype_prop_wounk.df %>%
  group_by(sample, source, celltype) %>%
  count(celltype) %>%
  ungroup(celltype) %>%
  mutate(value = round(n/sum(n)*1700))

require(ggplot2)
p1 = ggplot(celltype_prop_nounk.plotting_data) +
  geom_col(aes(x = sample, y=value,fill=celltype), position = position_fill(reverse = TRUE)) +
  coord_flip() + 
  scale_fill_manual(values = celltype_pe_color[seq(4)], name = "Cell Type")+
  facet_grid(source ~ ., scales = "free_y", space = "free_y", switch = "y")+
  theme(strip.placement = "outside",
        # strip.switch.pad.grid = unit(1, "cm"),
        axis.title.y = element_blank()) #element_text(vjust = -15))
p1
ggsave("./03_annotation/figures/celltype_nounk_pe_prop_by_source_sample_bar.jpeg", width = 6, height = 5)

# require(scales)
# hex3 <- hue_pal()(3)
# celltype_color <- c(hex3, "lightgrey")
# names(celltype_color) <- c("EPI_ICM", "PE", "TE", "unknown")
# 
celltype_prop.plotting_data = celltype_prop.df %>%
  group_by(sample, source, celltype) %>%
  count(celltype)
# ungroup(celltype) %>%
# mutate(value = round(n/sum(n)*1700))
# celltype_prop.plotting_data$celltype = as.character(celltype_prop.plotting_data$celltype)

require(ggplot2)
p2 = ggplot(celltype_prop.plotting_data) +
  geom_col(aes(x = sample, y=n, fill=celltype), position = position_fill(reverse = TRUE)) +
  coord_flip() + 
  scale_fill_manual(values = celltype_pe_color, name = "Cell Type") +
  facet_grid(source ~ ., scales = "free_y", space = "free_y", switch = "y")+
  theme(strip.placement = "outside",
        # strip.switch.pad.grid = unit(1, "cm"),
        axis.title.y = element_blank()) #element_text(vjust = -15))
p2
ggsave("./03_annotation/figures/celltype_pe_prop_by_source_sample_bar.jpeg", width = 6, height = 5)
p1/p2
ggsave(plot = p1/p2, "./03_annotation/figures/celltype_pe_prop_by_source_sample_bar_both.jpeg", width = 8, height = 8)
       
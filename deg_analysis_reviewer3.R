## This script is to identify GO for cell types and unknown cells
##

# load packages and functions
source('./script/Initiation.R', chdir = TRUE)

# load data
seurat.list = readRDS("./03_annotation/results/annotated_seurat_list_revision_2.rds")
split_obj = readRDS('04_integration/results/annotated_integrated_blastoids_revision2')
plotting.data = readRDS('results/plotting_data_revision_2.RDS')

split_obj$celltype = split_obj$annot_major
split_obj$clusters = split_obj$rpca_clusters


# only labeled
split_obj_lab = subset(split_obj,  celltype != "unknown")

# only unk
split_obj_unk = subset(split_obj,  celltype == "unknown")


split_obj_lab =  RunUMAP(split_obj_lab, reduction = "integrated.rpca", dims = 1:20)
DimPlot(split_obj_lab, group.by = "celltype", reduction = "umap.rpca")

split_obj_unk =  RunUMAP(split_obj_unk, reduction = "integrated.rpca", dims = 1:20)
DimPlot(split_obj_unk, group.by = "celltype", reduction = "umap.rpca")

################################################################################
# What are the identities of the labeled cells?



# transfer cluster labels to each object (from RPCA integration clusters)
# i=1
seurat.list = lapply(seurat.list, function(seurat.obj) {
  # seurat.obj = seurat.list[[i]]
  plotting.data_sub = plotting.data[plotting.data$sample == unique(seurat.obj$sample), ]
  seurat.obj$clusters = plotting.data_sub$rpca_clusters
  return(seurat.obj)
})

# get DEGs for each dataset using integrated.clusters  (RPCA)
topmarkers.all.ls = lapply(seurat.list, function(seurat.obj) {
  Idents(seurat.obj) = seurat.obj$clusters
  top_markers = FindAllMarkers(seurat.obj)
  return(top_markers)
})

# saveRDS(topmarkers.all.ls, file = "./06_deg_analysis/results/wilcox_deg_list_rpcaclust_unfiltered_all.rds")
topmarkers.all.ls = readRDS("./06_deg_analysis/results/wilcox_deg_list_rpcaclust_unfiltered_all.rds")

# 0.25LFC_1e2apval

clusters = Reduce(union, lapply(topmarkers.all.ls, function(x) x$cluster))

# label each df by dataset
topmarkers.all.ls = lapply(seq_along(topmarkers.all.ls), function(i) {
  topmarkers.all.ls[[i]]$sample = names(topmarkers.all.ls)[i]
  return(topmarkers.all.ls[[i]])
})
names(topmarkers.all.ls) = sapply(topmarkers.all.ls, function(x) unique(x$sample))

# keep genes appearing in at least 2 datasets
keep_genes.ls = lapply(seq_along(clusters), function(i)
{
  cluster_genes = lapply(topmarkers.all.ls, function(topmarkers) {
    gene_names = topmarkers$gene[topmarkers$cluster == i]
    return(gene_names)
  })
  cluster_genes = do.call("c", cluster_genes)
  keep_genes = cluster_genes[duplicated(cluster_genes)]
  return(keep_genes)
})

# filter DEG tables to only retain duplicated genes
topmarkers.all_dupfilt.ls = lapply(topmarkers.all.ls, function(topmarkers.all)
{
  topmarkers.all_sub.ls = lapply(seq_along(clusters), function(i) {
    topmarkers = topmarkers.all[topmarkers.all$cluster == (i - 1), ]
    topmarkers = topmarkers[topmarkers$gene %in% keep_genes.ls[[i]], ]
    return(topmarkers)
  })
  
  topmarkers.all_sub = do.call("rbind", topmarkers.all_sub.ls)
  return(topmarkers.all_sub)
})

sapply(topmarkers.all_dupfilt.ls, dim)
sapply(topmarkers.all_dupfilt.ls, function(x)
  table(x$cluster))

# filter out low significance genes
topmarkers.all_filt.ls = lapply(topmarkers.all_dupfilt.ls, function(topmarkers){
  topmarkers = topmarkers %>%
    filter((p_val_adj < 0.01))
  return(topmarkers)
})

sapply(topmarkers.all_filt.ls, dim)
sapply(topmarkers.all_filt.ls, function(x)
  table(x$cluster))

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
  combined_df <- filled_dfs %>%
    reduce(function(x, y)
      full_join(x, y, by = "Var1", suffix = c(".x", ".y"))) %>%
    replace(is.na(.), 0)
  
  # Rename columns for clarity (optional)
  colnames(combined_df) <- c("Var1", names(ndegs_per_clust.ls))
  
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

# rank genes by LFC for each cluster
topmarkers.all_ranked.ls = lapply(topmarkers.all_filt.ls, function(topmarkers) {
  topmarkers = topmarkers %>%
    group_by(cluster) %>%
    mutate(rank = rank(-avg_log2FC))
  
  return(topmarkers)
})

# aggregate cluster rank across datasets
# define functions
CheckGeneDist = function(topmarkers.clust.df) {
  pre_m = topmarkers.clust.df %>%
    # select(sample, gene) %>%
    split(.$sample) %>%
    map( ~ .x %>% select(-sample))
  pre_m = sapply(pre_m, function(x)
    x$gene)
  m = make_comb_mat(pre_m)
  UpSet(m)
}

# define variables
topmarkers.all_ranked.df = do.call("rbind", topmarkers.all_ranked.ls)
clusters = unique(topmarkers.all_ranked.df$cluster)

require(RobustRankAggreg)
aggregated_ranks.ls = lapply(seq_along(clusters), function(i) {
  print(i)
  topmarkers.clust.df = topmarkers.all_ranked.df[topmarkers.all_ranked.df$cluster == clusters[i], ]
  pre_ranked_list = topmarkers.clust.df %>%
    dplyr::select(sample, gene, rank, p_val_adj, cluster) %>%
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
  sum(x$Score < 0.1)) #number of ranked genes w/ p-value>0.01

# saveRDS(aggregated_ranks.ls, "./06_deg_analysis/results/aggregated_rank_gene_rpca_clust_list.rds")
aggregated_ranks.ls = readRDS("./06_deg_analysis/results/aggregated_rank_gene_rpca_clust_list.rds")

## -----------------------------------------------------------------------------
## Pathway Enrichment Analysis

# Next, we run a pathway enrichment analysis for the aggregated ranks to check for cluster enrichment.
# we run 4 times - 1 using all ranked genes, 3 using genes filtered by p-values 0.01, 0.05, 0.1

aggregated_ranks.ls

# load library
require(enrichR)

# Perform enrichment analysis using enrichR
dbs <-
  c(
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023"
  )

cluster_pathways_enrichment.ls = lapply(seq_along(aggregated_ranks.ls), function(i) {
  genes = aggregated_ranks.ls[[i]]$Name
  enrichr_results <- enrichr(genes, dbs)
  return(enrichr_results)
})

saveRDS(
  cluster_pathways_enrichment.ls,
  "./06_deg_analysis/results/pathway_enrichment_dup_rpca_clust_list.rds"
)

names(cluster_pathways_enrichment.ls) = paste0("cluster_", seq(31) - 1)
# define function for sanity check
CheckTermDist = function(cluster_terms.ls) {
  m = make_comb_mat(cluster_terms.ls)
  UpSet(m)
}

# plotEnrich(cluster_pathways_enrichment.ls[[1]]$KEGG_2021_Human)

## ----KEGG pathway
cluster_kegg_pathways_enrichment.ls = lapply(cluster_pathways_enrichment.ls, function(enrich_list) {
  return(enrich_list$KEGG_2021_Human)
})
# get only terms
cluster_kegg_pathways_terms.ls = lapply(cluster_kegg_pathways_enrichment.ls, function(enrich_list) {
  enrich_list = enrich_list[enrich_list$Adjusted.P.value < 0.01, ]
  return(enrich_list$Term)
})
# remove clusters with no terms
cluster_kegg_pathways_terms.ls = cluster_kegg_pathways_terms.ls[sapply(cluster_kegg_pathways_terms.ls, length) >
                                                                  0]
# sanity check
CheckTermDist(cluster_kegg_pathways_terms.ls)


## ----GO BP pathway
cluster_GO_BP_pathways_enrichment.ls = lapply(cluster_pathways_enrichment.ls, function(enrich_list) {
  return(enrich_list$GO_Biological_Process_2023)
})
# get only terms
cluster_GO_BP_pathways_terms.ls = lapply(cluster_GO_BP_pathways_enrichment.ls, function(enrich_list) {
  enrich_list = enrich_list[enrich_list$Adjusted.P.value < 0.01, ]
  return(enrich_list$Term)
})
# remove clusters with no terms
cluster_GO_BP_pathways_terms.ls = cluster_GO_BP_pathways_terms.ls[sapply(cluster_GO_BP_pathways_terms.ls, length) >
                                                                    1]
# sanity check
CheckTermDist(cluster_GO_BP_pathways_terms.ls)


## ----GO MF pathway
cluster_GO_MF_pathways_enrichment.ls = lapply(cluster_pathways_enrichment.ls, function(enrich_list) {
  return(enrich_list$GO_Molecular_Function_2023)
})
# get only terms
cluster_GO_MF_pathways_terms.ls = lapply(cluster_GO_MF_pathways_enrichment.ls, function(enrich_list) {
  enrich_list = enrich_list[enrich_list$Adjusted.P.value < 0.01, ]
  return(enrich_list$Term)
})
# remove clusters with no terms
cluster_GO_MF_pathways_terms.ls = cluster_GO_MF_pathways_terms.ls[sapply(cluster_GO_MF_pathways_terms.ls, length) >
                                                                    0]
# sanity check
CheckTermDist(cluster_GO_MF_pathways_terms.ls)


## ----GO CC pathway
cluster_GO_CC_pathways_enrichment.ls = lapply(cluster_pathways_enrichment.ls, function(enrich_list) {
  return(enrich_list$GO_Cellular_Component_2023)
})
# get only terms
cluster_GO_CC_pathways_terms.ls = lapply(cluster_GO_CC_pathways_enrichment.ls, function(enrich_list) {
  # enrich_list = cluster_GO_CC_pathways_enrichment.ls[[1]]
  enrich_list = enrich_list[enrich_list$Adjusted.P.value < 0.01, ]
  return(enrich_list$Term)
})
# remove clusters with no terms
cluster_GO_CC_pathways_terms.ls = cluster_GO_CC_pathways_terms.ls[sapply(cluster_GO_CC_pathways_terms.ls, length) >
                                                                    0]
# sanity check
CheckTermDist(cluster_GO_CC_pathways_terms.ls)

######### ----------------------------------------------------------------------
## Compare the terms

## EPI_ICM clusters - 12,5,-16
## PE I clusters - 1, 25
## PE II clusters - 4, -13
## TE clusters - 19, 7, 17, -26

## --- KEGG
cluster_KEGG_pathways_enrichment_primary.ls = cluster_kegg_pathways_enrichment.ls[c(12, 5, 16, 1, 25, 4, 13, 26, 19, 7, 17) +
                                                                                    1]
cluster_KEGG_pathways_terms_primary.ls = cluster_kegg_pathways_terms.ls[c(12, 5, 16, 1, 25, 4, 13, 26, 19, 7, 17) +
                                                                          1]
# sanity check
CheckTermDist(cluster_KEGG_pathways_terms_primary.ls)

# actually check the terms
plotEnrich(cluster_KEGG_pathways_enrichment_primary.ls$cluster_12)
plotEnrich(cluster_KEGG_pathways_enrichment_primary.ls$cluster_1)
plotEnrich(cluster_KEGG_pathways_enrichment_primary.ls$cluster_13)
plotEnrich(cluster_KEGG_pathways_enrichment_primary.ls$cluster_19)



## -- GO BP
cluster_GO_BP_pathways_enrichment_primary.ls = cluster_GO_BP_pathways_enrichment.ls[c(12, 5, 16, 1, 25, 4, 13, 26, 19, 7, 17) +
                                                                                      1]
cluster_GO_BP_pathways_terms_primary.ls = cluster_GO_BP_pathways_terms.ls[c(12, 5, 16, 1, 25, 4, 13, 26, 19, 7, 17) +
                                                                            1]

## -- GO CC
cluster_GO_CC_pathways_enrichment_primary.ls = cluster_GO_CC_pathways_enrichment.ls[c(12, 5, 16, 1, 25, 4, 13, 26, 19, 7, 17) +
                                                                                      1]
cluster_GO_CC_pathways_terms_primary.ls = cluster_GO_CC_pathways_terms.ls[c(12, 5, 16, 1, 25, 4, 13, 26, 19, 7, 17) +
                                                                            1]

## -- GO MF
cluster_GO_MF_pathways_enrichment_primary.ls = cluster_GO_MF_pathways_enrichment.ls[c(12, 5, 16, 1, 25, 4, 13, 26, 19, 7, 17) +
                                                                                      1]
cluster_GO_MF_pathways_terms_primary.ls = cluster_GO_MF_pathways_terms.ls[c(12, 5, 16, 1, 25, 4, 13, 26, 19, 7, 17) +
                                                                            1]


# sanity check
CheckTermDist(cluster_GO_BP_pathways_terms_primary.ls)
CheckTermDist(cluster_GO_CC_pathways_terms_primary.ls)
CheckTermDist(cluster_GO_MF_pathways_terms_primary.ls)

PlotClusterSummary = function(cluster_pathways_enrichment_primary.ls){
  require(enrichR)
  p1 = plotEnrich(cluster_pathways_enrichment_primary.ls$cluster_12, showTerms = 5, title = "EPI_ICM")
  p2 = plotEnrich(cluster_pathways_enrichment_primary.ls$cluster_1, showTerms = 5, title = "PE I")
  p3 = plotEnrich(cluster_pathways_enrichment_primary.ls$cluster_4, showTerms = 5, title = "PE II")
  p4 = plotEnrich(cluster_pathways_enrichment_primary.ls$cluster_19, showTerms = 5, title = "TE")
  p_all = wrap_plots(p1, p2, p3, p4, ncol = 1)
  return(p_all)
}

# actually check the terms
PlotClusterSummary(cluster_KEGG_pathways_enrichment_primary.ls)
PlotClusterSummary(cluster_GO_BP_pathways_enrichment_primary.ls)
PlotClusterSummary(cluster_GO_CC_pathways_enrichment_primary.ls)
PlotClusterSummary(cluster_GO_MF_pathways_enrichment_primary.ls)


p1 = plotEnrich(cluster_GO_BP_pathways_enrichment_primary.ls$cluster_12, showTerms = 5)
p2 = plotEnrich(cluster_GO_BP_pathways_enrichment_primary.ls$cluster_1, showTerms = 5)
p3 = plotEnrich(cluster_GO_BP_pathways_enrichment_primary.ls$cluster_13, showTerms = 5)
p4 = plotEnrich(cluster_GO_BP_pathways_enrichment_primary.ls$cluster_19, showTerms = 5)
wrap_plots(p1, p2, p3, p4, ncol = 1)

p1 = plotEnrich(cluster_GO_CC_pathways_enrichment_primary.ls$cluster_12)
p2 = plotEnrich(cluster_GO_CC_pathways_enrichment_primary.ls$cluster_1)
p3 = plotEnrich(cluster_GO_CC_pathways_enrichment_primary.ls$cluster_13)
p4 = plotEnrich(cluster_GO_CC_pathways_enrichment_primary.ls$cluster_19)
wrap_plots(p1, p2, p3, p4)

p1 = plotEnrich(cluster_GO_MF_pathways_enrichment_primary.ls$cluster_12)
p2 = plotEnrich(cluster_GO_MF_pathways_enrichment_primary.ls$cluster_1)
p3 = plotEnrich(cluster_GO_MF_pathways_enrichment_primary.ls$cluster_13)
p4 = plotEnrich(cluster_GO_MF_pathways_enrichment_primary.ls$cluster_19)
wrap_plots(p1, p2, p3, p4)

## Show only exclusive terms
GetTermsFromEnrich = function(enrich_list){
  terms_list = lapply(enrich_list, function(x) x$Term)
  return(terms_list)
}
GetExclusiveTerms = function(terms_list){
  all_terms = unlist(terms_list)
  intersect_terms = unique(all_terms[duplicated(all_terms)])
  exclusive_terms = all_terms[!(all_terms %in% intersect_terms)]
  terms_exc_list = lapply(terms_list, function(terms){
    terms = terms[terms %in% exclusive_terms]
    return(terms)
  })
  return(terms_exc_list)
}
GetEnrich_ExclusiveTerms = function(enrich_list, terms_exc_list){
  enrich_exc_list = lapply(seq_along(enrich_list), function(i) {
    x = enrich_list[[i]][enrich_list[[i]]$Term %in% terms_exc_list[[i]],]
    return(x)
    })
  names(enrich_exc_list) = names(enrich_list)
  return(enrich_exc_list)
}

## EPI_ICM clusters - 12,5,-16
## PE I clusters - 1, 25
## PE II clusters - 4, -13
## TE clusters - 19, 7, 17, -26

CheckTermDist(cluster_GO_BP_pathways_terms_primary.ls)
cluster_GO_BP_pathways_enrichment_primary_exc.ls = GetEnrich_ExclusiveTerms(cluster_GO_BP_pathways_enrichment_primary.ls, GetExclusiveTerms(GetTermsFromEnrich(cluster_GO_BP_pathways_enrichment_primary.ls)))
CheckTermDist(cluster_GO_BP_pathways_enrichment_primary_exc.ls)

p1 = plotEnrich(cluster_GO_BP_pathways_enrichment_primary_exc.ls$cluster_12)
p2 = plotEnrich(cluster_GO_BP_pathways_enrichment_primary_exc.ls$cluster_1)
p3 = plotEnrich(cluster_GO_BP_pathways_enrichment_primary_exc.ls$cluster_4)
p4 = plotEnrich(cluster_GO_BP_pathways_enrichment_primary_exc.ls$cluster_19)
wrap_plots(p1, p2, p3, p4)



library(simplifyEnrichment)

go_id = sapply(strsplit(cluster_GO_BP_pathways_enrichment_primary.ls$cluster_12$Term, " "), function(x) x[length(x)])
go_id = gsub("\\(","",go_id)
go_id = gsub(")","",go_id)
mat = GO_similarity(go_id)
pdf("./06_deg_analysis/results/simpleGO_BP_terms_clust12_dupDEGs.pdf")
simplifyGO(mat)
dev.off()


## sanity checks
DimPlot(split_obj,
        reduction = "umap.rpca",
        group.by = "tech",
        label = T)
DimPlot(split_obj,
        reduction = "umap.rpca",
        group.by = "harmony_clusters",
        label = T)
DimPlot(split_obj,
        reduction = "umap.rpca",
        group.by = "annot_major",
        label = T) + scale_color_manual(values = celltype_color)

DimPlot(split_obj,
        reduction = "umap.harmony",
        group.by = "tech",
        label = T)
DimPlot(split_obj,
        reduction = "umap.harmony",
        group.by = "harmony_clusters",
        label = T)
DimPlot(split_obj,
        reduction = "umap.harmony",
        group.by = "annot_major",
        label = T) + scale_color_manual(values = celltype_color)
DimPlot(split_obj,
        reduction = "umap.harmony",
        group.by = "sample",
        label = T)


################################################################################
# pathways inhibited by Kagawa for enhanced implantation
# - Hippo pathway
# - TGF-beta
# - ERK
library(KEGGREST)
library(org.Hs.eg.db)
kegg <- org.Hs.egPATH2EG
mapped <- mappedkeys(kegg)
kegg2 <- as.list(kegg[mapped])

KEGGREST::listDatabases()
org <- keggList("organism")
head(org)
# Pull all pathways for AT  
pathways.list <- keggList("pathway", "hsa")
head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
head(genes.by.pathway)

geneList <- DE.table$P.Value
names(geneList) <- DE.table$Gene
head(geneList)

# Wilcoxon test for each pathway
pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                             function(pathway) {
                               pathway.genes <- genes.by.pathway[[pathway]]
                               list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                               list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                               scores.in.pathway <- geneList[list.genes.in.pathway]
                               scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                               if (length(scores.in.pathway) > 0){
                                 p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                               } else{
                                 p.value <- NA
                               }
                               return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
                             }
))

# Assemble output table
outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- pathways.list[outdat$pathway.code]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat <- outdat[order(outdat$p.value),]
head(outdat)
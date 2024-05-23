

library(SeuratDisk)

print("Loading data")
seurat.list = readRDS("annotat")



## 2. Merge dataset and preprocess ---------------------------------------------

print("Merging data")
obj = merge(seurat.list[[1]], seurat.list[seq(2,length(seurat.list))])

print("Splitting data")
DefaultAssay(obj) = "RNA"

obj = JoinLayers(obj)

obj = split(obj, assay = "RNA", f = obj$sample)
obj@assays$RNA@layers[grep("scale.data",names(obj@assays$RNA@layers))] = NULL

print("Processing data")
plan("multicore", workers = 32)
obj <- FindVariableFeatures(obj)
obj$ERCC.perc[is.na(obj$ERCC.perc)] = 0
obj = ScaleData(obj, vars.to.regress = c("ERCC.perc","S.Score", "G2M.Score", "nFeature_RNA", "nCount_RNA"))
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
split_obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
# split_obj = split_obj[rowSums(split_obj)>0,]

save(split_obj, file = "preintegrated_seurat_obj.rda")


## 3. Integrate datasets -------------------------------------------------------

# A. CCA
split_obj = IntegrateLayers(object = split_obj, method = CCAIntegration, assay = "RNA",
                              orig.reduction = "pca", new.reduction = "integrated.cca", 
                              verbose = FALSE
)

# B. CCA w reference
split_obj = IntegrateLayers(object = split_obj, method = CCAIntegration, assay = "RNA",
                              orig.reduction = "pca", new.reduction = "integrated.cca_ref", 
                              reference = seq(3),
                              verbose = FALSE
)

# C. RPCA
split_obj <- IntegrateLayers(
  object = split_obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  normalization.method = "LogNormalize",
  k.weight = 80,
  verbose = TRUE
)

# D. RPCA w reference
split_obj <- IntegrateLayers(
  object = split_obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca_ref",
  normalization.method = "LogNormalize", reference = seq(3),
  k.weight = 80,
  verbose = TRUE
)

# E. MNN
split_obj <- IntegrateLayers(
  object = split_obj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = TRUE
)

# F. Harmony
split_obj <- IntegrateLayers(
  object = split_obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  normalization.method = "LogNormalize",
  k.weight = 80,
  verbose = TRUE
)

## 4. Generate embeddings and clusters ----------------------------------------
split_obj = JoinLayers(split_obj)

vec_int = c("cca", "cca_ref", "rpca", "rpca_ref", "mnn", "harmony")
for(i in seq_along(vec_int)){
  split_obj = FindNeighbors(split_obj, reduction = paste0("integrated.",vec_int[i])) %>% 
    FindClusters(resolution = 2, cluster.name = paste0(vec_int[i],"_clusters")) %>% 
    RunUMAP(reduction = paste0("integrated.",vec_int[i]), dims = 1:30, reduction.name = paste0("umap.",vec_int[i]))
}

## 5. Save Object  -------------------------------------------------------------
saveRDS(split_obj, "annotated_integrated_blastoids_revision2")

require(SeuratDisk)
SaveH5Seurat(split_obj, "annotated_integrated_blastoids_revision2")

## 5. Visualize and Export different integrations

umap_vec = c("umap.cca", "umap.cca_ref", "umap.rpca", "umap.rpca_ref", "umap.mnn", "umap.harmony")

for(i in seq_along(vec_int)){
  int_method = vec_int[i]
  umap_int = paste0("umap.", int_method)
  DimPlot(split_obj, group.by = c("sample", "tech", "ref", "annot_major"), reduction = umap_int, pt.size = 0.5, shuffle = TRUE)
  ggsave(paste0("04_integration/umap_emb/",int_method,"_umap.png"), width = 10, height = 6)
}

setwd('/home/balubao/Documents/Research/blastoid_comparison/')

source('script/Initiation.R')
# source('../GetGeneMetadata.R')

## Fan object assembly
counts = Read10X("../Data/Fan2021/filtered_feature_bc_matrix/", strip.suffix = T)
# counts = Read10X("../Data/Fan2021/raw_feature_bc_matrix/", strip.suffix = T)

library(SoupX)

# get all file paths
filepaths = system("ls Fan2021_SRR140*", intern = TRUE)
filepaths = strsplit(filepaths, " ")

# for loop counts
require(SoupX)
fan.list = list()
for(i in seq_along(filepaths)){
  
  #load filtered
  toc = Read10X(file.path(filepaths[i],"filtered_feature_bc_matrix"), strip.suffix = T)
  
  #load raw
  tod = Read10X(file.path(filepaths[i],"raw_feature_bc_matrix"), strip.suffix = T)
  
  #run soup
  sc = SoupChannel(tod, toc)
  sc = autoEstCont(sc)
  counts = adjustCounts(sc)
  
  # Add experiment name to cell barcodes to make them unique
  colnames(counts) = paste0(filepaths[i], "_", colnames(counts))
  
  # Store the result
  fan.list[[filepaths[i]]] = counts
  

}

# Combine all count matricies into one matrix
counts = do.call(cbind, fan.list)
seurat.obj = CreateSeuratObject(counts, 
                                project = "Fan2021", 
                                min.cells = 5,
                                min.features = 200)
seurat.obj$sample = sapply(strsplit(colnames(seurat.obj), split = "_"),"[[",1)

saveRDS(seurat.obj, "Fan2021_D6_strained.rds")

#cluster qc

#run doubletfinder
require(DoubletFinder)





Fan2021.seurat = CreateSeuratObject(counts, 
                                    project = "Fan2021", 
                                    min.cells = 5,
                                    min.features = 200)
outlist = QC_10x(list(Fan2021.seurat), Sourcelabel = "Fan2021", limit = 1000)
Fan2021.seurat = outlist[[1]][[1]]
Fan2021.seurat = GenerateQCMetrics(Fan2021.seurat)
PlotQCMetrics(Fan2021.seurat@meta.data, nCount_thresh = 6e5)
outlist[[2]][[1]]

Fan2021.seurat$sample = "day6_blastoids"

Fan2021.seurat = subset(Fan2021.seurat, (nCount_RNA<60000) & (nCount_RNA>1e4) & (mitoRatio<0.2)) 



#As in script
Fan2021.seurat = CreateSeuratObject(counts, 
                                    project = "Fan2021", 
                                    min.cells = 10,
                                    min.features = 100)
Fan2021.seurat = subset(Fan2021.seurat, (nCount_RNA<60000) & (nCount_RNA>7500) & (mitoRatio<0.2))


Fan2021.seurat = SCT_regress_cc(Fan2021.seurat)
Fan2021.seurat = RunPCA(Fan2021.seurat)
ElbowPlot(Fan2021.seurat)
# JackStrawPlot(Fan2021.seurat)
Fan2021.seurat = FindNeighbors(Fan2021.seurat, dims = 1:20) %>% FindClusters(resolution = 2, algorithm = 4, method = "igraph")
Fan2021.seurat = RunUMAP(Fan2021.seurat, dims = 1:15)

Fan2021.seurat = GetCellTypes(Fan2021.seurat, tissue="Embryo5")

mod.ls = GetAllModules(Fan2021.seurat, assay = "SCT")
QC.score = mod.ls$QCScores
Fan2021.seurat = mod.ls$seurat.object

DimPlot(Fan2021.seurat, reduction = "umap", group.by = c("seurat_clusters", "customclassif"))
FeaturePlot(Fan2021.seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio", "riboRatio"), order = T)
FeaturePlot(Fan2021.seurat, features = GetScoreNames(AllMarkers), order = T)

plts = DimPlot(Fan2021.seurat, reduction = "umap", group.by = c("seurat_clusters", "customclassif"), combine = FALSE)
plts1 = FeaturePlot(Fan2021.seurat, features = GetScoreNames(AllMarkers), order = T, combine = FALSE)
wrap_plots(c(plts, plts1))

### Plot heatymap of cell type markers
df2 = do.call("cbind",GetSCTypeMarkers("Embryo6"))
df2 = melt(df2)
df2 = df2[!duplicated(df2$value),]

Fan2021.seurat = GetResidual(Fan2021.seurat, df2$value)
plot_heatmap_mod(dataset = Fan2021.seurat,
                 # n = 6,
                 markers = df2$value,
                 sort_var = c("customclassif","EPI1","PE1","TE1"),
                 anno_var = c("customclassif","EPI1","PE1","TE1"),
                 anno_colors = list("Set1","Reds","Greens","Purples"),
                 hm_limit = c(-1,0,1),
                 hm_colors = c("purple","black","yellow"), 
                 row_split = df2$Var2)

saveRDS(Fan.seurat, '/home/balubao/Documents/Research/Data/Fan2021/Fan2021_uni_pp.rds')


## Fan object assembly

## Load data
fan.list = readRDS("../Data/Fan2021/Fan2021_list.rds")

fan.list = lapply(fan.list, GenerateQCMetrics)
fan.list = lapply(fan.list, function(x){
  x = RenameCells(x, paste0(unique(x$sample),'_',colnames(x)))
  return(x)})
tmp = lapply(fan.list, function(x){
  tmp=x@meta.data
  return(tmp)})
tmp = rbindlist(tmp, fill = T)



PlotQCMetrics(tmp,
              nCount_thresh = 5e5,
              nFeature_thresh = 7e3,
              mitoRatio_thresh = 0.2,
              riboRatio_thresh = NULL,
              log10FeatPerCount_thresh = 0.8, assay="RNA")


### Process
Fan.seurat = fan.list[[1]]
for(i in seq(2,4)){Fan.seurat=merge(Fan.seurat, fan.list[[i]])}

PlotQCMetrics(Fan.seurat@meta.data,
              nCount_thresh = 1e4,
              nFeature_thresh = 3e3,
              mitoRatio_thresh = 0.2,
              riboRatio_thresh = NULL,
              log10FeatPerCount_thresh = 0.8, assay="RNA")

Fan.seurat = subset(Fan.seurat, (mitoRatio < 0.2) & (nFeature_RNA > 3e3))

# ## Normalize
Fan.seurat <- SCT_regress_cc(Fan.seurat)

## QC
# revise GenerateQC to include Assay + slot
PlotQCMetrics(Fan.seurat@meta.data,
              nCount_thresh = 3e3,
              nFeature_thresh = 1e3,
              mitoRatio_thresh = 0.2,
              riboRatio_thresh = 0.4,
              log10FeatPerCount_thresh = 0.8, assay="SCT")

## Compute Modules
module.list = GetAllModules(Fan.seurat, 'SCT')
QC.table = module.list$QCScores
Fan.seurat = module.list$seurat.object

# module.list = GetBlastoidModules(Fan.seurat, 'SCT')
# QC.table = module.list$QCScores
# Fan.seurat = module.list$seurat.object
# module.list = GetPotencyModules(Fan.seurat, 'SCT')
# QC.table = c(QC.table, module.list$QCScores)
# Fan.seurat = module.list$seurat.object
# module.list = GetAmnionModules(Fan.seurat, 'SCT')
# QC.table = c(QC.table, module.list$QCScores)
# Fan.seurat = module.list$seurat.object

write.table(QC.table, "../Data/Fan2021/ScoreQC.txt")


## Downstream processing
Fan.seurat = RunPCA(Fan.seurat) 
ElbowPlot(Fan.seurat)
Fan.seurat = FindNeighbors(Fan.seurat) %>% FindClusters()
Fan.seurat = RunUMAP(Fan.seurat, dims = 1:20)

# saveRDS(Fan.seurat, '/home/balubao/Documents/Research/Data/Fan2021/Fan2021_v2_pp.rds')
# saveRDS(Fan.seurat, '/home/balubao/Documents/Research/Data/Fan2021/Fan2021_uni_pp.rds')
Fan.seurat = readRDS('/home/balubao/Documents/Research/Data/Fan2021/Fan2021_v2_pp.rds')

####################################
## Fan21 annotation and validation 
####################################

Fan.seurat = GetCellTypes(Fan.seurat, tissue = "Embryo7")
Fan.seurat = GetCellType(Fan.seurat)

# table(cbind(epi.idx, pe.idx, te.idx),rbind(epi.idx, pe.idx, te.idx))
DimPlot(Fan.seurat, group.by = c("seurat_clusters", "celltype"))
p1=FracPlot(Fan.seurat)

Fan.seurat$celltype_frac = Fan.seurat$seurat_clusters
levels(Fan.seurat$celltype_frac) = c("ELC","ELC","NC","NC","NC","TLC","NC","ELC","NC","ELC","ELC","TLC","ELC","ELC","HLC","TLC")
Fan.seurat$celltype_sfrac = paste0(Fan.seurat$celltype_frac,"_",Fan.seurat$seurat_clusters)

p2=DimPlot(Fan.seurat, group.by = c("celltype","seurat_clusters", "celltype_frac", "customclassif"))
wrap_plots(list(p1,p2))
####################################

## Sanity Checks
spl = c(DimPlot(Fan.seurat, reduction = 'umap', group.by = c('Phase'), combine = F),
        FeaturePlot(Fan.seurat, 
                    features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
                    order = T, combine = F))
wrap_plots(spl, ncol=2)

DimPlot(Fan.seurat, reduction = 'umap')
DimPlot(Fan.seurat, reduction = 'umap', group.by = c('sample','Phase','SCT_snn_res.0.8'))

FeaturePlot(Fan.seurat, 
            features = c('EPI1','PE1','TE1','ICM1',"AMN_E1","AMN_L1","STB1","CTB1","EVT1"), 
            order = T)
# FeaturePlot(Fan.seurat, 
            # features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
            # order = T)






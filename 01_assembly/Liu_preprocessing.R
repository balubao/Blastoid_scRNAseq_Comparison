setwd('/home/balubao/Documents/Research/blastoid_comparison/')

source('script/Initiation.R')
# source('../GetGeneMetadata.R')


## Liu object assembly

## Load data
ACCPATH = "../Data/Liu2021/PRJNA658478_Lui2021.txt"
acc = read.table(ACCPATH, header = T)

liu.list = readRDS("../Data/Liu2021/Lui2021_list.rds")
names(liu.list) = c("SRR12491161","SRR12491162")
for(i in seq_along(liu.list)){
  liu.list[[i]]$sample = paste0('Liu2021_',names(liu.list)[i])
  liu.list[[i]]$source = 'Liu2021'}

liu.list = lapply(liu.list, GenerateQCMetrics)
liu.list = lapply(liu.list, function(x){
  x = RenameCells(x, paste0(unique(x$sample),'_',colnames(x)))
  return(x)})
tmp = lapply(liu.list, function(x){
  tmp=x@meta.data
  return(tmp)})
# tmp = do.call('rbind', tmp)
tmp = as.data.frame(rbindlist(tmp, fill = T))
# RecreateSeuratObject
PlotQCMetrics(tmp,
              nCount_thresh = 3e3, #mad_thresh(tmp$nCount_RNA, direction = 'lower'),
              nFeature_thresh = 1e3, #mad_thresh(tmp$nFeature_RNA, direction = 'lower'),
              mitoRatio_thresh = 0.2,
              riboRatio_thresh = NULL,
              log10FeatPerCount_thresh = 0.8, assay="RNA")

liu.list = lapply(liu.list, function(x){subset(x, (mitoRatio < 0.2) & (nFeature_RNA > 1e3) & (nCount_RNA > 3e3))})

## Normalize
# liu.list = lapply(liu.list, SCT_regress_cc)

### Process
# Liu.seurat = merge(liu.list[[1]],liu.list[[2]])
Liu.seurat = liu.list[[1]]
Liu.seurat = SCT_regress_cc(Liu.seurat)
# for(i in seq(2,4)){Liu.seurat=merge(Liu.seurat, liu.list[[i]])}


## QC
# revise GenerateQC to include Assay + slot
PlotQCMetrics(Liu.seurat@meta.data,
              nCount_thresh = 3e3,
              nFeature_thresh = 1e3,
              mitoRatio_thresh = 0.2,
              riboRatio_thresh = 0.4,
              log10FeatPerCount_thresh = 0.8, assay="RNA")

## Compute Modules
module.list = GetAllModules(Liu.seurat, 'SCT')
QC.table = module.list$QCScores
Liu.seurat = module.list$seurat.object
# module.list = GetBlastoidModules(Liu.seurat, 'SCT')
# QC.table = module.list$QCScores
# Liu.seurat = module.list$seurat.object
# module.list = GetPotencyModules(Liu.seurat, 'SCT')
# QC.table = c(QC.table, module.list$QCScores)
# Liu.seurat = module.list$seurat.object
# module.list = GetAmnionModules(Liu.seurat, 'SCT')
# QC.table = c(QC.table, module.list$QCScores)
# Liu.seurat = module.list$seurat.object

write.table(QC.table, "../Data/Liu2021/ScoreQC.txt")


## Downstream processing
Liu.seurat = RunPCA(Liu.seurat) 
ElbowPlot(Liu.seurat)
Liu.seurat = FindNeighbors(Liu.seurat) %>% FindClusters()
Liu.seurat = RunUMAP(Liu.seurat, dims = 1:20)

# saveRDS(Liu.seurat, '/home/balubao/Documents/Research/Data/Liu2021/Liu2021_v2_pp.rds')
Liu.seurat = readRDS('/home/balubao/Documents/Research/Data/Liu2021/Liu2021_v2_pp.rds')

####################################
## Liu21 annotation and validation 
####################################

Liu.seurat = GetCellTypes(Liu.seurat, tissue = "Embryo7")
Liu.seurat = GetCellType(Liu.seurat)

# table(cbind(epi.idx, pe.idx, te.idx),rbind(epi.idx, pe.idx, te.idx))
DimPlot(Liu.seurat, group.by = c("seurat_clusters", "celltype"))
p1=FracPlot(Liu.seurat)

Liu.seurat$celltype_frac = Liu.seurat$seurat_clusters
levels(Liu.seurat$celltype_frac) = c("HLC", "NC", "ELC", "TLC", "TLC", "TLC", "TLC", "IM", "HLC", "ELC", "HLC", "NC", "IM", "ELC")
Liu.seurat$celltype_sfrac = paste0(Liu.seurat$celltype_frac,"_",Liu.seurat$seurat_clusters)

p2=DimPlot(Liu.seurat, group.by = c("celltype","seurat_clusters", "celltype_frac", "customclassif"))
wrap_plots(list(p1,p2))
####################################


## Sanity Checks
spl = c(DimPlot(Liu.seurat, reduction = 'umap', group.by = c('Phase','sample'), combine = F),
        FeaturePlot(Liu.seurat, 
                    features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
                    order = T, combine = F))
wrap_plots(spl, ncol=2)

DimPlot(Liu.seurat, reduction = 'umap')
DimPlot(Liu.seurat, reduction = 'umap', group.by = c('sample','Phase','SCT_snn_res.0.8'))

FeaturePlot(Liu.seurat, 
            features = c('EPI1','PE1','ICM1','Naive1','Primed1','TE1','AMN1',"AMN_E1","AMN_L1","CTB1","EVT1","STB1"), 
            order = T)
FeaturePlot(Liu.seurat, 
            features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
            order = T)






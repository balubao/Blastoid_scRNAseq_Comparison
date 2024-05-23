setwd('/home/balubao/Documents/Research/blastoid_comparison/')

source('script/Initiation.R')
# source('../GetGeneMetadata.R')


## Liu object assembly

## Load data
ACCPATH = "../Data/Yu2021/PRJNA686089_Yu2021.txt"
acc = read.table(ACCPATH, header = T)

yu.list = readRDS("../Data/Yu2021/Yu2021_list.rds")

names(yu.list) = paste0("SRR132682",seq(12,27))
for(i in seq_along(yu.list)){
  yu.list[[i]]$sample = paste0('Yu2021_',names(yu.list)[i])
  yu.list[[i]]$source = 'Yu2021'}


yu.list = lapply(yu.list, GenerateQCMetrics)
yu.list = lapply(yu.list, function(x){
  x = RenameCells(x, paste0(unique(x$sample),'_',colnames(x)))
  return(x)})
tmp = lapply(yu.list, function(x){
  tmp=x@meta.data
  return(tmp)})
tmp = as.data.frame(rbindlist(tmp, fill = T))
#tmp = rbindlist(tmp, fill = T, use.names = T)

PlotQCMetrics(tmp,
              nCount_thresh = mad_thresh(tmp$nCount_RNA, direction = 'lower'),
              nFeature_thresh = 1e3, #mad_thresh(tmp$nFeature_RNA, direction = 'lower'),
              mitoRatio_thresh = 0.2,
              riboRatio_thresh = NULL,
              log10FeatPerCount_thresh = 0.8, assay="RNA")

yu.list = lapply(yu.list, function(x){subset(x, (mitoRatio < 0.2) & (nFeature_RNA > 1e3))})

# run each object through the pp pipeline
yu.list = lapply(yu.list, function(Yu.seurat){
  #normalize
  Yu.seurat = SCT_regress_cc(Yu.seurat)
  
  ## Compute Modules
  module.list = GetAllModules(Yu.seurat, 'SCT')
  QC.table = module.list$QCScores
  Yu.seurat = module.list$seurat.object
  
  ## Downstream processing
  Yu.seurat = FindVariableFeatures(Yu.seurat, assay = 'SCT', nfeatures = 3000)
  Yu.seurat = RunPCA(Yu.seurat) 
  ElbowPlot(Yu.seurat)
  Yu.seurat = FindNeighbors(Yu.seurat) %>% FindClusters()
  Yu.seurat = RunUMAP(Yu.seurat, dims = 1:20)
  
  return(Yu.seurat)
})

### Process
Yu.seurat = yu.list[[1]]
for(i in seq(2,length(yu.list))){Yu.seurat=merge(Yu.seurat, yu.list[[i]])}

table(Yu.seurat$sample)
Yu.seurat$sample = factor(Yu.seurat$sample)
Yu.seurat$acc = Yu.seurat$sample
# levels(Yu.seurat$sample) = paste0("Yu2021_",acc$sample_title)


## QC
# revise GenerateQC to include Assay + slot
PlotQCMetrics(Yu.seurat@meta.data,
              nCount_thresh = 3e3,
              nFeature_thresh = 1e3,
              mitoRatio_thresh = 0.2,
              riboRatio_thresh = 0.4,
              log10FeatPerCount_thresh = 0.8, assay="RNA")

Yu.seurat = RecreateSeuratObject(Yu.seurat)
Yu.seurat = SCT_regress_cc(Yu.seurat)

PlotQCMetrics(Yu.seurat@meta.data,
              nCount_thresh = 3e3,
              nFeature_thresh = 1e3,
              mitoRatio_thresh = 0.2,
              riboRatio_thresh = 0.4,
              log10FeatPerCount_thresh = 0.8, assay="SCT")

## Compute Modules
module.list = GetAllModules(Yu.seurat, 'SCT')
QC.table = module.list$QCScores
Yu.seurat = module.list$seurat.object

# module.list = GetBlastoidModules(Yu.seurat, 'SCT')
# QC.table = module.list$QCScores
# Yu.seurat = module.list$seurat.object
# module.list = GetPotencyModules(Yu.seurat, 'SCT')
# QC.table = c(QC.table, module.list$QCScores)
# Yu.seurat = module.list$seurat.object
# module.list = GetAmnionModules(Yu.seurat, 'SCT')
# QC.table = c(QC.table, module.list$QCScores)
# Yu.seurat = module.list$seurat.object

write.table(QC.table, "../Data/Yu2021/ScoreQC.txt")


## Downstream processing
Yu.seurat = FindVariableFeatures(Yu.seurat, assay = 'SCT', nfeatures = 3000)
Yu.seurat = RunPCA(Yu.seurat) 
ElbowPlot(Yu.seurat)
Yu.seurat = FindNeighbors(Yu.seurat) %>% FindClusters()
Yu.seurat = RunUMAP(Yu.seurat, dims = 1:20)

# metaexp = data.frame(InternalID = c("LW36","LW58","LW59","LW60","LW61"),
#                      Method=c("TH","HT","HT","HT","HT"),
#                      Stage=c("D8","D3","D6","D9","D9"),
#                      Media=c("5iLA","5iLA","5iLA","5iLA","PXGL"))
# 
# keep.c = Yu2021.seurat$dataset %in% metaexp$InternalID
# Yu2021.seurat = subset(Yu2021.seurat,cells = colnames(Yu2021.seurat)[keep.c])
# keep.c = Yu2021.seurat$dataset %in% c("LW36","LW60","LW61")
# Yu2021.seurat = subset(Yu2021.seurat,cells = colnames(Yu2021.seurat)[keep.c])

#mergred Yu datasets
saveRDS(Yu.seurat, '/home/balubao/Documents/Research/Data/Yu2021/Yu2021_v2_pp.rds')

## Sanity Checks
spl_ls = lapply(yu.list, function(Yu.seurat){
  spl = c(DimPlot(Yu.seurat, reduction = 'umap', group.by = c('sample','Phase'), combine = F),
          FeaturePlot(Yu.seurat, 
                      features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
                      order = T, combine = F))
  wrap_plots(spl, ncol=2)
})
wrap_plots(spl_ls, ncol=2)
spl_ls = lapply(yu.list, function(Yu.seurat){
  spl = FeaturePlot(Yu.seurat, 
                features = c('EPI1','PE1','TE1','ICM1',"AMN_E1","AMN_L1","STB1","CTB1","EVT1"), 
                order = T)
})
wrap_plots(spl_ls, ncol=2)

spl = c(DimPlot(Yu.seurat, reduction = 'umap', group.by = c('sample','Phase'), combine = F),
        FeaturePlot(Yu.seurat, 
                    features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
                    order = T, combine = F))
wrap_plots(spl, ncol=2)

DimPlot(Yu.seurat, reduction = 'umap')
DimPlot(Yu.seurat, reduction = 'umap', group.by = c('sample','Phase','SCT_snn_res.0.8'))

FeaturePlot(Yu.seurat, 
            features = c('EPI1','PE1','TE1','ICM1',"AMN_E1","AMN_L1","STB1","CTB1","EVT1"), 
            order = T)
# FeaturePlot(Yu.seurat, 
#             features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
#             order = T)

####################################
Yu.seurat = read_rds('/home/balubao/Documents/Research/Data/Yu2021/Yu2021_v2_pp.rds')

####################################
## Yu21 annotation and validation 
####################################

Yu.seurat = GetCellTypes(Yu.seurat, tissue = "Embryo7")
Yu.seurat = GetCellType(Yu.seurat)

# table(cbind(epi.idx, pe.idx, te.idx),rbind(epi.idx, pe.idx, te.idx))
DimPlot(Yu.seurat, group.by = c("seurat_clusters", "celltype"))
p1=FracPlot(Yu.seurat)

Yu.seurat$celltype_frac = Yu.seurat$seurat_clusters
levels(Yu.seurat$celltype_frac) = c("ELC", "TLC", "ELC", "TLC", "ELC", "TLC", "IM", "TLC", "TLC", "ELC", "NC", "IM", "HLC")
Yu.seurat$celltype_sfrac = paste0(Yu.seurat$celltype_frac,"_",Yu.seurat$seurat_clusters)

p2=DimPlot(Yu.seurat, group.by = c("celltype","seurat_clusters", "celltype_frac", "customclassif"))
wrap_plots(list(p1,p2))

####################################


ggplot(Yanagida.seurat@meta.data, aes(x=AMN_L1, y=AMN_E1))+
  geom_point()

ggplot(Kagawa.seurat@meta.data, aes(x=TE1, y=AMN_E1))+
  geom_point()


ggplot(Petropoulos.seurat@meta.data, aes(x=AMN_E1, y=STB1))+
  geom_density_2d_filled()


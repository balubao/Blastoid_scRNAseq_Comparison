#############################################
## Natural Reference

Sozen_nat.seurat = readRDS('/home/balubao/Documents/Research/Data/Sozen2021/Sozen2021_nat_r2f.rds')

Sozen_nat.seurat$orig.ident = "Sozen2021_nat"

## normalize
Sozen_nat.seurat = SCT_regress_cc(Sozen_nat.seurat)

## Compute Modules
module.list = GetBlastoidModules(Sozen_nat.seurat, 'SCT')
QC.table = module.list$QCScores
Sozen_nat.seurat = module.list$seurat.object
module.list = GetPotencyModules(Sozen_nat.seurat, 'SCT')
QC.table = c(QC.table, module.list$QCScores)
Sozen_nat.seurat = module.list$seurat.object
module.list = GetAmnionModules(Sozen_nat.seurat, 'SCT')
QC.table = c(QC.table, module.list$QCScores)
Sozen_nat.seurat = module.list$seurat.object

write.table(QC.table, "../Data/Sozen2021/ScoreQC_nat.txt")


## Downstream processing
DefaultAssay(Sozen_nat.seurat) = "SCT"
Sozen_nat.seurat = RunPCA(Sozen_nat.seurat) 
ElbowPlot(Sozen_nat.seurat)
Sozen_nat.seurat = FindNeighbors(Sozen_nat.seurat) %>% FindClusters()
Sozen_nat.seurat = RunUMAP(Sozen_nat.seurat, dims = 1:10)

saveRDS(Sozen_nat.seurat, '/home/balubao/Documents/Research/Data/Sozen2021/Sozen_nat_pp.rds')


## Sanity Checks
spl = c(DimPlot(Sozen_nat.seurat, reduction = 'umap', group.by = c('Day','Phase'), combine = F),
        FeaturePlot(Sozen_nat.seurat, 
                    features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
                    order = T, combine = F))
wrap_plots(spl, ncol=2)

DimPlot(Sozen_nat.seurat, reduction = 'umap')
DimPlot(Sozen_nat.seurat, reduction = 'umap', group.by = c('Day','Phase','SCT_snn_res.0.8'))

FeaturePlot(Sozen_nat.seurat, 
            features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
            order = T)

#####################################
sozen.list = readRDS("../Data/Sozen2021/Sozen2021_list.rds")

sozen.list = lapply(sozen.list, GenerateQCMetrics)
sozen.list = lapply(sozen.list, function(x){
  x = RenameCells(x, paste0(unique(x$sample),'_',colnames(x)))
  return(x)})
Sozen.seurat = sozen.list[[1]]
# tmp = lapply(sozen.list, function(x){
  tmp=x@meta.data
  return(tmp)})
# tmp = rbindlist(tmp, fill = T)

#####################################
## Blastoid Model
library(SeuratDisk)
# Sozen.seurat = LoadH5Seurat('/home/balubao/Documents/Research/Data/Sozen2021/Sozen2021_r2f.h5seurat')


Sozen.seurat$orig.ident = "Sozen2021"
Sozen.seurat$sample = "Sozen2021_blast"
# "hEPSCs/D5_hEP-structures/D6_hEP-structures_multiplexed_samples"

## QC
# revise GenerateQC to include Assay + slot
Sozen.seurat = GenerateQCMetrics(Sozen.seurat, assay="RNA")
PlotQCMetrics(Sozen.seurat@meta.data, assay="RNA",
              nCount_thresh  = mad_thresh(Sozen.seurat$nCount_RNA),
              nFeature_thresh = mad_thresh(Sozen.seurat$nFeature_RNA), 
              mitoRatio_thresh = mad_thresh(Sozen.seurat$mitoRatio, direction = 'upper'),
              riboRatio_thresh = mad_thresh(Sozen.seurat$riboRatio, direction = 'upper'))
Sozen.seurat = subset(Sozen.seurat, (nFeature_RNA > 1e3) & (mitoRatio < 0.1))
# Sozen.seurat = subset(Sozen.seurat, (nFeature_TPM > mad_thresh(nFeature_TPM)) & (mitoRatio < mad_thresh(mitoRatio, direction = "upper")))

## normalize
Sozen.seurat = SCT_regress_cc(Sozen.seurat)

## Compute Modules
module.list = GetAllModules(Sozen.seurat, 'SCT')
QC.table = module.list$QCScores
Sozen.seurat = module.list$seurat.object

# module.list = GetBlastoidModules(Sozen.seurat, 'SCT')
# QC.table = module.list$QCScores
# Sozen.seurat = module.list$seurat.object
# module.list = GetPotencyModules(Sozen.seurat, 'SCT')
# QC.table = c(QC.table, module.list$QCScores)
# Sozen.seurat = module.list$seurat.object
# module.list = GetAmnionModules(Sozen.seurat, 'SCT')
# QC.table = c(QC.table, module.list$QCScores)
# Sozen.seurat = module.list$seurat.object

write.table(QC.table, "../Data/Sozen2021/ScoreQC.txt")


## Downstream processing
DefaultAssay(Sozen.seurat) = "SCT"
Sozen.seurat = RunPCA(Sozen.seurat) 
ElbowPlot(Sozen.seurat)
Sozen.seurat = FindNeighbors(Sozen.seurat) %>% FindClusters()
Sozen.seurat = RunUMAP(Sozen.seurat, dims = 1:10)

# saveRDS(Sozen.seurat, '/home/balubao/Documents/Research/Data/Sozen2021/Sozen_v2_pp.rds')
Sozen.seurat = readRDS('/home/balubao/Documents/Research/Data/Sozen2021/Sozen_v2_pp.rds')

####################################
## Sozen21 annotation and validation 
####################################

Sozen.seurat = GetCellTypes(Sozen.seurat, tissue = "Embryo7")
Sozen.seurat = GetCellType(Sozen.seurat)

# table(cbind(epi.idx, pe.idx, te.idx),rbind(epi.idx, pe.idx, te.idx))
DimPlot(Sozen.seurat, group.by = c("seurat_clusters", "celltype"))
p1=FracPlot(Sozen.seurat)

Sozen.seurat$celltype_frac = Sozen.seurat$seurat_clusters
levels(Sozen.seurat$celltype_frac) = c("HLC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "ELC", "TLC", "TLC")
Sozen.seurat$celltype_sfrac = paste0(Sozen.seurat$celltype_frac,"_",Sozen.seurat$seurat_clusters)

p2=DimPlot(Sozen.seurat, group.by = c("celltype","seurat_clusters", "celltype_frac", "customclassif"))
wrap_plots(list(p1,p2))
####################################


## Sanity Checks
spl = c(DimPlot(Sozen.seurat, reduction = 'umap', group.by = c('Day','Phase'), combine = F),
        FeaturePlot(Sozen.seurat, 
                    features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
                    order = T, combine = F))
wrap_plots(spl, ncol=2)

DimPlot(Sozen.seurat, reduction = 'umap')
DimPlot(Sozen.seurat, reduction = 'umap', group.by = c('Day','Phase','SCT_snn_res.0.8'))

FeaturePlot(Sozen.seurat, 
            features = c('EPI1','PE1','TE1','ICM1',"AMN_E1","AMN_L1","STB1","CTB1","EVT1"), 
            order = T)
FeaturePlot(Sozen.seurat, 
            features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
            order = T)

Sozen_nat.seurat = Sozen_nat.seurat[[1]]
Sozen_nat.seurat = GetCellTypes(Sozen_nat.seurat, tissue = "Embryo6")
DimPlot(Sozen_nat.seurat, group.by = c("Day","customclassif", "seurat_clusters"))




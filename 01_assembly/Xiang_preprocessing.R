setwd('/home/balubao/Documents/Research/blastoid_comparison/')

source('script/Initiation.R')
source('../GetGeneMetadata.R')


## Petropoulos object assembly

## Load data
CTSPATH = "../Data/Xiang2020/PRJNA562548_Xiang2020_counts.txt.gz"
METAPATH = "../Data/Xiang2020/PRJNA562548_Xiang2020_counts.txt.summary.gz"
ACCPATH = "../Data/Xiang2020/PRJNA562548_Xiang2020.txt"
ISMPATH = "../Data/Xiang2020/Xiang2020_CollatedInsertMetrics.txt"
getcellmeta_Xiang = function(acc){
  tmp = sapply(strsplit(acc$sample_title, split = "_"), '[[', 2)
  idx = str_locate(tmp, 'A')[,1]
  idx[is.na(idx)] = str_locate(tmp[is.na(idx)], 'N')[,1]
  
  Day = sapply(str_split(tmp, 'A'), '[[', 1) 
  Day = sapply(str_split(Day, 'N'), '[[', 1)
  
  label1 = rep(NA,length(Day))
  idx = str_detect(tmp, 'A')
  label1[idx] = sapply(str_split(tmp[idx], 'A'), function(x){paste0('A',x[length(x)])}) 
  idx = str_detect(tmp, 'N')
  label1[idx] = sapply(str_split(tmp[idx], 'N'), function(x){paste0('N',x[length(x)])}) 
  
  cellmeta = data.frame(Day=Day,
                        label1=label1,
                        row.names = paste0(acc$run_accession,'_',acc$sample_title))
  return(cellmeta)
}
Xiang.seurat = BuildSMARTseq2(CTSPATH = CTSPATH,
                                    METAPATH = METAPATH,
                                    ACCPATH = ACCPATH,
                                    getcellmeta = getcellmeta_Xiang,
                                    project = 'Xiang2020')
saveRDS(Xiang.seurat, '/home/balubao/Documents/Research/Data/Xiang2020/Xiang2020_raw.rds')
Xiang.seurat = readRDS('/home/balubao/Documents/Research/Data/Xiang2020/Xiang2020_raw.rds')
  
### Process
tpm_est = GetTPM_noeff(seurat.obj = Xiang.seurat,
                 CTSPATH = CTSPATH, 
                 ACCPATH = ACCPATH)

tpm_assay = CreateAssayObject(counts = tpm_est)
Xiang.seurat[['TPM']] = tpm_assay
DefaultAssay(Xiang.seurat) = "TPM"

# library(quminorm)
# qumi_assay = CreateAssayObject(counts = quminorm(tpm_est, mc.cores = 4))
# Xiang.seurat[['qUMI']] = qumi_assay
# DefaultAssay(Xiang.seurat) = "qUMI"

## QC
Xiang.seurat = GenerateQCMetrics(Xiang.seurat, assay="TPM")
PlotQCMetrics(Xiang.seurat@meta.data, 
              nCount_thresh  = mad_thresh(Xiang.seurat$nCount_TPM),
              nFeature_thresh = mad_thresh(Xiang.seurat$nFeature_TPM), 
              mitoRatio_thresh = mad_thresh(Xiang.seurat$mitoRatio, direction = "upper"))
# Xiang.seurat = subset(Xiang.seurat, nUMI > mad_thresh(nUMI))

## SCTransform
Xiang.seurat = SCT_SM2_regress_cc(Xiang.seurat, assay = "TPM")

## Compute Modules
module.list = GetBlastoidModules(Xiang.seurat, 'SCT_TPM')
QC.table = module.list$QCScores
Xiang.seurat = module.list$seurat.object
module.list = GetPotencyModules(Xiang.seurat, 'SCT_TPM')
QC.table = c(QC.table, module.list$QCScores)
Xiang.seurat = module.list$seurat.object
module.list = GetAmnionModules(Xiang.seurat, 'SCT_TPM')
QC.table = c(QC.table, module.list$QCScores)
Xiang.seurat = module.list$seurat.object

write.table(QC.table, "../Data/Xiang2020/ScoreQC.txt")

## Downstream processing
Xiang.seurat = RunPCA(Xiang.seurat) 
ElbowPlot(Xiang.seurat)
Xiang.seurat = FindNeighbors(Xiang.seurat) %>% FindClusters()
Xiang.seurat = RunUMAP(Xiang.seurat, dims = 1:20)

# saveRDS(Xiang.seurat, '/home/balubao/Documents/Research/Data/Xiang2020/Xiang2020_pp.rds')

## Sanity Checks
spl = c(DimPlot(Xiang.seurat, reduction = 'umap', group.by = c('Day','Phase'), combine = F),
        FeaturePlot(Xiang.seurat, 
                    features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
                    order = T, combine = F))
wrap_plots(spl, ncol=2)

DimPlot(Xiang.seurat, reduction = 'umap')
DimPlot(Xiang.seurat, reduction = 'umap', group.by = c('celltype_Man','celltype_Ont','Phase','SCT_snn_res.0.8'))

FeaturePlot(Xiang.seurat, 
            features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
            order = T)

###########################
## lets assume the relevant components of the data come from D6-8
Xiang.seurat = readRDS('/home/balubao/Documents/Research/Data/Xiang2020/Xiang2020_pp.rds')

relevant_days = c("D6","D7")
Xiang.seurat_sub = subset(Xiang.seurat, Day %in% relevant_days)


DimPlot(Xiang.seurat_sub, group.by = "Day")

Xiang.seurat_sub = FindVariableFeatures(Xiang.seurat_sub)
Xiang.seurat_sub = RunPCA(Xiang.seurat_sub) 
ElbowPlot(Xiang.seurat_sub)
Xiang.seurat_sub = FindNeighbors(Xiang.seurat_sub, dims = 1:10) %>% FindClusters()
Xiang.seurat_sub = RunUMAP(Xiang.seurat_sub, dims = 1:10)
DimPlot(Xiang.seurat_sub, group.by = c("seurat_clusters","Day"))

## Test with Xiang.seurat data
Xiang.seurat_sub = GetCellTypes(Xiang.seurat_sub, assay = "SCT_TPM", tissue = "Embryo6")
Xiang.seurat = GetCellTypes(Xiang.seurat, assay = "SCT_TPM", tissue = "Embryo6")

p1 = DimPlot(Xiang.seurat_sub, reduction = "umap", label = TRUE, group.by = 'customclassif')
p2 = DimPlot(Xiang.seurat, reduction = "umap", label = TRUE, group.by = 'customclassif')

p1+p2


DimPlot(Xiang.seurat, reduction = "umap", label = TRUE, group.by = c("seurat_clusters","Day",'customclassif'))

# Xiang.seurat_epi = subset(Xiang.seurat, seurat_clusters %in% c(3,4,5,8))
Xiang.seurat_epi = subset(Xiang.seurat_sub, seurat_clusters %in% c(0,4))
Xiang.seurat_epi = FindVariableFeatures(Xiang.seurat_epi)
Xiang.seurat_epi = RunPCA(Xiang.seurat_epi) 
ElbowPlot(Xiang.seurat_epi)
Xiang.seurat_epi = FindNeighbors(Xiang.seurat_epi, dims = 1:10) %>% FindClusters(resolution = 3)
Xiang.seurat_epi = RunUMAP(Xiang.seurat_epi, dims = 1:10)
Xiang.seurat_epi = GetCellTypes(Xiang.seurat_epi, assay = "SCT_TPM", tissue = "Embryo6")
DimPlot(Xiang.seurat_epi, reduction = "umap", label = TRUE, group.by = c("seurat_clusters","Day",'customclassif'))

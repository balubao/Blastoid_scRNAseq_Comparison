setwd('/home/balubao/Documents/Research/blastoid_comparison/')

source('script/Initiation.R')
source('../GetGeneMetadata.R')


## Kagawa object assembly

## Load data
CTSPATH = "../Data/Kagawa2022/PRJNA737134_Kagawa2022_counts.txt.gz"
METAPATH = "../Data/Kagawa2022/PRJNA737134_Kagawa2022_counts.txt.summary.gz"
ACCPATH = "../Data/Kagawa2022/PRJNA737134_Kagawa2022.txt"
ISMPATH = "../Data/Kagawa2022/Kagawa2022_CollatedInsertMetrics.txt"
getcellmeta_kagawa = function(acc){
  tmp = strsplit(acc$sample_title, ' ')
  tmp = lapply(tmp, function(x){unlist(strsplit(x, '-'))})
  cellmeta = data.frame(Origin=sapply(tmp, "[[", 1),
                        Day=sapply(tmp, "[[", 2),
                        row.names = paste0(acc$run_accession,'_',acc$sample_title))
  return(cellmeta)
}
Kagawa.seurat = BuildSMARTseq2(CTSPATH = CTSPATH,
                              METAPATH = METAPATH,
                              ACCPATH = ACCPATH,
                              getcellmeta = getcellmeta_kagawa,
                              project = 'Kagawa2022')
# saveRDS(Kagawa.seurat, '/home/balubao/Documents/Research/Data/Kagawa2022/Kagawa2022_raw.rds')
Kagawa.seurat = readRDS('/home/balubao/Documents/Research/Data/Kagawa2022/Kagawa2022_raw.rds')
Kagawa.seurat = subset(Kagawa.seurat, Origin == "blastoid")

### Process
tpm_est = GetTPM_noeff(seurat.obj = Kagawa.seurat,
                 CTSPATH = CTSPATH, 
                 ACCPATH = ACCPATH)

tpm_assay = CreateAssayObject(counts = tpm_est)
Kagawa.seurat[['TPM']] = tpm_assay
DefaultAssay(Kagawa.seurat) = "TPM"
 
# library(quminorm)
# qumi_assay = CreateAssayObject(counts = quminorm(tpm_est, shape = 1.5, mc.cores = 4))
# Kagawa.seurat[['qUMI']] = qumi_assay
# DefaultAssay(Kagawa.seurat) = "qUMI"
# 
# rpk_est = GetRPK(seurat.obj = Kagawa.seurat,
#                        CTSPATH = CTSPATH, 
#                        ACCPATH = ACCPATH)
# 
# rpk_assay = CreateAssayObject(counts = rpk_est)
# Kagawa.seurat[['RPK']] = rpk_assay
# 
# library(quminorm)
# qumi_assay = CreateAssayObject(counts = quminorm(rpk_est, shape = 1.5, mc.cores = 4))
# Kagawa.seurat[['qUMI_rpk']] = qumi_assay
# DefaultAssay(Kagawa.seurat) = "qUMI_rpk"



## QC
# revise GenerateQC to include Assay + slot
Kagawa.seurat = GenerateQCMetrics(Kagawa.seurat, assay="TPM")
PlotQCMetrics(Kagawa.seurat@meta.data, assay="TPM",
              nCount_thresh  = mad_thresh(Kagawa.seurat$nCount_TPM),
              nFeature_thresh = mad_thresh(Kagawa.seurat$nFeature_TPM), 
              mitoRatio_thresh = 0.2)
Kagawa.seurat = subset(Kagawa.seurat, nFeature_TPM > mad_thresh(nFeature_TPM))

## Normalize
Kagawa.seurat <- SCT_SM2_regress_cc(seurat.object = Kagawa.seurat, assay = "TPM") #SCT with relaxed parameters for SM2 case

## Compute Modules
module.list = GetBlastoidModules(Kagawa.seurat, 'SCT_TPM')
QC.table = module.list$QCScores
Kagawa.seurat = module.list$seurat.object
module.list = GetPotencyModules(Kagawa.seurat, 'SCT_TPM')
QC.table = c(QC.table, module.list$QCScores)
Kagawa.seurat = module.list$seurat.object
module.list = GetAmnionModules(Kagawa.seurat, 'SCT_TPM')
QC.table = c(QC.table, module.list$QCScores)
Kagawa.seurat = module.list$seurat.object

write.table(QC.table, "../Data/Kagawa2022/ScoreQC_nat.txt")


## Downstream processing
Kagawa.seurat = RunPCA(Kagawa.seurat) 
ElbowPlot(Kagawa.seurat)
Kagawa.seurat = FindNeighbors(Kagawa.seurat) %>% FindClusters()
Kagawa.seurat = RunUMAP(Kagawa.seurat, dims = 1:20)

# saveRDS(Kagawa.seurat, '/home/balubao/Documents/Research/Data/Kagawa2022/Kagawa2022_pp.rds')
Kagawa.seurat = readRDS('/home/balubao/Documents/Research/Data/Kagawa2022/Kagawa2022_pp.rds')

## Sanity Checks
spl = c(DimPlot(Kagawa.seurat, reduction = 'umap', group.by = c('Day','Phase'), combine = F),
        FeaturePlot(Kagawa.seurat, 
            features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
            order = T, combine = F))
wrap_plots(spl, ncol=2)

DimPlot(Kagawa.seurat, reduction = 'umap')
DimPlot(Kagawa.seurat, reduction = 'umap', group.by = c('Day','Phase','SCT_snn_res.0.8'))

FeaturePlot(Kagawa.seurat, 
            features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
            order = T)






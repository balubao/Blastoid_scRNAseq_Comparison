setwd('/home/balubao/Documents/Research/blastoid_comparison/')

source('script/Initiation.R')
source('../GetGeneMetadata.R')


## Yanagida object assembly

## Load data
CTSPATH = "../Data/Yanagida2021/PRJNA720968_Yanagida2021_counts.txt"
METAPATH = "../Data/Yanagida2021/PRJNA720968_Yanagida2021_counts.txt.summary"
ACCPATH = "../Data/Yanagida2021/PRJNA720968_Yanagida2021.txt"
getcellmeta_Yana = function(acc){
  tmp = strsplit(acc$sample_title, split = "[.]")
  cellmeta = data.frame(Day=sapply(tmp, '[[', 1),
                        group=sapply(tmp, '[[', 2),
                        section = sapply(tmp, '[[', 3),
                        row.names = paste0(acc$run_accession,'_',acc$sample_title))
  return(cellmeta)
}
Yanagida.seurat = BuildSMARTseq2(CTSPATH = CTSPATH,
               METAPATH = METAPATH,
               ACCPATH = ACCPATH,
               getcellmeta = getcellmeta_Yana,
               project = 'Yanagida2021')
saveRDS(Yanagida.seurat, '/home/balubao/Documents/Research/Data/Yanagida2021/Yanagida2021_raw.rds')
Yanagida.seurat = readRDS('/home/balubao/Documents/Research/Data/Yanagida2021/Yanagida2021_raw.rds')

## Break into respective groups
Yanagida_nat.seurat = subset(Yanagida.seurat, Day %in% c("Day7","Day6","Day5"))
Yanagida.seurat = subset(Yanagida.seurat, Day %in% c("Gata3_Day3","Gata3_Day4"))

# ##################### Natural
# 
# # rsem-prepare-reference --gtf Homo_sapiens.GRCh38.105.ERCC.gtf Homo_sapiens.GRCh38.dna.primary_assembly.ERCC.fa rsem_ERCC
# # 
# # ../software/RSEM-1.2.25/rsem-calculate-expression -p 8 --paired-end \
# # --star --star-path star \
# # --estimate-rspd \
# # --append-names \
# # --no-bam-output \
# # --calc-ci --single-cell-prior \
# # LPS_6h_sim_2M_1.fq LPS_6h_sim_2M_2.fq \
# # ../ref/mouse_ref LPS_6h_sim_2M
# # 
# # rsem-calculate-expression -p 1 \
# # --star --star-path /sw/csi/star/2.6.1d/el7.5_gnu6.4.0/STAR-2.6.1d \
# # --estimate-rspd \
# # --append-names \
# # --no-bam-output \
# # --calc-ci --single-cell-prior \
# # ../data/fastq/ERR1041403.fastq.gz \
# # ../../genome_ref/homo_sapiens/rsem_ERCC ERR1041403
# 
# Bam files are truncated, and can not be used reliably. We advise to either remap all datasets - this takes time. 
# Thus, we recommend postponing the process till after the conference and proceeding with the intial tentative results.
# 




##################### Blastoid
## Preprocessing

### Process
tpm_est = GetTPM_noeff(seurat.obj = Yanagida.seurat,
                       CTSPATH = CTSPATH, 
                       ACCPATH = ACCPATH)

tpm_assay = CreateAssayObject(counts = tpm_est)
Yanagida.seurat[['TPM']] = tpm_assay
DefaultAssay(Yanagida.seurat) = "TPM"



## QC
# revise GenerateQC to include Assay + slot
Yanagida.seurat = GenerateQCMetrics(Yanagida.seurat, assay="TPM")
PlotQCMetrics(Yanagida.seurat@meta.data, 
              nCount_thresh  = mad_thresh(Yanagida.seurat$nCount_TPM),
              nFeature_thresh = mad_thresh(Yanagida.seurat$nFeature_TPM), 
              mitoRatio_thresh = mad_thresh(Yanagida.seurat$mitoRatio, direction = "upper"))
# Yanagida.seurat = subset(Yanagida.seurat, nUMI > mad_thresh(nUMI))

## SCTransform
Yanagida.seurat = SCT_SM2_regress_cc(Yanagida.seurat, assay = "TPM")

## Compute Modules
module.list = GetBlastoidModules(Yanagida.seurat, 'SCT_TPM')
QC.table = module.list$QCScores
Yanagida.seurat = module.list$seurat.object
module.list = GetPotencyModules(Yanagida.seurat, 'SCT_TPM')
QC.table = c(QC.table, module.list$QCScores)
Yanagida.seurat = module.list$seurat.object
module.list = GetAmnionModules(Yanagida.seurat, 'SCT_TPM')
QC.table = c(QC.table, module.list$QCScores)
Yanagida.seurat = module.list$seurat.object

write.table(QC.table, "../Data/Yanagida2021/ScoreQC.txt")

## Downstream processing
Yanagida.seurat = RunPCA(Yanagida.seurat) 
ElbowPlot(Yanagida.seurat)
Yanagida.seurat = FindNeighbors(Yanagida.seurat) %>% FindClusters()
Yanagida.seurat = RunUMAP(Yanagida.seurat, dims = 1:20)

saveRDS(Yanagida.seurat, '/home/balubao/Documents/Research/Data/Yanagida2021/Yanagida2021_pp.rds')
Yanagida.seurat=readRDS('/home/balubao/Documents/Research/Data/Yanagida2021/Yanagida2021_pp.rds')

## Sanity Checks
spl = c(DimPlot(Yanagida.seurat, reduction = 'umap', group.by = c('Day','Phase'), combine = F),
        FeaturePlot(Yanagida.seurat, 
                    features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
                    order = T, combine = F))
wrap_plots(spl, ncol=2)

DimPlot(Yanagida.seurat, reduction = 'umap')
DimPlot(Yanagida.seurat, reduction = 'umap', group.by = c('celltype_Man','celltype_Ont','Phase','SCT_snn_res.0.8'))

FeaturePlot(Yanagida.seurat, 
            features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
            order = T)
# FeaturePlot(Yanagida.seurat, 
#             features = c('EPI1','PE1','TE1','AMN1','ICM1','Naive1', 'Primed1'), 
#             order = T)





## Checkpoint
# saveRDS(Yanagida_nat.seurat, '/home/balubao/Documents/Research/Data/Yanagida2021/Yanagida2021_r2f.rds')
# saveRDS(Yanagida_nat.seurat, '/home/balubao/Documents/Research/Data/Yanagida2021/Yanagida2021_nat_r2f.rds')


################
#### Natural Component
Yanagida.seurat = readRDS('/home/balubao/Documents/Research/Data/Yanagida2021/Yanagida2021_raw.rds')

## Break into respective groups
Yanagida_nat.seurat = subset(Yanagida.seurat, Day %in% c("Day7","Day6","Day5"))

### Process
tpm_est = GetTPM_noeff(seurat.obj = Yanagida_nat.seurat,
                       CTSPATH = CTSPATH, 
                       ACCPATH = ACCPATH)
tpm_assay = CreateAssayObject(counts = tpm_est)
Yanagida_nat.seurat[['TPM']] = tpm_assay
DefaultAssay(Yanagida_nat.seurat) = "TPM"

## QC
# revise GenerateQC to include Assay + slot
Yanagida_nat.seurat = GenerateQCMetrics(Yanagida_nat.seurat, assay="TPM")
PlotQCMetrics(Yanagida_nat.seurat@meta.data, 
              nCount_thresh  = mad_thresh(Yanagida_nat.seurat$nCount_TPM),
              nFeature_thresh = mad_thresh(Yanagida_nat.seurat$nFeature_TPM), 
              mitoRatio_thresh = mad_thresh(Yanagida_nat.seurat$mitoRatio, direction = "upper"))
# Yanagida_nat.seurat = subset(Yanagida_nat.seurat, nUMI > mad_thresh(nUMI))

## SCTransform
Yanagida_nat.seurat = SCT_SM2_regress_cc(Yanagida_nat.seurat, assay = "TPM")

## Compute Modules
module.list = GetBlastoidModules(Yanagida_nat.seurat, 'SCT_TPM')
QC.table = module.list$QCScores
Yanagida_nat.seurat = module.list$seurat.object
module.list = GetPotencyModules(Yanagida_nat.seurat, 'SCT_TPM')
QC.table = c(QC.table, module.list$QCScores)
Yanagida_nat.seurat = module.list$seurat.object
module.list = GetAmnionModules(Yanagida_nat.seurat, 'SCT_TPM')
QC.table = c(QC.table, module.list$QCScores)
Yanagida_nat.seurat = module.list$seurat.object

write.table(QC.table, "../Data/Yanagida2021/ScoreQC_nat.txt")

## Downstream processing
Yanagida_nat.seurat = RunPCA(Yanagida_nat.seurat) 
ElbowPlot(Yanagida_nat.seurat)
Yanagida_nat.seurat = FindNeighbors(Yanagida_nat.seurat) %>% FindClusters()
Yanagida_nat.seurat = RunUMAP(Yanagida_nat.seurat, dims = 1:10)

saveRDS(Yanagida_nat.seurat, '/home/balubao/Documents/Research/Data/Yanagida2021/Yanagida2021_nat_pp.rds')

## Sanity Checks
spl = c(DimPlot(Yanagida_nat.seurat, reduction = 'umap', group.by = c('Day','Phase'), combine = F),
        FeaturePlot(Yanagida_nat.seurat, 
                    features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
                    order = T, combine = F))
wrap_plots(spl, ncol=2)
DimPlot(Yanagida_nat.seurat, reduction = 'umap')
DimPlot(Yanagida_nat.seurat, reduction = 'umap', group.by = c('celltype_Man','celltype_Ont','Phase','SCT_snn_res.0.8'))
FeaturePlot(Yanagida_nat.seurat, 
            features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
            order = T)


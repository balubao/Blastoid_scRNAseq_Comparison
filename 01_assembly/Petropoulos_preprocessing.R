setwd('/home/balubao/Documents/Research/blastoid_comparison/')

source('script/Initiation.R')



## Petropoulos object assembly

## Load data
CTSPATH = "../Data/Petropoulos2016/PRJEB11202_Petropolous2016_counts.txt.gz"
METAPATH = "../Data/Petropoulos2016/PRJEB11202_Petropolous2016_counts.txt.summary.gz"
ACCPATH = "../Data/Petropoulos2016/PRJEB11202_Petropolous2016.txt"
getcellmeta_Petro = function(acc){
  tmp = strsplit(acc$sample_title, split = "[.]")
  cellmeta = data.frame(Day=sapply(tmp,'[[',1),
                        Day2=sapply(lapply(tmp, function(x){x[seq(length(x)-2)]}),function(x){paste(x,collapse="_")}),
                        label1=sapply(tmp, function(x){x[length(x)-1]}),
                        label2 = sapply(tmp, function(x){x[length(x)]}),
                        row.names = paste0(acc$run_accession,'_',acc$sample_title))
  return(cellmeta)
}
Petropoulos.seurat = BuildSMARTseq2(CTSPATH = gzfile(CTSPATH),
                                 METAPATH = gzfile(METAPATH),
                                 ACCPATH = ACCPATH,
                                 getcellmeta = getcellmeta_Petro,
                                 project = 'Petropoulos2016')
# saveRDS(Petropoulos.seurat, '/home/balubao/Documents/Research/Data/Petropoulos2016/Petropoulos2016_raw.rds')
Petropoulos.seurat = readRDS('/home/balubao/Documents/Research/Data/Petropoulos2016/Petropoulos2016_raw.rds')

### Process
### Process
tpm_est = GetTPM_noeff(seurat.obj = Petropoulos.seurat,
                       CTSPATH = CTSPATH, 
                       ACCPATH = ACCPATH)

tpm_assay = CreateAssayObject(counts = tpm_est)
Petropoulos.seurat[['TPM']] = tpm_assay
DefaultAssay(Petropoulos.seurat) = "TPM"

# library(quminorm)
# qumi_assay = CreateAssa0yObject(counts = quminorm(tpm_est, shape = 1.5, mc.cores = 4))
# Petropoulos.seurat[['TPM']] = qumi_assay
# DefaultAssay(Petropoulos.seurat) = "qUMI"
# 
# rpk_est = GetRPK(seurat.obj = Petropoulos.seurat,
#                        CTSPATH = CTSPATH, 
#                        ACCPATH = ACCPATH)
# 
# rpk_assay = CreateAssayObject(counts = rpk_est)
# Petropoulos.seurat[['RPK']] = rpk_assay
# 
# library(quminorm)
# qumi_assay = CreateAssayObject(counts = quminorm(rpk_est, shape = 1.5, mc.cores = 4))
# Petropoulos.seurat[['qUMI_rpk']] = qumi_assay
# DefaultAssay(Petropoulos.seurat) = "qUMI_rpk"



## QC
# revise GenerateQC to include Assay + slot
Petropoulos.seurat = GenerateQCMetrics(Petropoulos.seurat, assay="TPM")
PlotQCMetrics(Petropoulos.seurat@meta.data, assay="TPM",
              nCount_thresh  = mad_thresh(Petropoulos.seurat$nCount_TPM),
              nFeature_thresh = mad_thresh(Petropoulos.seurat$nFeature_TPM), 
              mitoRatio_thresh = mad_thresh(Petropoulos.seurat$mitoRatio, direction = 'upper'))
Petropoulos.seurat = subset(Petropoulos.seurat, (nFeature_TPM > mad_thresh(nFeature_TPM)) & (mitoRatio < mad_thresh(mitoRatio, direction = "upper")))

# ## Normalize
# seurat.list = SplitObject(Petropoulos.seurat, split.by = 'Day')
# seurat.list = lapply(seurat.list, SCT_SM2_regress_cc, assay="TPM")
# seurat.list = lapply(seurat.list,function(seurat.object){
#   module.list = GetBlastoidModules(seurat.object, 'SCT_TPM')
#   QC.table = module.list$QCScores
#   seurat.object = module.list$seurat.object
#   module.list = GetPotencyModules(seurat.object, 'SCT_TPM')
#   QC.table = c(QC.table, module.list$QCScores)
#   seurat.object = module.list$seurat.object
#   module.list = GetAmnionModules(seurat.object, 'SCT_TPM')
#   QC.table = c(QC.table, module.list$QCScores)
#   seurat.object = module.list$seurat.object
# })
# Petropoulos.seurat = seurat.list[[1]]
# for(i in seq(2,length(seurat.list))){Petropoulos.seurat = merge(Petropoulos.seurat, seurat.list[[i]])}
# Petropoulos.seurat = FindVariableFeatures(Petropoulos.seurat, assay = 'TPM') %>% ScaleData()
# integration.features = SelectIntegrationFeatures(seurat.list, assay = rep("SCT_TPM", length(seurat.list)))
# seurat.list = PrepSCTIntegration(seurat.list, anchor.features = integration.features, assay = rep("SCT_TPM", length(seurat.list)))
# anchors = FindIntegrationAnchors(seurat.list, assay = rep("SCT_TPM", length(seurat.list)))
# 
# # integrate data
# integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 80)
# Petropoulos.seurat = integrated

Petropoulos.seurat <- SCT_SM2_regress_cc(Petropoulos.seurat, assay = "TPM") #SCT with relaxed parameters for SM2 case
 
# Petropoulos.seurat <- SM2_regress_cc(Petropoulos.seurat, assay = "TPM") #SCT with relaxed parameters for SM2 case

 
## Compute Modules
# module.list = GetBlastoidModules(Petropoulos.seurat, 'SCT_TPM')
# QC.table = module.list$QCScores
# Petropoulos.seurat = module.list$seurat.object
# module.list = GetPotencyModules(Petropoulos.seurat, 'SCT_TPM')
# QC.table = c(QC.table, module.list$QCScores)
# Petropoulos.seurat = module.list$seurat.object
# module.list = GetAmnionModules(Petropoulos.seurat, 'SCT_TPM')
# QC.table = c(QC.table, module.list$QCScores)
# Petropoulos.seurat = module.list$seurat.object

module.list = GetAllModules(Petropoulos.seurat, 'SCT_TPM')
QC.table = module.list$QCScores
Petropoulos.seurat = module.list$seurat.object

write.table(QC.table, "../Data/Petropoulos2016/ScoreQC.txt")


## Downstream processing
DefaultAssay(Petropoulos.seurat) = "SCT_TPM"
Petropoulos.seurat = RunPCA(Petropoulos.seurat, assay = "SCT_TPM") 
ElbowPlot(Petropoulos.seurat)
Petropoulos.seurat = FindNeighbors(Petropoulos.seurat) %>% FindClusters()
Petropoulos.seurat = RunUMAP(Petropoulos.seurat, dims = 1:10)

# saveRDS(Petropoulos.seurat, '/home/balubao/Documents/Research/Data/Petropoulos2016/Petropoulos_pp.rds')
# Petropoulos.seurat=readRDS('/home/balubao/Documents/Research/Data/Petropoulos2016/Petropoulos_pp.rds')

## Sanity Checks
spl = c(DimPlot(Petropoulos.seurat, reduction = 'umap', group.by = c('Day','Phase'), combine = F),
        FeaturePlot(Petropoulos.seurat, 
                    features = c('mitoRatio','riboRatio','S.Score','G2M.Score','CC.Score'), 
                    order = T, combine = F))
wrap_plots(spl, ncol=2)

DimPlot(Petropoulos.seurat, reduction = 'umap')
DimPlot(Petropoulos.seurat, reduction = 'umap', group.by = c('Day','Phase','SCT_snn_res.0.8'))

FeaturePlot(Petropoulos.seurat, 
            features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
            order = T)



plotting.data = Petropoulos.seurat@meta.data
tmp = reshape2::melt(plotting.data[])
tmp2= tmp[tmp$variable %in% c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"),]
ggplot(tmp2,
       aes(x=Day, y=value))+
  geom_violin()+
  facet_grid(rows = vars(variable))+
  theme_bw()

ggplot(Petropoulos.seurat@meta.data, aes(x=Day,y=Primed1))+geom_violin()


#####################
## filter for reference use (E5-7)
Petropoulos.seurat = subset(Petropoulos.seurat, Day %in% c("E6")) #c("E5","E6","E7"))
# Petropoulos.seurat = CreateSeuratObject(counts = Petropoulos.seurat@assays$RNA@counts, 
#                                         meta.data = Petropoulos.seurat@meta.data,
#                                         min.features = 100, 
#                                         min.cells = 5)
AddAz
# normalize
Petropoulos.seurat <- SCT_SM2_regress_cc(Petropoulos.seurat, assay = "TPM") #SCT with relaxed parameters for SM2 case
## Downstream processing
DefaultAssay(Petropoulos.seurat) = "TPM"
Petropoulos.seurat = SM2_regress_cc(Petropoulos.seurat, assay = "TPM", block.size = 3000)
Petropoulos.seurat = RunPCA(Petropoulos.seurat) 
ElbowPlot(Petropoulos.seurat)
Petropoulos.seurat = FindNeighbors(Petropoulos.seurat) %>% FindClusters()
Petropoulos.seurat = RunUMAP(Petropoulos.seurat, dims = 1:10)
# saveRDS(Petropoulos.seurat, '/home/balubao/Documents/Research/Data/Petropoulos2016/Petropoulos_567_pp.rds')


Petropoulos.seurat = GetAllModules(Petropoulos.seurat, assay = "SCT_TPM")
Petropoulos.seurat = Petropoulos.seurat$seurat.object

DimPlot(Petropoulos.seurat, group.by = c('sample','Day'))
FeaturePlot(Petropoulos.seurat, features = c('EPI1','PE1','ICM1',"TE1","AMN_E1"))
FeaturePlot(Petropoulos.seurat, 
            features = c('EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1"), 
            order = T)
M = cor(Embeddings(Petropoulos.seurat, reduction = 'pca'),Petropoulos.seurat@meta.data[,c('Day','EPI1','PE1','ICM1','Naive1','Primed1',"PSC1",'TE1','AMN1',"AMN_E1","AMN_L1","STB1")])
ggcorrplot(M)


ggplot(Petropoulos.seurat@meta.data, aes(x=Day,y=TE1))+geom_violin()


## Test for results that agree with expected - 
Petropoulos.seurat = ScaleData(Petropoulos.seurat, assay = "RNA")
Petropoulos.seurat = ScaleData(Petropoulos.seurat, assay = "TPM")
Petropoulos.seurat = GetCellTypes(Petropoulos.seurat, tissue = "Embryo2", assay = "SCT_TPM", labname = "sclab_shared")
Petropoulos.seurat = GetCellTypes(Petropoulos.seurat, tissue = "Embryo3", assay = "SCT_TPM", labname = "sclab_gerri")
Petropoulos.seurat = GetCellTypes(Petropoulos.seurat, tissue = "Embryo4", assay = "SCT_TPM", labname = "sclab_petro")
Petropoulos.seurat = GetCellTypes(Petropoulos.seurat, tissue = "Embryo5", assay = "SCT_TPM", labname = "sclab_stri")
Petropoulos.seurat = GetCellTypes(Petropoulos.seurat, tissue = "Embryo6", assay = "SCT_TPM", labname = "sclab_all")
DimPlot(Petropoulos.seurat, group.by = "customclassif", split.by = "Day")


DimPlot(Petropoulos.seurat, group.by = c("Day","sclab_shared","sclab_gerri","sclab_petro","sclab_stri"))
ggplot(Petropoulos.seurat@meta.data, aes(x=Day, fill=customclassif))+
  geom_bar(position = "fill")

p1=ggplot(Petropoulos.seurat@meta.data, aes(x=Day, fill=sclab_shared))+
  geom_bar(position = "fill")
p2=ggplot(Petropoulos.seurat@meta.data, aes(x=Day, fill=sclab_gerri))+
  geom_bar(position = "fill")
p3=ggplot(Petropoulos.seurat@meta.data, aes(x=Day, fill=sclab_petro))+
  geom_bar(position = "fill")
p4=ggplot(Petropoulos.seurat@meta.data, aes(x=Day, fill=sclab_stri))+
  geom_bar(position = "fill")
p1|p2|p3

p1=ggplot(Petropoulos.seurat@meta.data, aes(x=Day, fill=sclab_shared))+
  geom_bar()
p2=ggplot(Petropoulos.seurat@meta.data, aes(x=Day, fill=sclab_gerri))+
  geom_bar()
p3=ggplot(Petropoulos.seurat@meta.data, aes(x=Day, fill=sclab_petro))+
  geom_bar()
p4=ggplot(Petropoulos.seurat@meta.data, aes(x=Day, fill=sclab_stri))+
  geom_bar(position = "fill")
p1|p2|p3

Petropoulos.seurat$NANOG2PDGFRA = GetAssayData(Petropoulos.seurat, assay = "TPM")["NANOG",]/
  GetAssayData(Petropoulos.seurat,assay = "TPM")["PDGFRA",]

Petropoulos.seurat$NANOG2PDGFRA[is.nan(Petropoulos.seurat$NANOG2PDGFRA)] = NA #0/0 results
Petropoulos.seurat$NANOG2PDGFRA[is.infinite(Petropoulos.seurat$NANOG2PDGFRA)] = 9.170156e+02 #1/0 results
use.c = !is.na(Petropoulos.seurat$NANOG2PDGFRA)
h = hclust(d = dist(Petropoulos.seurat$NANOG2PDGFRA[use.c]))
plot(h)




Petropoulos.seurat
df2 = do.call("cbind",AllMarkers)
df2 = melt(df2)
df2 = df2[!duplicated(df2$value),]
DoHeatmap(Petropoulos.seurat, features = df2$value, group.by = "sclab_stri") + NoLegend()
DoMultiBarHeatmap(Petropoulos.seurat, features = df2$value, additional.group.by = 'Day', group.by = "sclab_stri")




Petropoulos.seurat = GetResidual(Petropoulos.seurat, df2$value)
tmp = GetAssayData(Petropoulos.seurat, assay = "SCT_TPM")
use.g =  match(df2$value,rownames(tmp))
use.g = use.g[!is.na(use.g)]
tmp = tmp[use.g,]
library(Scillus)
plot_heatmap(dataset = Petropoulos.seurat,
             # n = 6,
             markers = df2$value,
             sort_var = c("sclab_all","Day","EPI1","PE1","TE1"),
             anno_var = c("sclab_all","Day","EPI1","PE1","TE1"),
             anno_colors = list("Set1","Accent","Reds","Greens","Purples"),
             hm_limit = c(-1,0,1),
             hm_colors = c("purple","black","yellow"))
# plot_heatmap_mod(dataset = Petropoulos.seurat,
#                  # n = 6,
#                  markers = df2$value,
#                  row_split = df2$Var2[df2$value %in% rownames(tmp)],
#                  sort_var = c("sclab_stri","Day"),
#                  anno_var = c("sclab_stri","Day"),
#                  anno_colors = list("Accent","Set1"),
#                  hm_limit = c(-1,0,1),
#                  hm_colors = c("purple","black","yellow"))


maxexp = rowMaxs(tmp)
# tmp = sweep(tmp, 1, maxexp, "/")
Heatmap(as.matrix(tmp), 
        row_split = df2$Var2[df2$value %in% rownames(tmp)], 
        column_split = Petropoulos.seurat$sclab_stri,
        show_column_names = FALSE)

 
## recap Fig4A in Striparo
df.ICM = data.frame(sclab = Petropoulos.seurat$sclab_stri,
                    POU5F1 = GetAssayData(Petropoulos.seurat, assay = "SCT_TPM")["POU5F1",])
ggplot(df.ICM, aes(x=rank(POU5F1), y=log10(POU5F1)))+geom_point()
df.ICM$NANOG = GetAssayData(Petropoulos.seurat, assay = "SCT_TPM")["NANOG",]
df.ICM$PDGFRA = GetAssayData(Petropoulos.seurat, assay = "SCT_TPM")["PDGFRA",]
df.ICM = df.ICM[log10(df.ICM$POU5F1)>0.7,]

df.ICM$NANOG2PDGFRA = df.ICM$NANOG/df.ICM$PDGFRA
df.ICM$NANOG2PDGFRA[is.nan(df.ICM$NANOG2PDGFRA)] = NA #0/0 results
df.ICM$NANOG2PDGFRA[is.infinite(df.ICM$NANOG2PDGFRA)] = 20 #1/0 results
use.c = !is.na(df.ICM$NANOG2PDGFRA)
h = hclust(d = dist(df.ICM$NANOG2PDGFRA[use.c]))
plot(h)

ggplot(df.ICM[use.c,])+
  geom_point(aes(x=h$order, 
                 color=sclab,
                 y=NANOG))
ggplot(df.ICM[use.c,])+
  geom_point(aes(x=h$order, 
                 color=sclab,
                 y=PDGFRA))

ggplot()+
  geom_point(aes(x=h$order, 
                 color=Petropoulos.seurat$sclab_stri[use.c],
                 y=GetAssayData(Petropoulos.seurat)["PDGFRA",use.c]))

ggplot()+
  geom_point(aes(x=Petropoulos.seurat$NANOG2PDGFRA, 
                 color=Petropoulos.seurat$sclab_stri,
                 y=GetAssayData(Petropoulos.seurat)["NANOG",]))
ggplot()+
  geom_point(aes(x=Petropoulos.seurat$NANOG2PDGFRA, 
                 color=Petropoulos.seurat$sclab_stri,
                 y=GetAssayData(Petropoulos.seurat)["PDGFRA",]))


#####################################################
## Genes are rapidly activiated and inhibited in early embryogenesis, therefore, we define the top markers along the trajectory of the different key lineages

Petropoulos.seurat_sub = subset(Petropoulos.seurat, Day == "E5")
Petropoulos.seurat_sub = DietSeurat(Petropoulos.seurat_sub)
# Petropoulos.seurat_sub = SCT_regress_cc(Petropoulos.seurat_sub, assay = "RNA")
Petropoulos.seurat_sub = NormalizeData(Petropoulos.seurat_sub) %>% FindVariableFeatures()
Petropoulos.seurat_sub = RunPCA(Petropoulos.seurat_sub) %>% FindNeighbors() %>% FindClusters()
Petropoulos.seurat_sub = RunUMAP(Petropoulos.seurat_sub, dims = 1:15)
ElbowPlot(Petropoulos.seurat_sub)
DimPlot(Petropoulos.seurat_sub, group.by = "customclassif")
DimPlot(Petropoulos.seurat_sub, group.by = "label1")
FeaturePlot(Petropoulos.seurat_sub, features = c("NANOG", "GATA6", "GATA3"))

seurat.list = SplitObject(Petropoulos.seurat_sub, split.by = "label1")
seurat.list = lapply(seurat.list, function(seurat.obj){
  seurat.obj = RecreateSeuratObject(seurat.obj)
  seurat.obj = NormalizeData(seurat.obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
})

FeaturePlot(Petropoulos.seurat, features = c("POU5F1", "CDH1", "NANOG"))
DimPlot(Petropoulos.seurat, group.by = "customclassif")


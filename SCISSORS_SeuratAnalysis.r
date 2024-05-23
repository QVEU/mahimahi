
################################################################################
### Load Libraries ###
################################################################################
remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

library(Seurat)
library(SeuratObject)
library(clustree)
library(tidyverse)

################################################################################
### Define File Paths ###
################################################################################

### Create Seurat Object. Reads from these samples were mapped to file path /data/lvd_qve/QVEU_Code/sequencing/template_fastas/refdata-gex-GRCh38-2020-A_PV_GFP_mRuby/GRCh38-2020-A_PV_GTF_mRuby/

dataset_loc <- "~/lab_share/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/"
ids <- c("RFP_C109S_PV", "WT_GFP_PV", "WT_GFP_RFP_Y88P_PV", "WT_GFP_RFP_D177A_PV", "RFP_Y88P_PV", "RFP_D177A_PV", "WT_GFP_RFP_C109S_PV", "WT_IRES_GFP_PV","Del_IRES_mRuby3_PV","WT_IRES_GFP_Del_IRES_mRuby3_PV","WT_IRES_mRuby3_MutPol_PV","Del_IRES_mRuby3_MutPol_PV","WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV", "Mock_5h_PV")
d10x.data <- sapply(ids, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x})
experiment.data <- do.call("cbind", d10x.data)

merged_seurat_object <- CreateSeuratObject(
  experiment.data,
  project = "Mutants",
  min.cells = 3,
  min.features = 10,
  names.field = 2,
  names.delim = "\\-")



################################################################################
### Add % mito gene expression to seurat object metadata
################################################################################

merged_seurat_object[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_object, pattern = "^MT-")
merged_seurat_object[["percent.virus"]] <- PercentageFeatureSet(merged_seurat_object, features = "PV")
merged_seurat_object[["mRuby"]] <- PercentageFeatureSet(merged_seurat_object, features = "mRuby")
merged_seurat_object[["GFP"]] <- PercentageFeatureSet(merged_seurat_object, features = "GFP")
merged_seurat_object[["percent.ribo"]] <- PercentageFeatureSet(merged_seurat_object, pattern = "^RP[SL]")

################################################################################
### Perform integration analysis ###
################################################################################

# split the dataset into a list
seurat_object_merged <- SplitObject(merged_seurat_object, split.by = "orig.ident")
#Remove the ridiculously huge chunk of code and instead make a for loop that find the infected threshold for virus_percentage and labels them according to low, high or not infected.
for (i in ids){
  if (i == "Mock_5h_PV"){
    seurat_object_merged[[i]]$InfectedStatus <- "Not_Infected"
    break
  }
  halfmin_virus = min(seurat_object_merged[[i]]$percent.virus[seurat_object_merged[[i]]$percent.virus>0])
  
  percent.virus_log10 = log10(seurat_object_merged[[i]]$percent.virus + halfmin_virus)
  
  viruspercentage = mixtools::normalmixEM(seurat_object_merged[[i]]$percent.virus[percent.virus_log10 > -2.5],k = 2)
  #call infected any cell above the mu + 2x sigma =  ----
  
  seurat_object_merged[[i]]$InfectedStatus <- "Not_Infected"
  seurat_object_merged[[i]]$InfectedStatus[seurat_object_merged[[i]]$percent.virus > (viruspercentage$mu[1]+viruspercentage$sigma[1]*2)] <- "Infected"
  
  table(seurat_object_merged[[i]]$InfectedStatus)
  
  seurat_object_merged[[i]]$InfectedStatus_groups <- "Not_Infected"
  seurat_object_merged[[i]]$InfectedStatus_groups[seurat_object_merged[[i]]$percent.virus > viruspercentage$mu[2]] <- "High"
  seurat_object_merged[[i]]$InfectedStatus_groups[seurat_object_merged[[i]]$percent.virus < viruspercentage$mu[2]] <- "Low"
  seurat_object_merged[[i]]$InfectedStatus_groups[seurat_object_merged[[i]]$percent.virus < (viruspercentage$mu[1]+viruspercentage$sigma[1]*2)]<- "Not_Infected" #annotate the threshold in the meta.data object
}

for (i in ids){
  assign(value = seurat_object_merged[[i]], x= paste(i))
}

for (i in ids)
{
  assign(value = read.table(file.path(dataset_loc,i,"outs/filtered_feature_bc_matrix/Doublet_scores.tsv"), header = F), x= paste0("doublets",i))
  
}
colnames(doubletsMock_5h_PV) <- c("doublet_scores","predicted_doublets")
Mock_5h_PV <- AddMetaData(Mock_5h_PV, metadata = doubletsMock_5h_PV$doublet_scores, col.name = "doublet_scores")
Mock_5h_PV <- AddMetaData(Mock_5h_PV, metadata = doubletsMock_5h_PV$predicted_doublets, col.name = "predicted_doublets")
head(Mock_5h_PV[[]])
table(Mock_5h_PV$predicted_doublets)

colnames(doubletsDel_IRES_mRuby3_MutPol_PV) <- c("doublet_scores","predicted_doublets")
Del_IRES_mRuby3_MutPol_PV <- AddMetaData(Del_IRES_mRuby3_MutPol_PV, metadata = doubletsDel_IRES_mRuby3_MutPol_PV$doublet_scores, col.name = "doublet_scores")
Del_IRES_mRuby3_MutPol_PV <- AddMetaData(Del_IRES_mRuby3_MutPol_PV, metadata = doubletsDel_IRES_mRuby3_MutPol_PV$predicted_doublets, col.name = "predicted_doublets")
head(Del_IRES_mRuby3_MutPol_PV[[]])
table(Del_IRES_mRuby3_MutPol_PV$predicted_doublets)

colnames(doubletsDel_IRES_mRuby3_PV) <- c("doublet_scores","predicted_doublets")
Del_IRES_mRuby3_PV <- AddMetaData(Del_IRES_mRuby3_PV, metadata = doubletsDel_IRES_mRuby3_PV$doublet_scores, col.name = "doublet_scores")
Del_IRES_mRuby3_PV <- AddMetaData(Del_IRES_mRuby3_PV, metadata = doubletsDel_IRES_mRuby3_PV$predicted_doublets, col.name = "predicted_doublets")
head(Del_IRES_mRuby3_PV[[]])
table(Del_IRES_mRuby3_PV$predicted_doublets)

colnames(doubletsRFP_C109S_PV) <- c("doublet_scores","predicted_doublets")
RFP_C109S_PV <- AddMetaData(RFP_C109S_PV, metadata = doubletsRFP_C109S_PV$doublet_scores, col.name = "doublet_scores")
RFP_C109S_PV <- AddMetaData(RFP_C109S_PV, metadata = doubletsRFP_C109S_PV$predicted_doublets, col.name = "predicted_doublets")
head(RFP_C109S_PV[[]])
table(RFP_C109S_PV$predicted_doublets)

colnames(doubletsRFP_D177A_PV) <- c("doublet_scores","predicted_doublets")
RFP_D177A_PV <- AddMetaData(RFP_D177A_PV, metadata = doubletsRFP_D177A_PV$doublet_scores, col.name = "doublet_scores")
RFP_D177A_PV <- AddMetaData(RFP_D177A_PV, metadata = doubletsRFP_D177A_PV$predicted_doublets, col.name = "predicted_doublets")
head(RFP_D177A_PV[[]])
table(RFP_D177A_PV$predicted_doublets)

colnames(doubletsRFP_Y88P_PV) <- c("doublet_scores","predicted_doublets")
RFP_Y88P_PV <- AddMetaData(RFP_Y88P_PV, metadata = doubletsRFP_Y88P_PV$doublet_scores, col.name = "doublet_scores")
RFP_Y88P_PV <- AddMetaData(RFP_Y88P_PV, metadata = doubletsRFP_Y88P_PV$predicted_doublets, col.name = "predicted_doublets")
head(RFP_Y88P_PV[[]])
table(RFP_Y88P_PV$predicted_doublets)

colnames(doubletsWT_GFP_PV) <- c("doublet_scores","predicted_doublets")
WT_GFP_PV <- AddMetaData(WT_GFP_PV, metadata = doubletsWT_GFP_PV$doublet_scores, col.name = "doublet_scores")
WT_GFP_PV <- AddMetaData(WT_GFP_PV, metadata = doubletsWT_GFP_PV$predicted_doublets, col.name = "predicted_doublets")
head(WT_GFP_PV[[]])
table(WT_GFP_PV$predicted_doublets)

colnames(doubletsWT_GFP_RFP_C109S_PV) <- c("doublet_scores","predicted_doublets")
WT_GFP_RFP_C109S_PV <- AddMetaData(WT_GFP_RFP_C109S_PV, metadata = doubletsWT_GFP_RFP_C109S_PV$doublet_scores, col.name = "doublet_scores")
WT_GFP_RFP_C109S_PV <- AddMetaData(WT_GFP_RFP_C109S_PV, metadata = doubletsWT_GFP_RFP_C109S_PV$predicted_doublets, col.name = "predicted_doublets")
head(WT_GFP_RFP_C109S_PV[[]])
table(WT_GFP_RFP_C109S_PV$predicted_doublets)

colnames(doubletsWT_GFP_RFP_D177A_PV) <- c("doublet_scores","predicted_doublets")
WT_GFP_RFP_D177A_PV <- AddMetaData(WT_GFP_RFP_D177A_PV, metadata = doubletsWT_GFP_RFP_D177A_PV$doublet_scores, col.name = "doublet_scores")
WT_GFP_RFP_D177A_PV <- AddMetaData(WT_GFP_RFP_D177A_PV, metadata = doubletsWT_GFP_RFP_D177A_PV$predicted_doublets, col.name = "predicted_doublets")
head(WT_GFP_RFP_D177A_PV[[]])
table(WT_GFP_RFP_D177A_PV$predicted_doublets)

colnames(doubletsWT_GFP_RFP_Y88P_PV) <- c("doublet_scores","predicted_doublets")
WT_GFP_RFP_Y88P_PV <- AddMetaData(WT_GFP_RFP_Y88P_PV, metadata = doubletsWT_GFP_RFP_Y88P_PV$doublet_scores, col.name = "doublet_scores")
WT_GFP_RFP_Y88P_PV <- AddMetaData(WT_GFP_RFP_Y88P_PV, metadata = doubletsWT_GFP_RFP_Y88P_PV$predicted_doublets, col.name = "predicted_doublets")
head(WT_GFP_RFP_Y88P_PV[[]])
table(WT_GFP_RFP_Y88P_PV$predicted_doublets)

colnames(doubletsWT_IRES_GFP_Del_IRES_mRuby3_PV) <- c("doublet_scores","predicted_doublets")
WT_IRES_GFP_Del_IRES_mRuby3_PV <- AddMetaData(WT_IRES_GFP_Del_IRES_mRuby3_PV, metadata = doubletsWT_IRES_GFP_Del_IRES_mRuby3_PV$doublet_scores, col.name = "doublet_scores")
WT_IRES_GFP_Del_IRES_mRuby3_PV <- AddMetaData(WT_IRES_GFP_Del_IRES_mRuby3_PV, metadata = doubletsWT_IRES_GFP_Del_IRES_mRuby3_PV$predicted_doublets, col.name = "predicted_doublets")
head(WT_IRES_GFP_Del_IRES_mRuby3_PV[[]])
table(WT_IRES_GFP_Del_IRES_mRuby3_PV$predicted_doublets)

colnames(doubletsWT_IRES_GFP_PV) <- c("doublet_scores","predicted_doublets")
WT_IRES_GFP_PV <- AddMetaData(WT_IRES_GFP_PV, metadata = doubletsWT_IRES_GFP_PV$doublet_scores, col.name = "doublet_scores")
WT_IRES_GFP_PV <- AddMetaData(WT_IRES_GFP_PV, metadata = doubletsWT_IRES_GFP_PV$predicted_doublets, col.name = "predicted_doublets")
head(WT_IRES_GFP_PV[[]])
table(WT_IRES_GFP_PV$predicted_doublets)

colnames(doubletsWT_IRES_mRuby3_MutPol_PV) <- c("doublet_scores","predicted_doublets")
WT_IRES_mRuby3_MutPol_PV <- AddMetaData(WT_IRES_mRuby3_MutPol_PV, metadata = doubletsWT_IRES_mRuby3_MutPol_PV$doublet_scores, col.name = "doublet_scores")
WT_IRES_mRuby3_MutPol_PV <- AddMetaData(WT_IRES_mRuby3_MutPol_PV, metadata = doubletsWT_IRES_mRuby3_MutPol_PV$predicted_doublets, col.name = "predicted_doublets")
head(WT_IRES_mRuby3_MutPol_PV[[]])
table(WT_IRES_mRuby3_MutPol_PV$predicted_doublets)

colnames(doubletsWT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV) <- c("doublet_scores","predicted_doublets")
WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV <- AddMetaData(WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV, metadata = doubletsWT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV$doublet_scores, col.name = "doublet_scores")
WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV <- AddMetaData(WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV, metadata = doubletsWT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV$predicted_doublets, col.name = "predicted_doublets")
head(WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV[[]])
table(WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV$predicted_doublets)
################################################################################
###Subset by sample according to VlnPlots and doublet scores
################################################################################
Mock_5h_PV = subset(Mock_5h_PV, nCount_RNA > 5000 & nCount_RNA < 40000 & nFeature_RNA >3000 & percent.ribo < 30 & percent.ribo > 10& percent.mt <20 & percent.mt >.01 & doublet_scores<0.57)

RFP_C109S_PV = subset(RFP_C109S_PV, nCount_RNA > 5000 & nFeature_RNA >2000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01 & doublet_scores<0.51)

WT_GFP_PV = subset(WT_GFP_PV, nCount_RNA > 5000 & nFeature_RNA >2000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.54)

WT_GFP_RFP_Y88P_PV = subset(WT_GFP_RFP_Y88P_PV, nCount_RNA > 4000 & nFeature_RNA >1500 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.55)

WT_GFP_RFP_D177A_PV = subset(WT_GFP_RFP_D177A_PV, nCount_RNA > 4000 & nFeature_RNA >2000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.50)

RFP_Y88P_PV = subset(RFP_Y88P_PV, nCount_RNA > 3000 & nFeature_RNA >2000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.47)

RFP_D177A_PV = subset(RFP_D177A_PV, nCount_RNA > 2000 & nFeature_RNA >1000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.52)

WT_GFP_RFP_C109S_PV = subset(WT_GFP_RFP_C109S_PV, nCount_RNA > 2000 & nFeature_RNA >2000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.46)

WT_IRES_GFP_PV = subset(WT_IRES_GFP_PV, nCount_RNA > 3000 & nFeature_RNA >2000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.45)

Del_IRES_mRuby3_PV = subset(Del_IRES_mRuby3_PV, nCount_RNA > 3000 & nFeature_RNA >2000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.45)

WT_IRES_GFP_Del_IRES_mRuby3_PV = subset(WT_IRES_GFP_Del_IRES_mRuby3_PV, nCount_RNA > 3000 & nFeature_RNA >2000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.53)

WT_IRES_mRuby3_MutPol_PV = subset(WT_IRES_mRuby3_MutPol_PV, nCount_RNA > 2000 & nFeature_RNA >1000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.57)

Del_IRES_mRuby3_MutPol_PV = subset(Del_IRES_mRuby3_MutPol_PV, nCount_RNA > 2000 & nFeature_RNA >1000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.57)

WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV = subset(WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV, nCount_RNA > 5000 & nFeature_RNA >2000 & percent.ribo < 30 & percent.ribo > 5& percent.mt <20 & percent.mt >.01& doublet_scores<0.54)

seurat_object_merged_filtered_list <- merge(RFP_C109S_PV, y = c(WT_GFP_PV, WT_GFP_RFP_Y88P_PV,
                                                                WT_GFP_RFP_D177A_PV,
                                                                RFP_Y88P_PV, RFP_D177A_PV, 
                                                                WT_GFP_RFP_C109S_PV, 
                                                                WT_IRES_GFP_PV, 
                                                                Del_IRES_mRuby3_PV,
                                                                WT_IRES_GFP_Del_IRES_mRuby3_PV,
                                                                WT_IRES_mRuby3_MutPol_PV,
                                                                Del_IRES_mRuby3_MutPol_PV,
                                                                WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV,
                                                                Mock_5h_PV))

table(seurat_object_merged_filtered_list@meta.data$predicted_doublets, seurat_object_merged_filtered_list$orig.ident)
################################################################################
#Evaluate batch effects. Mock infected was run with CVB3 and EVA71 treatments. IRES mutants and MutPol mutants with cotransfections and WT IRES were run together.
# Cotransfections of transcriptional mutants were run together with WT GFP.
################################################################################
exp.mat <- read.table(file = "~/lab_share/Projects/CM_kb/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
                      as.is = TRUE, row.names = 1)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


################################################################################
#### Rename WT per batch for possible correction
################################################################################


table(Idents(seurat_object_merged_filtered_list))

x <- NormalizeData(seurat_object_merged_filtered_list, verbose = FALSE)
x <- FindVariableFeatures(x, verbose = FALSE)
all.genes.x <- rownames(x)
x <- ScaleData(x, features = all.genes.x, verbose = FALSE)
x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes)
x <- RunPCA(x, features = c(s.genes, g2m.genes), verbose = FALSE)
PCAPlot(x)
merged_seurat <- ScaleData(x, features= all.genes.x, verbose = T, vars.to.regress = c("S.Score", "G2M.Score"))
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(merged_seurat), nfeatures.print = 10)

ElbowPlot(merged_seurat, ndims = 50)

merged_seurat <- FindNeighbors(merged_seurat, dims = 1:23, verbose = FALSE)
################################################################################
### Choose resolution
################################################################################

seurat_object <- FindClusters(merged_seurat, resolution = seq(0.1,1, by=0.1))

clustree(seurat_object)
DimPlot(seurat_object,group.by = 'RNA_snn_res.0.4', split.by = 'orig.ident')
DimPlot(seurat_object,group.by = 'RNA_snn_res.0.5', split.by = 'orig.ident')
DimPlot(seurat_object,group.by = 'RNA_snn_res.0.4')
Idents(seurat_object) <- seurat_object$RNA_snn_res.0.5

merged_seurat_0.5 <- RunUMAP(seurat_object, dims = 1:23, verbose = FALSE)
DimPlot(merged_seurat_0.5,group.by = 'RNA_snn_res.0.4', split.by = 'orig.ident', reduction = 'umap')
DimPlot(merged_seurat_0.5,group.by = 'RNA_snn_res.0.5', split.by = 'orig.ident', reduction = 'umap')
DimPlot(merged_seurat_0.5,group.by = 'RNA_snn_res.0.7', reduction = 'umap')
#####################################
### Visualize potential batch effects from replicate WT samples.
####################################
merged_seurat_0.5.cells <- Cells(merged_seurat_0.5)
Idents(object = merged_seurat_0.5, merged_seurat_0.5.cells[1561:3258]) <- 'Batch1 WT GFP'
Idents(object = merged_seurat_0.5, merged_seurat_0.5.cells[10583:12093]) <- 'Batch2 WT GFP'
WTs <- subset(merged_seurat_0.5, idents  = c('Batch2 WT GFP', 'Batch1 WT GFP'))
DimPlot(WTs, reduction = 'umap')
#####################################
### No significant batch effect found. Reset the Idents in the Suerat object.
####################################
Idents(merged_seurat_0.5) <- merged_seurat_0.5$RNA_snn_res.0.5
table(merged_seurat_0.5$orig.ident, merged_seurat_0.5$RNA_snn_res.0.5)
#####################################
### Save the object and metadata.
####################################
DimPlot(merged_seurat_0.5, reduction = 'umap', split.by = 'orig.ident')

saveRDS(object = merged_seurat_0.5, file = "/data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/mutsFinal.rds")

scissors.all_meta <- merged_seurat_0.5@meta.data
write.csv(scissors.all_meta, file = "/data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/scissors.metadata0.5.csv")
        
joined_seurat.markers <- FindAllMarkers(merged_seurat_0.5, only.pos = FALSE)
joined_seurat.markers %>%
group_by(cluster)
write.csv(joined_seurat.markers, file = "/data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/markers.0.5.csv")

#####################################
### Make a volcano plot per cluster.
####################################
library(ggrepel)
joined_seurat.markers$delabel <- NA
joined_seurat.markers$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
        
joined_seurat.markers$diffexpressed[joined_seurat.markers$avg_log2FC > 1 & joined_seurat.markers$p_val_adj < 0.05] <- "UP"
        
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
joined_seurat.markers$diffexpressed[joined_seurat.markers$avg_log2FC < -1 & joined_seurat.markers$p_val_adj < 0.05] <- "DOWN"
        
joined_seurat.markers$delabel[joined_seurat.markers$diffexpressed != "NO"] <- joined_seurat.markers$gene[joined_seurat.markers$diffexpressed != "NO"]
        
#Create a column for the shape manual
joined_seurat.markers$clusters <- factor(joined_seurat.markers$cluster)
        
p1 <- ggplot(joined_seurat.markers, aes(avg_log2FC, -log(p_val_adj,10), shape = factor(cluster))) + # -log10 conversion 
          geom_point(size=1) +
          xlab(expression("log"[2]*"(Fold Change)")) + 
          ylab(expression("-log"[10]*"Pvalue"))+
          geom_vline(xintercept=c(-1, 1), col="grey", linetype = 'dashed') +
          geom_hline(yintercept=-log10(0.05), col="grey", linetype = 'dashed')+
          coord_cartesian(ylim = c(0, 30000))+
          geom_point(data=subset(joined_seurat.markers, diffexpressed == "UP"), color="darkred")+
          geom_point(data=subset(joined_seurat.markers, diffexpressed == "DOWN"), color="darkblue")+
          #scale_color_brewer(palette="Paired") + theme(
          #legend.position = c(1, 0),
          #legend.justification = c("right", "bottom"),
          #legend.box.just = "right",
          #legend.margin = margin(6, 6, 6, 6)
          #)+
          labs(shape = 'Clusters') + geom_point(data=subset(joined_seurat.markers, diffexpressed == "NO"), color="grey")+
          scale_shape_manual(values=1:nlevels(joined_seurat.markers$clusters))+
          facet_wrap(~factor(cluster)) + labs(title = "Differential Expression across Clusters") +NoLegend()
        
p1 + geom_label_repel(size = 2, aes(label = delabel), max.overlaps = 1000)
    

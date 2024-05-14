################################################################################
### Load Libraries ###
################################################################################

library(Seurat)
library(tidyverse)

################################################################################
### Define File Paths ###
################################################################################

################################################################################
### Get list of CellRanger files
################################################################################


################################################################################
### Read in data
################################################################################

### Create Seurat Object

dataset_loc <- "~/lab_share/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/"
ids <- c("RFP_C109S_PV", "WT_GFP_PV", "WT_GFP_RFP_Y88P_PV", "WT_GFP_RFP_D177A_PV", "RFP_Y88P_PV", "RFP_D177A_PV", "WT_GFP_RFP_C109S_PV", "WT_IRES_GFP_PV","Del_IRES_mRuby3_PV","WT_IRES_GFP_Del_IRES_mRuby3_PV","WT_IRES_mRuby3_MutPol_PV","Del_IRES_mRuby3_MutPol_PV","WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV")
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
### Filter cells for % mito gene expression <=20 and # detected genes >= 1000 ###
################################################################################

seurat_object_merged_filtered <- subset(merged_seurat_object,
                                        subset = `percent.mt` <= 20 &
                                          nFeature_RNA >= 1000)

################################################################################
### Perform integration analysis ###
################################################################################


seurat_object_merged_filtered_list <- SplitObject(seurat_object_merged_filtered, split.by = "orig.ident")
#Remove the ridiculously huge chunk of code and instead make a for loop that find the infected threshold for virus_percentage and labels them according to low, high or not infected.
for (i in ids){
  
  halfmin_virus = min(seurat_object_merged_filtered_list[[i]]$percent.virus[seurat_object_merged_filtered_list[[i]]$percent.virus>0])
  
  percent.virus_log10 = log10(seurat_object_merged_filtered_list[[i]]$percent.virus + halfmin_virus)
  
  viruspercentage = mixtools::normalmixEM(seurat_object_merged_filtered_list[[i]]$percent.virus[percent.virus_log10 > -2.5],k = 2)
  #call infected any cell above the mu + 2x sigma =  ----
  
  seurat_object_merged_filtered_list[[i]]$InfectedStatus <- "Not_Infected"
  seurat_object_merged_filtered_list[[i]]$InfectedStatus[seurat_object_merged_filtered_list[[i]]$percent.virus > (viruspercentage$mu[1]+viruspercentage$sigma[1]*2)] <- "Infected"
  
  table(seurat_object_merged_filtered_list[[i]]$InfectedStatus)
  
  seurat_object_merged_filtered_list[[i]]$InfectedStatus_groups <- "Not_Infected"
  seurat_object_merged_filtered_list[[i]]$InfectedStatus_groups[seurat_object_merged_filtered_list[[i]]$percent.virus > viruspercentage$mu[2]] <- "High"
  seurat_object_merged_filtered_list[[i]]$InfectedStatus_groups[seurat_object_merged_filtered_list[[i]]$percent.virus < viruspercentage$mu[2]] <- "Low"
  seurat_object_merged_filtered_list[[i]]$InfectedStatus_groups[seurat_object_merged_filtered_list[[i]]$percent.virus < (viruspercentage$mu[1]+viruspercentage$sigma[1]*2)]<- "Not_Infected" #annotate the threshold in the meta.data object
}
################################################################################
### Split the object for individual QC and read in doublets files
################################################################################
for (i in ids){
  assign(value = seurat_object_merged_filtered_list[[i]], x= paste(i))
}

for (i in ids)
{
  assign(value = read.table(file.path(dataset_loc,i,"outs/filtered_feature_bc_matrix/Doublet_scores.tsv"), header = F), x= paste0("doublets",i))
  
}
################################################################################
#### Add doublet scores metadata from scrublet output
################################################################################
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
### Subset samples by VlnPlots and doublet scores
################################################################################
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
                                         WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV))
table(seurat_object_merged_filtered_list@meta.data$predicted_doublets, seurat_object_merged_filtered_list$orig.ident)
seurat_object_merged_filtered_list <- SplitObject(seurat_object_merged_filtered_list, split.by = 'orig.ident')

#### normalize each dataset individually and find 2000 variable features ####
seurat_object_merged_filtered_list <- lapply(X = seurat_object_merged_filtered_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
################################################################################
### Perform Integration ###
#### Select integration features ####
################################################################################
features <- SelectIntegrationFeatures(object.list = seurat_object_merged_filtered_list)
seurat_object_merged_filtered_list <- lapply(X = seurat_object_merged_filtered_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

##### Find integration anchors #### 
anchors <- FindIntegrationAnchors(object.list = seurat_object_merged_filtered_list,
                                  anchor.features = features,
                                  reduction = "rpca",
                                  dims = 1:50)

#### Integrate datasets ####
seurat_object_integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

################################################################################
### Cell cycle scoring, dimensionality reduction and clustering ###
################################################################################
exp.mat <- read.table(file = "~/lab_share/Projects/CM_kb/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
                      as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_object_integrated <- CellCycleScoring(seurat_object_integrated, s.features = s.genes, g2m.features = g2m.genes)
DefaultAssay(seurat_object_integrated) <- "integrated"
seurat_object_integrated <- ScaleData(seurat_object_integrated, verbose = FALSE)
seurat_object_integrated <- RunPCA(seurat_object_integrated, npcs = 30, verbose = FALSE)
seurat_object_integrated <- RunUMAP(seurat_object_integrated, reduction = "pca", dims = 1:30)
seurat_object_integrated <- FindNeighbors(seurat_object_integrated, reduction = "pca", dims = 1:30)

################################################################################
### cluster cells from 0.1 to 1 resolution ###
################################################################################
for(i in 1:length(seq(0.1, 1, 0.1))){
  seurat_object_integrated <- FindClusters(seurat_object_integrated, resolution = seq(0.1, 1, 0.05)[i])
}
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.1', label =TRUE)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.15', label =TRUE)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.2', label =TRUE)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.25', label =TRUE)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.3', label =TRUE)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.35', label =TRUE)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.4', label =TRUE)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.45', label =TRUE)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.5', label =TRUE)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.55', label =TRUE)
clustree(seurat_object_integrated, assay = 'integrated')

################################################################################
### Dimensionality reduction and clustering of RNA assay ###
################################################################################
Idents(seurat_object_integrated) <- 'integrated_snn_res.0.4'
DefaultAssay(seurat_object_integrated) <- 'RNA'
joined_seurat <- JoinLayers(seurat_object_integrated)

mutants.all_meta <- joined_seurat@meta.data
write.csv(mutants.all_meta, file = "/data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/ScissorsMetadata0.4.csv")
joined_seurat.markers <- FindAllMarkers(joined_seurat, only.pos = FALSE, logfc.threshold = 0.5)
joined_seurat.markers %>%
  group_by(cluster)
write.csv(joined_seurat.markers, file = "/data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/markers.0.4.csv")

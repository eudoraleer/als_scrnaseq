#########################################################################################################
# Title: T cell responses at diagnosis of amyotrophic lateral sclerosis predict disease progression
# Data: 10X 5' scRNA-Sequencing and immune profiling
# Creator: Lu Pan, lu.pan@ki.se
# Date last modified: 2020-04-05
# Content: data processing steps for Fig. 3,4,5, S3
#########################################################################################################

#########################################################################################################
# 1. Data cleaning and preprocessing for scRNA-Seq
library(Seurat)
library(harmony)
library(DoubletFinder)

setwd("./")
m1 <- readRDS("ALS_Meta_scRNA_Seq.RDS")
data <- NULL

for(j in 1:nrow(m1)){
  x <- Read10X(m1$Dir[j])
  x <- CreateSeuratObject(x, min.cells = 3, min.features = 200)
  x <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- x(x, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 5)
  x$GROUP <- m1$GROUP[j]
  x$SUBGROUP <- m1$SUBGROUP[j]
  x$SAMPLE_ID <- m1$SAMPLE_ID[j]
  x$CELL <- row.names(x@meta.data)
  x$BATCH <- m1$BATCH[j]
  data[[j]] <- x
}

data <- merge(data[1], data[c(2:length(data))])
data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data, features = VariableFeatures(data))
data <- RunPCA(data, features = VariableFeatures(data))
data <- RunHarmony(data, "BATCH")
data <- FindNeighbors(data, reduction = "harmony", dims = 1:30)
data <- FindClusters(data, resolution = 3.5)
degs <- FindAllMarkers(data, min.pct = 0.25, logfc.threshold = 0.25)
degs <- degs[which(degs$avg_log2FC > 0 & degs$p.adjust < 0.05),]
saveRDS(data, "data1.RDS")
saveRDS(degs, "degs.RDS")
#########################################################################################################
# 2. Post-manual annotation processing
data <- readRDS("data1.RDS")
annot <- read.table("ALS_Annotations.txt", sep = "\t", header = T)
data$CELL_TYPE <- annot[match(data$seurat_clusters, annot$CLUSTER),"CELL_TYPE"]
data$CELL_TYPE_LEVEL1 <- annot[match(data$seurat_clusters, annot$CLUSTER),"CELL_TYPE_LEVEL1"]
homotypic.prop <- modelHomotypic(data$CELL_TYPE)
nExp_poi <- round(0.075*nrow(data@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_2165", sct = FALSE)
data <- subset(data, subset = DF.classifications_0.25_0.09_2165 == "Singlet")
data <- RunTSNE(data, reduction = "harmony", dims = 1:30)
saveRDS(data, "data2.RDS")

#########################################################################################################
# 3. Data cleaning and processing for scImmune profiling
library(scRepertoire)
data <- readRDS("data2.RDS")
m2 <- readRDS("ALS_Meta_scImmune.RDS")

imm_data <- NULL
for(i in 1:nrow(m2)){
  x <- read.csv(m2$FILE[i], header = T)
  contig_list[[i]] <- x
}

names(contig_list) <- m2$SAMPLE_ID

tab <- combineTCR(contig_list, 
                  samples = names(contig_list), 
                  ID = m2[match(names(contig_list), m2$SAMPLE_ID),"GROUP"], cells ="T-AB")
tab <- addVariable(tab, name = "GROUP",
                   variables = m2[match(names(contig_list), m2$SAMPLE_ID),"GROUP"])
tab <- addVariable(tab, name = "SUBGROUP", 
                   variables = m2[match(names(contig_list), m2$SAMPLE_ID),"SUBGROUP"])

data$OID <- row.names(data@meta.data)
data$NID <- paste(data$SAMPLE_ID,data$GROUP,data$CELL, sep = "_")
row.names(data@meta.data) <- data$NID
colnames(data@assays$RNA@counts) <- data$NID
data <- combineExpression(tab, data, cloneCall="gene+nt",  proportion = FALSE,
                          clonecall_types=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
slot(data, "meta.data")$cloneType <- factor(slot(data, "meta.data")$cloneType, 
                                            levels = c("Hyperexpanded (100 < X <= 500)", 
                                                       "Large (20 < X <= 100)", 
                                                       "Medium (5 < X <= 20)", 
                                                       "Small (1 < X <= 5)", 
                                                       "Single (0 < X <= 1)", NA))
row.names(data@meta.data) <- data$OID
colnames(data@assays$RNA@counts) <- data$OID
saveRDS(data, "data3.RDS")

#########################################################################################################
# 4. Detectable TCR only
data$TCR <- ifelse(!is.na(data$`TCR_gene+nt_ClonalType`), "TCR", "NOT_DETECTED/NON-TCR")
tcr_data <- subset(data, subset = TCR == "TCR")
tcr_data <- RunHarmony(tcr_data, "BATCH")
tcr_data <- RunTSNE(tcr_data, reduction = "harmony", dims = 1:30)
saveRDS(tcr_data, "tcr_data.RDS")

#########################################################################################################


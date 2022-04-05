#########################################################################################################
# Title: T cell responses at diagnosis of amyotrophic lateral sclerosis predict disease progression
# Data: 10X 5' scRNA-Sequencing and immune profiling
# Creator: Lu Pan, lu.pan@ki.se
# Date last modified: 2020-04-05
# Content: Plots / Data for Fig. 3,4,5, S3
#########################################################################################################

#########################################################################################################
# Fig. 3
library(Seurat)
library(ComplexHeatmap)
setwd("./")
data <- readRDS("data3.RDS")

# Fig. 3A tSNE
fig3a <- DimPlot(data, reduction = "tsne", split.by = "GROUP", group.by = "CELL_TYPE", label.size = 3, repel = T)

# Fig. 3B Proportion Data
cellcounts <- data.frame(table(data@meta.data[,c("CELL_TYPE","SAMPLE_ID")]))
cellcounts <- dcast(cellcounts, CELL_TYPE ~ SAMPLE_ID, value.var = "Freq")
row.names(cellcounts) <- cellcounts$CELL_TYPE
cellcounts <- cellcounts[,grep("CELL_TYPE", colnames(cellcounts), ignore.case = T, invert = T)]
groups <- data@meta.data[match(colnames(cellcounts), data$SAMPLE_ID),"GROUP"]
subgroups <- data@meta.data[match(colnames(cellcounts), data$SAMPLE_ID),"SUBGROUP"]

N <- table(data$GROUP)
plotx <- NULL
fig3b <- NULL
for(j in 1:nrow(cellcounts)){
  g1 <- sum(cellcounts[j,which(groups == "ALS")])
  g2 <- sum(cellcounts[j,which(groups != "ALS")])
  if(g1 > 5 & g2 > 5){
    p <- prop.test(x = c(g1, g2), n = c(N[which(names(N) == "ALS")], N[which(names(N) != "ALS")]), alternative = "two.sided")
    fc <- mean(as.numeric(as.character(cellcounts[j,groups == "ALS"]/as.numeric(colSums(cellcounts))[groups == "ALS"])))/mean(as.numeric(as.character(cellcounts[j,groups != "ALS"]/as.numeric(colSums(cellcounts))[groups != "ALS"])))
    fig3b <- rbind(fig3b, data.frame(COMPARISON = "ALS vs Non-Inflammatory Control\n",
                                         P = p$p.value,
                                         FC = fc,
                                         CELL_TYPE = row.names(cellcounts)[j],
                                         LOG10FC = log10(fc),
                                         NEG_LOG10P = -log10(p$p.value),
                                         COMPARE1 = g2,
                                         COMPARE2 = g1))
  }
}

# Fig. 3C Median Expressions
Idents(data) <- "GROUP"
celltypes <- unique(data$CELL_TYPE)

fig3c_data <- NULL
for(i in 1:length(celltypes)){
  current <- subset(data, subset = CELL_TYPE == celltypes[i])
  current <- data.frame(log1p(AverageExpression(current, verbose = FALSE)$RNA))
  data_median <- plyr::rbind.fill(data_median, data.frame(CELL_TYPE = celltypes[i], GENE = row.names(current), current))
}

#########################################################################################################
# Fig. 4
tcr_data <- readRDS("tcr_data.RDS")
tcr_data <- subset(tcr_data, subset = CELL_TYPE_LEVEL1 == "T") # kept only concordant cells with both TCR and T-cell types annotation
tcr_data$CELL_TYPE <- factor(tcr_data$CELL_TYPE, levels = sort(unique(tcr_data$CELL_TYPE)))

# Fig. 4A clonal expansion overlay
fig4a <- clonalOverlay(tcr_data, reduction = "tsne", freq.cutpoint = 30, bins = 20, facet = "GROUP")

# Fig. 4B Proportions data
groups <- unique(tcr_data$SUBGROUP)
Idents(tcr_data) <- "CELL_TYPE"
fig4b_data <- NULL
for(i in 1:length(groups)){
  plotx <- subset(tcr_data, subset = SUBGROUP == groups[i])
  temp <- data.frame(table(plotx$CELL_TYPE))
  plotx$CELL_TYPE <- paste(plotx$CELL_TYPE, " (", temp[match(plotx$CELL_TYPE, temp$Var1),"Freq"], ")", sep = "")
  
  fig4b_data <- rbind(fig4b_data, data.frame(GROUP = groups[i], occupiedscRepertoire(plotx, x.axis = "cluster", proportion = T, exportTable = TRUE)))
}

# Fig. 4C Proportions data
celltypes <- unique(fig4b_data$cluster)
fig4c_data <- NULL

for(i in 1:length(celltypes)){
  current <- summary[summary$cluster == celltypes[i],]
  current <- dcast(current, cloneType ~ GROUP, value.var = "value")
  current[is.na(current)] <- 0
  row.names(current) <- current$cloneType
  current <- current[,grep("cloneType", colnames(current), ignore.case = T, invert = T)]
  chisq <- chisq.test(current, simulate.p.value = T)
  fig4c_data <- rbind(fig4c_data, data.frame(CELL_TYPE = celltypes[i], P = chisq$p.value))
}

#########################################################################################################
# Fig. 5
selected <- c("GATA3","EOMES","TBX21","RORC")
fig5a <- FeaturePlot(tcr_data, features = selected, keep.scale = "all", split.by = "SUBGROUP")

#########################################################################################################
# Fig. S3
data <- readRDS("data3.RDS")
selected <- c("CD3E","CD4","CD8A","TRDC","FOXP3",
              "CCR7","CCL5","GNLY","PRF1","GZMA",
              "FCGR3A","CD14","CD1C","XCR1","CLEC4C",
              "NCAM1","GNG11","CD19","MKI67","IGHM",
              "IGHD","IGHG1","IGHG2","IGHG3")
              
figs3 <- FeaturePlot(data, features = selected, keep.scale = "all")

#########################################################################################################




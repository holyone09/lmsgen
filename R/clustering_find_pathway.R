GSM2257302_expr_matrix<-readRDS(gzfile("./dataset/GSM2257302_exp_matrix.rds"))

library(Seurat)
library(dplyr)
library(Matrix)
GSM2257302_pbmc<- CreateSeuratObject(raw.data = GSM2257302_expr_matrix)
GSM2257302_pbmc <- NormalizeData(object = GSM2257302_pbmc)
GSM2257302_pbmc <- ScaleData(object = GSM2257302_pbmc)
GSM2257302_pbmc <- FindVariableGenes(object = GSM2257302_pbmc, do.plot = T)
GSM2257302_pbmc <- RunPCA(object = GSM2257302_pbmc, pc.genes = GSM2257302_pbmc@var.genes, do.print = TRUE, pcs.print = 1:5,genes.print = 5)
##find clusters
GSM2257302_pbmc <- FindClusters(object = GSM2257302_pbmc, reduction.type = "pca", dims.use = 1:12, resolution = 0.6, print.output = 0, save.SNN = T,force.recalc=TRUE)
GSM2257302_pbmc <- RunTSNE(object =GSM2257302_pbmc, dims.use = 1:12, do.fast = T)

GSM2257302_cur.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
GSM2257302_new.cluster.ids <- c("APS+MPS", "Somitomere", "PXM", "cDM", "Sclerotome","LatM","ESC","Earlysomite")
GSM2257302_pbmc@ident <- plyr::mapvalues(GSM2257302_pbmc@ident, from =GSM2257302_cur.cluster.ids, to = GSM2257302_new.cluster.ids)

TSNEPlot(GSM2257302_pbmc,do.label = F)
##gene analysis
GSM2257302_pbmc.markers_rex_all <- FindAllMarkers(GSM2257302_pbmc, only.pos = TRUE, min.pct = 0, thresh.use = 0)
GSM2257302_all_rex_all<-GSM2257302_pbmc.markers_rex_all %>% group_by(cluster)
GSM2257302_top100 <- GSM2257302_pbmc.markers_rex_all %>% group_by(cluster) %>% top_n(100, avg_logFC)
pdf("result/GSM2257302_top100.pdf")
DoHeatmap(object = GSM2257302_pbmc, genes.use =GSM2257302_top100$gene, slim.col.label = TRUE, remove.key = TRUE,cex.row = 1,cex.col=1,group.cex = 5.5,
          ,group.order = c("ESC","APS+MPS","PXM","LatM","Earlysomite","cDM","Somitomere","Sclerotome"))
dev.off()

#path finder
library(pathfindR)
GSM2257302_top100_input<-as.data.frame(GSM2257302_top100[,c(7,2,5)])
GSM2257302_top100_output <- run_pathfindR(GSM2257302_top100_input)

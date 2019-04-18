source("./cig_tools.R")
GSM2257302_data<-read.delim("dataset/GSM2257302_All_samples_sc_tpm.txt",header =T ,sep="\t",stringsAsFactors =F)
GSM2257302_data<-subset(GSM2257302_data,GSM2257302_data[,1]!="")
GSM2257302_data.uni <- unique( GSM2257302_data[ , 1 ] )
GSM2257302_data.info<-GSM2257302_data[,1:2]
GSM2257302_indx<-vector()

for(i in 1:length(GSM2257302_data.uni)){
  tmp<-as.character(subset(GSM2257302_data.info,GSM2257302_data$geneSymbol==GSM2257302_data.uni[i])$geneID)
  GSM2257302_indx<-c(GSM2257302_indx,tmp[1])
}
rownames(GSM2257302_data)<-GSM2257302_data[,2]

GSM2257302_expr_matrix<-as.matrix(GSM2257302_data[,3:ncol(GSM2257302_data)])
#rownames(GSM2257302_expr_matrix)<-GSM2257302_data[,2]
GSM2257302_expr_matrix<-GSM2257302_expr_matrix[GSM2257302_indx,]
GSM2257302_genes<-GSM2257302_data[GSM2257302_indx,1:2]
GSM2257302_samples<-colnames(GSM2257302_data)[c(-1,-2)]
GSM2257302_pseudo_times<-vector()
for(i in 1:length(GSM2257302_samples)){
  GSM2257302_pseudo_times<-cbind(GSM2257302_pseudo_times,strsplit(GSM2257302_samples,"[.]")[[i]][1])
}
GSM2257302_pseudo_times<-as.vector(GSM2257302_pseudo_times)

hmcol = colorRampPalette(c("blue","white","red"))(n = 50)

#candi_scl<-GSM2257302_scl$gene
#candi_scl<-c("COL3A1","COL2A1","TBX15","FAM110B","RGMB","FBN2","LUM","PTN","PAM")
#candi_scl<-c("PTN","ALK","PTPRB","CTNNB1","SRC","CTTN","ACTB","SDC3","NCAN","DCN","BCAN","SDC1","PTPRZ1")
#candi_scl<-as.character(read.delim("GO_0043589_skin_morphogenesis.txt",header=F)[,1])
candi_scl<-unique(gene_info_ra$SYMBOL)
#candi_scl<-c("COL3A1","CYP1B1","DCN","IGFBP5","SERPINF1","PTN","SFRP1","SEMA3A","LIMCH1","MCTP1")
#candi_scl<-c("FOXC2","FN1","HAS2","IGFBP5","LAMB1","PDGFRB","TGFB2","SEMA3A","FOXP1","PHPT1","MALAT1")
#candi_scl<-c("COL3A1","IGFBP5","NKX3-1","PTCH1","SFRP1","SOX4","SOX11","TGFB2","SEMA3A","FOXP1","HHIP")
#candi_scl<-c("RUNX1T1","COL3A1","EFNA5","IGFBP5","PTCH1","PTN","PTPRS","SFRP1","SRSF6","SOX11","TGFB2","SEMA3A")
candi_scl<-c("COL2A1","FBN2","FN1","LAMB1","PTCH1","SFRP1","SOX4","SOX11","TBX15","TGFB2","LRIG3")
scl_gene<-candi_scl[candi_scl %in% GSM2257302_genes$geneSymbol]

GSM2257302_scl_indx<-get_index(GSM2257302_genes,scl_gene)
GSM2257302_scl_gene_list<-list(p_gene=GSM2257302_scl_indx)
GSM2257302_state<-c("H7hESC","APS","DLL1PXM","D2_25somitomere","Earlysomite","cDM","Sclerotome")

GSM2257302_scl_sep<-gen_matrix(GSM2257302_expr_matrix,GSM2257302_scl_gene_list,GSM2257302_state,GSM2257302_pseudo_times)

GSM2257302_esc_scl.dist<-gen_interacting_mat(GSM2257302_scl_sep[[1]],GSM2257302_genes)
GSM2257302_aps_scl.dist<-gen_interacting_mat(GSM2257302_scl_sep[[2]],GSM2257302_genes)
GSM2257302_pxm_scl.dist<-gen_interacting_mat(GSM2257302_scl_sep[[3]],GSM2257302_genes)
GSM2257302_som_scl.dist<-gen_interacting_mat(GSM2257302_scl_sep[[4]],GSM2257302_genes)
GSM2257302_esom_scl.dist<-gen_interacting_mat(GSM2257302_scl_sep[[5]],GSM2257302_genes)
GSM2257302_dom_scl.dist<-gen_interacting_mat(GSM2257302_scl_sep[[6]],GSM2257302_genes)
GSM2257302_scl_scl.dist<-gen_interacting_mat(GSM2257302_scl_sep[[7]],GSM2257302_genes)
num_gene<-c(nrow(GSM2257302_esc_scl.dist$dist_mat),nrow(GSM2257302_aps_scl.dist$dist_mat),
            nrow(GSM2257302_pxm_scl.dist$dist_mat),nrow(GSM2257302_som_scl.dist$dist_mat),
            nrow(GSM2257302_esom_scl.dist$dist_mat),nrow(GSM2257302_dom_scl.dist$dist_mat),
            nrow(GSM2257302_scl_scl.dist$dist_mat))
names(num_gene)<-c("GSM2257302_esc_scl.dist","GSM2257302_aps_scl.dist",
                   "GSM2257302_pxm_scl.dist","GSM2257302_som_scl.dist",
                   "GSM2257302_esom_scl.dist","GSM2257302_dom_scl.dist",
                   "GSM2257302_scl_scl.dist")
idx_grp<-names(num_gene[num_gene==max(num_gene)])
smt_indx<-rownames(GSM2257302_scl_scl.dist$dist_mat)[GSM2257302_scl_scl.dist$clust_info$order]
#esc
library(pheatmap)
#esc_scl_dist_mat<-GSM2257302_esc_scl.dist$dist_mat
#esc_scl_genes<-smt_indx[smt_indx %in% rownames(esc_scl_dist_mat)]
#esc_scl<-pheatmap(esc_scl_dist_mat[esc_scl_genes,esc_scl_genes],color  = rev(hmcol),
esc_scl_dist_mat<-intep_mat(smt_indx,GSM2257302_esc_scl.dist$dist_mat)
esc_scl<-pheatmap(esc_scl_dist_mat[smt_indx,smt_indx],color  = rev(hmcol),
                   main="embryonic morphogenesis in ESC",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#aps
library(pheatmap)
#aps_scl_dist_mat<-GSM2257302_aps_scl.dist$dist_mat
#aps_scl_genes<-smt_indx[smt_indx %in% rownames(aps_scl_dist_mat)]
#aps_scl<-pheatmap(aps_scl_dist_mat[aps_scl_genes,aps_scl_genes],color  = rev(hmcol),
aps_scl_dist_mat<-intep_mat(smt_indx,GSM2257302_aps_scl.dist$dist_mat)
aps_scl<-pheatmap(aps_scl_dist_mat[smt_indx,smt_indx],color  = rev(hmcol),
                   main="embryonic morphogenesis in APS",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#pxm
library(pheatmap)
#pxm_scl_dist_mat<-GSM2257302_pxm_scl.dist$dist_mat
#pxm_scl_genes<-smt_indx[smt_indx %in% rownames(pxm_scl_dist_mat)]
#pxm_scl<-pheatmap(pxm_scl_dist_mat[pxm_scl_genes,pxm_scl_genes],color  = rev(hmcol),
pxm_scl_dist_mat<-intep_mat(smt_indx,GSM2257302_pxm_scl.dist$dist_mat)
pxm_scl<-pheatmap(pxm_scl_dist_mat[smt_indx,smt_indx],color  = rev(hmcol),
                   main="embryonic morphogenesis in PXM",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#somitomere
library(pheatmap)
som_scl_dist_mat<-GSM2257302_som_scl.dist$dist_mat
som_scl_genes<-smt_indx[smt_indx %in% rownames(som_scl_dist_mat)]
som_scl<-pheatmap(som_scl_dist_mat[som_scl_genes,som_scl_genes],color  = rev(hmcol),
                   main="embryonic morphogenesis in Somitomere",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#Earlysomite
library(pheatmap)
#esom_scl_dist_mat<-GSM2257302_esom_scl.dist$dist_mat
#esom_scl_genes<-smt_indx[smt_indx %in% rownames(esom_scl_dist_mat)]
#esom_scl<-pheatmap(esom_scl_dist_mat[esom_scl_genes,esom_scl_genes],color  = rev(hmcol),
esom_scl_dist_mat<-intep_mat(smt_indx,GSM2257302_esom_scl.dist$dist_mat)
esom_scl<-pheatmap(esom_scl_dist_mat[smt_indx,smt_indx],color  = rev(hmcol),
                    main="embryonic morphogenesis in Earlysomite",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)

#Dermomyotome
library(pheatmap)
#dom_scl_dist_mat<-GSM2257302_dom_scl.dist$dist_mat
#dom_scl_genes<-smt_indx[smt_indx %in% rownames(dom_scl_dist_mat)]
#dom_scl<-pheatmap(dom_scl_dist_mat[dom_scl_genes,dom_scl_genes],color  = rev(hmcol),
dom_scl_dist_mat<-intep_mat(smt_indx,GSM2257302_dom_scl.dist$dist_mat)
dom_scl<-pheatmap(dom_scl_dist_mat[smt_indx,smt_indx],color  = rev(hmcol),
                   main="embryonic morphogenesis in Dermomyotome",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)

#Sclerotome
library(pheatmap)
scl_scl_dist_mat<-GSM2257302_scl_scl.dist$dist_mat
scl_scl_genes<-smt_indx[smt_indx %in% rownames(scl_scl_dist_mat)]
scl_scl<-pheatmap(scl_scl_dist_mat[scl_scl_genes,scl_scl_genes],color  = rev(hmcol),
                   main="embryonic morphogenesis in Sclerotome",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)

library(gridExtra)
#grid.arrange(esc_scl[[4]],aps_scl[[4]],pxm_scl[[4]],som_scl[[4]],esom_scl[[4]],dom_scl[[4]],scl_scl[[4]],nrow=1)
pdf("/Volumes/Home-IRP$/publish/new_pub/final_codes/heatmaps/GO_0048598_embryonic_morphogenesis.pdf")
grid.arrange(esc_scl[[4]],aps_scl[[4]],nrow=1)
grid.arrange(pxm_scl[[4]],som_scl[[4]],nrow=1)
grid.arrange(esom_scl[[4]],dom_scl[[4]],nrow=1)
grid.arrange(scl_scl[[4]],nrow=1)
#grid.arrange(esom_scl[[4]],scl_scl[[4]],nrow=1)
dev.off()

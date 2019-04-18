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
#candi_wnt<-as.character(read.delim("wnt_geneset.txt",header=F)[,1])
candi_wnt<-c("RHOA","CCND1","CAMK2D","CCND2","CSNK1A1","CSNK2A1","CSNK2A2","CTBP1","CTNNB1","SERPINF1","PPP3CA",
             "SFRP1","SFRP2","TCF7L2","TP53","WNT5A","WNT8A","CER1","DKK1","DAAM1","BAMBI","CACYBP","DKK4","LEF1",
             "CHD8","SENP2","TBL1XR1","TCF7L1","PRICKLE2")
wnt_gene<-candi_wnt[candi_wnt %in% GSM2257302_genes$geneSymbol]

GSM2257302_wnt_indx<-get_index(GSM2257302_genes,wnt_gene)
GSM2257302_gene_list<-list(p_gene=GSM2257302_wnt_indx)
#GSM2257302_state<-c("H7hESC","APS","DLL1PXM","D2_25somitomere","Earlysomite","Sclerotome")
GSM2257302_state<-c("H7hESC","APS","DLL1PXM","D2_25somitomere","Earlysomite","cDM","Sclerotome")

GSM2257302_wnt_sep<-gen_matrix(GSM2257302_expr_matrix,GSM2257302_gene_list,GSM2257302_state,GSM2257302_pseudo_times)

#GSM2257302_esc_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[1]],GSM2257302_genes)
#GSM2257302_aps_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[2]],GSM2257302_genes)
#GSM2257302_pxm_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[3]],GSM2257302_genes)
#GSM2257302_som_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[4]],GSM2257302_genes)
#GSM2257302_esom_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[5]],GSM2257302_genes)
#GSM2257302_dom_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[6]],GSM2257302_genes)
#GSM2257302_scl_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[7]],GSM2257302_genes)

#num_gene<-c(nrow(GSM2257302_esc_wnt.dist$dist_mat),nrow(GSM2257302_aps_wnt.dist$dist_mat),
#            nrow(GSM2257302_pxm_wnt.dist$dist_mat),nrow(GSM2257302_som_wnt.dist$dist_mat),
#            nrow(GSM2257302_esom_wnt.dist$dist_mat),nrow(GSM2257302_dom_wnt.dist$dist_mat),
#            nrow(GSM2257302_scl_wnt.dist$dist_mat))
#names(num_gene)<-c("GSM2257302_esc_wnt.dist","GSM2257302_aps_wnt.dist",
#                   "GSM2257302_pxm_wnt.dist","GSM2257302_som_wnt.dist",
#                   "GSM2257302_esom_wnt.dist","GSM2257302_dom_wnt.dist",
#                   "GSM2257302_scl_scl.dist")
#idx_grp<-names(num_gene[num_gene==max(num_gene)])
#smt_indx<-rownames(GSM2257302_esc_wnt.dist$dist_mat)[GSM2257302_esc_wnt.dist$clust_info$order]


#esc
GSM2257302_esc_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[1]],GSM2257302_genes)
library(pheatmap)
esc_wnt_order<-GSM2257302_esc_wnt.dist$clust_info$order
esc_wnt_dist_mat<-GSM2257302_esc_wnt.dist$dist_mat
esc_wnt_indx<-rownames(esc_wnt_dist_mat)[esc_wnt_order]
esc_wnt<-pheatmap(esc_wnt_dist_mat[esc_wnt_order,esc_wnt_order],color  = rev(hmcol),
                   main="WNT in ESC",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#aps
GSM2257302_aps_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[2]],GSM2257302_genes)
library(pheatmap)
aps_wnt_dist_mat<-GSM2257302_aps_wnt.dist$dist_mat
aps_wnt_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(aps_wnt_dist_mat)]

aps_wnt<-pheatmap(aps_wnt_dist_mat[aps_wnt_genes,aps_wnt_genes],color  = rev(hmcol),
                   main="WNT in APS",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#pxm
GSM2257302_pxm_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[3]],GSM2257302_genes)
library(pheatmap)
pxm_wnt_dist_mat<-GSM2257302_pxm_wnt.dist$dist_mat
pxm_wnt_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(pxm_wnt_dist_mat)]
pxm_wnt<-pheatmap(pxm_wnt_dist_mat[pxm_wnt_genes,pxm_wnt_genes],color  = rev(hmcol),
                   main="WNT in PXM",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#somitomere
GSM2257302_som_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[4]],GSM2257302_genes)
library(pheatmap)
som_wnt_dist_mat<-GSM2257302_som_wnt.dist$dist_mat
som_wnt_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(som_wnt_dist_mat)]
som_wnt<-pheatmap(som_wnt_dist_mat[som_wnt_genes,som_wnt_genes],color  = rev(hmcol),
                   main="WNT in Somitomere",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#Earlysomite
GSM2257302_esom_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[5]],GSM2257302_genes)
library(pheatmap)
esom_wnt_dist_mat<-GSM2257302_esom_wnt.dist$dist_mat
esom_wnt_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(esom_wnt_dist_mat)]
esom_wnt<-pheatmap(esom_wnt_dist_mat[esom_wnt_genes,esom_wnt_genes],color  = rev(hmcol),
                    main="WNT in Earlysomite",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)

#Dermomyotome
GSM2257302_dom_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[6]],GSM2257302_genes)
library(pheatmap)
dom_wnt_dist_mat<-GSM2257302_dom_wnt.dist$dist_mat
dom_wnt_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(dom_wnt_dist_mat)]
dom_wnt<-pheatmap(dom_wnt_dist_mat[dom_wnt_genes,dom_wnt_genes],color  = rev(hmcol),
                   main="WNT in Dermomyotome",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)

#Sclerotome
GSM2257302_scl_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[7]],GSM2257302_genes)
library(pheatmap)
scl_wnt_dist_mat<-GSM2257302_scl_wnt.dist$dist_mat
scl_wnt_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(scl_wnt_dist_mat)]
scl_wnt<-pheatmap(scl_wnt_dist_mat[scl_wnt_genes,scl_wnt_genes],color  = rev(hmcol),
                   main="WNT in Sclerotome",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)

library(gridExtra)
#grid.arrange(esc_scl[[4]],aps_scl[[4]],pxm_scl[[4]],som_scl[[4]],esom_scl[[4]],dom_scl[[4]],scl_scl[[4]],nrow=1)
pdf("/Volumes/Home-IRP$/publish/new_pub/final_codes/heatmaps/wnt_signalingpathway_final.pdf")
grid.arrange(esc_wnt[[4]],aps_wnt[[4]],nrow=1)
grid.arrange(pxm_wnt[[4]],som_wnt[[4]],nrow=1)
grid.arrange(esom_wnt[[4]],dom_wnt[[4]],nrow=1)
#grid.arrange(esom_wnt[[4]],scl_wnt[[4]],nrow=1)
grid.arrange(scl_wnt[[4]],nrow=1)
dev.off()

##tgfb
candi_tgfb<-as.character(read.delim("TGF_beta_geneset.txt",header=F)[,1])
tgfb_gene<-candi_tgfb[candi_tgfb %in% GSM2257302_genes$geneSymbol]

GSM2257302_tgfb_indx<-get_index(GSM2257302_genes,tgfb_gene)
GSM2257302_tgfb_gene_list<-list(p_gene=GSM2257302_tgfb_indx)
GSM2257302_state<-c("H7hESC","APS","DLL1PXM","D2_25somitomere","Earlysomite","cDM","Sclerotome")

GSM2257302_tgfb_sep<-gen_matrix(GSM2257302_expr_matrix,GSM2257302_tgfb_gene_list,GSM2257302_state,GSM2257302_pseudo_times)

#esc
GSM2257302_esc_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[1]],GSM2257302_genes)
library(pheatmap)
esc_tgfb_order<-GSM2257302_esc_tgfb.dist$clust_info$order
esc_tgfb_dist_mat<-GSM2257302_esc_tgfb.dist$dist_mat
esc_tgfb_indx<-rownames(esc_tgfb_dist_mat)[esc_tgfb_order]
esc_tgfb<-pheatmap(esc_tgfb_dist_mat[esc_tgfb_order,esc_tgfb_order],color  = rev(hmcol),
                   main="TGF-B in ESC",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#aps
GSM2257302_aps_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[2]],GSM2257302_genes)
library(pheatmap)
aps_tgfb_dist_mat<-GSM2257302_aps_tgfb.dist$dist_mat
aps_tgfb_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(aps_tgfb_dist_mat)]

aps_tgfb<-pheatmap(aps_tgfb_dist_mat[aps_tgfb_genes,aps_tgfb_genes],color  = rev(hmcol),
                   main="TGF-B in APS",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#pxm
GSM2257302_pxm_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[3]],GSM2257302_genes)
library(pheatmap)
pxm_tgfb_dist_mat<-GSM2257302_pxm_tgfb.dist$dist_mat
pxm_tgfb_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(pxm_tgfb_dist_mat)]
pxm_tgfb<-pheatmap(pxm_tgfb_dist_mat[pxm_tgfb_genes,pxm_tgfb_genes],color  = rev(hmcol),
                   main="TGF-B in PXM",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#somitomere
GSM2257302_som_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[4]],GSM2257302_genes)
library(pheatmap)
som_tgfb_dist_mat<-GSM2257302_som_tgfb.dist$dist_mat
som_tgfb_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(som_tgfb_dist_mat)]
som_tgfb<-pheatmap(som_tgfb_dist_mat[som_tgfb_genes,som_tgfb_genes],color  = rev(hmcol),
                   main="TGF-B in Somitomere",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
#Earlysomite
GSM2257302_esom_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[5]],GSM2257302_genes)
library(pheatmap)
esom_tgfb_dist_mat<-GSM2257302_esom_tgfb.dist$dist_mat
esom_tgfb_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(esom_tgfb_dist_mat)]
esom_tgfb<-pheatmap(esom_tgfb_dist_mat[esom_tgfb_genes,esom_tgfb_genes],color  = rev(hmcol),
                    main="TGF-B in Earlysomite",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)

#Dermomyotome
GSM2257302_dom_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[6]],GSM2257302_genes)
library(pheatmap)
dom_tgfb_dist_mat<-GSM2257302_dom_tgfb.dist$dist_mat
dom_tgfb_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(dom_tgfb_dist_mat)]
dom_tgfb<-pheatmap(dom_tgfb_dist_mat[dom_tgfb_genes,dom_tgfb_genes],color  = rev(hmcol),
                    main="TGF-B in Dermomyotome",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)

#Sclerotome
GSM2257302_scl_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[7]],GSM2257302_genes)
library(pheatmap)
scl_tgfb_dist_mat<-GSM2257302_scl_tgfb.dist$dist_mat
scl_tgfb_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(scl_tgfb_dist_mat)]
scl_tgfb<-pheatmap(scl_tgfb_dist_mat[scl_tgfb_genes,scl_tgfb_genes],color  = rev(hmcol),
                   main="TGF-B in Sclerotome",cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)

library(gridExtra)
#grid.arrange(esc_wnt[[4]],aps_wnt[[4]],pxm_wnt[[4]], som_wnt[[4]],esom_wnt[[4]],dom_wnt[[4]],scl_wnt[[4]],
#             esc_tgfb[[4]],aps_tgfb[[4]],pxm_tgfb[[4]],som_tgfb[[4]],esom_tgfb[[4]],dom_tgfb[[4]],scl_tgfb[[4]],nrow=2)
pdf("/Volumes/Home-IRP$/publish/new_pub/final_codes/heatmaps/tgfb_signalingpathway_final.pdf")
grid.arrange(esc_tgfb[[4]],aps_tgfb[[4]],nrow=1)
grid.arrange(pxm_tgfb[[4]],som_tgfb[[4]],nrow=1)
grid.arrange(esom_tgfb[[4]],dom_tgfb[[4]],nrow=1)
grid.arrange(scl_tgfb[[4]],nrow=1)
#grid.arrange(esom_tgfb[[4]],scl_tgfb[[4]],nrow=1)
dev.off()

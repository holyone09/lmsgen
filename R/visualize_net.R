source("./cig_tools.R")
GSM2257302_data<-read.delim("/Volumes/Home-IRP$/development/new_method/GSM2257302_All_samples_sc_tpm.txt",header =T ,sep="\t",stringsAsFactors =F)
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
GSM2257302_state<-c("H7hESC","APS","DLL1PXM","D2_25somitomere","Earlysomite","cDM","Sclerotome")
##wnt
candi_wnt<-as.character(read.delim("wnt_geneset.txt",header=F)[,1])
wnt_gene<-candi_wnt[candi_wnt %in% GSM2257302_genes$geneSymbol]

GSM2257302_wnt_indx<-get_index(GSM2257302_genes,wnt_gene)
GSM2257302_wnt_gene_list<-list(p_gene=GSM2257302_wnt_indx)
GSM2257302_wnt_sep<-gen_matrix(GSM2257302_expr_matrix,GSM2257302_wnt_gene_list,GSM2257302_state,GSM2257302_pseudo_times)

#dist matrix
GSM2257302_esc_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[1]],GSM2257302_genes)
GSM2257302_aps_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[2]],GSM2257302_genes)
GSM2257302_pxm_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[3]],GSM2257302_genes)
GSM2257302_som_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[4]],GSM2257302_genes)
GSM2257302_esom_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[5]],GSM2257302_genes)
GSM2257302_dom_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[6]],GSM2257302_genes)
GSM2257302_scl_wnt.dist<-gen_interacting_mat(GSM2257302_wnt_sep[[7]],GSM2257302_genes)
#gen indices
esc_wnt_indx<-rownames(GSM2257302_esc_wnt.dist$dist_mat)[GSM2257302_esc_wnt.dist$clust_info$order]
aps_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(GSM2257302_aps_wnt.dist$dist_mat)]
pxm_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(GSM2257302_pxm_wnt.dist$dist_mat)]
som_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(GSM2257302_som_wnt.dist$dist_mat)]
esom_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(GSM2257302_esom_wnt.dist$dist_mat)]
dom_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(GSM2257302_dom_wnt.dist$dist_mat)]
scl_genes<-esc_wnt_indx[esc_wnt_indx %in% rownames(GSM2257302_scl_wnt.dist$dist_mat)]
#adjacent matrix
thred=0.31
esc_wnt_dave<-colMeans(GSM2257302_esc_wnt.dist$dist_mat)
esc_wnt_dave<-subset(esc_wnt_dave,esc_wnt_dave<thred)

aps_wnt_dave<-colMeans(GSM2257302_aps_wnt.dist$dist_mat[aps_genes,aps_genes])
aps_wnt_dave<-subset(aps_wnt_dave,aps_wnt_dave<thred)

pxm_wnt_dave<-colMeans(GSM2257302_pxm_wnt.dist$dist_mat[pxm_genes,pxm_genes])
pxm_wnt_dave<-subset(pxm_wnt_dave,pxm_wnt_dave<thred)

som_wnt_dave<-colMeans(GSM2257302_som_wnt.dist$dist_mat[som_genes,som_genes])
som_wnt_dave<-subset(som_wnt_dave,som_wnt_dave<thred)

esom_wnt_dave<-colMeans(GSM2257302_esom_wnt.dist$dist_mat[esom_genes,esom_genes])
esom_wnt_dave<-subset(esom_wnt_dave,esom_wnt_dave<0.35)

dom_wnt_dave<-colMeans(GSM2257302_dom_wnt.dist$dist_mat[dom_genes,dom_genes])
dom_wnt_dave<-subset(dom_wnt_dave,dom_wnt_dave<0.38)

scl_wnt_dave<-colMeans(GSM2257302_scl_wnt.dist$dist_mat[scl_genes,scl_genes])
scl_wnt_dave<-subset(scl_wnt_dave,scl_wnt_dave<0.38)
library(WGCNA)
h7hesc_wnt_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_esc_wnt.dist$dist_mat[names(esc_wnt_dave),names(esc_wnt_dave)], power=6, type='signed')
aps_wnt_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_aps_wnt.dist$dist_mat[names(aps_wnt_dave),names(aps_wnt_dave)], power=6, type='signed')
pxm_wnt_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_pxm_wnt.dist$dist_mat[names(pxm_wnt_dave),names(pxm_wnt_dave)], power=6, type='signed')
som_wnt_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_som_wnt.dist$dist_mat[names(som_wnt_dave),names(som_wnt_dave)], power=6, type='signed')
esom_wnt_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_esom_wnt.dist$dist_mat[names(esom_wnt_dave),names(esom_wnt_dave)], power=6, type='signed')
dom_wnt_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_dom_wnt.dist$dist_mat[names(dom_wnt_dave),names(dom_wnt_dave)], power=6, type='signed')
scl_wnt_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_scl_wnt.dist$dist_mat[names(scl_wnt_dave),names(scl_wnt_dave)], power=6, type='signed')
#network matrix
h7hesc_wnt_net<-draw_net_new(h7hesc_wnt_adj_matrix,0)
aps_wnt_net<-draw_net_new(aps_wnt_adj_matrix,0)
pxm_wnt_net<-draw_net_new(pxm_wnt_adj_matrix,0)
som_wnt_net<-draw_net_new(som_wnt_adj_matrix,0)
esom_wnt_net<-draw_net_new(esom_wnt_adj_matrix,0)
dom_wnt_net<-draw_net_new(dom_wnt_adj_matrix,0)
scl_wnt_net<-draw_net_new(scl_wnt_adj_matrix,0)
#vis. network
par(mfrow = c(3, 3))
par(font =2)

par(mar = c(1, 1, 1, 1), oma = c(4, 4, 0.5, 0.5))
plot(h7hesc_wnt_net[[1]],edge.curved=.2,layout=h7hesc_wnt_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4, main="ESC")
plot(aps_wnt_net[[1]],edge.curved=.2,layout=aps_wnt_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="APS")
plot(pxm_wnt_net[[1]],edge.curved=.2,layout=pxm_wnt_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="PXM")
plot(som_wnt_net[[1]],edge.curved=.2,layout=som_wnt_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Somitomere")
plot(esom_wnt_net[[1]],edge.curved=.2,layout=esom_wnt_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Earlysomite")
plot(dom_wnt_net[[1]],edge.curved=.2,layout=dom_wnt_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Dermomyotome")
plot(scl_wnt_net[[1]],edge.curved=.2,layout=scl_wnt_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Sclerotome")

##tgfb
candi_tgfb<-as.character(read.delim("TGF_beta_geneset.txt",header=F)[,1])
tgfb_gene<-candi_tgfb[candi_tgfb %in% GSM2257302_genes$geneSymbol]

GSM2257302_tgfb_indx<-get_index(GSM2257302_genes,tgfb_gene)
GSM2257302_tgfb_gene_list<-list(p_gene=GSM2257302_tgfb_indx)
GSM2257302_tgfb_sep<-gen_matrix(GSM2257302_expr_matrix,GSM2257302_tgfb_gene_list,GSM2257302_state,GSM2257302_pseudo_times)

#dist matrix
GSM2257302_esc_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[1]],GSM2257302_genes)
GSM2257302_aps_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[2]],GSM2257302_genes)
GSM2257302_pxm_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[3]],GSM2257302_genes)
GSM2257302_som_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[4]],GSM2257302_genes)
GSM2257302_esom_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[5]],GSM2257302_genes)
GSM2257302_dom_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[6]],GSM2257302_genes)
GSM2257302_scl_tgfb.dist<-gen_interacting_mat(GSM2257302_tgfb_sep[[7]],GSM2257302_genes)
#gen indices
esc_tgfb_indx<-rownames(GSM2257302_esc_tgfb.dist$dist_mat)[GSM2257302_esc_tgfb.dist$clust_info$order]
aps_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(GSM2257302_aps_tgfb.dist$dist_mat)]
pxm_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(GSM2257302_pxm_tgfb.dist$dist_mat)]
som_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(GSM2257302_som_tgfb.dist$dist_mat)]
esom_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(GSM2257302_esom_tgfb.dist$dist_mat)]
dom_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(GSM2257302_dom_tgfb.dist$dist_mat)]
scl_genes<-esc_tgfb_indx[esc_tgfb_indx %in% rownames(GSM2257302_scl_tgfb.dist$dist_mat)]
#adjacent matrix
thred=0.34
esc_tgfb_dave<-colMeans(GSM2257302_esc_tgfb.dist$dist_mat)
esc_tgfb_dave<-subset(esc_tgfb_dave,esc_tgfb_dave<thred)

aps_tgfb_dave<-colMeans(GSM2257302_aps_tgfb.dist$dist_mat[aps_genes,aps_genes])
aps_tgfb_dave<-subset(aps_tgfb_dave,aps_tgfb_dave<thred)

pxm_tgfb_dave<-colMeans(GSM2257302_pxm_tgfb.dist$dist_mat[pxm_genes,pxm_genes])
pxm_tgfb_dave<-subset(pxm_tgfb_dave,pxm_tgfb_dave<thred)

som_tgfb_dave<-colMeans(GSM2257302_som_tgfb.dist$dist_mat[som_genes,som_genes])
som_tgfb_dave<-subset(som_tgfb_dave,som_tgfb_dave<thred)

esom_tgfb_dave<-colMeans(GSM2257302_esom_tgfb.dist$dist_mat[esom_genes,esom_genes])
esom_tgfb_dave<-subset(esom_tgfb_dave,esom_tgfb_dave<0.4)

dom_tgfb_dave<-colMeans(GSM2257302_dom_tgfb.dist$dist_mat[dom_genes,dom_genes])
dom_tgfb_dave<-subset(dom_tgfb_dave,dom_tgfb_dave<0.44)

scl_tgfb_dave<-colMeans(GSM2257302_scl_tgfb.dist$dist_mat[scl_genes,scl_genes])
scl_tgfb_dave<-subset(scl_tgfb_dave,scl_tgfb_dave<0.44)
library(WGCNA)
h7hesc_tgfb_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_esc_tgfb.dist$dist_mat[names(esc_tgfb_dave),names(esc_tgfb_dave)], power=6, type='signed')
aps_tgfb_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_aps_tgfb.dist$dist_mat[names(aps_tgfb_dave),names(aps_tgfb_dave)], power=6, type='signed')
pxm_tgfb_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_pxm_tgfb.dist$dist_mat[names(pxm_tgfb_dave),names(pxm_tgfb_dave)], power=6, type='signed')
som_tgfb_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_som_tgfb.dist$dist_mat[names(som_tgfb_dave),names(som_tgfb_dave)], power=6, type='signed')
esom_tgfb_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_esom_tgfb.dist$dist_mat[names(esom_tgfb_dave),names(esom_tgfb_dave)], power=6, type='signed')
dom_tgfb_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_dom_tgfb.dist$dist_mat[names(dom_tgfb_dave),names(dom_tgfb_dave)], power=6, type='signed')
scl_tgfb_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_scl_tgfb.dist$dist_mat[names(scl_tgfb_dave),names(scl_tgfb_dave)], power=6, type='signed')
#network matrix
h7hesc_tgfb_net<-draw_net_new(h7hesc_tgfb_adj_matrix,0)
aps_tgfb_net<-draw_net_new(aps_tgfb_adj_matrix,0)
pxm_tgfb_net<-draw_net_new(pxm_tgfb_adj_matrix,0)
som_tgfb_net<-draw_net_new(som_tgfb_adj_matrix,0)
esom_tgfb_net<-draw_net_new(esom_tgfb_adj_matrix,0)
dom_tgfb_net<-draw_net_new(dom_tgfb_adj_matrix,0)
scl_tgfb_net<-draw_net_new(scl_tgfb_adj_matrix,0)
#vis. network
par(mfrow = c(3, 3))
par(font =2)

par(mar = c(1, 1, 1, 1), oma = c(4, 4, 0.5, 0.5))
plot(h7hesc_tgfb_net[[1]],edge.curved=.2,layout=h7hesc_tgfb_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4, main="ESC")
plot(aps_tgfb_net[[1]],edge.curved=.2,layout=aps_tgfb_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="APS")
plot(pxm_tgfb_net[[1]],edge.curved=.2,layout=pxm_tgfb_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="PXM")
plot(som_tgfb_net[[1]],edge.curved=.2,layout=som_tgfb_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Somitomere")
plot(esom_tgfb_net[[1]],edge.curved=.2,layout=esom_tgfb_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Earlysomite")
plot(dom_tgfb_net[[1]],edge.curved=.2,layout=dom_tgfb_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Dermomyotome")
plot(scl_tgfb_net[[1]],edge.curved=.2,layout=scl_tgfb_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Sclerotome")

##glu
candi_glu<-as.character(read.delim("glu_geneset.txt",header=F)[,1])
glu_gene<-candi_glu[candi_glu %in% GSM2257302_genes$geneSymbol]

GSM2257302_glu_indx<-get_index(GSM2257302_genes,glu_gene)
GSM2257302_glu_gene_list<-list(p_gene=GSM2257302_glu_indx)
GSM2257302_glu_sep<-gen_matrix(GSM2257302_expr_matrix,GSM2257302_glu_gene_list,GSM2257302_state,GSM2257302_pseudo_times)

#dist matrix
GSM2257302_esc_glu.dist<-gen_interacting_mat(GSM2257302_glu_sep[[1]],GSM2257302_genes)
GSM2257302_aps_glu.dist<-gen_interacting_mat(GSM2257302_glu_sep[[2]],GSM2257302_genes)
GSM2257302_pxm_glu.dist<-gen_interacting_mat(GSM2257302_glu_sep[[3]],GSM2257302_genes)
GSM2257302_som_glu.dist<-gen_interacting_mat(GSM2257302_glu_sep[[4]],GSM2257302_genes)
GSM2257302_esom_glu.dist<-gen_interacting_mat(GSM2257302_glu_sep[[5]],GSM2257302_genes)
GSM2257302_dom_glu.dist<-gen_interacting_mat(GSM2257302_glu_sep[[6]],GSM2257302_genes)
GSM2257302_scl_glu.dist<-gen_interacting_mat(GSM2257302_glu_sep[[7]],GSM2257302_genes)
#gen indices
esc_glu_indx<-rownames(GSM2257302_esc_glu.dist$dist_mat)[GSM2257302_esc_glu.dist$clust_info$order]
aps_genes<-esc_glu_indx[esc_glu_indx %in% rownames(GSM2257302_aps_glu.dist$dist_mat)]
pxm_genes<-esc_glu_indx[esc_glu_indx %in% rownames(GSM2257302_pxm_glu.dist$dist_mat)]
som_genes<-esc_glu_indx[esc_glu_indx %in% rownames(GSM2257302_som_glu.dist$dist_mat)]
esom_genes<-esc_glu_indx[esc_glu_indx %in% rownames(GSM2257302_esom_glu.dist$dist_mat)]
dom_genes<-esc_glu_indx[esc_glu_indx %in% rownames(GSM2257302_dom_glu.dist$dist_mat)]
scl_genes<-esc_glu_indx[esc_glu_indx %in% rownames(GSM2257302_scl_glu.dist$dist_mat)]
#adjacent matrix
thred=0.35
esc_glu_dave<-colMeans(GSM2257302_esc_glu.dist$dist_mat)
esc_glu_dave<-subset(esc_glu_dave,esc_glu_dave<thred)

aps_glu_dave<-colMeans(GSM2257302_aps_glu.dist$dist_mat[aps_genes,aps_genes])
aps_glu_dave<-subset(aps_glu_dave,aps_glu_dave<thred)

pxm_glu_dave<-colMeans(GSM2257302_pxm_glu.dist$dist_mat[pxm_genes,pxm_genes])
pxm_glu_dave<-subset(pxm_glu_dave,pxm_glu_dave<0.4)

som_glu_dave<-colMeans(GSM2257302_som_glu.dist$dist_mat[som_genes,som_genes])
som_glu_dave<-subset(som_glu_dave,som_glu_dave<0.38)

esom_glu_dave<-colMeans(GSM2257302_esom_glu.dist$dist_mat[esom_genes,esom_genes])
esom_glu_dave<-subset(esom_glu_dave,esom_glu_dave<0.44)

dom_glu_dave<-colMeans(GSM2257302_dom_glu.dist$dist_mat[dom_genes,dom_genes])
dom_glu_dave<-subset(dom_glu_dave,dom_glu_dave<0.44)

scl_glu_dave<-colMeans(GSM2257302_scl_glu.dist$dist_mat[scl_genes,scl_genes])
scl_glu_dave<-subset(scl_glu_dave,scl_glu_dave<0.50)
library(WGCNA)
h7hesc_glu_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_esc_glu.dist$dist_mat[names(esc_glu_dave),names(esc_glu_dave)], power=6, type='signed')
aps_glu_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_aps_glu.dist$dist_mat[names(aps_glu_dave),names(aps_glu_dave)], power=6, type='signed')
pxm_glu_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_pxm_glu.dist$dist_mat[names(pxm_glu_dave),names(pxm_glu_dave)], power=6, type='signed')
som_glu_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_som_glu.dist$dist_mat[names(som_glu_dave),names(som_glu_dave)], power=6, type='signed')
esom_glu_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_esom_glu.dist$dist_mat[names(esom_glu_dave),names(esom_glu_dave)], power=6, type='signed')
dom_glu_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_dom_glu.dist$dist_mat[names(dom_glu_dave),names(dom_glu_dave)], power=6, type='signed')
scl_glu_adj_matrix <- adjacency.fromSimilarity(1-GSM2257302_scl_glu.dist$dist_mat[names(scl_glu_dave),names(scl_glu_dave)], power=6, type='signed')
#network matrix
h7hesc_glu_net<-draw_net_new(h7hesc_glu_adj_matrix,0)
aps_glu_net<-draw_net_new(aps_glu_adj_matrix,0)
pxm_glu_net<-draw_net_new(pxm_glu_adj_matrix,0)
som_glu_net<-draw_net_new(som_glu_adj_matrix,0)
esom_glu_net<-draw_net_new(esom_glu_adj_matrix,0)
dom_glu_net<-draw_net_new(dom_glu_adj_matrix,0)
scl_glu_net<-draw_net_new(scl_glu_adj_matrix,0)
#vis. network
par(mfrow = c(3, 3))
par(font =2)

par(mar = c(1, 1, 1, 1), oma = c(4, 4, 0.5, 0.5))
plot(h7hesc_glu_net[[1]],edge.curved=.2,layout=h7hesc_glu_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4, main="ESC")
plot(aps_glu_net[[1]],edge.curved=.2,layout=aps_glu_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="APS")
plot(pxm_glu_net[[1]],edge.curved=.2,layout=pxm_glu_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="PXM")
plot(som_glu_net[[1]],edge.curved=.2,layout=som_glu_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Somitomere")
plot(esom_glu_net[[1]],edge.curved=.2,layout=esom_glu_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Earlysomite")
plot(dom_glu_net[[1]],edge.curved=.2,layout=dom_glu_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Dermomyotome")
plot(scl_glu_net[[1]],edge.curved=.2,layout=scl_glu_net[[2]],vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black",edge.arrow.size=.4,main="Sclerotome")

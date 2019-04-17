##pre-processing
source("utility/cig_tools.R")
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
candi_genes<-list()
candi_genes[[1]]<-as.character(read.delim("dataset/wnt_geneset.txt",header=F)[,1])
candi_genes[[2]]<-as.character(read.delim("dataset/TGF_beta_geneset.txt",header=F)[,1])
candi_genes[[3]]<-as.character(read.delim("dataset/emb_geneset.txt",header=F)[,1])

candi_gene=candi_genes[[1]]
ref_gene=GSM2257302_genes
pseudo_times=GSM2257302_pseudo_times
expr_matrix=GSM2257302_expr_matrix

hmcol = colorRampPalette(c("blue","white","red"))(n = 50)
genes<-candi_gene[candi_gene %in% ref_gene$geneSymbol]
gene_indx<-get_index(ref_gene,genes)
gene_list<-list(p_gene=gene_indx)
cell_state<-c("H7hESC","APS","DLL1PXM","D2_25somitomere","Earlysomite","cDM","Sclerotome")

state_sep<-gen_matrix(expr_matrix,gene_list,cell_state,pseudo_times)
gene.dist<-list()
for(i in 1:length(state_sep)){
  gene.dist[[i]]<-gen_interacting_mat(state_sep[[1]],GSM2257302_genes)
}

library(pheatmap)
hmcol = colorRampPalette(c("blue","white","red"))(n = 50)
dist_mat<-gene.dist[[1]]$dist_mat
gene_order<-gene.dist[[1]]$clust_info$order
pheatmap(dist_mat[gene_order,gene_order],color  = rev(hmcol),
         cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)

##function
gen_exp_patterns<-function(candi_gene,ref_gene,pseudo_times,expr_matrix){
  hmcol = colorRampPalette(c("blue","white","red"))(n = 50)
  genes<-candi_gene[candi_gene %in% ref_gene$geneSymbol]
  gene_indx<-get_index(ref_gene,genes)
  gene_list<-list(p_gene=gene_indx)
  cell_state<-c("H7hESC","APS","DLL1PXM","D2_25somitomere","Earlysomite","cDM","Sclerotome")
  
  state_sep<-gen_matrix(expr_matrix,gene_list,cell_state,pseudo_times)
  gene.dist<-list()
  for(i in 1:length(state_sep)){
    gene.dist[[i]]<-gen_interacting_mat(state_sep[[1]],GSM2257302_genes)
  }
  
  fig.list<-list()
  for(i in 1:length(gene.dist)){
    fig.list[[i]]<-draw_matrix(gene.dist[[i]],gene.dist[[i]]$clust_info$order)
  }
  
  for(i in 1:length(fig.list)){
    fig.list[[i]]
  }
}
draw_matrix<-function(g.dist,gene_order){
  library(pheatmap)
  hmcol = colorRampPalette(c("blue","white","red"))(n = 50)
  dist_mat<-g.dist$dist_mat
  pheatmap(dist_mat[gene_order,gene_order],color  = rev(hmcol),
                    cluster_rows=F,cluster_cols=F,border_color = NA,fontsize_col =6,fontsize_row =6)
  
}

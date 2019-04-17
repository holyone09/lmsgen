GSM2257302_data<-read.delim("dataset/GSM2257302_All_samples_sc_tpm.txt",header =T ,sep="\t",stringsAsFactors =F)
GSM2257302_data<-subset(GSM2257302_data,GSM2257302_data[,1]!="")
GSM2257302_data.uni <- unique( GSM2257302_data$geneSymbol )
which(duplicated(GSM2257302_data.uni))
GSM2257302_data.info<-GSM2257302_data[,1:2]
GSM2257302_indx<-data.frame(geneSymbol=vector(),geneID=vector())

for(i in 1:length(GSM2257302_data.uni)){
  tmp<-subset(GSM2257302_data.info,GSM2257302_data$geneSymbol==GSM2257302_data.uni[i])
  GSM2257302_indx<-rbind(GSM2257302_indx,tmp[1,])
}
rownames(GSM2257302_data)<-GSM2257302_data[,2]

GSM2257302_expr_matrix<-as.matrix(GSM2257302_data[,3:ncol(GSM2257302_data)])
#rownames(GSM2257302_expr_matrix)<-GSM2257302_data[,2]
GSM2257302_genes<-GSM2257302_data[GSM2257302_indx[,2],1:2]
#which(duplicated(GSM2257302_indx$geneSymbol))
GSM2257302_expr_matrix<-GSM2257302_expr_matrix[GSM2257302_indx[,2],]
rownames(GSM2257302_expr_matrix)<-GSM2257302_indx[,1]

saveRDS(GSM2257302_expr_matrix,file="./dataset/GSM2257302_exp_matrix.rds")

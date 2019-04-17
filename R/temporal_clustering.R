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

c10<-read.table(file="result/cl0.txt",header=F,stringsAsFactors=F)
c11<-read.table(file="result/cl1.txt",header=F,stringsAsFactors=F)
c12<-read.table(file="result/cl2.txt",header=F,stringsAsFactors=F)
c13<-read.table(file="result/cl3.txt",header=F,stringsAsFactors=F)
c14<-read.table(file="result/cl4.txt",header=F,stringsAsFactors=F)
#c15<-read.table(file="result/cl5.txt",header=F,stringsAsFactors=F)
c16<-read.table(file="result/cl6.txt",header=F,stringsAsFactors=F)
c17<-read.table(file="result/cl7.txt",header=F,stringsAsFactors=F)
c10_sel<-vector()
for(i in 1:nrow(c10)){
  if(strsplit(c10[,1] ,"[.]")[[i]][1]=="APS")
  {
    c10_sel<-cbind(c10_sel,c10[i,1])
  }
    
}
  
###selected status
GSM2257302_expr_matrix_sel<-cbind(
  GSM2257302_expr_matrix[,c10_sel],
  GSM2257302_expr_matrix[,c11[,1]],
  GSM2257302_expr_matrix[,c12[,1]],
  GSM2257302_expr_matrix[,c13[,1]],
  GSM2257302_expr_matrix[,c14[,1]],
  #GSM2257302_expr_matrix[,c15[,1]],
  GSM2257302_expr_matrix[,c16[,1]],
  GSM2257302_expr_matrix[,c17[,1]]
)
library(stringr)
GSM2257302_samples_sel<-colnames(GSM2257302_expr_matrix_sel)
GSM2257302_pseudo_times_sel<-c(rep("APS",len=length(c10_sel)),rep("Somitomere",len=length(c11[,1])),
                               rep("PXM",len=length(c12[,1])),rep("Dermomyotome",len=length(c13[,1])),
                               rep("Sclerotome",len=length(c14[,1])),
                               #rep("LatM",len=length(c15[,1])),
                               rep("ESC",len=length(c16[,1])),rep("Earlysomite",len=length(c17[,1]))
)

proc_GSM2257302 <- TSCAN::preprocess(GSM2257302_expr_matrix_sel)
colnames(proc_GSM2257302) <- GSM2257302_samples_sel
GSM2257302_clust <- exprmclust_imsgen(proc_GSM2257302, clusternum = 7,reduce=T)

##
GSM2257302_orderTSCAN <- TSCAN::TSCANorder(GSM2257302_clust, orderonly = FALSE,flip = T)
GSM2257302_pseudotime_order_tscan <- as.character(GSM2257302_orderTSCAN$sample_name)
GSM2257302_orderTSCAN_cluster<-GSM2257302_orderTSCAN[GSM2257302_samples_sel,]
GSM2257302_cell_pseudotime<-as.data.frame(cbind(GSM2257302_orderTSCAN_cluster$sample_name,GSM2257302_pseudo_times_sel))
colnames(GSM2257302_cell_pseudotime)<-c("sample_name","cell_type")
GSM2257302_cell_pseudotime$pseudotime_order_tscan <- NA
GSM2257302_cell_pseudotime$pseudotime_order_tscan <- GSM2257302_orderTSCAN_cluster$Pseudotime

library(ggbeeswarm)
library(ggplot2)
library(ggthemes)
ggplot(GSM2257302_cell_pseudotime, 
       aes(x = pseudotime_order_tscan, 
           y = cell_type, colour = cell_type)) + 
  scale_y_discrete(limits=c("ESC","APS","PXM","Somitomere","Earlysomite","Dermomyotome","Sclerotome"))+
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Pseudotime") + ylab("Timepoint") 
#+ggtitle("Cells ordered by TSCAN pseudotime")


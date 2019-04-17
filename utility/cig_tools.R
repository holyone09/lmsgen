exprmclust_imsgen <- function (data, clusternum = 2:9, modelNames = "VVV", reduce = T, cluster = NULL,pca_num=1:20) {
  library(mclust)
  library(igraph)
  set.seed(12345)
  if (reduce) {
    sdev <- prcomp(t(data), scale = T)$sdev[pca_num]
    x <- pca_num
    optpoint <- which.min(sapply(2:10, function(i) {
      x2 <- pmax(0, x - i)
      sum(lm(sdev ~ x + x2)$residuals^2)
    }))
    pcadim = optpoint + 1
    tmpdata <- t(apply(data, 1, scale))
    colnames(tmpdata) <- colnames(data)
    tmppc <- prcomp(t(tmpdata), scale = T)
    pcareduceres <- t(tmpdata) %*% tmppc$rotation[, 1:pcadim]
  }
  else {
    pcareduceres <- t(data)
  }
  if (is.null(cluster)) {   
    clusternum <- clusternum[clusternum > 1]
    res <- suppressWarnings(Mclust(pcareduceres, G = clusternum, modelNames = modelNames))
    clusterid <- apply(res$z, 1, which.max)
    clunum <- res$G
  } else {
    clunum <- length(unique(cluster))
    clusterid <- cluster
  }
  clucenter <- matrix(0, ncol = ncol(pcareduceres), nrow = clunum)
  for (cid in 1:clunum) {
    clucenter[cid, ] <- colMeans(pcareduceres[names(clusterid[clusterid == cid]), , drop = F])
  }
  dp <- as.matrix(dist(clucenter))
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  list(pcareduceres = pcareduceres, MSTtree = dp_mst, clusterid = clusterid, clucenter = clucenter)
}
draw_net_new<-function(adjMatrix,th=0){
  source("./GENIE3_R/genie3.R")
  library(igraph)
  library(RCy3)
  link_tbx6.list <- get.link.list(adjMatrix, report.max=250)
  edge_listsi_tbx6 <- link_tbx6.list[!duplicated(link_tbx6.list),]
  Gsi_tbx6 <- graph.data.frame(edge_listsi_tbx6,directed = T)
  
  Asi_tbx6 <- get.adjacency(Gsi_tbx6,sparse = F,attr = "im",type = "both")
  
  g_arasi_tbx6 <- graph.adjacency(Asi_tbx6,mode = "directed",weighted = T)
  dll_gD_tbx6 <- igraph::simplify(g_arasi_tbx6)
  degAll_tbx6 <- igraph::degree(dll_gD_tbx6, v = igraph::V(dll_gD_tbx6), mode = "all")
  betAll_tbx6 <- igraph::betweenness(dll_gD_tbx6, v = igraph::V(dll_gD_tbx6), directed = T) / (((igraph::vcount(dll_gD_tbx6) - 1) * (igraph::vcount(dll_gD_tbx6)-2)) / 2)
  betAll_tbx6.norm <- (betAll_tbx6 - min(betAll_tbx6))/(max(betAll_tbx6) - min(betAll_tbx6))
  rm(betAll_tbx6)
  dsAll_tbx6 <- igraph::similarity.dice(dll_gD_tbx6, vids = igraph::V(dll_gD_tbx6), mode = "all")
  
  dll_gD_tbx6 <- igraph::set.vertex.attribute(dll_gD_tbx6, "degree", index = igraph::V(dll_gD_tbx6), value = degAll_tbx6)
  dll_gD_tbx6 <- igraph::set.vertex.attribute(dll_gD_tbx6, "betweenness", index = igraph::V(dll_gD_tbx6), value = betAll_tbx6.norm)
  
  F1 <- function(x) {data.frame(V4 = dsAll_tbx6[which(igraph::V(dll_gD_tbx6)$name == as.character(x$from.gene)), which(igraph::V(dll_gD_tbx6)$name == as.character(x$to.gene))])}
  dataSet_tbx6.ext <- plyr::ddply(edge_listsi_tbx6, .variables=c("from.gene", "to.gene", "im"), function(x) data.frame(F1(x)))
  
  
  dll_gD_tbx6 <- igraph::set.edge.attribute(dll_gD_tbx6, "weight", index = igraph::E(dll_gD_tbx6), value = 0)
  dll_gD_tbx6 <- igraph::set.edge.attribute(dll_gD_tbx6, "similarity", index = igraph::E(dll_gD_tbx6), value = 0)
  
  for (i in 1:nrow(dataSet_tbx6.ext))
  {
    igraph::E(dll_gD_tbx6)[as.character(dataSet_tbx6.ext$from.gene) %--% as.character(dataSet_tbx6.ext$to.gene)]$weight <- as.numeric(dataSet_tbx6.ext$im)
    igraph::E(dll_gD_tbx6)[as.character(dataSet_tbx6.ext$from.gene) %--% as.character(dataSet_tbx6.ext$to.gene)]$similarity <- as.numeric(dataSet_tbx6.ext$V4)
  }
  #E(dll_gD)$width <- 1+E(dll_gD)$weight/12
  #V(dll_gD_tbx6)$size <- degAll_tbx6*3
  #cat(max(E(dll_gD_tbx6)$weight))
  if(th<=0.0){
    th=mean(E(dll_gD_tbx6)$weight)
  }
  tryCatch(E(dll_gD_tbx6)[ weight > th ]$color <- "red",error=function(e) {print("error")})
  #E(dll_gD_tbx6)[ weight <= mean(E(dll_gD_tbx6)$weight) ]$color <- "grey"
  tryCatch(E(dll_gD_tbx6)[ (weight <= th)&E(dll_gD_tbx6)$weight!=0.0 ]$color <- "grey",error=function(e) {print("error")})
  tryCatch(E(dll_gD_tbx6)[ (weight <= th)&E(dll_gD_tbx6)$weight!=0.0 ]$shape<- "dashed",error=function(e) {print("error")})
  #E(dll_gD_tbx6)$width <-0.5
  #E(dll_gD_tbx6)$width <-1+E(dll_gD_tbx6)$weight/12
  #l_es_tbx6<- layout.circle(dll_gD_tbx6)
  l_es_tbx6<- layout.fruchterman.reingold(dll_gD_tbx6)
  #plot(dll_gD_tbx6,edge.curved=.2,layout=l_es_tbx6,vertex.label.cex=0.7,vertex.label.font=2,vertex.size=20,vertex.label.color="black")
  rn_list<-list(dll_gD_tbx6,l_es_tbx6)
  #l_tbx6<- layout_with_fr(dll_gD_tbx6)
  #plot(dll_gD_tbx6,edge.curved=.2,layout=l_tbx6)
  rn_list
}

intep_mat<-function(indx,dist_mat){
  idx_diff<-setdiff(indx,rownames(dist_mat))
  dist_mat_f<-dist_mat
  for(i in 1:length(idx_diff)){
    dist_mat_f<-rbind(dist_mat_f,rep(1,times=length(colnames(dist_mat_f))))
    dist_mat_f<-cbind(dist_mat_f,rep(1,times=length(rownames(dist_mat_f))))
  }
  diag(dist_mat_f)<-0
  rownames(dist_mat_f)<-c(rownames(dist_mat),idx_diff)
  colnames(dist_mat_f)<-c(colnames(dist_mat),idx_diff)
  
  dist_mat_f
}

gen_interacting_mat_dist<-function(exp_mat,g_list){
  
  mat_act.dist = log10(exp_mat + 1)
  rownames(mat_act.dist)<-g_list
  mat_act.dist = as.matrix(dist(mat_act.dist))
  mat_act.dist = mat_act.dist/max(mat_act.dist)
  
  mat_act.dist
  
}
#wnt,tgfb,glu,plu=0.75,gollagen=0.7,skin=0.1,neg reg cell motil=0.09,appanse dev=0.4,embryonic morphogenesis=0.7
gen_interacting_mat<-function(exp_mat,gene_ref,thr=0.75){
  pseudoCount_act = log10(exp_mat + 1)
  mat_act.dist = pseudoCount_act[rowMeans(pseudoCount_act)>thr,]
  rownames(mat_act.dist)<-gene_ref[rownames(mat_act.dist),]$geneSymbol
  mat_act.dist = as.matrix(dist(mat_act.dist))
  mat_act.dist = mat_act.dist/max(mat_act.dist)
  mat_act.dist.hclust<- hclust(dist(mat_act.dist))
  results<-list(dist_mat=mat_act.dist,clust_info=mat_act.dist.hclust)
  results
}

get_index<-function(genesDB,genes){
  tmp<-vector()
  gene_indx<-vector()
  for(i in 1:length(genes)){
    tmp<-as.character(subset(genesDB,genesDB$geneSymbol==genes[i])$geneID)
    gene_indx<-c(gene_indx,tmp)
  }
  gene_indx
}

gen_matrix<-function(exp_mat,g_list,exp_state,time_indx){
  
  results<-list()
  if(length(g_list)==1){
    for(j in 1:length(exp_state)){
      results[[j]]<-exp_mat[c(g_list[[1]]),time_indx==exp_state[j]]
    }
  }
  else{
    
  for(i in 1:(length(g_list)-1)){
    
    tmp<-list()
   # cat(i,c(g_list[[1]],g_list[[i+1]]))
    for(j in 1:length(exp_state)){
      tmp[[j]]<-exp_mat[c(g_list[[1]],g_list[[i+1]]),time_indx==exp_state[j]]
    }
    
    results[[i]]<-tmp
  }
  }
  results
}
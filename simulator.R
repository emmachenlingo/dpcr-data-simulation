## simulation pipeline functions
library(MASS)
library(sn)
library(ggplot2)
library(flowCore)
library(FlowSOM)
library(flowPeaks)
library(flowClust)
library(ddPCRclust)
library(SamSPECTRAL)
library(flowMerge)
# library(fossil)
library(aricode)
library(e1071)
library(caret)
library(DepthProc)
library(dplyr)
library(patchwork)
library(ggpubr)
library(dbscan)
library(mgcv)
library(spatstat)
library(sf)
library(sp)
library(maptools)
library(raster)
library(rstan)
library(tidyverse)
library(cowplot)
library(splines)
library(cluster)
library(clue)
library(nnet)


## plotting function
plotting_func<-function(df_data,plot_name,sim){
  png(file=paste0(plot_name,sim,'.png'),width = 500,height = 300)
  p<-ggplot(data=df_data, aes(channel1, channel2, colour = factor(group)))+ 
    geom_point(size=4,show.legend = FALSE) +labs(x = "Green Channel",y='Red Channel')+theme(text = element_text(size = 50),
                                                                                            panel.grid.major = element_blank(),
                                                                                              panel.grid.minor = element_blank(),
                                                                                              panel.background = element_blank(),
                                                                                              axis.line = element_line(colour = "black"),
                                                                                              plot.margin = margin(t = 20,  # Top margin
                                                                                                                   r = 20,  # Right margin
                                                                                                                   b = 20,  # Bottom margin
                                                                                                                   l = 20))
  print(p)
  dev.off()
}


eval_func<-function(class_group,true_group){
  adj_randindex<-ARI(class_group,true_group)
  
  return(adj_randindex)
}


## gridding function
gridding_density<-function(data,s1=15,s2=15){
  range_G<-seq(min(data[,1])-abs(min(data[,1])/10),max(data[,1])+abs(min(data[,1])/10),length.out = s1)
  range_R<-seq(min(data[,2])-abs(min(data[,2])/10),max(data[,2])+abs(min(data[,2])/10),length.out = s2)
  
  intensity_mat<-matrix(0,ncol=s1,nrow=s2)
  for (i in 1:(s1-1)){
    for (j in 1:(s2-1)){
      intensity_mat[j,i]<-sum((data[,1]>=range_G[i]& data[,1]<range_G[i+1]) & (data[,2]>=range_R[j]& data[,2]<range_R[j+1])) 
    }
  }
  return (list(intensity_mat,range_G,range_R))
}


## search for each central points in each cluster
# central_points_func<-function(data,method='neighbourhood',s1=15,s2=15,perc=0.1,plotname=NULL,plot=FALSE){
#   unique_group<-unique(data$group)
#   data_index<-1:nrow(data)
#   outliers<-NULL
#   outliers_index<-NULL
#   non_outliers<-NULL
#   non_outliers_index<-NULL
#   
#   if (method=='neighbourhood'){
#     for (i in unique_group) {
#       data_group<-data[data$group==i,]
#       data_group_n<-nrow(data_group)
#       data_group_index<-data_index[data$group==i]
#       intensity_cal<-gridding_density(data_group,s1=s1,s2=s1)
#       range_G<-intensity_cal[[2]]
#       range_R<-intensity_cal[[3]]
#       unit_G<-range_G[2]-range_G[1]
#       unit_R<-range_R[2]-range_R[1]
#       neighbour_count<-NULL
#       
#       for (j in 1:data_group_n){
#         neighbour_count<-c(neighbour_count,sum((data_group[,1]>=(data_group[j,1]-unit_R))&(data_group[,1]<(data_group[j,1]+unit_R))&(data_group[,2]>=(data_group[j,2]-unit_G))&(data_group[,2]<(data_group[j,2]+unit_G))))
#       }
#       
#       data_group_reordered<-data_group[order(neighbour_count),]
#       data_group_index_reordered<-data_group_index[order(neighbour_count)]
#       
#       outliers<-rbind(outliers,data_group_reordered[1:round(data_group_n*perc),])
#       outliers_index<-c(outliers_index,data_group_index_reordered[1:round(data_group_n*perc)])
#       non_outliers<-rbind(non_outliers,data_group_reordered[(round(data_group_n*perc)+1):data_group_n,])
#       non_outliers_index<-c(non_outliers_index,data_group_index_reordered[(round(data_group_n*perc)+1):data_group_n])
#       
#       if(plot){
#         png(file=paste0(plotname,method,'.png'),width = 500,height = 300)
#         plot(outliers,col='red',xlim=c(min(data[,1])-10,max(data[,1])+10),ylim=c(min(data[,2])-10,max(data[,2])+10),xlab='Channel2', ylab='Channel3')
#         points(non_outliers)
#         dev.off()
#       }
#     }
#     
#   }
#   
#   return(list(outliers,outliers_index,non_outliers,non_outliers_index))
# }

central_points_func<-function(data,perc=0.1,plotname=NULL,plot=FALSE){
  unique_group<-unique(data$group)
  data_index<-1:nrow(data)
  outliers<-NULL
  outliers_index<-NULL
  non_outliers<-NULL
  non_outliers_index<-NULL
  
  data_sil<-silhouette_coef(data,data$group)[[1]]
  
  for (i in unique_group) {
    data_group<-data[data$group==i,]
    data_group_n<-nrow(data_group)
    data_group_index<-data_index[data$group==i]
    
    group_sil<-data_sil[data_sil$cluster==i,'width']
    data_group_reordered<-data_group[order(group_sil),]
    data_group_index_reordered<-data_group_index[order(group_sil)]
      
    outliers<-rbind(outliers,data_group_reordered[1:round(data_group_n*perc),])
    outliers_index<-c(outliers_index,data_group_index_reordered[1:round(data_group_n*perc)])
    non_outliers<-rbind(non_outliers,data_group_reordered[(round(data_group_n*perc)+1):data_group_n,])
    non_outliers_index<-c(non_outliers_index,data_group_index_reordered[(round(data_group_n*perc)+1):data_group_n])
      
    if(plot){
      png(file=paste0(plotname,method,'.png'),width = 500,height = 300)
      plot(outliers,col='red',xlim=c(min(data[,1])-10,max(data[,1])+10),ylim=c(min(data[,2])-10,max(data[,2])+10),xlab='Channel2', ylab='Channel3')
      points(non_outliers)
      dev.off()
    }
  }
    
  return(list(outliers,outliers_index,non_outliers,non_outliers_index))
}


euclidean.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

misclass_index<-function(data,class_group,true_centers){
  ## wrongly classified data points
  unique_data_group<-unique(class_group)
  unique_data_group<-unique_data_group[unique_data_group!=0]
  dist_mat<-NULL
  for (i in unique_data_group){
    if (!is.null(nrow(data[class_group==i,]))){
      group_mean<-apply(data[class_group==i,],2,mean)
    } else {group_mean<-data[class_group==i,]}
    dist<-NULL
    for (j in 1:nrow(true_centers)){
      dist<-c(dist,euclidean.dist(group_mean,true_centers[j,]))
    }
    dist_mat<-rbind(dist_mat,dist)
    
  }
  
  data_group_new<-class_group+100
  
  if (length(unique_data_group)<=nrow(true_centers)){
    assign_res<-solve_LSAP(dist_mat)
    
    for (k in 1:length(unique_data_group)){
       data_group_new[class_group==unique_data_group[k]]<-assign_res[k]
    }
    
  } else {
    assign_res<-solve_LSAP(t(dist_mat))

    for (k in 1:length(assign_res)){
      data_group_new[class_group==unique_data_group[assign_res[k]]]<-k
    }
  }
  
  return(data_group_new)
}


# silhouette_coef<-function(data,clustering,dist_method='euclidean',plot=FALSE,plot_name='orig',sim='orig'){
silhouette_coef<-function(data,clustering,plot=FALSE,plot_name='orig',sim='orig'){
  
  # si <- silhouette(clustering, dist(data, "euclidean"))
  ndim<-ncol(data)
  si <- approxSilhouette(data[,-ndim], clustering)
  
  sil_tab<-as.data.frame(si) %>% group_by(cluster) %>% summarise(mean_sil=round(mean(width),2))
  
  # find the mean of each cluster
  
  clust_pos<- data %>% group_by(group) %>% summarise(mean_x=mean(channel1),mean_y=mean(channel2))
  colnames(clust_pos)[1]<-'cluster'
  df_merged <- merge(x=sil_tab,y=clust_pos, by="cluster", all.x=TRUE)
  
  if(plot){
    png(file=paste0(plot_name,sim,'.png'),width = 500,height = 300)
    p<-ggplot(data=data, aes(channel1, channel2, colour = factor(group)))+ 
      geom_point(size=0.7,show.legend = FALSE) +labs(x = "Green Channel",y='Red Channel')+theme(panel.grid.major = element_blank(),
                                                                                                panel.grid.minor = element_blank(),
                                                                                                panel.background = element_blank(),
                                                                                                axis.line = element_line(colour = "black"),
                                                                                                plot.margin = margin(t = 20,  # Top margin
                                                                                                                     r = 20,  # Right margin
                                                                                                                     b = 20,  # Bottom margin
                                                                                                                     l = 20)) + annotate("text", x=df_merged$mean_x, y=df_merged$mean_y, label= df_merged$mean_sil)
    
    print(p)
    dev.off()
  }
  
  return(list(si,df_merged))
  
}

mahal_dist<-function(data){
  index_sel<-NULL
  index_all<-1:nrow(data)
  group_num<-length(unique(data$group))
  clust_pos<- data %>% group_by(group) %>% summarise(mean_x=mean(channel1),mean_y=mean(channel2))
  # Cutoff value for distances from Chi-Sqaure Dist. 
  cutoff <- qchisq(p = 0.99, df = ncol(data))
  
  for (i in 1:group_num){
    index_group<-index_all[data$group==i]
    data_group<-data[data$group==i,c(1,2)]
    data_group_cov<-cov(data_group)
    distances <- mahalanobis(x = data_group, center = as.numeric(clust_pos[clust_pos$group==i,c(2,3)]), cov = as.matrix(data_group_cov))
    
    ## Display observation whose distance greater than cutoff value
    index_sel<-c(index_sel,index_group[distances <= cutoff])
  } 
  return(index_sel)
}

skew_t_dist_fit<-function(data){
  X <- rep(1,nrow(data))
  fit <-  sn::mst.mple(x=X,y=cbind(data$x,data$y), opt.method="nlminb")
  return(fit)
}



## skew-t distribution as offset in the poisson point process
skew_t_dist_pp_sample<-function(data,fit,dense=1){
  data_log_dense<-function(x,y) {
    dense<-sn::dmst(cbind(x,y), xi=as.numeric(fit$dp$beta),Omega=fit$dp$Omega,alpha = fit$dp$alpha,nu=fit$dp$nu,log=TRUE) 
    return(dense)}
  
  wind0<-c(min(data[,1])-abs(min(data[,1])/10),max(data[,1])+abs(min(data[,1])/10))
  wind1<-c(min(data[,2])-abs(min(data[,2])/10),max(data[,2])+abs(min(data[,2])/10))
  
  pp <- as.ppp(data,c(wind0,wind1))
  x<-data$x
  y<-data$y
  pp_model<-ppm(pp~x+y+offset(data_log_dense))
  
  lamb<-function(x,y,beta0,beta1,beta2) {dense*exp(data_log_dense(x,y)+beta0+beta1*x+beta2*y)}
  pp_sim <- rpoispp(lamb, beta0=coef(pp_model)[1],beta1=coef(pp_model)[2],beta2=coef(pp_model)[3],win=owin(wind0,wind1))
  
  pp_data<-data.frame(cbind(x=pp_sim$x,y=pp_sim$y))
  return(list(pp_data,pp_model))
}

rains_pp_sample<-function(data,fit,fit2,range=NULL,dense=1){
  
  data_log_dense<-function(x,y) {
    dense<-sn::dmst(cbind(x,y), xi=as.numeric(fit$dp$beta),Omega=fit$dp$Omega,alpha = fit$dp$alpha,nu=fit$dp$nu,log=TRUE) 
    return(dense)}
  
  if (is.null(range)){
    wind0<-c(min(data[,1])-abs(min(data[,1])/10),max(data[,1])+abs(min(data[,1])/10))
    wind1<-c(min(data[,2])-abs(min(data[,2])/10),max(data[,2])+abs(min(data[,2])/10))
  } else {
    wind0=range[1,]
    wind1=range[2,]
  }
  
  lamb<-function(x,y,beta0,beta1,beta2) {dense*exp(data_log_dense(x,y)+beta0+beta1*x+beta2*y)}
  pp_sim <- rpoispp(lamb, beta0=coef(fit2)[1],beta1=coef(fit2)[2],beta2=coef(fit2)[3],win=owin(wind0,wind1))
  
  pp_data<-data.frame(cbind(x=pp_sim$x,y=pp_sim$y))
  return(pp_data)
}

skew_t_dist_pp_sample_constraints<-function(data,fit,pars,dense=1){
  data_log_dense<-function(x,y) {
    dense<-sn::dmst(cbind(x,y), xi=as.numeric(fit$dp$beta),Omega=fit$dp$Omega,alpha = fit$dp$alpha,nu=fit$dp$nu,log=TRUE) 
    return(dense)}
  
  constraint_func1<-function(x,y,par=pars[1,]) {
    return((y-par[2]-par[1]*x)^2)
  }
  
  constraint_func2<-function(x,y,par=pars[2,]) {
    return((y-par[2]-par[1]*x)^2)
  }
  
  constraint_func3<-function(x,y,par=pars[3,]) {
    return((y-par[2]-par[1]*x)^2)
  }
  
  wind0<-c(min(data[,1])-abs(min(data[,1])/10),max(data[,1])+abs(min(data[,1])/10))
  wind1<-c(min(data[,2])-abs(min(data[,2])/10),max(data[,2])+abs(min(data[,2])/10))
  
  pp <- as.ppp(data,c(wind0,wind1))
  x<-data$x
  y<-data$y
  pp_model<-ppm(pp~x+y+offset(data_log_dense)+constraint_func1+constraint_func2+constraint_func3)
  
  lamb<-function(x,y,beta0,beta1,beta2,beta3,beta4,beta5) {dense*exp(data_log_dense(x,y)+beta0+beta1*x+beta2*y+beta3*constraint_func1(x,y)+beta4*constraint_func2(x,y)+beta5*constraint_func3(x,y))}
  pp_sim <- rpoispp(lamb, beta0=coef(pp_model)[1],beta1=coef(pp_model)[2],beta2=coef(pp_model)[3],beta3=coef(pp_model)[4],beta4=coef(pp_model)[5],beta5=coef(pp_model)[6],win=owin(wind0,wind1))
  
  pp_data<-data.frame(cbind(x=pp_sim$x,y=pp_sim$y))
  return(list(pp_data,pp_model))
}




simulation_plotcheck<-function(data,plot_name,sim){
  #png(file=paste0(plot_name,sim,'.png'),width = 1000,height = 600)
  p1 <- ggplot(data[[1]], aes(x,y)) +
    geom_point() + 
    #geom_density_2d() +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
    theme_classic() 
  
  p2 <- ggplot(data[[2]], aes(x,y)) +
    geom_point() + 
    #geom_density_2d() +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
    theme_classic() 
  
  p3 <- ggplot(data[[3]], aes(x,y)) +
    geom_point() + 
    #geom_density_2d() +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
    theme_classic() 
  
  p4 <- ggplot(data[[4]], aes(x,y)) +
    geom_point() + 
    #geom_density_2d() +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
    theme_classic()
  
  prow <- plot_grid (p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",plot.margin = margin(l=0,r = 0)),
                     p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     align = 'vh',
                     labels = c("A","B","C","D"),
                     hjust = 0,
                     nrow = 2
  )
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, ncol = 1, rel_heights = c(1.5, 0.2))
  #print(p)
  #dev.off()
  ggsave(paste0(plot_name,sim,".png"),plot = p,width = 8,
         height = 8, dpi = 350)
}


simulation_plotcheck_depth<-function(data1,data2,plot_name,sim,method='Projection'){
  cluster1_mahadist<-mahalanobis(data1[[1]], colMeans(data1[[1]]), cov(data1[[1]]))
  cluster1_depth1<-DepthProc::depth(data1[[1]],data1[[1]],method = method)
  cluster1_depth2<-DepthProc::depth(data1[[1]],data2[[1]],method = method)
  cluster1_depth<-data.frame(OrigDepth=cluster1_depth1,SimDepth=cluster1_depth2,MahaDist=cluster1_mahadist)
  
  cluster2_mahadist<-mahalanobis(data1[[2]], colMeans(data1[[2]]), cov(data1[[2]]))
  cluster2_depth1<-DepthProc::depth(data1[[2]],data1[[2]],method = method)
  cluster2_depth2<-DepthProc::depth(data1[[2]],data2[[2]],method = method)
  cluster2_depth<-data.frame(OrigDepth=cluster2_depth1,SimDepth=cluster2_depth2,MahaDist=cluster2_mahadist)
  
  cluster3_mahadist<-mahalanobis(data1[[3]], colMeans(data1[[3]]), cov(data1[[3]]))
  cluster3_depth1<-DepthProc::depth(data1[[3]],data1[[3]],method = method)
  cluster3_depth2<-DepthProc::depth(data1[[3]],data2[[3]],method = method)
  cluster3_depth<-data.frame(OrigDepth=cluster3_depth1,SimDepth=cluster3_depth2,MahaDist=cluster3_mahadist)
  
  cluster4_mahadist<-mahalanobis(data1[[4]], colMeans(data1[[4]]), cov(data1[[4]]))
  cluster4_depth1<-DepthProc::depth(data1[[4]],data1[[4]],method = method)
  cluster4_depth2<-DepthProc::depth(data1[[4]],data2[[4]],method = method)
  cluster4_depth<-data.frame(OrigDepth=cluster4_depth1,SimDepth=cluster4_depth2,MahaDist=cluster4_mahadist)
  
  p1 <- ggplot(cluster1_depth, aes(x = OrigDepth, y = SimDepth, colour = MahaDist)) +
    geom_point()+ geom_abline(intercept = 0, slope = 1)+
    theme_classic() 
  
  p2 <- ggplot(cluster2_depth, aes(x = OrigDepth, y = SimDepth, colour = MahaDist)) +
    geom_point()+ geom_abline(intercept = 0, slope = 1)+
    theme_classic() 
  
  p3 <- ggplot(cluster3_depth, aes(x = OrigDepth, y = SimDepth, colour = MahaDist)) +
    geom_point()+ geom_abline(intercept = 0, slope = 1)+
    theme_classic() 
  
  p4 <- ggplot(cluster4_depth, aes(x = OrigDepth, y = SimDepth, colour = MahaDist)) +
    geom_point()+ geom_abline(intercept = 0, slope = 1)+
    theme_classic() 
  
  prow <- plot_grid (p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.margin = margin(l=0,r = 0)),
                     p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     align = 'vh',
                     labels = c("A","B","C","D"),
                     hjust = 0,
                     nrow = 2
  )
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, ncol = 1, rel_heights = c(1.5, 0.2))
  
  ggsave(paste0(plot_name,sim,".png"),plot = p,width = 10,
         height = 8, dpi = 350)
  
}


approxSilhouette <- function(x, clusters) {
  x <- as.matrix(x)
  uclust <- sort(unique(clusters))
  averaged <- list(length(uclust))
  clust.var <- numeric(length(uclust))
  
  for (i in seq_along(uclust)) {
    current <- uclust[i]==clusters
    xcurrent <- x[current,,drop=FALSE]
    centroid <- colMeans(xcurrent)
    averaged[[i]] <- centroid
    clust.var[i] <- sum(colMeans(sweep(xcurrent, 2, centroid)^2))
  }
  
  self.dist <- other.dist <- rep(Inf, nrow(x))
  other.clust <- integer(nrow(x))
  tx <- t(x)
  
  for (i in seq_along(uclust)) {
    D <- sqrt(colSums((tx - averaged[[i]])^2) + clust.var[i])
    
    is.self <- uclust[i]==clusters
    self.dist[is.self] <- D[is.self]
    
    is.other <- !is.self
    other.D <- D[is.other]
    better <- other.D < other.dist[is.other]
    other.dist[is.other][better] <- other.D[better]
    other.clust[is.other][better] <- i
  }
  
  result<-data.frame(
    cluster=clusters,
    other=uclust[other.clust],
    width=(other.dist - self.dist)/pmax(other.dist, self.dist),
    row.names=rownames(x)
  )
  
  return(result)
}


## silhouette coefficient plots
silhouette_plot<-function(sil_data,plot_name,sim){
  cluster_sil<-list()
  for (i in 1:4){
    cluster_sil_tmp<-NULL
    for (j in 1:100){
      cluster_sil_tmp<-c(cluster_sil_tmp,sil_data[[j]][i,2])
    }
    cluster_sil[[i]]<-data.frame(silhouette_coef=cluster_sil_tmp)
  }
  
  
  png(file=paste0(plot_name,sim,'.png'),width = 1000,height = 600)
  p1 <- ggplot(cluster_sil[[4]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  
  p2 <- ggplot(cluster_sil[[3]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  p3 <- ggplot(cluster_sil[[1]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  p4 <- ggplot(cluster_sil[[2]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  prow <- plot_grid (p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",plot.margin = margin(l=0,r = 0)),
                     p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     align = 'vh',
                     labels = c("A","B","C","D"),
                     hjust = 0,
                     nrow = 2
  )
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, ncol = 1, rel_heights = c(1.5, 0.2))
  print(p)
  dev.off()
}

silhouette_single_plot<-function(sil_data,plot_name,sim){
  cluster_sil<-list()
  for (i in 1:4){
    cluster_sil[[i]]<-data.frame(silhouette_coef=sil_data[sil_data[,1]==i,3])
  }
  
  png(file=paste0(plot_name,sim,'.png'),width = 1000,height = 600)
  p1 <- ggplot(cluster_sil[[4]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  
  p2 <- ggplot(cluster_sil[[3]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  p3 <- ggplot(cluster_sil[[1]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  p4 <- ggplot(cluster_sil[[2]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  prow <- plot_grid (p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",plot.margin = margin(l=0,r = 0)),
                     p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     align = 'vh',
                     labels = c("A","B","C","D"),
                     hjust = 0,
                     nrow = 2
  )
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, ncol = 1, rel_heights = c(1.5, 0.2))
  print(p)
  dev.off()
}

eval_plot<-function(result_list,plot_name,sim){
  ### argument: 
  ### result_list should be a list of evaluation results, such as rand index
  
  methods_name<-c('kmeans','kmeans_initials','cmeans','cmeans_initials','dbscan', 'fsom','flowpeaks','flowclust',
                  'flowclust_auto','flowmerge_moreclust','samspectral','samspectral_auto')
  
  result_df<-NULL
  for (i in 1:length(methods_name)){
    result_df<-rbind(result_df,cbind(result_list[[methods_name[i]]],methods_name[i]))
  }
  
  result_df<-data.frame(result_df)
  colnames(result_df)<-c('ri','adj_ri','ami','methods')
  result_df$adj_ri<-as.numeric(result_df$adj_ri)
  
  png(file=paste0(plot_name,sim,'.png'),width = 1000,height = 400)
  
  p<-ggplot(result_df, aes(x=methods, y=adj_ri)) +
    geom_boxplot()+theme_classic()+labs(x = "Methods",y='Adjusted Rand Index')
  
  print(p)
  dev.off()
}


sens_plot<-function(result_list,true_lambda,plot_name,sim){
  ### argument: 
  ### result_list should be a list of 
  
  methods_name<-c('kmeans','kmeans_initials','cmeans','cmeans_initials','dbscan', 'fsom','flowpeaks','flowclust',
                  'flowclust_auto','flowmerge_moreclust','samspectral','samspectral_auto')
  
  result_df<-NULL
  for (i in 1:length(methods_name)){
    result_df<-rbind(result_df,cbind((result_list[[methods_name[i]]]-true_lambda)/true_lambda,methods_name[i]))
  }
  
  result_df<-data.frame(result_df)
  colnames(result_df)<-c('rel_bias','methods')
  result_df$rel_bias<-as.numeric(result_df$rel_bias)
  result_df<-result_df[!is.na(result_df$rel_bias),]
  
  png(file=paste0(plot_name,sim,'.png'),width = 1000,height = 400)
  
  p<-ggplot(result_df, aes(x=methods, y=rel_bias)) +
    geom_boxplot()+theme_classic()+labs(x = "methods",y='relative bias')+ylim(-1, 1)+
    geom_hline(yintercept=0, linetype="dashed", color = "red")
  
  print(p)
  dev.off()
}

##calico
calico<-function(ddpcr_data,cluster_num){
  colnames(ddpcr_data)<-c('Ch1','Ch2')
  ddpcr_data_subsample <- data.frame()
  
  for (i in 1:4) {
    set.seed(i)
    tmp <- ddpcr_data[sample(nrow(ddpcr_data), floor(0.9 * nrow(ddpcr_data))), ]
    tmp$Rep <- i
    ddpcr_data_subsample <- rbind(ddpcr_data_subsample, tmp)
  }
  
  
  ddpcr_data_subsample$Rep <- factor(ddpcr_data_subsample$Rep,levels=seq(1,max(ddpcr_data_subsample$Rep)))
  max_ch1 <- min(by(ddpcr_data_subsample,ddpcr_data_subsample$Rep,function(X) {X[which.max(X$Ch1),1]},simplify=T))
  max_ch2 <- min(by(ddpcr_data_subsample,ddpcr_data_subsample$Rep,function(X) {X[which.max(X$Ch2),2]},simplify=T))
  min_ch1 <- max(by(ddpcr_data_subsample,ddpcr_data_subsample$Rep,function(X) {X[which.min(X$Ch1),1]},simplify=T))
  min_ch2 <- max(by(ddpcr_data_subsample,ddpcr_data_subsample$Rep,function(X) {X[which.min(X$Ch2),2]},simplify=T))
  
  
  ddpcr_data_subset_index <- (ddpcr_data$Ch1 > min_ch1 & ddpcr_data$Ch1 < max_ch1 & ddpcr_data$Ch2 > min_ch2 & ddpcr_data$Ch2 < max_ch2)
  ddpcr_data_subset<-ddpcr_data[ddpcr_data_subset_index,]
  #convert to grid, use sliding window to average in nearby cells on grid
  #http://gis.stackexchange.com/questions/24588/converting-point-data-into-gridded-dataframe-for-histogram-analysis-using-r
  
  ddpcr_copy <- ddpcr_data_subset
  
  coordinates(ddpcr_copy) <- ~Ch2+Ch1
  ddpcr_raster <- raster(ncols=600,nrows=600)
  extent(ddpcr_raster) <- extent(ddpcr_copy)
  ddpcr_raster <- rasterize(ddpcr_copy, ddpcr_raster, 1, background = 0, fun = function(X,...) {
    if (length(X) > 0) {
      1
    }
    else { 0 }
  })
  ddpcr_raster <- flip(ddpcr_raster,direction='y')
  
  ddpcr_raster_nonzero_focal <- raster::focal(ddpcr_raster,w=matrix(1/9,nc=9,nr=9))
  ddpcr_raster_nonzero <- as.data.frame(which(as.matrix(ddpcr_raster_nonzero_focal) > 0, arr.ind = T))
  ddpcr_raster_kmeans <- kmeans(ddpcr_raster_nonzero,cluster_num,nstart = 25)
  
  #get extents of raster
  #xmin xmax ymin ymax
  raster_extents <- extent(ddpcr_raster)
  raster_extents <- c(raster_extents[1],raster_extents[2],raster_extents[3],raster_extents[4])
  
  #get real coordinates
  raster_centers <- ddpcr_raster_kmeans$centers
  raster_centers_locs <- as.matrix(t(apply(raster_centers,1,function(X) {
    tmp_y <- raster_extents[1] + ((raster_extents[2] - raster_extents[1]) * X[2]/480)
    tmp_x <- raster_extents[3] + ((raster_extents[4] - raster_extents[3]) * X[1]/480)
    c(tmp_x,tmp_y)
  })))
  
  #perform second k-means round to get real variances and centers
  set.seed(12345)
  ddpcr_kc <- kmeans(ddpcr_data_subset,centers=raster_centers_locs,trace=T,nstart=25)
  
  ddpcr_data$group<-0
  ddpcr_data$group[ddpcr_data_subset_index]<-ddpcr_kc$cluster
  
  return(ddpcr_data$group)

}


dpcp <- function(df_data, cluster_num) {
  
  notnoise.refdata <- subset(df_data,
                             df_data$group > 0)
  
  #Calculate position of reference sample centers
  cent.ref <- lapply(unique(notnoise.refdata$group), function(y) {
    vic.centers <- mean(subset(notnoise.refdata[, 1], notnoise.refdata$group == y))
    fam.centers <- mean(subset(notnoise.refdata[, 2], notnoise.refdata$group == y))
    
    cbind.data.frame("channel1" = vic.centers, "channel2" = fam.centers)
  })
  centersDB <- do.call(rbind.data.frame, cent.ref)
  
  centersDB$dist <- raster::pointDistance(centersDB, c(0, 0), lonlat = FALSE)
  
  empty <- subset(centersDB, centersDB$dist == min(centersDB$dist))
  row.names(empty) <- 1
  noEmpty <- subset(centersDB, centersDB$dist != min(centersDB$dist))
  row.names(noEmpty) <- 1:nrow(noEmpty)
  
  targetFam <- which.min(noEmpty$channel1)
  targetVic <- which.min(noEmpty$channel2)
  
  centMat <- rbind(empty[, c(1, 2)], noEmpty[targetFam, c(1, 2)],
                   noEmpty[targetVic, c(1, 2)])
  
  row.names(centMat) <- c("Empty", "Target1", "Target2")
  
  if (cluster_num==4){
    doublePos <- centMat[2, ] + centMat[3, ] - centMat[1, ]
    
    allcenters <- rbind(centMat, doublePos)
    row.names(allcenters) <- c(row.names(centMat),paste(row.names(centMat)[2], row.names(centMat)[3], sep = " + "))
  } 
  
  if (cluster_num==3){
    allcenters <- centMat
    row.names(allcenters) <- row.names(centMat)
  }
  
  
  cm1 <- e1071::cmeans(df_data[,c(1,2)], centers = allcenters, iter.max = 3,
                       dist = "euclidean", method = "ufcl")
  
  #Use row names of centers table as cluster name
  clusname <- rownames(allcenters)[sort(unique(cm1$cluster))]
  cmclus <- factor(cm1$cluster)
  levels(cmclus) <- clusname
  
  cm.first <- cbind.data.frame(df_data[,c(1,2)], "cluster" = cmclus)
  
  sep.data <- split.data.frame(cm.first[,c(1,2)], cm.first[,3])
  
  clus.cent <- lapply(sep.data, colMeans)
  
  clus.cent <- do.call(rbind, clus.cent)
  
  noclus.cent <- subset(allcenters,
                        !(rownames(allcenters) %in% row.names(clus.cent)))
  
  new.cent <- rbind.data.frame(clus.cent, noclus.cent)
  
  new.cent <- new.cent[order(match(rownames(new.cent),
                                   rownames(allcenters))), ]
  
  cm <- e1071::cmeans(df_data[,c(1,2)], centers = new.cent, iter.max = 3,
                      dist = "euclidean", method = "ufcl")
  
  #Use row names of centers table as cluster name
  clusname <- rownames(allcenters)[sort(unique(cm$cluster))]
  cmclus <- factor(cm$cluster)
  levels(cmclus) <- clusname
  
  #Use row names of centers table as membership column names
  colnames(cm$membership) <- rownames(allcenters)
  
  clusters <- cbind.data.frame(df_data[,c(1,2)], "cluster" = cmclus)
  
  clus <- table(clusters$cluster)
  trueclus <- names(clus)[clus > 4]
  
  if (length(trueclus) <= 4) {
    probability <- 0.5
  } else if (length(trueclus) >= 5 & length(trueclus) <= 7) {
    probability <- 0.5
  } else {
    probability <- 0.5
  }
  
  max.mem <- apply(cm$membership, 1, max)
  
  rain <- list("rain" = subset(clusters, max.mem < probability),
               "norain" = subset(clusters, max.mem >= probability))
  
  mahala.dis <- if (nrow(rain$rain) < 2) {
    clusters
  } else {
    sapply(rownames(allcenters), function(y) {
      norain.sub <- subset(rain$norain[, c(1, 2)], rain$norain[, 3] == y)
      
      if (nrow(norain.sub) < 3 | rcond(cov(norain.sub))  <= 1e-5) {
        eucl.dist <- raster::pointDistance(
          rain$rain[, c(1, 2)], allcenters[which(rownames(allcenters) == y), ],
          lonlat = FALSE)
      } else {
        mahalanobis(rain$rain[, c(1, 2)], colMeans(norain.sub),
                    cov(norain.sub), tol = 1e-5)
      }
    })
  }
  
  reClus <- if (!is.null(colnames(mahala.dis)) &
                all(colnames(mahala.dis) %in%  c("channel1", "channel2", "cluster"))) {
    mahala.dis
  } else {
    rainClus <- cbind(
      rain$rain[, c(1, 2)],
      "cluster" = colnames(mahala.dis)[apply(mahala.dis, 1, which.min)])
    rbind(rain$norain, rainClus)
  }
  
  return(reClus)
  
}

## perform classification
clustering_method<-function(directory,data,sim,true_group,initials,cluster_num=4,noncentrals=NULL,plot=FALSE) {
  ## argument:
  ### data: should be a matrix or data frame
  setwd(directory)
  ## kmeans as a reference
  kmeans.gr <- kmeans(data, centers = cluster_num, iter.max = 50,nstart = 100)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=kmeans.gr$cluster)
  kmeans_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if (plot) {
    plotting_func(df_data,'kmeans_cluster_',sim)
  }
  
  kmeans_adjri_noncentral<-eval_func(df_data$group[noncentrals],true_group[noncentrals])
  
  rematch_kmeans<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  kmeans_lambda1_sens<-(-log(1-(sum(rematch_kmeans==2)+sum(rematch_kmeans==3))/length(rematch_kmeans)))
  kmeans_lambda2_sens<-(-log(1-(sum(rematch_kmeans==4)+sum(rematch_kmeans==3))/length(rematch_kmeans)))
  
  kmeans_cluster1_sens<-sum(rematch_kmeans==1)
  
  df_data$rematch<-rematch_kmeans
  
  kmeans_cluster_center<-cbind(df_data %>%
    group_by(rematch) %>%
    summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='kmeans')
  
  ## kmeans with centers input as initial values
  kmeans_initials.gr <- kmeans(data, centers = initials, iter.max = 50,nstart = 100)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=kmeans_initials.gr$cluster)
  kmeans_initials_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if (plot) {
    plotting_func(df_data,'kmeans_initials_cluster_',sim)
  }
  
  kmeans_initials_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_kmeans_initials<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  kmeans_initials_lambda1_sens<-(-log(1-(sum(rematch_kmeans_initials==2)+sum(rematch_kmeans_initials==3))/length(rematch_kmeans_initials)))
  kmeans_initials_lambda2_sens<-(-log(1-(sum(rematch_kmeans_initials==4)+sum(rematch_kmeans_initials==3))/length(rematch_kmeans_initials)))
  
  kmeans_initials_cluster1_sens<-sum(rematch_kmeans_initials==1)
  
  df_data$rematch<-rematch_kmeans_initials
  
  kmeans_initials_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                          summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='kmeans_initials')
  
  ## fuzzy c means
  cmeans.gr <- cmeans(data, centers = cluster_num, iter.max = 50)
  cmeans_cluster<-apply(cmeans.gr$membership,1,which.max)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=cmeans_cluster)
  cmeans_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if (plot) {
    plotting_func(df_data,'cmeans_cluster_',sim)
  }
  
  cmeans_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_cmeans<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  cmeans_lambda1_sens<-(-log(1-(sum(rematch_cmeans==2)+sum(rematch_cmeans==3))/length(rematch_cmeans)))
  cmeans_lambda2_sens<-(-log(1-(sum(rematch_cmeans==4)+sum(rematch_cmeans==3))/length(rematch_cmeans)))
  
  cmeans_cluster1_sens<-sum(rematch_cmeans==1)
  
  df_data$rematch<-rematch_cmeans
  
  cmeans_cluster_center<-cbind(df_data %>%
                                 group_by(rematch) %>%
                                 summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='cmeans')
  
  
  ## fuzzy c means with centers input as initial values
  cmeans_initials.gr <- cmeans(data, centers = initials, iter.max = 50)
  cmeans_cluster_initials<-apply(cmeans_initials.gr$membership,1,which.max)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=cmeans_cluster_initials)
  cmeans_initials_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if (plot) {
    plotting_func(df_data,'cmeans_initials_cluster_',sim)
  }
  
  cmeans_initials_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_cmeans_initials<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  cmeans_initials_lambda1_sens<-(-log(1-(sum(rematch_cmeans_initials==2)+sum(rematch_cmeans_initials==3))/length(rematch_cmeans_initials)))
  cmeans_initials_lambda2_sens<-(-log(1-(sum(rematch_cmeans_initials==4)+sum(rematch_cmeans_initials==3))/length(rematch_cmeans_initials)))
  
  cmeans_initials_cluster1_sens<-sum(rematch_cmeans_initials==1)
  
  df_data$rematch<-rematch_cmeans_initials
  
  cmeans_initials_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                          summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='cmeans_initials')
  
  ## dbscan
  dbscan_res <- dbscan(data,eps = 0.15, minPts = 5)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=dbscan_res$cluster)
  dbscan_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  ## F_score
  if (plot) {
    plotting_func(df_data,'dbscan_cluster_',sim)
  }
  
  dbscan_adjri_noncentral<-eval_func(df_data$group[noncentrals],true_group[noncentrals])
  
  rematch_dbscan<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  dbscan_lambda1_sens<-(-log(1-(sum(rematch_dbscan==2)+sum(rematch_dbscan==3))/length(rematch_dbscan)))
  dbscan_lambda2_sens<-(-log(1-(sum(rematch_dbscan==4)+sum(rematch_dbscan==3))/length(rematch_dbscan)))
  
  dbscan_cluster1_sens<-sum(rematch_dbscan==1)
  dbscan_clusternum<-length(unique(rematch_dbscan[rematch_dbscan!=0]))
  
  df_data$rematch<-rematch_dbscan
  
  dbscan_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                          summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='dbscan')
  
  
  ## dbscan + cmeans
  tryCatch( { 
    dpcp_group <- dpcp(df_data=df_data, cluster_num=cluster_num)
    dpcp_group$cluster_ind<-ifelse(dpcp_group$cluster=='Empty',1,ifelse(dpcp_group$cluster=='Target2',2,ifelse(dpcp_group$cluster=='Target1 + Target2',3,4)))
    dpcp_group<-dpcp_group[order(as.numeric(rownames(dpcp_group))),]
    df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=dpcp_group$cluster_ind)
    dpcp_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
    
    ## F_score
    if (plot) {
      plotting_func(df_data,'dpcp_cluster_',sim)
    }
    
    dpcp_adjri_noncentral<-eval_func(df_data$group[noncentrals],true_group[noncentrals])
    
    df_data$rematch<-dpcp_group$cluster_ind
    
    dpcp_lambda1_sens<-(-log(1-(sum(df_data$group==2)+sum(df_data$group==3))/length(df_data$group)))
    dpcp_lambda2_sens<-(-log(1-(sum(df_data$group==4)+sum(df_data$group==3))/length(df_data$group)))
    
    dpcp_cluster1_sens<-sum(df_data$group==1)
    
    dpcp_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                   summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='dpcp')
  },error = function(e) {
    dpcp_adjri<<-NA
    dpcp_adjri_noncentral<<-NA
    dpcp_lambda1_sens<<-NA
    dpcp_lambda2_sens<<-NA
    dpcp_cluster1_sens<<-NA
    dpcp_cluster_center<<-cbind(rematch=NA,channel1_center=NA,channel2_center=NA,method='dpcp')
  
    })
    
  ## flowsom
  fsom<-list()
  fsom$map<-SOM(as.matrix(data))
  fsom<-BuildMST(fsom)
  ## hc (hierarchical clustering is used?) in meta-clustering
  metaClustering <- as.character(metaClustering_consensus(fsom$map$codes,k = cluster_num))
  final_cluster<-metaClustering[fsom$map$mapping[, 1]]
  
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=as.numeric(final_cluster))
  fsom_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if(plot){
    plotting_func(df_data,'flowsom_cluster_',sim)
  }
  
  fsom_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_fsom<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  fsom_lambda1_sens<-(-log(1-(sum(rematch_fsom==2)+sum(rematch_fsom==3))/length(rematch_fsom)))
  fsom_lambda2_sens<-(-log(1-(sum(rematch_fsom==4)+sum(rematch_fsom==3))/length(rematch_fsom)))
  
  fsom_cluster1_sens<-sum(rematch_fsom==1)
  
  df_data$rematch<-rematch_fsom
  
  fsom_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                 summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='fsom')
  
  ## flowpeaks
  fp<-flowPeaks(df_data[,c(1,2)])
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=fp$peaks.cluster)
  flowpeaks_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if(plot){
    plotting_func(df_data,'flowpeaks_',sim)
  }
  
  flowpeaks_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_flowpeaks<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  flowpeaks_lambda1_sens<-(-log(1-(sum(rematch_flowpeaks==2)+sum(rematch_flowpeaks==3))/length(rematch_flowpeaks)))
  flowpeaks_lambda2_sens<-(-log(1-(sum(rematch_flowpeaks==4)+sum(rematch_flowpeaks==3))/length(rematch_flowpeaks)))
  
  flowpeaks_cluster1_sens<-sum(rematch_flowpeaks==1)
  flowpeaks_clusternum<-length(unique(rematch_flowpeaks))
  
  df_data$rematch<-rematch_flowpeaks
  
  flowpeaks_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                               summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='flowpeaks')
  
  ## flowclust
  ## cluster number defined beforehand
  df_data<-data.frame(channel1=data[,1],channel2=data[,2])
  res <- flowClust(df_data, varNames=c("channel1", "channel2"), K=cluster_num)
  
  df_data<-data.frame(channel1=data[,1],channel2=data[,2], group=res@label)
  flowclust_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if(plot){
    plotting_func(df_data,'flowclust_',sim)
  }
  
  flowclust_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_flowclust<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  flowclust_lambda1_sens<-(-log(1-(sum(rematch_flowclust==2)+sum(rematch_flowclust==3))/length(rematch_flowclust)))
  flowclust_lambda2_sens<-(-log(1-(sum(rematch_flowclust==4)+sum(rematch_flowclust==3))/length(rematch_flowclust)))
  
  flowclust_cluster1_sens<-sum(rematch_flowclust==1)
  
  df_data$rematch<-rematch_flowclust
  
  flowclust_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                    summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='flowclust')
  
  ## flowclust
  ### cluster number determined automatically
  df_data<-data.frame(channel1=data[,1],channel2=data[,2])
  res <- flowClust(df_data, varNames=c("channel1", "channel2"), K=1:cluster_num)
  
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=res[[which.max(BIC(res))]]@label)
  flowclust_auto_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if(plot){
    plotting_func(df_data,'flowclust_auto_',sim)
  }
  
  flowclust_auto_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_flowclust_auto<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  flowclust_auto_lambda1_sens<-(-log(1-(sum(rematch_flowclust_auto==2)+sum(rematch_flowclust_auto==3))/length(rematch_flowclust_auto)))
  flowclust_auto_lambda2_sens<-(-log(1-(sum(rematch_flowclust_auto==4)+sum(rematch_flowclust_auto==3))/length(rematch_flowclust_auto)))
  
  flowclust_auto_cluster1_sens<-sum(rematch_flowclust_auto==1)
  flowclust_auto_clusternum<-length(unique(rematch_flowclust_auto))
  
  df_data$rematch<-rematch_flowclust_auto
  
  flowclust_auto_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                    summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='flowclust_auto')
  
  
  ## flowmerge with more clusters
  data("rituximab")
  exprs(rituximab)<-rbind(exprs(rituximab),matrix(0,nrow=nrow(data)-nrow(exprs(rituximab)),ncol=8))
  rituximab@exprs[,c(1,2)]<-as.matrix(data)
  flowClust.res <- flowClust(rituximab, varNames=c(colnames(rituximab)[1:2]), K=1:(cluster_num+2),level=1)
  flowClust.maxBIC<-flowClust.res[[which.max(BIC(flowClust.res))]];
  
  ## ----stage3-------------------------------------------------------------------
  flowClust.flowobj<-flowObj(flowClust.maxBIC,rituximab);
  flowClust.merge<-merge(flowClust.flowobj,metric="entropy");
  i<-fitPiecewiseLinreg(flowClust.merge);
  flowClust.mergeopt<-flowClust.merge[[i]];
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=flowClust.mergeopt@label)
  flowmerge_moreclust_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  
  if(plot){
    plotting_func(df_data,'flowmerge_moreclust_',sim)
  }
  
  flowmerge_moreclust_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_flowmerge_moreclust<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  flowmerge_moreclust_lambda1_sens<-(-log(1-(sum(rematch_flowmerge_moreclust==2)+sum(rematch_flowmerge_moreclust==3))/length(rematch_flowmerge_moreclust)))
  flowmerge_moreclust_lambda2_sens<-(-log(1-(sum(rematch_flowmerge_moreclust==4)+sum(rematch_flowmerge_moreclust==3))/length(rematch_flowmerge_moreclust)))
  
  flowmerge_moreclust_cluster1_sens<-sum(rematch_flowmerge_moreclust==1)
  flowmerge_moreclust_clusternum<-length(unique(rematch_flowmerge_moreclust))
  
  df_data$rematch<-rematch_flowmerge_moreclust
  
  flowmerge_moreclust_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                         summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='flowmerge_moreclust')
  
  ## SamSpectral
  L <- SamSPECTRAL(data.points=as.matrix(data),dimensions=c(1,2),number.of.clusters=cluster_num,normal.sigma = 200,separation.factor = 0.7)
  df_data<-cbind(df_data[,c(1,2)],group=L)
  samspectral_adjri<-eval_func(class_group=df_data$group[!is.na(L)],true_group=true_group[!is.na(L)])
  
  if(plot){
    plotting_func(df_data,'samspectralclust_',sim)
  }
  
  samspectral_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals & !is.na(L)],true_group=true_group[noncentrals  & !is.na(L)])
  
  rematch_samspectral<-misclass_index(data=data[!is.na(L),c(1,2)],class_group=df_data$group[!is.na(L)],true_centers=initials)
  samspectral_lambda1_sens<-(-log(1-(sum(rematch_samspectral==2)+sum(rematch_samspectral==3))/length(rematch_samspectral)))
  samspectral_lambda2_sens<-(-log(1-(sum(rematch_samspectral==4)+sum(rematch_samspectral==3))/length(rematch_samspectral)))
  
  samspectral_cluster1_sens<-sum(rematch_samspectral==1)
  
  df_data<-df_data[!is.na(L),]
  df_data$rematch<-rematch_samspectral
  
  samspectral_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                              summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='samspectral')
  
  ## automatic one
  L1 <- SamSPECTRAL(data.points=as.matrix(data),dimensions=c(1,2),normal.sigma = 200,separation.factor = 0.7)
  df_data<-cbind(data[,c(1,2)],group=L1)
  
  samspectral_auto_adjri<-eval_func(class_group=df_data$group[!is.na(L1)],true_group=true_group[!is.na(L1)])
  
  
  if(plot){
    plotting_func(df_data,'samspectralclust_auto_',sim)
  }
  
  samspectral_auto_adjri_noncentral<-eval_func(class_group=df_data$group[!is.na(L1) & noncentrals],true_group=true_group[!is.na(L1) & noncentrals])
  
  rematch_samspectral_auto<-misclass_index(data=data[!is.na(L1),c(1,2)],class_group=df_data$group[!is.na(L1)],true_centers=initials)
  samspectral_auto_lambda1_sens<-(-log(1-(sum(rematch_samspectral_auto==2)+sum(rematch_samspectral_auto==3))/length(rematch_samspectral_auto)))
  samspectral_auto_lambda2_sens<-(-log(1-(sum(rematch_samspectral_auto==4)+sum(rematch_samspectral_auto==3))/length(rematch_samspectral_auto)))
  
  samspectral_auto_cluster1_sens<-sum(rematch_samspectral_auto==1)
  samspectral_auto_clusternum<-length(unique(rematch_samspectral_auto))
  
  df_data<-df_data[!is.na(L1),]
  df_data$rematch<-rematch_samspectral_auto
  
  samspectral_auto_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                      summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='samspectral_auto')
  
  ## dpcr 2d methods
  ## calico
  calico.gr <- calico(data, cluster_num = cluster_num)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=calico.gr)
  calico_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  
  if (plot) {
    plotting_func(df_data,'calico_cluster_',sim)
  }
  
  calico_adjri_noncentral<-eval_func(df_data$group[noncentrals],true_group[noncentrals])
  
  rematch_calico<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  calico_lambda1_sens<-(-log(1-(sum(rematch_calico==2)+sum(rematch_calico==3))/length(rematch_calico)))
  calico_lambda2_sens<-(-log(1-(sum(rematch_calico==4)+sum(rematch_calico==3))/length(rematch_calico)))
  
  calico_cluster1_sens<-sum(rematch_calico==1)
  
  df_data$rematch<-rematch_calico
  
  calico_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                           summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='calico')
  
  
  adj_ri<-data.frame(kmeans=kmeans_adjri,kmeans_initials=kmeans_initials_adjri,cmeans=cmeans_adjri,cmeans_initials=cmeans_initials_adjri,dbscan=dbscan_adjri,dpcp=dpcp_adjri,fsom=fsom_adjri,flowpeaks=flowpeaks_adjri,
                     flowclust=flowclust_adjri,flowclust_auto=flowclust_auto_adjri,flowmerge_moreclust=flowmerge_moreclust_adjri,samspectral=samspectral_adjri,samspectral_auto=samspectral_auto_adjri,
                     calico=calico_adjri)
  
  adj_ri_noncentral<-data.frame(kmeans=kmeans_adjri_noncentral,kmeans_initials=kmeans_initials_adjri_noncentral,cmeans=cmeans_adjri_noncentral,cmeans_initials=cmeans_initials_adjri_noncentral,dbscan=dbscan_adjri_noncentral,dpcp=dpcp_adjri_noncentral,
                                fsom=fsom_adjri_noncentral,flowpeaks=flowpeaks_adjri_noncentral,flowclust=flowclust_adjri_noncentral,flowclust_auto=flowclust_auto_adjri_noncentral,flowmerge_moreclust=flowmerge_moreclust_adjri_noncentral,
                                samspectral=samspectral_adjri_noncentral,samspectral_auto=samspectral_auto_adjri_noncentral,calico=calico_adjri_noncentral)
  
  lambda1_sens<-data.frame(kmeans=kmeans_lambda1_sens,kmeans_initials=kmeans_initials_lambda1_sens,cmeans=cmeans_lambda1_sens,cmeans_initials=cmeans_initials_lambda1_sens,dbscan=dbscan_lambda1_sens,dpcp=dpcp_lambda1_sens,
                           fsom=fsom_lambda1_sens,flowpeaks=flowpeaks_lambda1_sens,flowclust=flowclust_lambda1_sens,flowclust_auto=flowclust_auto_lambda1_sens,flowmerge_moreclust=flowmerge_moreclust_lambda1_sens,
                           samspectral=samspectral_lambda1_sens,samspectral_auto=samspectral_auto_lambda1_sens,calico=calico_lambda1_sens)
  
  lambda2_sens<-data.frame(kmeans=kmeans_lambda2_sens,kmeans_initials=kmeans_initials_lambda2_sens,cmeans=cmeans_lambda2_sens,cmeans_initials=cmeans_initials_lambda2_sens,dbscan=dbscan_lambda2_sens,dpcp=dpcp_lambda2_sens,
                           fsom=fsom_lambda2_sens,flowpeaks=flowpeaks_lambda2_sens,flowclust=flowclust_lambda2_sens,flowclust_auto=flowclust_auto_lambda2_sens,flowmerge_moreclust=flowmerge_moreclust_lambda2_sens,
                           samspectral=samspectral_lambda2_sens,samspectral_auto=samspectral_auto_lambda2_sens,calico=calico_lambda2_sens)
  
  cluster1_sens<-data.frame(kmeans=kmeans_cluster1_sens,kmeans_initials=kmeans_initials_cluster1_sens,cmeans=cmeans_cluster1_sens,cmeans_initials=cmeans_initials_cluster1_sens,dbscan=dbscan_cluster1_sens,dpcp=dpcp_cluster1_sens,
                            fsom=fsom_cluster1_sens,flowpeaks=flowpeaks_cluster1_sens,flowclust=flowclust_cluster1_sens,flowclust_auto=flowclust_auto_cluster1_sens,flowmerge_moreclust=flowmerge_moreclust_cluster1_sens,
                            samspectral=samspectral_cluster1_sens,samspectral_auto=samspectral_auto_cluster1_sens,calico=calico_cluster1_sens)
  
  clusternum_count<-data.frame(dbscan=dbscan_clusternum,flowpeaks=flowpeaks_clusternum,flowclust_auto=flowclust_auto_clusternum,flowmerge_moreclust=flowmerge_moreclust_clusternum,
                               samspectral_auto=samspectral_auto_clusternum)
  
  cluster_centers<-rbind(kmeans_cluster_center,kmeans_initials_cluster_center,cmeans_cluster_center,cmeans_initials_cluster_center,dbscan_cluster_center,dpcp_cluster_center,
                         fsom_cluster_center,flowpeaks_cluster_center,flowclust_cluster_center,flowclust_auto_cluster_center,flowmerge_moreclust_cluster_center,
                         samspectral_cluster_center,samspectral_auto_cluster_center,calico_cluster_center)
  
  return(list(adj_ri,adj_ri_noncentral,lambda1_sens,lambda2_sens,cluster1_sens,clusternum_count,cluster_centers))
  
}


## perform classification
clustering_method_3clusters<-function(directory,data,sim,true_group,initials,cluster_num=3,noncentrals=NULL,plot=FALSE) {
  ## argument:
  ### data: should be a matrix or data frame
  setwd(directory)
  
  ## kmeans as a reference
  kmeans.gr <- kmeans(data, centers = cluster_num, iter.max = 50,nstart = 100)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=kmeans.gr$cluster)
  kmeans_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  
  if (plot) {
    plotting_func(df_data,'kmeans_cluster_',sim)
  }
  
  kmeans_adjri_noncentral<-eval_func(df_data$group[noncentrals],true_group[noncentrals])
  
  rematch_kmeans<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  kmeans_lambda1_sens<-(-log(1-(sum(rematch_kmeans==2))/length(rematch_kmeans)))
  kmeans_lambda2_sens<-(-log(1-(sum(rematch_kmeans==3))/length(rematch_kmeans)))
  
  kmeans_cluster1_sens<-sum(rematch_kmeans==1)
  
  df_data$rematch<-rematch_kmeans
  
  kmeans_cluster_center<-cbind(df_data %>%
                                 group_by(rematch) %>%
                                 summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='kmeans')
  
  ## kmeans with centers input as initial values
  kmeans_initials.gr <- kmeans(data, centers = initials, iter.max = 50,nstart = 100)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=kmeans_initials.gr$cluster)
  kmeans_initials_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if (plot) {
    plotting_func(df_data,'kmeans_initials_cluster_',sim)
  }
  
  kmeans_initials_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_kmeans_initials<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  kmeans_initials_lambda1_sens<-(-log(1-(sum(rematch_kmeans_initials==2))/length(rematch_kmeans_initials)))
  kmeans_initials_lambda2_sens<-(-log(1-(sum(rematch_kmeans_initials==3))/length(rematch_kmeans_initials)))
  
  kmeans_initials_cluster1_sens<-sum(rematch_kmeans_initials==1)
  
  df_data$rematch<-rematch_kmeans_initials
  
  kmeans_initials_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                          summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='kmeans_initials')
  
  
  ## fuzzy c means
  cmeans.gr <- cmeans(data, centers = cluster_num, iter.max = 50)
  cmeans_cluster<-apply(cmeans.gr$membership,1,which.max)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=cmeans_cluster)
  cmeans_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if (plot) {
    plotting_func(df_data,'cmeans_cluster_',sim)
  }
  
  cmeans_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_cmeans<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  cmeans_lambda1_sens<-(-log(1-(sum(rematch_cmeans==2))/length(rematch_cmeans)))
  cmeans_lambda2_sens<-(-log(1-(sum(rematch_cmeans==3))/length(rematch_cmeans)))
  
  cmeans_cluster1_sens<-sum(rematch_cmeans==1)
  
  df_data$rematch<-rematch_cmeans
  
  cmeans_cluster_center<-cbind(df_data %>%
                                 group_by(rematch) %>%
                                 summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='cmeans')
  
  
  ## fuzzy c means with centers input as initial values
  cmeans_initials.gr <- cmeans(data, centers = initials, iter.max = 50)
  cmeans_cluster_initials<-apply(cmeans_initials.gr$membership,1,which.max)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=cmeans_cluster_initials)
  cmeans_initials_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if (plot) {
    plotting_func(df_data,'cmeans_initials_cluster_',sim)
  }
  
  cmeans_initials_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_cmeans_initials<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  cmeans_initials_lambda1_sens<-(-log(1-(sum(rematch_cmeans_initials==2))/length(rematch_cmeans_initials)))
  cmeans_initials_lambda2_sens<-(-log(1-(sum(rematch_cmeans_initials==3))/length(rematch_cmeans_initials)))
  
  cmeans_initials_cluster1_sens<-sum(rematch_cmeans_initials==1)
  
  df_data$rematch<-rematch_cmeans_initials
  
  cmeans_initials_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                          summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='cmeans_initials')
  
  ## dbscan
  dbscan_res <- dbscan(data,eps = 0.15, minPts = 5)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=dbscan_res$cluster)
  dbscan_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  ## F_score
  if (plot) {
    plotting_func(df_data,'dbscan_cluster_',sim)
  }
  
  dbscan_adjri_noncentral<-eval_func(df_data$group[noncentrals],true_group[noncentrals])
  
  rematch_dbscan<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  dbscan_lambda1_sens<-(-log(1-(sum(rematch_dbscan==2))/length(rematch_dbscan)))
  dbscan_lambda2_sens<-(-log(1-(sum(rematch_dbscan==3))/length(rematch_dbscan)))
  
  dbscan_cluster1_sens<-sum(rematch_dbscan==1)
  dbscan_clusternum<-length(unique(rematch_dbscan[rematch_dbscan!=0]))
  
  df_data$rematch<-rematch_dbscan
  
  dbscan_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                 summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='dbscan')
  
  
  ## dbscan + cmeans
  tryCatch( { 
    dpcp_group <- dpcp(df_data=df_data, cluster_num=cluster_num)
    dpcp_group$cluster_ind<-ifelse(dpcp_group$cluster=='Empty',1,ifelse(dpcp_group$cluster=='Target2',2,ifelse(dpcp_group$cluster=='Target1 + Target2',3,4)))
    dpcp_group<-dpcp_group[order(as.numeric(rownames(dpcp_group))),]
    df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=dpcp_group$cluster_ind)
    dpcp_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
    
    ## F_score
    if (plot) {
      plotting_func(df_data,'dpcp_cluster_',sim)
    }
    
    dpcp_adjri_noncentral<-eval_func(df_data$group[noncentrals],true_group[noncentrals])
    
    df_data$rematch<-dpcp_group$cluster_ind
    
    dpcp_lambda1_sens<-(-log(1-(sum(df_data$group==2))/length(df_data$group)))
    dpcp_lambda2_sens<-(-log(1-(sum(df_data$group==3))/length(df_data$group)))
    
    dpcp_cluster1_sens<-sum(df_data$group==1)
    
    dpcp_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                 summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='dpcp')
  },error = function(e) {
    dpcp_adjri<<-NA
    dpcp_adjri_noncentral<<-NA
    dpcp_lambda1_sens<<-NA
    dpcp_lambda2_sens<<-NA
    dpcp_cluster1_sens<<-NA
    dpcp_cluster_center<<-cbind(rematch=NA,channel1_center=NA,channel2_center=NA,method='dpcp')
    
  })
  
  ## flowsom
  fsom<-list()
  fsom$map<-SOM(as.matrix(data))
  fsom<-BuildMST(fsom)
  ## hc (hierarchical clustering is used?) in meta-clustering
  metaClustering <- as.character(metaClustering_consensus(fsom$map$codes,k = cluster_num))
  final_cluster<-metaClustering[fsom$map$mapping[, 1]]
  
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=as.numeric(final_cluster))
  fsom_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if(plot){
    plotting_func(df_data,'flowsom_cluster_',sim)
  }
  
  fsom_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_fsom<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  fsom_lambda1_sens<-(-log(1-(sum(rematch_fsom==2))/length(rematch_fsom)))
  fsom_lambda2_sens<-(-log(1-(sum(rematch_fsom==3))/length(rematch_fsom)))
  
  fsom_cluster1_sens<-sum(rematch_fsom==1)
  
  df_data$rematch<-rematch_fsom
  
  fsom_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                               summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='fsom')
  
  ## flowpeaks
  fp<-flowPeaks(df_data[,c(1,2)])
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=fp$peaks.cluster)
  flowpeaks_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if(plot){
    plotting_func(df_data,'flowpeaks_',sim)
  }
  
  flowpeaks_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_flowpeaks<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  flowpeaks_lambda1_sens<-(-log(1-(sum(rematch_flowpeaks==2))/length(rematch_flowpeaks)))
  flowpeaks_lambda2_sens<-(-log(1-(sum(rematch_flowpeaks==3))/length(rematch_flowpeaks)))
  
  flowpeaks_cluster1_sens<-sum(rematch_flowpeaks==1)
  flowpeaks_clusternum<-length(unique(rematch_flowpeaks))
  
  df_data$rematch<-rematch_flowpeaks
  
  flowpeaks_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                    summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='flowpeaks')
  
  ## flowclust
  ## cluster number defined beforehand
  df_data<-data.frame(channel1=data[,1],channel2=data[,2])
  res <- flowClust(df_data, varNames=c("channel1", "channel2"), K=cluster_num)
  
  df_data<-data.frame(channel1=data[,1],channel2=data[,2], group=res@label)
  flowclust_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if(plot){
    plotting_func(df_data,'flowclust_',sim)
  }
  
  flowclust_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_flowclust<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  flowclust_lambda1_sens<-(-log(1-(sum(rematch_flowclust==2))/length(rematch_flowclust)))
  flowclust_lambda2_sens<-(-log(1-(sum(rematch_flowclust==3))/length(rematch_flowclust)))
  
  flowclust_cluster1_sens<-sum(rematch_flowclust==1)
  
  df_data$rematch<-rematch_flowclust
  
  flowclust_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                    summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='flowclust')
  
  ## flowclust
  ### cluster number determined automatically
  df_data<-data.frame(channel1=data[,1],channel2=data[,2])
  res <- flowClust(df_data, varNames=c("channel1", "channel2"), K=1:cluster_num)
  
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=res[[which.max(BIC(res))]]@label)
  flowclust_auto_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  if(plot){
    plotting_func(df_data,'flowclust_auto_',sim)
  }
  
  flowclust_auto_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_flowclust_auto<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  flowclust_auto_lambda1_sens<-(-log(1-(sum(rematch_flowclust_auto==2))/length(rematch_flowclust_auto)))
  flowclust_auto_lambda2_sens<-(-log(1-(sum(rematch_flowclust_auto==3))/length(rematch_flowclust_auto)))
  
  flowclust_auto_cluster1_sens<-sum(rematch_flowclust_auto==1)
  flowclust_auto_clusternum<-length(unique(rematch_flowclust_auto))
  
  df_data$rematch<-rematch_flowclust_auto
  
  flowclust_auto_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                         summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='flowclust_auto')
  
  
  ## flowmerge with more clusters
  data("rituximab")
  exprs(rituximab)<-rbind(exprs(rituximab),matrix(0,nrow=nrow(data)-nrow(exprs(rituximab)),ncol=8))
  rituximab@exprs[,c(1,2)]<-as.matrix(data)
  flowClust.res <- flowClust(rituximab, varNames=c(colnames(rituximab)[1:2]), K=1:(cluster_num+2),level=1)
  flowClust.maxBIC<-flowClust.res[[which.max(BIC(flowClust.res))]];
  
  ## ----stage3-------------------------------------------------------------------
  flowClust.flowobj<-flowObj(flowClust.maxBIC,rituximab);
  flowClust.merge<-merge(flowClust.flowobj,metric="entropy");
  i<-fitPiecewiseLinreg(flowClust.merge);
  flowClust.mergeopt<-flowClust.merge[[i]];
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=flowClust.mergeopt@label)
  flowmerge_moreclust_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  
  if(plot){
    plotting_func(df_data,'flowmerge_moreclust_',sim)
  }
  
  flowmerge_moreclust_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals],true_group=true_group[noncentrals])
  
  rematch_flowmerge_moreclust<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  flowmerge_moreclust_lambda1_sens<-(-log(1-(sum(rematch_flowmerge_moreclust==2))/length(rematch_flowmerge_moreclust)))
  flowmerge_moreclust_lambda2_sens<-(-log(1-(sum(rematch_flowmerge_moreclust==3))/length(rematch_flowmerge_moreclust)))
  
  flowmerge_moreclust_cluster1_sens<-sum(rematch_flowmerge_moreclust==1)
  flowmerge_moreclust_clusternum<-length(unique(rematch_flowmerge_moreclust))
  
  df_data$rematch<-rematch_flowmerge_moreclust
  
  flowmerge_moreclust_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                              summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='flowmerge_moreclust')
  
  ## SamSpectral
  L <- SamSPECTRAL(data.points=as.matrix(data),dimensions=c(1,2),number.of.clusters=cluster_num,normal.sigma = 200,separation.factor = 0.7)
  df_data<-cbind(df_data[,c(1,2)],group=L)
  samspectral_adjri<-eval_func(class_group=df_data$group[!is.na(L)],true_group=true_group[!is.na(L)])
  
  if(plot){
    plotting_func(df_data,'samspectralclust_',sim)
  }
  
  samspectral_adjri_noncentral<-eval_func(class_group=df_data$group[noncentrals & !is.na(L)],true_group=true_group[noncentrals  & !is.na(L)])
  
  rematch_samspectral<-misclass_index(data=data[!is.na(L),c(1,2)],class_group=df_data$group[!is.na(L)],true_centers=initials)
  samspectral_lambda1_sens<-(-log(1-(sum(rematch_samspectral==2))/length(rematch_samspectral)))
  samspectral_lambda2_sens<-(-log(1-(sum(rematch_samspectral==3))/length(rematch_samspectral)))
  
  samspectral_cluster1_sens<-sum(rematch_samspectral==1)
  
  df_data$rematch<-rematch_samspectral
  
  samspectral_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                      summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='samspectral')
  
  ## automatic one
  L1 <- SamSPECTRAL(data.points=as.matrix(data),dimensions=c(1,2),normal.sigma = 200,separation.factor = 0.7)
  df_data<-cbind(df_data[,c(1,2)],group=L1)
  
  samspectral_auto_adjri<-eval_func(class_group=df_data$group[!is.na(L1)],true_group=true_group[!is.na(L1)])
  
  
  if(plot){
    plotting_func(df_data,'samspectralclust_auto_',sim)
  }
  
  samspectral_auto_adjri_noncentral<-eval_func(class_group=df_data$group[!is.na(L1) & noncentrals],true_group=true_group[!is.na(L1) & noncentrals])
  
  rematch_samspectral_auto<-misclass_index(data=data[!is.na(L1),c(1,2)],class_group=df_data$group[!is.na(L1)],true_centers=initials)
  samspectral_auto_lambda1_sens<-(-log(1-(sum(rematch_samspectral_auto==2))/length(rematch_samspectral_auto)))
  samspectral_auto_lambda2_sens<-(-log(1-(sum(rematch_samspectral_auto==3))/length(rematch_samspectral_auto)))
  
  samspectral_auto_cluster1_sens<-sum(rematch_samspectral_auto==1)
  samspectral_auto_clusternum<-length(unique(rematch_samspectral_auto))
  
  df_data$rematch<-rematch_samspectral_auto
  
  samspectral_auto_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                           summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='samspectral_auto')
  ## dpcr 2d methods
  ## calico
  calico.gr <- calico(data, cluster_num = cluster_num)
  df_data<-data.frame(channel1=data[,1],channel2=data[,2],group=calico.gr)
  calico_adjri<-eval_func(class_group=df_data$group,true_group=true_group)
  
  
  if (plot) {
    plotting_func(df_data,'calico_cluster_',sim)
  }
  
  calico_adjri_noncentral<-eval_func(df_data$group[noncentrals],true_group[noncentrals])
  
  rematch_calico<-misclass_index(data=data[,c(1,2)],class_group=df_data$group,true_centers=initials)
  calico_lambda1_sens<-(-log(1-(sum(rematch_calico==2))/length(rematch_calico)))
  calico_lambda2_sens<-(-log(1-(sum(rematch_calico==3))/length(rematch_calico)))
  
  calico_cluster1_sens<-sum(rematch_calico==1)
  
  df_data$rematch<-rematch_calico
  
  calico_cluster_center<-cbind(df_data %>%group_by(rematch) %>%
                                 summarise_at(vars(c(channel1,channel2)), list(center = mean)),method='calico')
  
  
  adj_ri<-data.frame(kmeans=kmeans_adjri,kmeans_initials=kmeans_initials_adjri,cmeans=cmeans_adjri,cmeans_initials=cmeans_initials_adjri,dbscan=dbscan_adjri,dpcp=dpcp_adjri,fsom=fsom_adjri,flowpeaks=flowpeaks_adjri,
                     flowclust=flowclust_adjri,flowclust_auto=flowclust_auto_adjri,flowmerge_moreclust=flowmerge_moreclust_adjri,samspectral=samspectral_adjri,samspectral_auto=samspectral_auto_adjri,
                     calico=calico_adjri)
  
  adj_ri_noncentral<-data.frame(kmeans=kmeans_adjri_noncentral,kmeans_initials=kmeans_initials_adjri_noncentral,cmeans=cmeans_adjri_noncentral,cmeans_initials=cmeans_initials_adjri_noncentral,dbscan=dbscan_adjri_noncentral,dpcp=dpcp_adjri_noncentral,
                                fsom=fsom_adjri_noncentral,flowpeaks=flowpeaks_adjri_noncentral,flowclust=flowclust_adjri_noncentral,flowclust_auto=flowclust_auto_adjri_noncentral,flowmerge_moreclust=flowmerge_moreclust_adjri_noncentral,
                                samspectral=samspectral_adjri_noncentral,samspectral_auto=samspectral_auto_adjri_noncentral,calico=calico_adjri_noncentral)
  
  lambda1_sens<-data.frame(kmeans=kmeans_lambda1_sens,kmeans_initials=kmeans_initials_lambda1_sens,cmeans=cmeans_lambda1_sens,cmeans_initials=cmeans_initials_lambda1_sens,dbscan=dbscan_lambda1_sens,dpcp=dpcp_lambda1_sens,
                           fsom=fsom_lambda1_sens,flowpeaks=flowpeaks_lambda1_sens,flowclust=flowclust_lambda1_sens,flowclust_auto=flowclust_auto_lambda1_sens,flowmerge_moreclust=flowmerge_moreclust_lambda1_sens,
                           samspectral=samspectral_lambda1_sens,samspectral_auto=samspectral_auto_lambda1_sens,calico=calico_lambda1_sens)
  
  lambda2_sens<-data.frame(kmeans=kmeans_lambda2_sens,kmeans_initials=kmeans_initials_lambda2_sens,cmeans=cmeans_lambda2_sens,cmeans_initials=cmeans_initials_lambda2_sens,dbscan=dbscan_lambda2_sens,dpcp=dpcp_lambda2_sens,
                           fsom=fsom_lambda2_sens,flowpeaks=flowpeaks_lambda2_sens,flowclust=flowclust_lambda2_sens,flowclust_auto=flowclust_auto_lambda2_sens,flowmerge_moreclust=flowmerge_moreclust_lambda2_sens,
                           samspectral=samspectral_lambda2_sens,samspectral_auto=samspectral_auto_lambda2_sens,calico=calico_lambda2_sens)
  
  cluster1_sens<-data.frame(kmeans=kmeans_cluster1_sens,kmeans_initials=kmeans_initials_cluster1_sens,cmeans=cmeans_cluster1_sens,cmeans_initials=cmeans_initials_cluster1_sens,dbscan=dbscan_cluster1_sens,dpcp=dpcp_cluster1_sens,
                            fsom=fsom_cluster1_sens,flowpeaks=flowpeaks_cluster1_sens,flowclust=flowclust_cluster1_sens,flowclust_auto=flowclust_auto_cluster1_sens,flowmerge_moreclust=flowmerge_moreclust_cluster1_sens,
                            samspectral=samspectral_cluster1_sens,samspectral_auto=samspectral_auto_cluster1_sens,calico=calico_cluster1_sens)
  
  clusternum_count<-data.frame(dbscan=dbscan_clusternum,flowpeaks=flowpeaks_clusternum,flowclust_auto=flowclust_auto_clusternum,flowmerge_moreclust=flowmerge_moreclust_clusternum,
                               samspectral_auto=samspectral_auto_clusternum)
  
  cluster_centers<-rbind(kmeans_cluster_center,kmeans_initials_cluster_center,cmeans_cluster_center,cmeans_initials_cluster_center,dbscan_cluster_center,dpcp_cluster_center,
                         fsom_cluster_center,flowpeaks_cluster_center,flowclust_cluster_center,flowclust_auto_cluster_center,flowmerge_moreclust_cluster_center,
                         samspectral_cluster_center,samspectral_auto_cluster_center,calico_cluster_center)
  
  return(list(adj_ri,adj_ri_noncentral,lambda1_sens,lambda2_sens,cluster1_sens,clusternum_count,cluster_centers))
  
}


data_sim_hpc<-function(data,sim_setup,pars_cluster,group_size,probs,sim,initials,ind=i,cluster_num=4,orthog=TRUE,singlemodal=TRUE,overlap=1,cv=0.1,rain_control=0.1,rain_perc=0.1,mu_drop=20000,sd_drop=1000,loops=10,directory=directory){
  
  ## list of lists
  true_lambda1s<-NULL
  true_lambda2s<-NULL
  true_parnum1<-NULL
  
  df_adj<-NULL
  df_adj_noncentral<-NULL
  df_lambda1_sens<-NULL
  df_lambda2_sens<-NULL
  df_cluster1_sens<-NULL
  df_clusternum<-NULL
  df_sil_orig<-NULL
  df_cluster_center<-NULL
  
  for (loop in 1:loops){
    set.seed(Sys.time()+loop+ind*10)
    p <- round(rnorm(1,mu_drop,sd_drop))
    sample_size<-rmultinom(1, p, probs)
    dim(sample_size)<-NULL
    tryCatch( { 
      ## for multi-modalities
      if (singlemodal) {
        
        sampled_data<-NULL
        
        sim_result<-skew_t_dist_pp_sample(data[[1]],sim_setup[[1]],dense=sample_size[1]/group_size[1])
        cluster_sim<-sim_result[[1]]
        cluster_group<-rep(1,nrow(cluster_sim))
        sampled_data<-rbind(sampled_data,cbind(cluster_sim,cluster_group))
        
        for (i in 2:cluster_num){
          sim_result<-skew_t_dist_pp_sample_constraints(data[[i]],sim_setup[[i]],pars_cluster[[i]],dense=sample_size[i]/group_size[i])
          cluster_sim<-sim_result[[1]]
          cluster_group<-rep(i,nrow(cluster_sim))
          sampled_data<-rbind(sampled_data,cbind(cluster_sim,cluster_group))
        }
        
        sampled_data<-data.frame(channel1=sampled_data[,1],channel2=sampled_data[,2],group=sampled_data[,3])
        
        sil_orig_all<-silhouette_coef(sampled_data,sampled_data$group,plot = FALSE,plot_name ='sim',sim='sim_sampleall')
        sil_orig_allpoints<-sil_orig_all[[1]]
        
        allpoints_index<-1:nrow(sil_orig_allpoints)
        high_sil_points_pos<-NULL
        
        for (g in 2:cluster_num){
          sil_orig_cluster_index<-sil_orig_allpoints[,1]==g
          sil_orig_cluster<-sil_orig_allpoints[sil_orig_cluster_index,3]
          sil_orig_cluster_pos<-allpoints_index[sil_orig_cluster_index]
          sil_order<-order(sil_orig_cluster)
          low_sil_num<-round(length(sil_orig_cluster)*rain_control)
          high_sil_points_pos<-c(high_sil_points_pos,sil_orig_cluster_pos[sil_order[(low_sil_num+1):length(sil_orig_cluster)]])
          sampled_data_new<-rbind(sampled_data[sampled_data$group==1,],sampled_data[high_sil_points_pos,])
        }
        
        if (rain_perc!=0){
          total_amount<-round(table(sampled_data$group)*(rain_perc/rain_control)*((1-rain_control)/(1-rain_perc)))
          sampled_data_morerain<-NULL
          for (j in 2:cluster_num){
            sim_result<-skew_t_dist_pp_sample_constraints(data[[j]],sim_setup[[j]],pars_cluster[[j]],dense=total_amount[j]/group_size[j])
            cluster_sim<-sim_result[[1]]
            cluster_group<-rep(j,nrow(cluster_sim))
            sampled_data_morerain<-rbind(sampled_data_morerain,cbind(cluster_sim,cluster_group))
          }
          
          sampled_data_morerain<-data.frame(channel1=sampled_data_morerain[,1],channel2=sampled_data_morerain[,2],group=sampled_data_morerain[,3])
          sampled_data_morerain<-rbind(sampled_data[sampled_data$group==1,],sampled_data_morerain)
          
          
          ## first generate more data, then remove the center
          sil_orig_all<-silhouette_coef(sampled_data_morerain,sampled_data_morerain$group,plot = FALSE,plot_name ='sim',sim='sim_sampleall')
          sil_orig_allpoints<-sil_orig_all[[1]]
          
          low_sil_points<-NULL
          for (l in 2:cluster_num) {
            sil_orig_allpoints_cluster<-sil_orig_allpoints[sil_orig_allpoints[,1]==l,3]
            sil_order<-order(sil_orig_allpoints_cluster)
            low_sil_num<-round(total_amount[l]*rain_control)
            low_sil_points<-rbind(low_sil_points,sampled_data_morerain[sil_orig_allpoints[,1]==l,][sil_order,][1:low_sil_num,])
          }
          
          sampled_data_new<-rbind(sampled_data_new,low_sil_points)
          
        }
        
      } else {
        ## the second distribution should have a bit of variation, let's assume the first distribution generates 80% of the data, so one big group and one small
        sampled_data<-NULL
        sim_setup2<-list()
        data2<-list()
        actual_perc<-NULL
        
        group_num<-rbinom(1,sample_size[1],prob=0.8)
        actual_perc<-c(actual_perc,group_num/sample_size[1])
        sim_result1<-skew_t_dist_pp_sample(data[[1]],sim_setup[[1]],dense=group_num/group_size[1])
        cluster_sim1<-sim_result1[[1]]
        cluster_group1<-rep(1,nrow(cluster_sim1))
        
        sim_setup_tmp<-sim_setup[[1]]
        
        orig_center<-sim_setup_tmp$dp$beta
        mah_least<-function(coord) {
          mah_dist<-mahalanobis(x = c(sim_setup_tmp$dp$beta[1],coord), center = orig_center, cov =sim_setup_tmp$dp$Omega)
          dev_overlap<-(mah_dist-overlap)^2
          return(dev_overlap)
        }
        
        new_coord<-optimize(mah_least,c(min(cluster_sim1[,2]),orig_center[2]))$minimum
        sim_setup_tmp$dp$beta<-c(sim_setup_tmp$dp$beta[1],new_coord)
        cov_mat_cv<-rnorm((ncol(cluster_sim1))*(ncol(cluster_sim1)+1)/2,1,cv)
        
        diag(sim_setup_tmp$dp$Omega)<-diag(sim_setup_tmp$dp$Omega)*cov_mat_cv[1:2]
        sim_setup_tmp$dp$Omega[row(sim_setup_tmp$dp$Omega)!=col(sim_setup_tmp$dp$Omega)]<-sim_setup_tmp$dp$Omega[row(sim_setup_tmp$dp$Omega)!=col(sim_setup_tmp$dp$Omega)]*cov_mat_cv[3]
        
        data_tmp<-data[[1]]
        data_tmp[,2]<-data_tmp[,2]-(orig_center[2]-new_coord)
        sim_result2<-skew_t_dist_pp_sample(data_tmp,sim_setup_tmp,dense=(sample_size[1]-group_num)/group_size[1])
        cluster_sim2<-sim_result2[[1]]
        cluster_group2<-rep(1,nrow(cluster_sim2))
        
        sim_setup2[[1]]<-sim_setup_tmp
        data2[[1]]<-data_tmp
        sampled_data<-rbind(sampled_data,cbind(cluster_sim1,g=cluster_group1),cbind(cluster_sim2,g=cluster_group2))
        
        for (i in 2:cluster_num){
          group_num<-rbinom(1,sample_size[i],prob=0.8)
          actual_perc<-c(actual_perc,group_num/sample_size[i])
          sim_result1<-skew_t_dist_pp_sample_constraints(data[[i]],sim_setup[[i]],pars_cluster[[i]],dense=group_num/group_size[i])
          cluster_sim1<-sim_result1[[1]]
          cluster_group1<-rep(i,nrow(cluster_sim1))
          
          sim_setup_tmp<-sim_setup[[i]]
          
          orig_center<-sim_setup_tmp$dp$beta
          mah_least<-function(coord) {
            mah_dist<-mahalanobis(x = c(sim_setup_tmp$dp$beta[1],coord), center = orig_center, cov =sim_setup_tmp$dp$Omega)
            dev_overlap<-(mah_dist-overlap)^2
            return(dev_overlap)
          }
          
          new_coord<-optimize(mah_least,c(min(cluster_sim1[,2]),orig_center[2]))$minimum
          sim_setup_tmp$dp$beta<-c(sim_setup_tmp$dp$beta[1],new_coord)
          cov_mat_cv<-rnorm((ncol(cluster_sim1))*(ncol(cluster_sim1)+1)/2,1,cv)
          
          diag(sim_setup_tmp$dp$Omega)<-diag(sim_setup_tmp$dp$Omega)*cov_mat_cv[1:2]
          sim_setup_tmp$dp$Omega[row(sim_setup_tmp$dp$Omega)!=col(sim_setup_tmp$dp$Omega)]<-sim_setup_tmp$dp$Omega[row(sim_setup_tmp$dp$Omega)!=col(sim_setup_tmp$dp$Omega)]*cov_mat_cv[3]
          
          data_tmp<-data[[i]]
          data_tmp[,2]<-data_tmp[,2]-(orig_center[2]-new_coord)
          sim_result2<-skew_t_dist_pp_sample(data_tmp,sim_setup_tmp,dense=(sample_size[i]-group_num)/group_size[i])
          cluster_sim2<-sim_result2[[1]]
          cluster_group2<-rep(i,nrow(cluster_sim2))
          
          sim_setup2[[i]]<-sim_setup_tmp
          data2[[i]]<-data_tmp
          sampled_data<-rbind(sampled_data,cbind(cluster_sim1,g=cluster_group1),cbind(cluster_sim2,g=cluster_group2))
        }
        
        sampled_data<-data.frame(channel1=sampled_data[,1],channel2=sampled_data[,2],group=sampled_data[,3])
        
        sil_orig_all<-silhouette_coef(sampled_data,sampled_data$group,plot = FALSE,plot_name ='sim',sim='sim_sampleall')
        sil_orig_allpoints<-sil_orig_all[[1]]
        
        allpoints_index<-1:nrow(sil_orig_allpoints)
        high_sil_points_pos<-NULL
        

        for (g in 2:cluster_num){
          sil_orig_cluster_index<-sil_orig_allpoints[,1]==g
          sil_orig_cluster<-sil_orig_allpoints[sil_orig_cluster_index,3]
          sil_orig_cluster_pos<-allpoints_index[sil_orig_cluster_index]
          sil_order<-order(sil_orig_cluster)
          low_sil_num<-round(length(sil_orig_cluster)*rain_control)
          high_sil_points_pos<-c(high_sil_points_pos,sil_orig_cluster_pos[sil_order[(low_sil_num+1):length(sil_orig_cluster)]])
        }
        sampled_data_new<-rbind(sampled_data[sampled_data$group==1,],sampled_data[high_sil_points_pos,])
        
        if (rain_perc!=0){
          total_amount<-round(table(sampled_data$group)*(rain_perc/rain_control)*((1-rain_control)/(1-rain_perc)))
          sampled_data_morerain<-NULL
          
          for (j in 2:cluster_num){
            sim_result1<-skew_t_dist_pp_sample_constraints(data[[j]],sim_setup[[j]],pars_cluster[[j]],dense=actual_perc[j]*total_amount[j]/group_size[j])
            cluster_sim1<-sim_result1[[1]]
            cluster_group1<-rep(j,nrow(cluster_sim1))
            
            sim_result2<-skew_t_dist_pp_sample(data2[[j]],sim_setup2[[j]],dense=(1-actual_perc[j])*total_amount[j]/group_size[j])
            cluster_sim2<-sim_result2[[1]]
            cluster_group2<-rep(j,nrow(cluster_sim2))
            sampled_data_morerain<-rbind(sampled_data_morerain,cbind(cluster_sim1,g=cluster_group1),cbind(cluster_sim2,g=cluster_group2))
          }
          
          sampled_data_morerain<-data.frame(channel1=sampled_data_morerain[,1],channel2=sampled_data_morerain[,2],group=sampled_data_morerain[,3])
          sampled_data_morerain<-rbind(sampled_data[sampled_data$group==1,],sampled_data_morerain)
            
          sil_orig_all<-silhouette_coef(sampled_data_morerain,sampled_data_morerain$group,plot = FALSE,plot_name ='sim',sim=sim)
          sil_orig_allpoints<-sil_orig_all[[1]]
          
          low_sil_points<-NULL
          for (l in 2:cluster_num) {
            sil_orig_allpoints_cluster<-sil_orig_allpoints[sil_orig_allpoints[,1]==l,3]
            sil_order<-order(sil_orig_allpoints_cluster)
            low_sil_num<-round(total_amount[l]*rain_control)
            low_sil_points<-rbind(low_sil_points,sampled_data_morerain[sil_orig_allpoints[,1]==l,][sil_order,][1:low_sil_num,])
          }
          
          sampled_data_new<-rbind(sampled_data_new,low_sil_points)
          
        }
      }
      
      ## for bad separation, just pull them towards each other
      
      ## for orthogonality
      if (!orthog){
        rotate_mat<-matrix(c(1,0.1,0.1,1),nrow=2,ncol=2)
        data_mat<-as.matrix(sampled_data_new[,c(1,2)])
        unorthog_mat<-data_mat%*%rotate_mat
        sampled_data_new[,c(1,2)]<-unorthog_mat
      }
      
      file_name<-paste0(sim,'.csv')
      write.csv(sampled_data_new,file=file_name)
      
      # sil_orig_all<-silhouette_coef(sampled_data_new,sampled_data_new$group,plot = FALSE,plot_name ='sim',sim=sim)
      # sil_orig_mean_tmp<-sil_orig_all[[2]][,2]
      
      ## looking for central points
      center_search<-central_points_func(sampled_data_new)
      non_central_index<-center_search[[2]]
      non_central_points_pos<- (1:nrow(sampled_data_new)) %in% non_central_index 
      sampled_non_central_points_pos<-non_central_points_pos

      lambda1<-(-log(1-(sum(sampled_data_new$group==2)+sum(sampled_data_new$group==3))/nrow(sampled_data_new))) # positive on x-axis
      lambda2<-(-log(1-(sum(sampled_data_new$group==4)+sum(sampled_data_new$group==3))/nrow(sampled_data_new))) # positive on y-axis
      true_lambda1s<-c(true_lambda1s,lambda1)
      true_lambda2s<-c(true_lambda2s,lambda2)
      
      parnum1<-sum(sampled_data_new$group==1)
      true_parnum1<-c(true_parnum1,parnum1)
      
      list_eval_tmp<-clustering_method(directory=directory,data=sampled_data_new[,c(1,2)],sim=sim,initials=initials,true_group=sampled_data_new$group,noncentrals=sampled_non_central_points_pos,plot=FALSE)
      df_tmp_adj<-list_eval_tmp[[1]]
      df_tmp_adj_noncentral<-list_eval_tmp[[2]]
      df_tmp_lambda1_sens<-list_eval_tmp[[3]]
      df_tmp_lambda2_sens<-list_eval_tmp[[4]]
      df_tmp_cluster1_sens<-list_eval_tmp[[5]]
      df_tmp_clusternum<-list_eval_tmp[[6]]
      df_tmp_cluster_center<-list_eval_tmp[[7]]
      
    }, error = function(e) {
      # sil_orig_mean_tmp<<-rep(NA,4)
      df_tmp_adj<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                              flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_adj_noncentral<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                         flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_lambda1_sens<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                       flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_lambda2_sens<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                       flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_cluster1_sens<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                        flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_clusternum<<-data.frame(dbscan=NA,flowpeaks=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral_auto=NA)
      df_tmp_cluster_center<<-data.frame(rematch=NA,channel1_center=NA,channel2_center=NA,method=NA)
    })
    
    # df_sil_orig<-rbind(df_sil_orig,sil_orig_mean_tmp)
    df_adj<-rbind(df_adj,df_tmp_adj)
    df_adj_noncentral<-rbind(df_adj_noncentral,df_tmp_adj_noncentral)
    df_lambda1_sens<-rbind(df_lambda1_sens,df_tmp_lambda1_sens)
    df_lambda2_sens<-rbind(df_lambda2_sens,df_tmp_lambda2_sens)
    df_cluster1_sens<-rbind(df_cluster1_sens,df_tmp_cluster1_sens)
    df_clusternum<-rbind(df_clusternum,df_tmp_clusternum)
    df_cluster_center<-rbind(df_cluster_center,df_tmp_cluster_center)
  }
  
  # file_name<-paste0('df_sil',ind,'.csv')
  # write.csv(df_sil_orig,file=file_name)
  
  df_true_values<-data.frame(true_lambda1s=true_lambda1s,true_lambda2s=true_lambda2s,true_parnum1=true_parnum1)
  file_name<-paste0('df_true_values',ind,'.csv')
  write.csv(df_true_values,file=file_name)
  
  
  file_name<-paste0('df_adj',ind,'.csv')
  write.csv(df_adj,file=file_name)
  
  file_name<-paste0('df_adj_noncentral',ind,'.csv')
  write.csv(df_adj_noncentral,file=file_name)
  
  file_name<-paste0('df_lambda1_sens',ind,'.csv')
  write.csv(df_lambda1_sens,file=file_name)
  
  file_name<-paste0('df_lambda2_sens',ind,'.csv')
  write.csv(df_lambda2_sens,file=file_name)
  
  file_name<-paste0('df_cluster1_sens',ind,'.csv')
  write.csv(df_cluster1_sens,file=file_name)
  
  file_name<-paste0('df_clusternum',ind,'.csv')
  write.csv(df_clusternum,file=file_name)
  
  file_name<-paste0('df_cluster_center',ind,'.csv')
  write.csv(df_cluster_center,file=file_name)
}

## change the initials, probs, data, sim_setup
## change the group label, 4 to 3
data_sim_hpc_3clusters<-function(data,sim_setup,group_size,probs,sim,initials,ind=i,cluster_num=3,orthog=TRUE,singlemodal=TRUE,overlap=1,cv=0.1,rain_control=0.1,rain_perc=0.1,mu_drop=20000,sd_drop=1000,loops=10,directory=directory){
  
  ## list of lists
  true_lambda1s<-NULL
  true_lambda2s<-NULL
  true_parnum1<-NULL
  
  df_adj<-NULL
  df_adj_noncentral<-NULL
  df_lambda1_sens<-NULL
  df_lambda2_sens<-NULL
  df_cluster1_sens<-NULL
  df_clusternum<-NULL
  df_sil_orig<-NULL
  df_cluster_center<-NULL
  
  for (loop in 1:loops){
    set.seed(Sys.time()+loop+ind*10)
    p <- round(rnorm(1,mu_drop,sd_drop))
    sample_size<-rmultinom(1, p, probs)
    dim(sample_size)<-NULL
    tryCatch( { 
      ## for multi-modalities
      if (singlemodal) {
        
        sampled_data<-NULL
        for (i in 1:cluster_num){
          sim_result<-skew_t_dist_pp_sample(data[[i]],sim_setup[[i]],dense=sample_size[i]/group_size[i])
          cluster_sim<-sim_result[[1]]
          cluster_group<-rep(i,nrow(cluster_sim))
          sampled_data<-rbind(sampled_data,cbind(cluster_sim,cluster_group))
        }
        
        sampled_data<-data.frame(channel1=sampled_data[,1],channel2=sampled_data[,2],group=sampled_data[,3])
        
        sil_orig_all<-silhouette_coef(sampled_data,sampled_data$group,plot = FALSE,plot_name ='sim',sim='sim_sampleall')
        sil_orig_allpoints<-sil_orig_all[[1]]
        
        allpoints_index<-1:nrow(sil_orig_allpoints)
        high_sil_points_pos<-NULL
        
        for (g in 2:cluster_num){
          sil_orig_cluster_index<-sil_orig_allpoints[,1]==g
          sil_orig_cluster<-sil_orig_allpoints[sil_orig_cluster_index,3]
          sil_orig_cluster_pos<-allpoints_index[sil_orig_cluster_index]
          sil_order<-order(sil_orig_cluster)
          low_sil_num<-round(length(sil_orig_cluster)*rain_control)
          high_sil_points_pos<-c(high_sil_points_pos,sil_orig_cluster_pos[sil_order[(low_sil_num+1):length(sil_orig_cluster)]])
          sampled_data_new<-rbind(sampled_data[sampled_data$group==1,],sampled_data[high_sil_points_pos,])
        }
        
        if (rain_perc!=0){
          total_amount<-round(table(sampled_data$group)*(rain_perc/rain_control)*((1-rain_control)/(1-rain_perc)))
          sampled_data_morerain<-NULL
          for (j in 2:cluster_num){
            sim_result<-skew_t_dist_pp_sample(data[[j]],sim_setup[[j]],dense=total_amount[j]/group_size[j])
            cluster_sim<-sim_result[[1]]
            cluster_group<-rep(j,nrow(cluster_sim))
            sampled_data_morerain<-rbind(sampled_data_morerain,cbind(cluster_sim,cluster_group))
          }
          
          sampled_data_morerain<-data.frame(channel1=sampled_data_morerain[,1],channel2=sampled_data_morerain[,2],group=sampled_data_morerain[,3])
          sampled_data_morerain<-rbind(sampled_data[sampled_data$group==1,],sampled_data_morerain)
          
          
          ## first generate more data, then remove the center
          sil_orig_all<-silhouette_coef(sampled_data_morerain,sampled_data_morerain$group,plot = FALSE,plot_name ='sim',sim='sim_sampleall')
          sil_orig_allpoints<-sil_orig_all[[1]]
          
          low_sil_points<-NULL
          for (l in 2:cluster_num) {
            sil_orig_allpoints_cluster<-sil_orig_allpoints[sil_orig_allpoints[,1]==l,3]
            sil_order<-order(sil_orig_allpoints_cluster)
            low_sil_num<-round(total_amount[l]*rain_control)
            low_sil_points<-rbind(low_sil_points,sampled_data_morerain[sil_orig_allpoints[,1]==l,][sil_order,][1:low_sil_num,])
          }
          
          sampled_data_new<-rbind(sampled_data_new,low_sil_points)
          
        }
        
      } else {
        ## the second distribution should have a bit of variation, let's assume the first distribution generates 80% of the data, so one big group and one small
        sampled_data<-NULL
        sim_setup2<-list()
        data2<-list()
        actual_perc<-NULL
        for (i in 1:cluster_num){
          group_num<-rbinom(1,sample_size[i],prob=0.8)
          actual_perc<-c(actual_perc,group_num/sample_size[i])
          sim_result1<-skew_t_dist_pp_sample(data[[i]],sim_setup[[i]],dense=group_num/group_size[i])
          cluster_sim1<-sim_result1[[1]]
          cluster_group1<-rep(i,nrow(cluster_sim1))
          
          sim_setup_tmp<-sim_setup[[i]]
          
          orig_center<-sim_setup_tmp$dp$beta
          mah_least<-function(coord) {
            mah_dist<-mahalanobis(x = c(sim_setup_tmp$dp$beta[1],coord), center = orig_center, cov =sim_setup_tmp$dp$Omega)
            dev_overlap<-(mah_dist-overlap)^2
            return(dev_overlap)
          }
          
          new_coord<-optimize(mah_least,c(min(cluster_sim1[,2]),orig_center[2]))$minimum
          sim_setup_tmp$dp$beta<-c(sim_setup_tmp$dp$beta[1],new_coord)
          cov_mat_cv<-rnorm((ncol(cluster_sim1))*(ncol(cluster_sim1)+1)/2,1,cv)
          
          diag(sim_setup_tmp$dp$Omega)<-diag(sim_setup_tmp$dp$Omega)*cov_mat_cv[1:2]
          sim_setup_tmp$dp$Omega[row(sim_setup_tmp$dp$Omega)!=col(sim_setup_tmp$dp$Omega)]<-sim_setup_tmp$dp$Omega[row(sim_setup_tmp$dp$Omega)!=col(sim_setup_tmp$dp$Omega)]*cov_mat_cv[3]
          
          data_tmp<-data[[i]]
          data_tmp[,2]<-data_tmp[,2]-(orig_center[2]-new_coord)
          sim_result2<-skew_t_dist_pp_sample(data_tmp,sim_setup_tmp,dense=(sample_size[i]-group_num)/group_size[i])
          cluster_sim2<-sim_result2[[1]]
          cluster_group2<-rep(i,nrow(cluster_sim2))
          
          sim_setup2[[i]]<-sim_setup_tmp
          data2[[i]]<-data_tmp
          sampled_data<-rbind(sampled_data,cbind(cluster_sim1,g=cluster_group1),cbind(cluster_sim2,g=cluster_group2))
        }
        
        sampled_data<-data.frame(channel1=sampled_data[,1],channel2=sampled_data[,2],group=sampled_data[,3])
        
        sil_orig_all<-silhouette_coef(sampled_data,sampled_data$group,plot = FALSE,plot_name ='sim',sim='sim_sampleall')
        sil_orig_allpoints<-sil_orig_all[[1]]
        
        allpoints_index<-1:nrow(sil_orig_allpoints)
        high_sil_points_pos<-NULL
        
        for (g in 2:cluster_num){
          sil_orig_cluster_index<-sil_orig_allpoints[,1]==g
          sil_orig_cluster<-sil_orig_allpoints[sil_orig_cluster_index,3]
          sil_orig_cluster_pos<-allpoints_index[sil_orig_cluster_index]
          sil_order<-order(sil_orig_cluster)
          low_sil_num<-round(length(sil_orig_cluster)*rain_control)
          high_sil_points_pos<-c(high_sil_points_pos,sil_orig_cluster_pos[sil_order[(low_sil_num+1):length(sil_orig_cluster)]])
          sampled_data_new<-rbind(sampled_data[sampled_data$group==1,],sampled_data[high_sil_points_pos,])
        }
        
        if (rain_perc!=0){
          total_amount<-round(table(sampled_data$group)*(rain_perc/rain_control)*((1-rain_control)/(1-rain_perc)))
          sampled_data_morerain<-NULL
          
          for (j in 2:cluster_num){
            sim_result1<-skew_t_dist_pp_sample(data[[j]],sim_setup[[j]],dense=actual_perc[j]*total_amount[j]/group_size[j])
            cluster_sim1<-sim_result1[[1]]
            cluster_group1<-rep(j,nrow(cluster_sim1))
            
            sim_result2<-skew_t_dist_pp_sample(data2[[j]],sim_setup2[[j]],dense=(1-actual_perc[j])*total_amount[j]/group_size[j])
            cluster_sim2<-sim_result2[[1]]
            cluster_group2<-rep(j,nrow(cluster_sim2))
            sampled_data_morerain<-rbind(sampled_data_morerain,cbind(cluster_sim1,g=cluster_group1),cbind(cluster_sim2,g=cluster_group2))
          }
          
          sampled_data_morerain<-data.frame(channel1=sampled_data_morerain[,1],channel2=sampled_data_morerain[,2],group=sampled_data_morerain[,3])
          sampled_data_morerain<-rbind(sampled_data[sampled_data$group==1,],sampled_data_morerain)
          
          sil_orig_all<-silhouette_coef(sampled_data_morerain,sampled_data_morerain$group,plot = FALSE,plot_name ='sim',sim=sim)
          sil_orig_allpoints<-sil_orig_all[[1]]
          
          low_sil_points<-NULL
          for (l in 2:cluster_num) {
            sil_orig_allpoints_cluster<-sil_orig_allpoints[sil_orig_allpoints[,1]==l,3]
            sil_order<-order(sil_orig_allpoints_cluster)
            low_sil_num<-round(total_amount[l]*rain_control)
            low_sil_points<-rbind(low_sil_points,sampled_data_morerain[sil_orig_allpoints[,1]==l,][sil_order,][1:low_sil_num,])
          }
          
          sampled_data_new<-rbind(sampled_data_new,low_sil_points)
          
        }
      }
      
      ## for bad separation, just pull them towards each other
      
      ## for orthogonality
      if (!orthog){
        rotate_mat<-matrix(c(1,0.1,0.1,1),nrow=2,ncol=2)
        data_mat<-as.matrix(sampled_data_new[,c(1,2)])
        unorthog_mat<-data_mat%*%rotate_mat
        sampled_data_new[,c(1,2)]<-unorthog_mat
      }
      
      file_name<-paste0(sim,'.csv')
      write.csv(sampled_data_new,file=file_name)
      
      # sil_orig_all<-silhouette_coef(sampled_data_new,sampled_data_new$group,plot = FALSE,plot_name ='sim',sim=sim)
      # sil_orig_mean_tmp<-sil_orig_all[[2]][,2]
      
      ## looking for central points
      center_search<-central_points_func(sampled_data_new)
      non_central_index<-center_search[[2]]
      non_central_points_pos<- (1:nrow(sampled_data_new)) %in% non_central_index 
      sampled_non_central_points_pos<-non_central_points_pos
      
      lambda1<-(-log(1-(sum(sampled_data_new$group==2))/nrow(sampled_data_new))) # positive on x-axis
      lambda2<-(-log(1-(sum(sampled_data_new$group==3))/nrow(sampled_data_new))) # positive on y-axis
      true_lambda1s<-c(true_lambda1s,lambda1)
      true_lambda2s<-c(true_lambda2s,lambda2)
      
      parnum1<-sum(sampled_data_new$group==1)
      true_parnum1<-c(true_parnum1,parnum1)
      
      list_eval_tmp<-clustering_method_3clusters(directory=directory,data=sampled_data_new[,c(1,2)],sim=sim,initials=initials,true_group=sampled_data_new$group,noncentrals=sampled_non_central_points_pos,plot=FALSE)
      df_tmp_adj<-list_eval_tmp[[1]]
      df_tmp_adj_noncentral<-list_eval_tmp[[2]]
      df_tmp_lambda1_sens<-list_eval_tmp[[3]]
      df_tmp_lambda2_sens<-list_eval_tmp[[4]]
      df_tmp_cluster1_sens<-list_eval_tmp[[5]]
      df_tmp_clusternum<-list_eval_tmp[[6]]
      df_tmp_cluster_center<-list_eval_tmp[[7]]
      
    }, error = function(e) {
      # sil_orig_mean_tmp<<-rep(NA,4)
      df_tmp_adj<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                              flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_adj_noncentral<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                         flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_lambda1_sens<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                       flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_lambda2_sens<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                       flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_cluster1_sens<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                        flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_clusternum<<-data.frame(dbscan=NA,flowpeaks=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral_auto=NA)
      df_tmp_cluster_center<<-data.frame(rematch=NA,channel1_center=NA,channel2_center=NA,method=NA)
    })
    
    # df_sil_orig<-rbind(df_sil_orig,sil_orig_mean_tmp)
    df_adj<-rbind(df_adj,df_tmp_adj)
    df_adj_noncentral<-rbind(df_adj_noncentral,df_tmp_adj_noncentral)
    df_lambda1_sens<-rbind(df_lambda1_sens,df_tmp_lambda1_sens)
    df_lambda2_sens<-rbind(df_lambda2_sens,df_tmp_lambda2_sens)
    df_cluster1_sens<-rbind(df_cluster1_sens,df_tmp_cluster1_sens)
    df_clusternum<-rbind(df_clusternum,df_tmp_clusternum)
    df_cluster_center<-rbind(df_cluster_center,df_tmp_cluster_center)
  }
  
  # file_name<-paste0('df_sil',ind,'.csv')
  # write.csv(df_sil_orig,file=file_name)
  
  df_true_values<-data.frame(true_lambda1s=true_lambda1s,true_lambda2s=true_lambda2s,true_parnum1=true_parnum1)
  file_name<-paste0('df_true_values',ind,'.csv')
  write.csv(df_true_values,file=file_name)
  
  
  file_name<-paste0('df_adj',ind,'.csv')
  write.csv(df_adj,file=file_name)
  
  file_name<-paste0('df_adj_noncentral',ind,'.csv')
  write.csv(df_adj_noncentral,file=file_name)
  
  file_name<-paste0('df_lambda1_sens',ind,'.csv')
  write.csv(df_lambda1_sens,file=file_name)
  
  file_name<-paste0('df_lambda2_sens',ind,'.csv')
  write.csv(df_lambda2_sens,file=file_name)
  
  file_name<-paste0('df_cluster1_sens',ind,'.csv')
  write.csv(df_cluster1_sens,file=file_name)
  
  file_name<-paste0('df_clusternum',ind,'.csv')
  write.csv(df_clusternum,file=file_name)
  
  file_name<-paste0('df_cluster_center',ind,'.csv')
  write.csv(df_cluster_center,file=file_name)
}

realdata_hpc<-function(orig_data,sample_size,sim,initials,ind,cluster_num=4,loops=10,directory=directory){
  
  ## list of lists
  true_lambda1s<-NULL
  true_lambda2s<-NULL
  true_parnum1<-NULL
  
  df_adj<-NULL
  df_adj_noncentral<-NULL
  df_lambda1_sens<-NULL
  df_lambda2_sens<-NULL
  df_cluster1_sens<-NULL
  df_clusternum<-NULL
  df_sil_orig<-NULL
  df_cluster_center<-NULL
  
  for (loop in 1:loops){
    set.seed(Sys.time()+loop+ind*10)
    
    tryCatch( { 
      ## for multi-modalities
      sampled_data_new<-orig_data[sample(1:nrow(orig_data),sample_size,replace=TRUE),]
      
      center_search<-central_points_func(sampled_data_new)
      non_central_index<-center_search[[2]]
      non_central_points_pos<- (1:nrow(sampled_data_new)) %in% non_central_index 
      sampled_non_central_points_pos<-non_central_points_pos
      
      lambda1<-(-log(1-(sum(sampled_data_new$group==2)+sum(sampled_data_new$group==3))/nrow(sampled_data_new))) # positive on x-axis
      lambda2<-(-log(1-(sum(sampled_data_new$group==4)+sum(sampled_data_new$group==3))/nrow(sampled_data_new))) # positive on y-axis
      true_lambda1s<-c(true_lambda1s,lambda1)
      true_lambda2s<-c(true_lambda2s,lambda2)
      
      parnum1<-sum(sampled_data_new$group==1)
      true_parnum1<-c(true_parnum1,parnum1)
      
      list_eval_tmp<-clustering_method(directory=directory,data=sampled_data_new[,c(1,2)],sim=sim,initials=initials,true_group=sampled_data_new$group,noncentrals=sampled_non_central_points_pos,plot=FALSE)
      df_tmp_adj<-list_eval_tmp[[1]]
      df_tmp_adj_noncentral<-list_eval_tmp[[2]]
      df_tmp_lambda1_sens<-list_eval_tmp[[3]]
      df_tmp_lambda2_sens<-list_eval_tmp[[4]]
      df_tmp_cluster1_sens<-list_eval_tmp[[5]]
      df_tmp_clusternum<-list_eval_tmp[[6]]
      df_tmp_cluster_center<-list_eval_tmp[[7]]
      
    }, error = function(e) {
      # sil_orig_mean_tmp<<-rep(NA,4)
      df_tmp_adj<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                              flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_adj_noncentral<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                         flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_lambda1_sens<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                       flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_lambda2_sens<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                       flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_cluster1_sens<<-data.frame(kmeans=NA,kmeans_initials=NA,cmeans=NA,cmeans_initials=NA,dbscan=NA,dpcp=NA,fsom=NA,flowpeaks=NA,
                                        flowclust=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral=NA,samspectral_auto=NA,calico=NA)
      df_tmp_clusternum<<-data.frame(dbscan=NA,flowpeaks=NA,flowclust_auto=NA,flowmerge_moreclust=NA,samspectral_auto=NA)
      df_tmp_cluster_center<<-data.frame(rematch=NA,channel1_center=NA,channel2_center=NA,method=NA)
    })
    
    # df_sil_orig<-rbind(df_sil_orig,sil_orig_mean_tmp)
    df_adj<-rbind(df_adj,df_tmp_adj)
    df_adj_noncentral<-rbind(df_adj_noncentral,df_tmp_adj_noncentral)
    df_lambda1_sens<-rbind(df_lambda1_sens,df_tmp_lambda1_sens)
    df_lambda2_sens<-rbind(df_lambda2_sens,df_tmp_lambda2_sens)
    df_cluster1_sens<-rbind(df_cluster1_sens,df_tmp_cluster1_sens)
    df_clusternum<-rbind(df_clusternum,df_tmp_clusternum)
    df_cluster_center<-rbind(df_cluster_center,df_tmp_cluster_center)
  }
  
  # file_name<-paste0('df_sil',ind,'.csv')
  # write.csv(df_sil_orig,file=file_name)
  
  df_true_values<-data.frame(true_lambda1s=true_lambda1s,true_lambda2s=true_lambda2s,true_parnum1=true_parnum1)
  file_name<-paste0('df_true_values',ind,'.csv')
  write.csv(df_true_values,file=file_name)
  
  
  file_name<-paste0('df_adj',ind,'.csv')
  write.csv(df_adj,file=file_name)
  
  file_name<-paste0('df_adj_noncentral',ind,'.csv')
  write.csv(df_adj_noncentral,file=file_name)
  
  file_name<-paste0('df_lambda1_sens',ind,'.csv')
  write.csv(df_lambda1_sens,file=file_name)
  
  file_name<-paste0('df_lambda2_sens',ind,'.csv')
  write.csv(df_lambda2_sens,file=file_name)
  
  file_name<-paste0('df_cluster1_sens',ind,'.csv')
  write.csv(df_cluster1_sens,file=file_name)
  
  file_name<-paste0('df_clusternum',ind,'.csv')
  write.csv(df_clusternum,file=file_name)
  
  file_name<-paste0('df_cluster_center',ind,'.csv')
  write.csv(df_cluster_center,file=file_name)
}

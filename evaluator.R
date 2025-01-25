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

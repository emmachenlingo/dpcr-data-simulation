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

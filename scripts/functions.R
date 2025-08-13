


## bin_regions: makes bins across regins (genes).Size should be same as window size from chromHMM. 
## states are needed to check that the features are within part of genome assigned to 
## states. 
## This function is called by other functions
bin_regions <- function(region,w_size,states){
  s_size=states%>%group_by(chrom)%>%summarize(min=min(start),max=max(start))
  rg_w_p = filter(region,strand=="+")%>%
    bed_makewindows(win_size = w_size)
  rg_w_m = filter(region,strand=="-")%>%
    bed_makewindows(win_size = w_size,reverse = TRUE)
  rg_w=rbind(rg_w_m,rg_w_p)
  outofchrom = rg_w%>%group_by(geneId)%>%
    summarise(start=min(start),end=max(end),chrom=unique(chrom))%>%
    full_join(s_size,by="chrom")%>%mutate(keep=ifelse(start>=min & end<=max,T,F))%>%
    filter(keep==F)
  if (nrow(outofchrom)>0){
    print("regions outside chrom")
    print(outofchrom%>%select(-keep))
    print("those are excluded")
    rg_w=rg_w%>%filter(!geneId%in%outofchrom$geneId)
  }
  rg_w
}

## id_to_genename: After aggr. clustering of apclusters, the geneIds in the cluster is
## messed up. This function gets the right geneids for each cluster.
## This function is called by other functions
id_to_genename <- function(clusterRes){
  cluster_df=tibble(geneindex=unlist(sapply(clusterRes@clusters,names)),
                    cluster=rep(1:length(clusterRes@clusters),
                                lapply(clusterRes@clusters,length)))
  names_df = tibble(geneId=colnames(clusterRes@sim),
                    id=as.character(1:ncol(clusterRes@sim)))
  tot_df = right_join(names_df,cluster_df,by=c("id"="geneindex"))%>%
    mutate(geneId=ifelse(is.na(geneId),id,geneId))%>%select(-id)%>%distinct()
  
}

## state_prob_per_feature: combines the emission matrix, the assigned states and 
## the gene positions data to get a tibble with the emission prob for each bin for
## each gene.  
state_prob_per_feature <- function(features,states,emission,w_size=200){
  rg_w=bin_regions(features,w_size,states) # generates bins across the gene regions that matches the state bins
  b_insect = bed_intersect(states,rg_w)%>%filter(.overlap==w_size)
  state_prob = full_join(b_insect ,emission,by=c("state.x"="state"))
  state_prob=state_prob%>%arrange(geneId.y,.win_id.y)
  keep=state_prob%>%select(-chrom,-starts_with(c("end","start",".","state","strand")))%>%ungroup() 
  keep
}

## cluster_to_df: combines the clustering results with emission probablilities per gene and bin.
## Also creates a combination factor that can be used to plot the genes in an honest way. 
## works also for ts clusters from the dtwclust package
cluster_to_df <- function(clust=c("ap","ts"),cluster,features,states,w_size=200){
  rg_w=bin_regions(features,w_size,states) # generates bins across the gene regions that matches the state bins
  b_insect = bed_intersect(states,rg_w)%>%filter(.overlap==w_size)
  if(clust=="ts"){
    cluster_df = tibble(geneId=names(cluster@datalist),cluster=cluster@cluster)
  } else if (clust=="ap"){
    if(length(setdiff(sapply(cluster@clusters,names),colnames(cluster@sim)))>0){
      cluster_df=id_to_genename(cluster)
    } else {
    cluster_df=tibble(geneId=unlist(sapply(cluster@clusters,names)),
                      cluster=rep(1:length(cluster@clusters),
                                  lapply(cluster@clusters,length)))
    }
  }
  else{
    stop("wrong or missing clust option")
    }
  df = full_join(cluster_df,b_insect,by=c("geneId"="geneId.y"))%>%dplyr::rename(bin=.win_id.y,state=state.x)
  df = df%>%group_by(bin,cluster)%>%mutate(n_genes=.overlap/w_size)%>%ungroup()
  df = df%>%arrange(geneId,bin)%>%
    group_by(geneId)%>%
    mutate(pat=(paste(state,collapse = "_")))%>%
    mutate(pat_l=nchar(pat))%>%group_by(pat)%>%mutate(n=n())%>%
    mutate(combination=paste(pat_l,n,pat,sep=":"))%>%ungroup()
  df = df%>%mutate(combination_factor=factor(combination,levels=sort(unique(df$combination))))
  df = df%>%select(-pat,-pat_l,-n,-combination)
  }


## plot_df_patternwise: plots the states across the genes, cluster wise
## plot_df_patternwise: plots the states across the genes, cluster wise
plot_df_patternwise <- function(df){
  # Calculate total genes per cluster for labels
  cluster_counts <- df %>%
    group_by(cluster) %>%
    summarise(total_genes = n_distinct(geneId), .groups = 'drop')
  
  # Create cluster labels with gene counts
  df <- df %>%
    left_join(cluster_counts, by = "cluster") %>%
    mutate(cluster_label = factor(paste0(cluster, ", n=", total_genes), 
                                  levels = paste0(sort(unique(cluster)), ", n=", total_genes[match(sort(unique(cluster)), cluster)])))
  
  pattern_plot=df%>%
    ggplot(aes(x=bin,y=n_genes,fill=state,col=state,group=combination_factor))+
    geom_col()+
    facet_wrap(~cluster_label,scales="free_y")+
    scale_fill_manual(values=cs)+
    scale_color_manual(values=cs)+
    scale_x_continuous(
      breaks = c(0, 5, 10, 15, 20)
    ) +
    xlab("")+
    guides(fill="none",col="none")+
    theme_minimal()+
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
  
  return(pattern_plot)
}

## plot_examplars
plot_examplars <- function(df){
  examplars <- df %>% filter(geneId%in%names(aps_results@exemplars)) %>% 
    ggplot(aes(x=bin,y=n_genes,fill=state,col=state))+
    geom_col()+
    facet_wrap(~cluster,scales="free_y")+
    scale_fill_manual(values=cs)+
    scale_color_manual(values=cs)+
    theme_minimal()+
    ylab("")+
    xlab("state")+
    guides(fill="none",col="none")+
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank())
}

## gene_pos_per_cluster: adds cluster identity to gene position, useful e.g. to 
## print bedfiles per cluster
gene_pos_per_cluster <- function(gene_pos,clusterRes,clust="ap"){
  if(clust!="ap"){stop("only works on apcluster functions")}
  if(length(setdiff(sapply(clusterRes@clusters,names),colnames(clusterRes@sim)))>0){
    cluster_df=id_to_genename(clusterRes)
  } else {
    cluster_df=cluster_df=tibble(geneId=unlist(sapply(cluster@clusters,names)),
                                 cluster=rep(1:length(cluster@clusters),
                                             lapply(cluster@clusters,length)))
  }
  gen_w_cluster = full_join(cluster_df,gene_pos,by=c("geneId"))
  
}





## not in use!
feature_to_mat <- function(feature,states,emission,w_size=200){
  bins = bin_regions(feature,w_size,states)
  b_int = bed_intersect(states,bins)%>%filter(.overlap==w_size)
  mat = emission[b_int$state.x,]
}

## not in use!
state_list_per_feature <- function(features,states,emission,w_size=200){
  feature_list = split(features,f=features$geneId)
  list_mat = lapply(feature_list,feature_to_mat,states)
}

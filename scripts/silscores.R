
clusters=ap_df%>%select(geneId,cluster)%>%unique()

distance_mat = proxy::dist(list_for_clust, method = "dtw_basic")
distance_df = as_tibble(distance_mat)%>%mutate(geneIdA=rownames(distance_mat))%>%pivot_longer(-geneIdA,names_to ="geneIdB")
distance_df = left_join(distance_df,clusters,by=c("geneIdA"="geneId"))%>%
  rename(clusterA=cluster)
distance_df = left_join(distance_df,clusters,by=c("geneIdB"="geneId"))%>%
  rename(clusterB=cluster)

silhouette_scores= function(distance_df){
## per gene, mean average to all clusters incl its own
  mean_dist =distance_df%>%group_by(geneIdA,clusterB)%>%
    filter(geneIdA!=geneIdB)%>%summarise(mean=mean(value),clusterA=unique(clusterA))
# a score (within cluster)
  a_scores =mean_dist%>%filter(clusterB==clusterA)%>%rename(A=mean)%>%select(-clusterB)
# b score (closest cluster)
  b_scores = mean_dist%>%filter(clusterB!=clusterA)%>%group_by(geneIdA)%>%summarise(B=min(mean),clusterB = clusterB[which.min(mean)]) %>%
  ungroup()
  s_scores = full_join(a_scores,b_scores)%>%mutate(sil=(B-A)/max(A,B))
#### Hurrah
}

sil = silhouette_scores(distance_df)

## test 
cluster_codes = clusters$cluster
names(cluster_codes)=clusters$geneId
sil_f=silhouette(cluster_codes[rownames(distance_mat)],distance_mat)


s_cluster_scores = sil%>%group_by(clusterA)%>%summarise(mean(sil))


s_scores%>%ggplot(aes(y=sil,x=as.factor(clusterA)))+geom_violin()

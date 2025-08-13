setwd("/groups/berger/user/vikas.shukla/States_AT/Analysis/Landscapes/AP_clustering/")
source("scripts/functions.R")
source("scripts/dtw_similarity_apcluster/dtw_similarity.R")
source("scripts/setup.R")
set.seed(997863)
#-------- Perform AP clustering and plot the clusters --------------------------
#Export all genes <400 bp and >4000 bp
gene_pos %>% filter(end-start<400 | end-start>4000) %>%
  mutate(score=".") %>% 
  select(chrom,start,end,geneId,score,strand) %>%
  write_tsv(file = "/groups/berger/user/vikas.shukla/States_AT/Analysis/Landscapes/AP_clustering/BED_files/st400_gt4000.bed",
            col_names = FALSE)

# Subsetting protein coding genes based on the gene length
gene_pos <- gene_pos %>% filter(end-start>=400 & end-start <=4000) 

# Subsetting to reduce size for tests is possible here!
gene_pos_100 <- gene_pos %>% sample_n(nrow(gene_pos))

# Include all bins overlapping, aligns with states bin
gene_pos_100 <- gene_pos_100 %>% mutate(start=start-start%%200,end=end+200-end%%200) 

# `state_prob_per_feature` combines the emission matrix, the states, and the gene positions data to 
# get a tibble with the emission prob for each bin for each gene 
for_clust <- state_prob_per_feature(gene_pos_100,states,emiss_mat,w_size=200)

# The dtw algorithm needs the data to be a list, this step needs some patience 
list_for_clust <- lapply(split(for_clust %>% select(-geneId.y), for_clust$geneId.y), as.matrix)
save(list_for_clust, file = file.path(dir_rdata,"dtw_matrix.Rdata"))

# This step takes the longest of time hence I wrote a sbatch script to do this and saved the results
if (file.exists("Rdata/aps_results.Rdata")) {
  load(file = "Rdata/aps_results.Rdata")
} else {
    aps_results = apcluster(s=dtw_similarity,x=list_for_clust,q=0)
    save(aps_results, file = file.path(dir_rdata,"aps_results.Rdata"))
}

# `cluster_to_df` combines the clustering results with the emissions prob for each bin and gene
ap_df = cluster_to_df(clust="ap",cluster=aps_results,features=gene_pos_100,states,w_size=200)
save(ap_df, file = file.path(dir_rdata,"ap_df.Rdata"))

# `plot_df_patternwise` generates the plots for each cluster
patterns <- plot_df_patternwise(ap_df)
svglite::svglite(filename = file.path(dir_plots,"all_clusters.svg"))
plot(patterns)
dev.off()
#------------------------ Plot the examplars -----------------------------------
# `plot_examplars` generates plots of the examplars
examplars <- plot_examplars(ap_df)
svglite::svglite(filename = file.path(dir_plots,"all_examplars.svg"))
plot(examplars)
dev.off()
#------------- Agglomerative clustering and Plotting the tree ------------------
# `aggExCluster` from apcluster performs agglomerative clustering of the clusters
if (file.exists("Rdata/aggres.Rdata")) {
  load(file = "Rdata/aggres.Rdata")
} else {
  aggres <- aggExCluster(x=aps_results)
  save(aggres, file = file.path(dir_rdata, "aggres.Rdata"))
}

svglite::svglite(filename = file.path(dir_plots,"Tree_v1.svg"))
plot(aggres, 
     nodePar=list(pch=NA, lab.cex=0.5), 
     main="Agglomerative clustering of clusters")
dev.off()
#------------- Plotting the clusters according to the tree order ---------------
cl_levels <- aggres@order
ap_df$cluster <- factor(ap_df$cluster, levels = cl_levels)
new_patterns <- plot_df_patternwise(ap_df)
svglite::svglite(filename = file.path(dir_plots,"all_clusters_tree_order.svg"))
plot(new_patterns)
dev.off()
#------------- Pick the number of clusters and plot the merged clusters --------
n_cl <- 28
agg_n <- aggres[[n_cl]]
agg_n_df <- cluster_to_df(clust="ap",cluster=agg_n,features=gene_pos_100,states,w_size=200)
agg_patters <- plot_df_patternwise(agg_n_df)
svglite::svglite(filename = file.path(dir_plots,"merged_clusters.svg"))
plot(agg_patters)
dev.off()
#------------- Plot merged clusters a/c to covar order -------------------------
load(file = "/groups/berger/user/vikas.shukla/States_AT/Analysis/Landscapes/AP_clustering/GO_merged/covar_levels_merged.Rdata")
agg_n_df$cluster <- factor(agg_n_df$cluster, levels = covar_levels)
agg_patters_covar <- plot_df_patternwise(agg_n_df)
plot(agg_patters_covar)
#------------- Rename the clusters and plot the merged clusters ----------------
rename_clusters28 <- c(
  "25" = "G1", "17" = "G2", "24" = "G3", "27" = "G4", "19" = "G5",
  "20" = "G6", "1" = "G7", "11" = "G8", "28" = "G9", "23" = "G10",
  "4" = "G11", "2" = "G12", "6" = "Flip", "13" = "P1", "22" = "G13",
  "21" = "P2", "16" = "P3", "12" = "P4", "9" = "P5", "18" = "P6",
  "8" = "P7", "26" = "P8", "14" = "R1", "7" = "P9", "5" = "PR1",
  "10" = "PR2", "3" = "R2", "15" = "R3") #The levels still follow covar order

# Rename clusters using the rename_clusters28 vector
agg_n_df_renamed <- agg_n_df %>%
  mutate(cluster = recode(cluster, !!!rename_clusters28))

# Set the levels of the 'cluster' column to preserve the order in rename_clusters28
agg_n_df_renamed$cluster <- factor(agg_n_df_renamed$cluster, levels = rename_clusters28)

# Plot after renaming with the preserved order
agg_patters_covar_renamed <- plot_df_patternwise(agg_n_df_renamed)
plot(agg_patters_covar_renamed)
#-------- Generate patternplots for the sub-clusters of the merged clusters ----
# Add the new cluster id to the cluster_df 
clust_with_subclust <- agg_n_df %>% 
  rename(tree_cluster=cluster) %>% 
  full_join(ap_df)

# for each new cluster (from aggres), plot the sub-clusters
sub_patterns=list()
for( i in unique(clust_with_subclust$tree_cluster)){
  df = clust_with_subclust%>%filter(tree_cluster==i)
  sub_patterns[[i]]=plot_df_patternwise(df)
}

#Add location to the sub_patterns
for (i in 1:n_cl){
  pdf(file.path(dir_subclusters, paste0("cluster_", i, ".pdf")))
  plot(sub_patterns[[i]])
  dev.off()
}

#---------------- Export bed files for the clusters ----------------------------
to_bed_with_cluster <- gene_pos_per_cluster(gene_pos_100,agg_n)
to_bed_with_cluster <- to_bed_with_cluster %>% 
  filter(cluster!="NA") %>%
  mutate(score=".") %>% 
  select(chrom, start, end, geneId, score, strand, cluster)

for(i in 1:max(to_bed_with_cluster$cluster)){
  tmp=filter(to_bed_with_cluster,cluster==i)%>%select(-cluster)
  write_delim(tmp,
              col_names=FALSE,
              file=file.path(dir_bed, paste0("cluster_", i, ".bed")), 
              delim = "\t")
}

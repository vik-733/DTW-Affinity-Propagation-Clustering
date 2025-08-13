## data to use:
source("setup_bj_data.R")
library(microbenchmark)

# pick
#gene_pos=gene_pos%>%filter(end-start>=600 & end-start<=3000)

## pick subset  
set.seed(997865)
genes_selected=gene_pos%>%sample_n(20000)
## include all bins overlapping, aligns with states bin
genes_selected=genes_selected%>%mutate(start=start-start%%200,end=end+200-end%%200)


# emission matrix (two formats)

emission=as.matrix(emiss_mat[,-1])
rownames(emission)=emiss_mat$state

### functions
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
state_prob_per_feature <- function(features,states,emission,w_size=200){
  rg_w=bin_regions(features,w_size,states) # generates bins across the gene regions that matches the state bins
  b_insect = bed_intersect(states,rg_w)%>%filter(.overlap==w_size)
  state_prob = full_join(b_insect ,emission,by=c("state.x"="state"))
  state_prob=state_prob%>%arrange(geneId.y,.win_id.y)
  keep=state_prob%>%select(-chrom,-starts_with(c("end","start",".","state","strand")))%>%ungroup() 
  keep
}


feature_to_mat <- function(feature,states,emission,w_size=200){
  bins = bin_regions(feature,w_size,states)
  b_int = bed_intersect(states,bins)%>%filter(.overlap==w_size)
  mat = emission[b_int$state.x,]
}

state_list_per_feature <- function(features,states,emission,w_size=200){
  feature_list = split(features,f=features$geneId)
  list_mat = lapply(feature_list,feature_to_mat,states,emission)
}

aa = function(list){
  list_for_clust= lapply(split(list,f=list$geneId.y),as.matrix)
}
alternative1 <- function(features,states,emission,w_width=200){
  for_clust=state_prob_per_feature(features,states,emiss_mat,w_size=200)
  list_for_clust = lapply(split(for_clust,f=for_clust$geneId.y),as.matrix)
}

alternative2 <- function(features,states,emission,w_width=200){
  for_clust=state_prob_per_feature(features,states,emiss_mat,w_size=200)
  ## split b
  X <- 76
  # Create a factor with X levels, each representing a segment of the data
  f <- cut(seq(nrow(for_clust)), breaks=X, labels=FALSE)
  # Use the factor to split the dataframe into a list of dataframes
  list_of_dfs <- split(for_clust, f)
  list_for_clust =mclapply(list_of_dfs,aa)
}


##test speed
#timing_result1 <- microbenchmark({
#  state_list_per_feature(genes_selected,states,emission=emission)
#}, times = 1)
#print(timing_result1)

## this is super slow - strange in fact - maybe check later why?
#timing_result2 <- microbenchmark({
#  alternative1(genes_selected,states,emiss_mat)
#},times=10)
#print(timing_result2)

timing_result3 <- microbenchmark({
  alternative2(genes_selected,states,emiss_mat)
},times=10)
print(timing_result3)
## alternative 2 with ncapply?


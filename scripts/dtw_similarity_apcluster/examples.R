
# This function is a small wraper script that allows us to use the dynamic time wraping (DTW) distance,
# implemented in the dtwclust package, when performing affinity propagation clustering using the apcluster.

source("dtw_similarity.R") 
library(apcluster)


# uni-variate case:
# in this toy example we include only four genomic features (here GeneA-D). We have one experiment ("H2A.X")
# Note that the number of bins can vary between genes! 

geneA = rnorm(12)
geneB = rnorm(9)
geneC = rnorm(20)
geneD = rnorm(15)

list_of_genes = list(geneA,geneB,geneC,geneD)

ap_dtw_uni = apcluster(s=dtw_similarity,x=list_of_genes,q=0)

## multi-variable case, in this small toy example we have three variables (histone variants: H2A.X","H3","H2A.W")


set.seed(71659)
geneA = matrix(rnorm(12),ncol=3,dimnames=list(c("bin1","bin2","bin3","bin4"),c("H2A.X","H3","H2A.W")))
geneB = matrix(rnorm(18),ncol=3,dimnames=list(c("bin1","bin2","bin3","bin4","bin5","bin6"),c("H2A.X","H3","H2A.W")))
geneC = matrix(rnorm(9),ncol=3,dimnames=list(c("bin1","bin2","bin3"),c("H2A.X","H3","H2A.W")))
geneD = matrix(rnorm(15),ncol=3,dimnames=list(c("bin1","bin2","bin3","bin4","bin5"),c("H2A.X","H3","H2A.W")))


list_of_matrices = list(geneA,geneB,geneC,geneD)



# cluster using apcluster (just as an example, it obviously does not make sense to cluster this data set!)
 
ap_dtw_multi <- apcluster(s=dtw_similarity,x=list_of_matrices,q=0)

# NOTE it often is a good idea to set use q=0 as results in less cluster. See also the apcluster manual.
# https://cran.r-project.org/web/packages/apcluster/vignettes/apcluster.pdf

# If the data set is large this method becomes slow. However one can use the leveraged version 
# (also explained in the apcluster manual).  It takes two extra arguments: 
# frac which is the fraction of the complete data set that should be used and 
# sweeps witch is the number of times the clustering should be done (the best result is reported). See also the apcluster manual.

apl_dtw_multi<- apclusterL(s=dtw_similarity,x=list_of_matrices,frac=0.75,sweeps=5,q=0)


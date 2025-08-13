require(dtwclust)

dtw_similarity<-function(x,sel=NA){
  if (any(is.na(sel))){
    dmat <- proxy::dist(x, method = "dtw_basic")
  }
  else{
    dmat <-proxy::dist(x=x,y=x[sel], method = "dtw_basic")
    }
  s <- exp(-as.matrix(dmat) / max(as.matrix(dmat)))
  as(s,"matrix")
}

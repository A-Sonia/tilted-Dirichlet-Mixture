#' Title
#'
#' @param par parameters: list containing two elements the eta vector of length K and the alpha matrix of size KxD
#' @param K  number of clusters
#' @param D dimension of the data
#' @param j number of cluster of interest
#' @param u0 hyperparameter
#' @param v0 hyperparameter
#' @param W observations
#' @param A allocation matrix
#'
#' @return returns the cluster adjactent to cluster j
#' @export
#'
#' @examples
adjacent=function(par,K,D,j,u0,v0,W,A){
  m=(par$alpha-v0)/(u0-v0)
  #modes=t(sapply(1:K,function(k) (m[k,]-1)/(sum(m[k,]-K))))
  modes=t(sapply(1:K,function(k) colMeans(rbind(W[A$C==k,],rep(1,D)))))
  dist=matrix(NA,K,1)
  for(i in (1:K)[-j])
  {
    dist[i,1]=sum((modes[i,]-modes[j,])^2)
  }
  which.min(dist)
}


#' Allocation function
#' Allocates observations to clusters
#'
#' @param L marginal likelihoods matrix (NxK)
#' @param N number of observations
#'
#' @return list containing two elements: a matrix Y of size NxK Y[i,k]=1 iif W_i is allocated to cluster k
#' and a vector C of size N given the number of the cluster associated with each observation
#' @export
#'
#' @examples
allocation=function(L,N){
  Y=c()
  C=c()
  for(i in 1:N)
  {
    proba=exp(L[i,])/(sum(exp(L[i,])))
    y=rmultinom(1,1,proba)#(proba==max(proba))#
    Y=rbind(Y,t(y))
    C=c(C,which(y==1))
  }
  list(Y=Y,C=C)
}

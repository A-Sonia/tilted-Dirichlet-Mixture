#' Title
#'
#' @param W angular observations: matrix of size NxD
#' @param par parameters: list containing two elements the eta vector of length K and the alpha matrix of size KxD
#' @param N number of observations
#' @param K number of clusters
#' @param D dimension of the observations
#'
#' @return a matrix of size NxK
#' @export
#'
#' @examples
Llk_dm=function(W,par,N,K,D){
  matrix(log(par$eta),ncol=K,nrow=nrow(sim$W),byrow=TRUE)+
    sapply(1:nrow(par$alpha),function(k){
      apply(sim$W,1,function(x) ddirichlet(x,par$alpha[k,],log=TRUE) )})
}



#' log-likelihood function
#'
#' @param W angular observations: matrix of size NxD
#' @param par parameters: list containing two elements the eta vector of length K and the alpha matrix of size KxD
#' @param N number of observations
#' @param K number of clusters
#' @param D dimension of the observations
#'
#' @return a matrix of size NxK
#' @export
#'
#' @examples
Llk=function(W,par,N,K,D){
  M=par$eta%*%(par$alpha/rowSums(par$alpha))
  total_L=t(apply(W,1,function(w){
    log(par$eta)+
      lgamma(rowSums(par$alpha))-
      rowSums(lgamma(par$alpha))-
      (1+rowSums(par$alpha))*log(M%*%w)[1,1]+
      as.vector(par$alpha%*%t(log(M)))+
      as.vector((par$alpha-1)%*%log(w))-
      log(D)
  }))
  total_L
}


#' Title
#'

#' @param A Allocation list (As returned by the allocation function)
#' @param L log-likelihood matrix of dimension NxK
#' @param N number of observations
#'
#' @return
#' @export
#'
#' @examples
augLlk=function(A,L,N){
  al=c()

  for(i in 1:N){
    al=c(al,L[i,A$C[i]])#-
  }#-log(sum(exp(L[i,]))))#-log(ncol(L)))
  al
}


#' simulate from a constrained dirichlet mixtures (for the simulation study)
#'
#' @param N number of observations to simulate
#' @param K number of clusters
#' @param D dimension ofthe observations
#' @param e0 hyperparameter for the distribution of eta (~dirichlet(e0))
#' @param u0 hyperparameter for the distribution of alphas
#' @param plt TRUE/FALSE plotting
#'
#' @return
#' @export
#'
#' @examples
simulation=function(N,K,D,e0,u0,plt=FALSE){
  cond=TRUE
  while(cond){
    eta=sort(as.numeric(rdirichlet(1,rep(e0,K))),decreasing=TRUE)
    alpha=diag(u0,ncol=D,nrow=K-1)+matrix(runif(D*(K-1),0,0.1),ncol=D,nrow=K-1)
    alpha=alpha/rowSums(alpha)
    alpha=rbind(alpha,(1/D-colSums(eta[-K]*alpha))/eta[K])
    cond=any(alpha<=0)
  }
  eta
  alpha
  eta%*%(alpha/rowSums(alpha))
  par=list(eta=eta,alpha=1+(u0-1)*alpha)
  W=c()
  Y=c()
  C=c()
  for(i in 1:N){
    y=rmultinom(1,1,par$eta)
    Y=rbind(Y,y)
    k=which(y==1)
    W=rbind(W,rdirichlet(1,par$alpha[k,]))
    C=c(C,k)
  }
  if(plt){
    pairs(W,col=C,xlim=c(0,1),ylim=c(0,1))
  }
  ind=order(W[,1])
  list(W=W[ind,],Y=Y[ind,],C=C[ind],par=par)
}

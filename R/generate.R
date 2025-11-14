
#' Generate n observations from a tDM
#'
#' @param n number of observations to generate
#' @param par parameter list
#' @param K number of clusters
#' @param D dimension of the data
#'
#' @return
#' @export
#'
#' @examples
generate=function(n,par,K,D){
  eta=as.numeric(par$eta)
  alpha=par$alpha
  m=colSums(eta*alpha/(rowSums(alpha)))
  X=P=W=U=Y=c()
  N=0
  while(N<n){
    y=t(rmultinom(1,size=1,prob=eta))
    a=alpha[which(y!=0),]
    x=rdirichlet(1,a)
    w=(x/m)/(sum(x/m))
    p=min(m/sum(m*w))
    u=runif(1)
    if(u<=p){N=N+1}
    P=c(P,p)
    U=c(U,u)
    X=rbind(X,x)
    W=rbind(W,w)
    Y=rbind(Y,y)
  }
  C=as.numeric((1:K)%*%t(Y))
  ind=order(W[,1])
  list(X=X[ind,],W=W[ind,],Y=Y[ind,],U=U[ind],P=P[ind],C=C[ind],I=(U[ind]<P[ind]))
}

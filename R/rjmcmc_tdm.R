#' MCMC algorithm to sample for the tDM parameters.
#'
#' @param sim simulated data: list with elements W (matrix of dimension N*K)
#' @param simu simulation number
#' @param S chain length
#' @param e0 hyperparameter
#' @param u0 hyperparameter
#' @param v0 hyperparameter
#' @param a0 hyperparameter
#' @param b0 hyperparameter
#' @param l0 hyperparameter
#' @param s0 hyperparameter
#' @param Kmin minimum number of clusters
#' @param Kmax maximum number of clusters
#' @param th thinning: save only on in "th" iteration
#' @param plt TRUE/FALSE plot the generated cluster each th^th chain
#' @param type tilted Dirichlet mixture (tdm) or unconstrained dirichelt mixture (dm)
#' @param save.path folder path to save the results
#'
#' @return
#' @export
#'
#' @examples
rjmcmc_tdm=function(sim,simu,
                    S=100,e0=2,
                    u0=50,v0=1,
                    a0=100,b0=100,
                    l0=1,s0=1,
                    Kmin=2,Kmax=10,
                    th=10,plt=FALSE,
                    type='tdm',
                    save.path=NULL){
  if(is.null(save.path)){
    message('folder path for saving the chains missing.')
  }
  beg=Sys.time()
  count=0
  W=as.matrix(sim$W)
  C=sim$C
  N=nrow(W)
  D=ncol(W)
  K=2
  par=list(eta=rep(1/K,K),
           alpha=v0+(u0-v0)*matrix(runif(K*D),ncol=D,nrow=K))

  if(type=='tdm'){
    L= Llk(W,par,N,K,D)
  }else{
    L= Llk_dm(W,par,N,K,D)
  }

  A=allocation(L,N)
  aL=augLlk(A,L,N)
  pos=sum(aL)+ddirichlet(par$eta,rep(e0,K),log=TRUE)
  C=A$C


  params=cbind(0,0,0,pos,K,cbind((par$eta),par$alpha))
  results=data.frame(params)
  colnames(results)=c('iter','move','posterior','ratio','K','eta',paste0('alpha',1:D))
  clustering=matrix(c(0,C),nrow=1)
  if(type=='tdm'){
    suff='.csv'
  }else{
    suff='_dm.csv'
  }
  write.table(results,
              file = file.path(save.path,paste0('Results',D,'_',simu,suff)),
              sep = ",", quote = FALSE,
              row.names = FALSE,
              col.names =TRUE)
  write.table(clustering,
              file = file.path(save.path,paste0('Allocations',D,'_',simu,suff)),
              sep = ",", quote = FALSE,
              row.names = FALSE,
              col.names =FALSE)


  for( s in 1:S){

    # First move --------------------------------------------------------------

    accept=0
    n=(colSums(A$Y)+1)/(N+1)
    eta1=rdirichlet(1,1+b0*n)
    m=t(sapply(1:K,function(k) colMeans(rbind(W[A$C==k,],rep(1,D)/D))))
    #alternatively:
    #m=(par$alpha-v0)/(u0-v0)
    m1=matrix(NA,ncol=D,nrow=K)
    for(k in 1:K){
      for(d in 1:D){
        m1[k,d]=rgamma(1,1+a0*m[k,d],scale=1)
      }
      m1[k,]=m1[k,]/sum(m1[k,])
    }



    par1=list(eta=as.numeric(eta1),
              alpha=v0+(u0-v0)*m1)

    P1=
      sum(sapply(1:K,function(k) ddirichlet(m1[k,],m[k,]*a0+1,log=TRUE)))+
      sum(ddirichlet(eta1,1+b0*n,log=TRUE))
    P=
      sum(sapply(1:K,function(k) ddirichlet(m[k,],m1[k,]*a0+1,log=TRUE)))+
      sum(ddirichlet(n,1+b0*eta1,log=TRUE))


    if(type=='tdm'){
      L1= Llk(W,par1,N,K,D)
    }else{
      L1= Llk_dm(W,par1,N,K,D)
    }
    A1=allocation(L1,N)
    aL1=as.numeric(augLlk(A1,L1,N))

    pos1=sum(aL1)+
      ddirichlet(par1$eta,rep(e0,K),log=TRUE)

    ratio=min(0,pos1-pos+P-P1)

    if(log(runif(1))<=ratio){
      par=par1
      A=A1
      pos=pos1
      L=L1
      accept=1
      accept_ratio=ratio
    }

    jump=runif(1,-1,1)
    if(K<=Kmin){jump=1}
    if(K>=Kmax){jump=-1}

    # # Split move --------------------------------------------------------------

    if(jump>=0){
      K1=K+1
      j=floor(runif(1,1,K1))
      nu=rbeta(1,s0,s0)
      eta1=par$eta
      eta1[j]=nu*par$eta[j]
      eta1=c(eta1,(1-nu)*par$eta[j])

      m=(par$alpha-v0)/(u0-v0)
      epsilon=c()
      P1=c()
      for(d in 1:D){
        lb=max(0,1/nu-(1-nu)/(nu*m[j,d]))
        ub=min(1/nu,1/m[j,d])
        x=rbeta(1,l0,l0)
        epsilon=c(epsilon,lb+(ub-lb)*x)
        P1=c(P1,dbeta(x,l0,l0,log=TRUE)-log(ub-lb))
      }
      m1=m
      m1[j,]=m[j,]*epsilon
      m1=rbind(m1,m[j,]*(1-nu*epsilon)/(1-nu))

      par1=list(eta=as.numeric(eta1),alpha=v0+(u0-v0)*m1)

      Jacobian=log(par$eta[j])+sum(log(m[j,]))-D*log(1-nu)
      prop=Jacobian-dbeta(nu,s0,s0,log=TRUE)-sum(P1)+log(K)

      if(type=='tdm'){
        L1= Llk(W,par1,N,K1,D)
      }else{
        L1= Llk_dm(W,par1,N,K1,D)
      }
      A1=allocation(L1,N)
      aL1=as.numeric(augLlk(A1,L1,N))

      pos1=sum(aL1)+ddirichlet(par1$eta,rep(e0,K1),log=TRUE)

      ratio=min(0,pos1-pos+prop)

      if((log(runif(1))<=ratio)*(K1==adjacent(par1,K1,D,j,u0,v0,W,A1))){
        K=K1
        par=par1
        A=A1
        pos=pos1
        L=L1
        accept=accept+2
        accept_ratio=ratio
      }
    }

    # # Merge move --------------------------------------------------------------

    if(jump<0){
      accept=-accept
      K1=K-1
      j=floor(runif(1,1,K+1))
      i=adjacent(par,K,D,j,u0,v0,W,A)

      eta1=par$eta
      eta1[j]=par$eta[j]+par$eta[i]
      eta1=eta1[-i]
      nu=par$eta[j]/(par$eta[j]+par$eta[i])

      m=(par$alpha-v0)/(u0-v0)
      m1=m
      m1[j,]=nu*m[j,]+(1-nu)*m[i,]
      m1=m1[-i,]
      epsilon=c()
      P1=c()
      for(d in 1:D){
        lb=max(0,1/nu-(1-nu)/(nu*(nu*m[j,d]+(1-nu)*m[i,d])))
        ub=min(1/nu,1/(nu*m[j,d]+(1-nu)*m[i,d]))
        e=m[j,d]/(nu*m[j,d]+(1-nu)*m[i,d])
        x=(e-lb)/(ub-lb)
        epsilon=c(epsilon,e)
        P1=c(P1,dbeta(x,l0,l0,log=TRUE)-log(ub-lb))
      }
      par1=list(eta=as.numeric(eta1),alpha=v0+(u0-v0)*m1)

      Jacobian=log(par$eta[j]+par$eta[i])+
        sum(sapply(1:D,function(d)
          log(par$eta[i]*m[i,d]+par$eta[j]*m[j,d])-log(par$eta[i])))
      prop=-sum(P1)-dbeta(nu,s0,s0,log=TRUE)-Jacobian+log(K)

      if(type=='tdm'){
        L1= Llk(W,par1,N,K1,D)
      }else{
        L1= Llk_dm(W,par1,N,K1,D)
      }
      A1=allocation(L1,N)
      aL1=as.numeric(augLlk(A1,L1,N))
      pos1=sum(aL1)+
        ddirichlet(par1$eta,rep(e0,K1),log=TRUE)

      ratio=min(0,pos1-pos+prop)

      if((log(runif(1))<=ratio)){
        K=K1
        par=par1
        A=A1
        pos=pos1
        L=L1
        accept=accept-2
        accept_ratio=ratio
      }
    }
    if(plt){
      if(accept!=0){
        print(paste0(s,'>>>',K,'>>>',length(unique(A$C)),'>>>',accept))
      }
      if(s%%100==0){
        print(s)
        pairs(W[,1:min(D,5)],col=A$C,xlim=c(0,1),ylim=c(0,1))
      }
    }
    #
    if(accept!=0){
      count=count+1
      if(type=='tdm'){
        suff='.csv'
      }else{
        suff='_dm.csv'
      }
      write.table(cbind(s,accept,pos,accept_ratio,K,cbind((par$eta),par$alpha)),
                  file = file.path(save.path,paste0('Results',D,'_',simu,suff)), sep = ",",
                  append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
      write.table(matrix(c(s,A$C),nrow=1),
                  file = file.path(save.path,paste0('Allocations',D,'_',simu,suff)), sep = ",",
                  append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }
  }
  print(paste0(count,'/',s))
  time=Sys.time()-beg
}

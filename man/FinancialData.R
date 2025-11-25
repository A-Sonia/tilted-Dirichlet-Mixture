#setwd("")
#save.path=""

devtools::load_all()

library(fGarch)
library(extraDistr)
library(evd)
library(zoo)
library(colorspace)
library(ggplot2)
library(mcclust)
library(tidyverse)
#devtools::install_github("sarawade/mcclust.ext")
library(mcclust.ext)
library(corrplot)

# Load the data  -------------------------------------------------------


data1 <- read.csv(file.path("~/Dropbox/tDM/data","data1.csv"))
data1=data1[,c('Date','S.P.Technology','S.P.Health.Care','S.P.Industrial','S.P.Energy','S.P.Financial','S.P.Materials')]
missing=apply(data1,1,function(x) any(is.na(x)))
D=ncol(data1)-1
Returns=data.frame(as.Date(data1[-1,1]))
for(d in 1:D)
{
  Returns=cbind(Returns,as.numeric(-diff(log(data1[,1+d]))))
}


colnames(Returns)=colnames(data1)


# Marginal analysis -------------------------------------------------------
N=nrow(Returns)
X0=U0=Z0=c()
par(mfrow=c(2,3),
    pty='s',mgp=c(2,1,0),
    omi=c(0,0,0,0)+0.1,
    mar=c(3,3,2,0),
    cex=1.2)
for(d in 1:D){
  fit3=  garchFit(~arma(1,0)+garch(1,1),Returns[,1+d],cond.dist='std')
  nu=coef(fit3)['shape']
  se=fit3@fit$se.coef['shape']
  plot(fit3,which=13,main='',pch=20)
  text(-2.5,4,colnames(Returns)[d+1])
  text(-2.5,3,bquote(hat(nu)==.(round(nu,2))~'('~.(round(se,2))~')'))
  X0=cbind(X0,residuals(fit3))
  U0=cbind(U0,pt(residuals(fit3,standardize=TRUE),nu))
  Z0=cbind(Z0,qfrechet(pt(residuals(fit3,standardize=TRUE),nu)))
}


R0=rowSums(Z0)
r0=quantile(R0,0.97)
W0=Z0/rowSums(Z0)
I0=(R0>=r0)
sum(I0)
pairs(W0[I0,],pch=20,xlim=c(0,1),ylim=c(0,1))
colMeans(W0[I0,])

par(mfrow=c(3,2),
    pty='m',
    mgp=c(1,0.5,0),
    omi=c(0,0,0.2,0)+0.1,mar=c(2,3,1,0),
    cex=1.2)
for( d in 1:D){
  plot(Returns[,1],log(Z0[,d]),
       col=alpha(8+I0,0.2+I0),
       type='h',xlab='',ylab='',main=colnames(Returns[1+d]))
  abline(h=0,lwd=0.4)
}


par(mfrow=c(2,1),
    pty='m',
    mgp=c(2,0.5,0),
    omi=c(0,0,0,0),
    mar=c(3,4,0.1,0.1),las=1)
plot(Returns[,1],Returns[,2],type='h',xlab='',ylab='Negative log-returns')
plot(Returns[,1],log(Z0[,2]),type='h',xlab='',ylab='GARCH-Rescaled negative log-returns (log)',col=alpha(I0+8,0.2+I0))


# Prepare angular data for the tDM model ----------------------------------


sim=list(W=W0[which(I0),],R=R0,Z=Z0,X=Returns,I=I0)
W=sim$W
dim(W)
D=ncol(W)
N=nrow(W)


# Parameters --------------------------------------------------------------


colors=c('#e31a1c','#33a02c',
                     'cyan',
                     'magenta',
                     '#ff7f00',
                     '#b2df8a',
                     'blue',
                     '#a65628',
                     '#6a3d9a',
                     '#cab2d6',
                     'yellow',
                     '#a6cee3',
                     '#fb9a99',
                     '#fdbf6f')


#Hyperparameters
e0=2
u0=20
v0=0.1
a0=1000
b0=100


l0=1
s0=1
Kmin=2
Kmax=15


N=nrow(sim$W)
D=ncol(sim$W)
thining.step=100
burnin=100
S=floor(1000*thining.step+burnin)


# Run the MCMC ------------------------------------------------------------
run=TRUE

simu=1
set.seed(1009+simu)
if(run==TRUE){
  time=rjmcmc_tdm(sim,simu,S,
                  e0,u0,v0,a0,b0,l0,s0,
                  Kmin,Kmax,
                  thining.step,
                  type='tdm',
                  plt=TRUE,
                  save.path=save.path)
  print(time)
}


# plot the results --------------------------------------------------------

Results <- read.csv(file.path(save.path,
                              paste0('Results',D,'_',simu,'.csv')))
Allocations <- read.csv(file.path(save.path,
                                  paste0('Allocations',D,'_',simu,'.csv')))

s=max(Results$iter)
thin=unique(Results$iter)
thin=thin[thin>0]
thinResults=Results[Results$iter%in%thin,]
uniqueResults=unique(thinResults[,c('iter','K','posterior')])

hist(uniqueResults$K)
plot(uniqueResults$iter,uniqueResults$posterior,col=uniqueResults$K,type='p',pch=20)


dpos=as.numeric((uniqueResults$posterior))
E=matrix(0,N,N)
cE=list()
C=list()
Cmean=matrix(0,N,N)
for(k in 1:length(thin)){
  s=thin[k]
  Ks=unique(thinResults[thinResults$iter==s,4])
  pars=list(eta=as.numeric(thinResults[thinResults$iter==s,5]),
            alpha=as.matrix(thinResults[thinResults$iter==s,5+(1:D)]))
  L=Llk(sim$W,pars,N,Ks,D)
  d_pos=dpos[uniqueResults$iter==s]
  proba=exp(L)/rowSums(exp(L))
  A=matrix(1,N,N)
  for( i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      A[i,j]=A[j,i]=(proba[i,]%*%proba[j,])
    }
  }
  E=E+A
  cE[[k]]=A
  C[[k]]=comp.psm((Allocations[k,-1]))
  Cmean=Cmean+C[[k]]
}
E=E/length(thin)
Cmean=Cmean/length(thin)
pE=E
diag(pE)=1
diag(Cmean)=1



uniqueResults$K_Star=sapply(1:length(thin),function(i) length(table(as.numeric(Allocations[i,-1]))))
uniqueResults$PEAR1=1-pear(cls=(Allocations[,-1]),psm=Cmean)
uniqueResults$PEAR=1-pear(cls=(Allocations[,-1]),psm=pE)
uniqueResults$VI1=VI.lb(cls=(Allocations[,-1]),psm=Cmean)
uniqueResults$VI=VI.lb(cls=(Allocations[,-1]),psm=pE)
uniqueResults$Hmg1=unlist(sapply(1:length(thin),function(i) sum(abs(C[[i]]-Cmean))/(2*N*(N-1))))
uniqueResults$Hmg=unlist(sapply(1:length(thin),function(i) sum(abs(C[[i]]-pE))/(2*N*(N-1))))

res=cbind(simu=simu,dim=D,uniqueResults)
choices=rbind(
  c( 'PEAR','1)Ahat',as.numeric(uniqueResults[which.min(uniqueResults$PEAR1),c('iter','K','K_Star','PEAR1')])),
  c( 'PEAR','2)Phat',as.numeric(uniqueResults[which.min(uniqueResults$PEAR),c('iter','K','K_Star','PEAR')])),

  c( 'VI','1)Ahat',as.numeric(uniqueResults[which.min(uniqueResults$VI1),c('iter','K','K_Star','VI1')])),
  c( 'VI','2)Phat',as.numeric(uniqueResults[which.min(uniqueResults$VI),c('iter','K','K_Star','VI')])),

  c( 'Binder', '1)Ahat',as.numeric(uniqueResults[which.min(uniqueResults$Hmg1),c('iter','K','K_Star','Hmg1')])),
  c( 'Binder', '2)Phat',as.numeric(uniqueResults[which.min(uniqueResults$Hmg),c('iter','K','K_Star','Hmg')]))
)

choices=cbind(simu,D,choices)
colnames(choices)=c('simu','dim','method','type','iter','K','K_Star','criterion')
choices=data.frame(choices)
choices$iter=as.numeric(  choices$iter)
choices$K=as.numeric(choices$K)
choices$criterion=as.numeric(  choices$criterion)
choices



crit=2

cl=choices[crit,'iter']
alloc=as.numeric(Allocations[which(as.numeric(Allocations[,1])==cl),-1])
colnames(sim$W)=substring(colnames(sim$X)[-1],5)
pdf(file.path(save.path,paste0('pairs_Fin_W_',simu,'_',crit,'.pdf')),12,12)
pairs(sim$W,col=colors[alloc],pch=20,lower.panel = NULL)
dev.off()

dataW=data.frame(W=sim$W,A=alloc,date=sim$X[sim$I,1])
freq=data.frame(table(dataW$A))
colnames(freq)=c('A','freq')
freq=freq[order(freq$freq,decreasing = TRUE),]
freq$rank=1:nrow(freq)
colnames(dataW)=c(colnames(sim$X)[-1],'A','date')
dataW=merge(dataW,freq)
dataW=dataW%>%pivot_longer(-c(freq,A,rank,date))%>%data.frame()



supp.labs <- paste0('Cluster num ',freq$rank,' (size=',freq$freq,')')
names(supp.labs) <- freq$rank
pdf(file.path(save.path,paste0('boxplot_fin_W_',simu,'_',crit,'.pdf')),width = 12,height = 10)
ggplot(dataW,aes(y=as.factor(name),x=value,fill=as.factor(A)))+
  geom_vline(xintercept=(0:D)/D,linetype=1,color=8,linewidth=0.1)+
  geom_vline(xintercept=1/D)+
  geom_boxplot()+
  facet_wrap(vars(rank),labeller =labeller(rank=supp.labs),ncol=3)+
  theme_bw()+theme(legend.position = 'none',text=element_text(size=15),panel.grid = element_blank())+
  scale_fill_manual(values=colors)+
  labs(x='Angular variables',y='')
dev.off()

pdf(file.path(save.path,paste0('time_fin_W_',simu,'_',crit,'.pdf')),width = 12,height = 8)
ggplot(dataW,aes(x=date,y=value,col=as.factor(A)))+
  geom_col(linewidth=1)+
  facet_wrap(vars(name),ncol=1)+
  geom_hline(yintercept = 1/D,linetype=2)+
  theme_bw()+theme(legend.position = 'none',text=element_text(size=15))+
  scale_color_manual(values=colors)+
  labs(x='Time',y='W')
dev.off()

dataX=data.frame(X=sim$Z[sim$I,],A=alloc,date=sim$X[sim$I,1])
freq=data.frame(table(dataX$A))
colnames(freq)=c('A','freq')
freq=freq[order(freq$freq,decreasing = TRUE),]
freq$rank=1:nrow(freq)

colnames(dataX)=c(colnames(sim$X)[-1],'A','date')
dataX=merge(dataX,freq)
dataX=dataX%>%pivot_longer(-c(freq,A,rank,date))%>%data.frame()

pdf(file.path(save.path,paste0('time_fin_Z_',simu,'_',crit,'.pdf')),width = 12,height = 8)
ggplot(dataX,aes(x=date,y=log(value),col=as.factor(A)))+
  geom_col(linewidth=1)+
  facet_wrap(vars(name),ncol=1)+
  theme_bw()+theme(legend.position = 'none',text=element_text(size=15))+
  scale_color_manual(values=colors)+
  geom_hline(yintercept = log(r0),linetype=2)+
  labs(x='Time',y='')
dev.off()

pdf(file.path(save.path,paste0('boxplot_fin_Z_',simu,'_',crit,'.pdf')),width = 12,height = 8)
ggplot(dataX,aes(y=as.factor(name),x=log(value),fill=as.factor(A)))+
  geom_vline(xintercept=log((1:10)*r0/10),linetype=1,linewidth=0.2,color=8)+
  geom_boxplot()+
  facet_wrap(vars(rank),labeller =labeller(rank=supp.labs),ncol=3)+
  theme_bw()+theme(legend.position = 'none',text=element_text(size=15),panel.grid  = element_blank())+
  scale_fill_manual(values=colors)+
  labs(x='Extreme rescaled returns (log)',y='')
dev.off()


colnames(pE)=rownames(pE)=as.character(sim$X[sim$I,'Date'])
H0=corrplot(pE,is.corr=FALSE,method='color',
            type = 'lower', order = 'hclust', tl.col = 'black',
            cl.ratio = 0.2, tl.srt = 0)

Hf=H0$corrPos
colnames(Hf)=c('date','datay','x','y','corr')
Hf$date=as.Date(Hf$date)
Hf=merge(Hf,dataX)

pdf(file.path(save.path,paste0('Adj_fin_',simu,'_',crit,'.pdf')),width = 12,height = 12)
ggplot(Hf,aes(x=x,y=y,fill=corr))+
  geom_tile()+
  geom_point(aes(y=max(y)-x+2,x=x,col=as.factor(A)),shape=15,size=3)+
  theme_bw()+
  guides(color='none')+
  theme(panel.grid = element_blank(),
        legend.position = 'top',
        legend.key.width = unit(2,'cm'),
        axis.text = element_blank())+
  scale_color_manual(values=colors)+
  scale_fill_binned_sequential('Grays',breaks=seq(0,1,0.1),name='')+
  labs(x='',y='')
dev.off()


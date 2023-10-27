library(MTS)
n=200
sim=300
#include xi_b1 xi_b2 and xi_s1 xi_s2 structure
rho=diag(c(0.5,0,0.5,0))
sig1=matrix(c(0.5,0,0,0.2),2,2)
sig2=matrix(c(0.2,0.15,0.15,0.1),2,2)
sig4=rbind(cbind(sig1,sig2),cbind(t(sig2),sig1))
sig=sig4-rho%*%sig4%*%rho

yt=array(0, c(n,4,sim)) 
cov_y=array(0, c(4,4,sim)) 
for (i in 1:sim){
  m1=VARMAsim(n,arlags=c(1),phi=rho,sigma=sig)
  yt[,,i]=m1$series
  cov_y[,,i]=cov(yt[,,i])
}
writeMat('xi_ybs200.mat',xi_yt=yt)

#get mean of cov_y by matlab
library(R.matlab)
Matlab$startServer()
isOpen <- open(matlab)
writeMat('xi_yt200.mat',xi_yt=yt)

setVariable(matlab, cov_ym=cov_y)
evaluate(matlab, "cov_ymean=mean(cov_ym,3)")
cov_ymeanr <- getVariable(matlab, "cov_ymean")
close(matlab)

#just include xi_b1 and xi_s1 structure
rho=0.5
p1=matrix(c(0.5,0,0,0.5),2,2)
sig=(1-rho^2)*matrix(c(0.5,0.2,0.2,0.5),2,2)

m1=VARMAsim(300,arlags=c(1),phi=p1,sigma=sig)
zt=m1$series
dim(zt)
cov(zt)

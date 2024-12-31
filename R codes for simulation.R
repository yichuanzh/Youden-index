# Constructing GP, Bayes, NB, PB, and BCB confidence intervals (CIs) for the Youden index
# and its cut-off point	in optimal cobminations of multivariate normal biomarkers with covariates
############################################################################

# Required packages to provide a parallel programming:
#       foreach
#	doParallel
#	mvtnorm
#	LaplacesDemon
#	expm
#	EnvStats
#############################################################################

library(foreach)
library(doParallel)

parallel::detectCores()

n.cores <- parallel::detectCores() - 3

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
  )

doParallel::registerDoParallel(cl = my.cluster)

foreach::getDoParRegistered()
foreach::getDoParWorkers()


M=1000       ##### Number of simulation to compute the GP, Bayes, NB, PB, and BCB bootstrap
p=2          ##### Number of biomarkers
rep=1000     ##### Number of replications to compute the coverage probability and average length 
alpha=.05    ##### Type-I error
vec=function(B){
dim=dim(B)
Bp=t(B)
B.vec=as.numeric(Bp)
B.vec
}

vec.inv=function(b,q){
matrix(b,nr=q,byrow=T)
}



z=c(1,1)     #### Covariate vector
q=length(z)+1

#Scenario
B1=matrix(c(.5,.8,1,.5,.2,1.2),nr=q)
B2=matrix(c(.9,1,1.4,.7,.6,1.8),nr=q)
Sigma1=matrix(c(1,.3,.3,1),nr=p)
Sigma2=matrix(c(1,.3,.3,1),nr=p)

n1=10; n2=10    #### sample sizes


# Function to compute the Youden Index and its cut-off point

parameter=function(B1,Sigma1,B2,Sigma2,z){
lambda=solve(Sigma1+Sigma2)%*%(t(B2)-t(B1))%*%as.matrix(c(1,t(z)))
mu1=c(t(lambda)%*%t(B1)%*%as.matrix(c(1,t(z))))
sigma1=sqrt(c(t(lambda)%*%Sigma1%*%lambda))
mu2=c(t(lambda)%*%t(B2)%*%as.matrix(c(1,t(z))))
sigma2=sqrt(c(t(lambda)%*%Sigma2%*%lambda))
mu=c(mu1,mu2)
sigma=c(sigma1,sigma2)
w2=which(mu==max(mu))
w1=which(mu==min(mu))
a=mu[w2]-mu[1]
b=sigma[w2]/sigma[w1]
if(b==1) c0=mean(mu) else c0=(mu[w1]*(b^2-1)-a+b*sqrt(a^2+(b^2-1)*sigma[w1]^2*
log(b^2)))/(b^2-1)
J=pnorm((mu[w2]-c0)/sigma[w2])+pnorm((c0-mu[w1])/sigma[w1])-1
c(c0,J)
}
par=parameter(B1,Sigma1,B2,Sigma2,z)

c0=par[1]; J=par[2]         #### Compute the true values of c0(z) and J(z)       

c0.hat=J.hat=c0.gv=J.gv=Le.c0.gv=Le.J.gv=Cov.c0.gv=Cov.J.gv=
c0.Bayes=J.Bayes=Le.c0.Bayes=Le.J.Bayes=Cov.c0.Bayes=Cov.J.Bayes=
Le.c0.nb=Le.J.nb=Cov.c0.nb=Cov.J.nb=
Le.c0.pb=Le.J.pb=Cov.c0.pb=Cov.J.pb=
Le.c0.bcb=Le.J.bcb=Cov.c0.bcb=Cov.J.bcb=c()


L=foreach(i=1:rep,.combine=rbind) %dopar% {
library(mvtnorm)
library(LaplacesDemon)
library(expm)
library(EnvStats)

# Simulation of the design matrices
Z1=matrix(c(rep(1,n1),runif(n1*(q-1),0,3)),nrow=n1)
Z2=matrix(c(rep(1,n2),runif(n2*(q-1),0,3)),nrow=n2)

# Simulation of the error matrices
E1=rmvnorm(n1, mean=rep(0,p), sigma=Sigma1)
E2=rmvnorm(n2, mean=rep(0,p), sigma=Sigma2)

# Simulation of the biomarkers results
X1=Z1%*%B1+E1
X2=Z2%*%B2+E2

# Estimation of the parameters B1, B2, Sigma1, and Sigma2
B1.hat=solve(t(Z1)%*%Z1)%*%t(Z1)%*%X1
B2.hat=solve(t(Z2)%*%Z2)%*%t(Z2)%*%X2
Q1=diag(n1)-Z1%*%solve(t(Z1)%*%Z1)%*%t(Z1)
Q2=diag(n2)-Z2%*%solve(t(Z2)%*%Z2)%*%t(Z2)
Sigma1.hat=round(1/n1*t(X1)%*%Q1%*%X1,6)
Sigma2.hat=round(1/n2*t(X2)%*%Q2%*%X2,6)

# MLE of c0(z) and J(z)
par.hat=parameter(B1.hat,Sigma1.hat,B2.hat,Sigma2.hat,z)
c0.hat[i]=par.hat[1]; J.hat[i]=par.hat[2]           ########  MLE

# Constructing the GP confidence interval 
G.c0=G.J=c()
for(j in 1:M){
W1=rwishart(n1-q, diag(p))
W2=rwishart(n2-q, diag(p))
G.Sigma1=n1*sqrtm(Sigma1.hat)%*%solve(W1)%*%sqrtm(Sigma1.hat)
G.Sigma2=n2*sqrtm(Sigma2.hat)%*%solve(W2)%*%sqrtm(Sigma2.hat)
N1.vec=c(rmvnorm(1,mean=rep(0,q*p),sigma=kronecker(solve(t(Z1)%*%Z1),G.Sigma1)))
N2.vec=c(rmvnorm(1,mean=rep(0,q*p),sigma=kronecker(solve(t(Z2)%*%Z2),G.Sigma2)))
N1=vec.inv(N1.vec,q)
N2=vec.inv(N2.vec,q)
G.B1=B1.hat-N1
G.B2=B2.hat-N2
G.par=parameter(G.B1,G.Sigma1,G.B2,G.Sigma2,z)
G.c0[j]=G.par[1]; G.J[j]=G.par[2]
}

# GP CI for c0(z)
c0.low.gv=quantile(G.c0,alpha/2,name=F); c0.up.gv=quantile(G.c0,1-alpha/2,name=F)
Le.c0.gv[i]=c0.up.gv-c0.low.gv                   ##### Length of GP CI for c0(z)
Cov.c0.gv[i]=c0.low.gv < c0 && c0 < c0.up.gv     ##### CP of GP CI for c0 (z) 

# GP CI for J(z)
J.low.gv=quantile(G.J,alpha/2,name=F); J.up.gv=quantile(G.J,1-alpha/2,name=F)
Le.J.gv[i]=J.up.gv-J.low.gv                      ##### Length of GP CI for J(z)
Cov.J.gv[i]=J.low.gv < J && J < J.up.gv          ##### CP of GP CI for J(z)


# Constructing the Bayesian credible interval
post.c0=post.J=c()

for(j in 1:M){
Sigma1.post.inverse=rwishart(n1-p-q+1,1/(n1-p-q+1)*round(solve(Sigma1.hat),6))
Sigma2.post.inverse=rwishart(n2-p-q+1,1/(n2-p-q+1)*round(solve(Sigma2.hat),6))
Sigma1.post=solve(Sigma1.post.inverse)
Sigma2.post=solve(Sigma2.post.inverse)
B1.vec=c(rmvnorm(1,mean=vec(B1.hat),sigma=kronecker(solve(t(Z1)%*%Z1),Sigma1.post)))
B2.vec=c(rmvnorm(1,mean=vec(B2.hat),sigma=kronecker(solve(t(Z2)%*%Z2),Sigma2.post)))
B1.post=vec.inv(B1.vec,q)
B2.post=vec.inv(B2.vec,q)
par.post=parameter(B1.post,Sigma1.post,B2.post,Sigma2.post,z)
post.c0[j]=par.post[1]; post.J[j]=par.post[2]
}

# Bayesian CI for c0(z)
c0.low.Bayes=quantile(post.c0,alpha/2,name=F); c0.up.Bayes=quantile(post.c0,1-alpha/2,name=F)
Le.c0.Bayes[i]=c0.up.Bayes-c0.low.Bayes               ##### Length of Bayesian CI for c0(z)
Cov.c0.Bayes[i]=c0.low.Bayes < c0 && c0 < c0.up.Bayes ##### CP of Bayesian CI for c0(z)  

# Bayesian CI for J(z)
J.low.Bayes=quantile(post.J,alpha/2,name=F); J.up.Bayes=quantile(post.J,1-alpha/2,name=F)
Le.J.Bayes[i]=J.up.Bayes-J.low.Bayes                  ##### Length of Bayesian CI for J(z)
Cov.J.Bayes[i]=J.low.Bayes < J && J < J.up.Bayes      ##### CP of Bayesian CI for J(z)

# Constructing the bootstrap confidence intervals 
boot.c0=boot.J=c()
for(j in 1:M){
W1=rwishart(n1-q, Sigma1.hat)
W2=rwishart(n2-q, Sigma2.hat)
N1.vec=c(rmvnorm(1,mean=vec(B1.hat),sigma=kronecker(solve(t(Z1)%*%Z1),Sigma1.hat)))
N2.vec=c(rmvnorm(1,mean=vec(B2.hat),sigma=kronecker(solve(t(Z2)%*%Z2),Sigma2.hat)))
N1=vec.inv(N1.vec,q)
N2=vec.inv(N2.vec,q)
boot.B1=N1
boot.B2=N2
boot.Sigma1=1/n1*W1
boot.Sigma2=1/n2*W2
boot.par=parameter(boot.B1,boot.Sigma1,boot.B2,boot.Sigma2,z)
boot.c0[j]=boot.par[1]; boot.J[j]=boot.par[2]
}
c0.bar=mean(boot.c0); J.bar=mean(boot.J)
c0.sd=sd(boot.c0); J.sd=sd(boot.J)

# NB CI for c0(z)
c0.low.nb=c0.bar-qnorm(1-alpha/2)*c0.sd; c0.up.nb=c0.bar+qnorm(1-alpha/2)*c0.sd
Le.c0.nb[i]=c0.up.nb-c0.low.nb                   ##### Length of NB CI for c0(z)
Cov.c0.nb[i]=c0.low.nb < c0 && c0 < c0.up.nb     ##### CP of NB CI for c0(z)  

# NB CI for J(z)
J.low.nb=max(J.bar-qnorm(1-alpha/2)*J.sd,0)
J.up.nb=min(J.bar+qnorm(1-alpha/2)*J.sd,1)
Le.J.nb[i]=J.up.nb-J.low.nb                     ##### Length of NB CI for J(z)
Cov.J.nb[i]=J.low.nb < J && J < J.up.nb         ##### CP of NB CI for J(z)  

# PB CI for c0(z)
c0.low.pb=quantile(boot.c0,alpha/2,name=F); c0.up.pb=quantile(boot.c0,1-alpha/2,name=F)
Le.c0.pb[i]=c0.up.pb-c0.low.pb                 ##### Length of PB CI for c0(z)
Cov.c0.pb[i]=c0.low.pb < c0 && c0 < c0.up.pb   ##### CP of PB CI for c0(z)  

# PB CI for J(z)
J.low.pb=quantile(boot.J,alpha/2,name=F); J.up.pb=quantile(boot.J,1-alpha/2,name=F)
Le.J.pb[i]=J.up.pb-J.low.pb                    ##### Length of PB CI for J(z)
Cov.J.pb[i]=J.low.pb < J && J < J.up.pb        ##### CP of PB CI for J(z)


k.c0=qnorm(pemp(c0.hat[i],boot.c0))
#BCB CI for c0(z)
c0.low.bcb=qemp(pnorm(2*k.c0+qnorm(alpha/2)),boot.c0)
c0.up.bcb=qemp(pnorm(2*k.c0+qnorm(1-alpha/2)),boot.c0)
Le.c0.bcb[i]=c0.up.bcb-c0.low.bcb                ##### Length of BCB CI for c0(z)
Cov.c0.bcb[i]=c0.low.bcb < c0 && c0 < c0.up.bcb  ##### CP of BCB CI for c0(z) 

k.J=qnorm(pemp(J.hat[i],boot.J))
#BCB CI for J(z)
J.low.bcb=qemp(pnorm(2*k.J+qnorm(alpha/2)),boot.J)
J.up.bcb=qemp(pnorm(2*k.J+qnorm(1-alpha/2)),boot.J)
Le.J.bcb[i]=J.up.bcb-J.low.bcb                ##### Length of BCB CI for J(z)
Cov.J.bcb[i]=J.low.bcb < J && J < J.up.bcb    ##### CP of BCB CI for J(z)
c(Cov.c0.gv[i],Cov.c0.Bayes[i],Cov.c0.nb[i],Cov.c0.pb[i],Cov.c0.bcb[i],
Cov.J.gv[i],Cov.J.Bayes[i],Cov.J.nb[i],Cov.J.pb[i],Cov.J.bcb[i],
Le.c0.gv[i],Le.c0.Bayes[i],Le.c0.nb[i],Le.c0.pb[i],Le.c0.bcb[i],
Le.J.gv[i],Le.J.Bayes[i],Le.J.nb[i],Le.J.pb[i],Le.J.bcb[i])
}

out=apply(L,2,mean)
CP.GP.c0=out[1]
CP.GP.J=out[6]
CP.Bayes.c0=out[2]
CP.Bayes.J=out[7]
CP.NB.c0=out[3]
CP.NB.J=out[8]
CP.PB.c0=out[4]
CP.PB.J=out[9]
CP.BcB.c0=out[5]
CP.BcB.J=out[10]

Le.GP.c0=out[11]
Le.GP.J=out[16]
Le.Bayes.c0=out[12]
Le.Bayes.J=out[17]
Le.NB.c0=out[13]
Le.NB.J=out[18]
Le.PB.c0=out[14]
Le.PB.J=out[19]
Le.BcB.c0=out[15]
Le.BcB.J=out[20]


Interval.c0=round(matrix(c(CP.GP.c0,CP.Bayes.c0,CP.NB.c0,CP.PB.c0,CP.BcB.c0,
Le.GP.c0,Le.Bayes.c0,Le.NB.c0,Le.PB.c0,Le.BcB.c0),byrow=T,nr=2),3)
dimnames(Interval.c0)=list(c("COV","LE"),c("GP","Bayes","NB","PB","BCB"))

Interval.J=round(matrix(c(CP.GP.J,CP.Bayes.J,CP.NB.J,CP.PB.J,CP.BcB.J,
Le.GP.J,Le.Bayes.J,Le.NB.J,Le.PB.J,Le.BcB.J),
byrow=T,nr=2),5)
dimnames(Interval.J)=list(c("COV","LE"),c("GP","Bayes","NB","PB","BCB"))


Interval.c0
Interval.J



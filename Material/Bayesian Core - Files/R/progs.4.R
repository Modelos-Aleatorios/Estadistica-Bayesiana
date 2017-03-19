# CHAPTER 4: R FUNCTIONS
# 04/02/2010

################################################################################

hm1=function(n,x0,sigma2)
{

# GAUSSIAN RANDOM WALK MH SAMPLER FOR THE STANDARD NORMAL DISTRIBUTION

x=rep(0,n)
x[1]=x0
for (i in 2:n)
{
y=rnorm(1,x[i-1],sqrt(sigma2))
if (runif(1)<=exp(-0.5*(y^2-x[i-1]^2))) x[i]=y else x[i]=x[i-1]
}
x
}

################################################################################

probitll=function(beta,y,X)
{

# LOG-LIKELIHOOD OF beta FOR THE PROBIT MODEL
# LOG-POSTERIOR OF beta FOR THE PROBIT MODEL UNDER FLAT PRIOR

lF1=pnorm(X%*%as.vector(beta),log.p=TRUE)
lF2=pnorm(-X%*%as.vector(beta),log.p=TRUE)
sum(y*lF1+(1-y)*lF2)

}

################################################################################

probitnoinflpost=function(beta,y,X)
{

# LOG-POSTERIOR OF beta FOR THE PROBIT MODEL UNDER NON-INFORMATIVE PRIOR

if (is.matrix(beta)==F) beta=as.matrix(t(beta))
n=dim(beta)[1]
k=dim(beta)[2]
pll=rep(0,n)
for (i in 1:n)
{
lF1=pnorm(X%*%beta[i,],log.p=TRUE)
lF2=pnorm(-X%*%beta[i,],log.p=TRUE)
pll[i]=sum(y*lF1+(1-y)*lF2)-(2*k-1)/4*log(t(beta[i,])%*%t(X)%*%X%*%beta[i,])+lgamma((2*k-1)/4)-(k/2)*log(pi)
}
pll
}

################################################################################


hmflatprobit=function(niter,y,X,scale)
{

# GAUSSIAN RANDOM WALK MH SAMPLER FOR THE PROBIT MODEL
# UNDER FLAT PRIOR

library(mnormt)
p=dim(X)[2]
mod=summary(glm(y~-1+X,family=binomial(link="probit")))
beta=matrix(0,niter,p)
beta[1,]=as.vector(mod$coeff[,1])
Sigma2=as.matrix(mod$cov.unscaled)
for (i in 2:niter)
{
tildebeta=rmnorm(1,beta[i-1,],scale*Sigma2)
llr=probitll(tildebeta,y,X)-probitll(beta[i-1,],y,X)
if (runif(1)<=exp(llr)) beta[i,]=tildebeta else beta[i,]=beta[i-1,]
}
beta
}

################################################################################

hmnoinfprobit=function(niter,y,X,scale)
{

# GAUSSIAN RANDOM WALK MH SAMPLER FOR THE PROBIT MODEL
# UNDER NON-INFORMATIVE PRIOR

library(mnormt)
p=dim(X)[2]
mod=summary(glm(y~-1+X,family=binomial(link="probit")))
beta=matrix(0,niter,p)
beta[1,]=as.vector(mod$coeff[,1])
# beta[1,]=rep(0.00001,p)
Sigma2=as.matrix(mod$cov.unscaled)
for (i in 2:niter)
{
tildebeta=rmnorm(1,beta[i-1,],scale*Sigma2)
llr=probitnoinflpost(tildebeta,y,X)-probitnoinflpost(beta[i-1,],y,X)
if (runif(1)<=exp(llr)) beta[i,]=tildebeta else beta[i,]=beta[i-1,]
}
beta
}

################################################################################

rntneg=function(n,mu,sigma2)
{

# PRODUCE n SAMPLES FROM THE SPECIFIED RIGHT-TRUNCATED TO 0
# NORMAL DISTRIBUTION

qnorm(runif(n)*pnorm(0,mu,sqrt(sigma2)),mu,sqrt(sigma2))
}

################################################################################

rntpos=function(n,mu,sigma2)
{

# PRODUCE n SAMPLES FROM THE SPECIFIED LEFT-TRUNCATED TO 0
# NORMAL DISTRIBUTION

-rntneg(n,-mu,sigma2)
}

################################################################################

gibbsprobit=function(niter,y,X)
{

# GIBBS SAMPLING FOR THE PROBIT MODEL UNDER FLAT PRIOR

library(mnormt)
p=dim(X)[2]
beta=matrix(0,niter,p)
z=rep(0,length(y))
mod=summary(glm(y~-1+X,family=binomial(link="probit")))
beta[1,]=as.vector(mod$coefficient[,1])
Sigma2=solve(t(X)%*%X)
for (i in 2:niter)
{
mean=X%*%beta[i-1,]
z[y==1]=rntpos(sum(y==1),mean[y==1],1)
z[y==0]=rntneg(sum(y==0),mean[y==0],1)
Mu=Sigma2%*%t(X)%*%z
beta[i,]=rmnorm(1,Mu,Sigma2)
}
beta
}

################################################################################

logitll=function(beta,y,X)
{

# LOG-LIKELIHOOD OF beta FOR THE LOGIT MODEL
# LOG-POSTERIOR OF beta FOR THE LOGIT MODEL UNDER NON-INFORMATIVE PRIOR

lF1=plogis(X%*%as.vector(beta),log.p=TRUE)
lF2=plogis(-X%*%as.vector(beta),log.p=TRUE)
sum(y*lF1+(1-y)*lF2)
}

################################################################################

logitnoinflpost=function(beta,y,X)
{

# LOG-POSTERIOR OF beta FOR THE LOGIT MODEL
# UNDER NON-INFORMATIVE PRIOR

if (is.matrix(beta)==F) beta=as.matrix(t(beta))
n=dim(beta)[1]
k=dim(beta)[2]
pll=rep(0,n)
for (i in 1:n)
{
lF1=plogis(X%*%beta[i,],log.p=TRUE)
lF2=plogis(-X%*%beta[i,],log.p=TRUE)
pll[i]=sum(y*lF1+(1-y)*lF2)-(2*k-1)/4*log(t(beta[i,])%*%t(X)%*%X%*%beta[i,]/2)+lgamma((2*k-1)/4)-(k/2)*log(2*pi)
}
pll
}

################################################################################

hmflatlogit=function(niter,y,X,scale)
{

# GAUSSIAN RANDOM WALK MH SAMPLER FOR THE LOGIT MODEL
# UNDER FLAT PRIOR

library(mnormt)
p=dim(X)[2]
mod=summary(glm(y~-1+X,family=binomial(link="logit")))
beta=matrix(0,niter,p)
beta[1,]=as.vector(mod$coeff[,1])
Sigma2=as.matrix(mod$cov.unscaled)
for (i in 2:niter)
{
tildebeta=rmnorm(1,beta[i-1,],scale*Sigma2)
llr=logitll(tildebeta,y,X)-logitll(beta[i-1,],y,X)
if (runif(1)<=exp(llr)) beta[i,]=tildebeta else beta[i,]=beta[i-1,]
}
beta
}

################################################################################

hmnoinflogit=function(niter,y,X,scale)
{

# GAUSSIAN RANDOM WALK MH SAMPLER FOR THE LOGIT MODEL
# UNDER NON-INFORMATIVE PRIOR

library(mnormt)
p=dim(X)[2]
mod=summary(glm(y~-1+X,family=binomial(link="logit")))
beta=matrix(0,niter,p)
beta[1,]=as.vector(mod$coeff[,1])
Sigma2=as.matrix(mod$cov.unscaled)
for (i in 2:niter)
{
tildebeta=rmnorm(1,beta[i-1,],scale*Sigma2)
llr=logitnoinflpost(tildebeta,y,X)-logitnoinflpost(beta[i-1,],y,X)
if (runif(1)<=exp(llr)) beta[i,]=tildebeta else beta[i,]=beta[i-1,]
}
beta
}

################################################################################

loglinll=function(beta,y,X)
{

# LOG-LIKELIHOOD OF beta FOR THE LOG-LINEAR MODEL
# LOG-POSTERIOR OF beta FOR THE LOG-LINEAR MODEL UNDER FLAT PRIOR

if (is.matrix(beta)==F) beta=as.matrix(t(beta))
n=dim(beta)[1]
pll=rep(0,n)
for (i in 1:n)
{
lF=exp(X%*%beta[i,])
pll[i]=sum(dpois(y,lF,log=TRUE))
}
pll
}

################################################################################

loglinnoinflpost=function(beta,y,X)
{

# LOG-POSTERIOR OF beta FOR THE LOG-LINEAR MODEL UNDER NON-INFORMATIVE PRIOR

if (is.matrix(beta)==F) beta=as.matrix(t(beta))
n=dim(beta)[1]
k=dim(beta)[2]
pll=rep(0,n)
for (i in 1:n)
{
lF=exp(X%*%beta[i,])
pll[i]=sum(dpois(y,lF,log=TRUE))-k/2*log(t(beta[i,])%*%t(X)%*%X%*%beta[i,]/2)+lgamma(k/2)-(k/2)*log(2*pi)
}
pll
}

################################################################################

hmflatloglin=function(niter,y,X,scale)
{

# GAUSSIAN RANDOM WALK MH SAMPLER FOR THE LOG-LINEAR MODEL
# UNDER FLAT PRIOR

library(mnormt)
p=dim(X)[2]
mod=summary(glm(y~-1+X,family=poisson()))
beta=matrix(0,niter,p)
beta[1,]=as.vector(mod$coeff[,1])
Sigma2=as.matrix(mod$cov.unscaled)
for (i in 2:niter)
{
tildebeta=rmnorm(1,beta[i-1,],scale*Sigma2)
llr=loglinll(tildebeta,y,X)-loglinll(beta[i-1,],y,X)
if (runif(1)<=exp(llr)) beta[i,]=tildebeta else beta[i,]=beta[i-1,]
}
beta
}

################################################################################

hmnoinfloglin=function(niter,y,X,scale)
{

# GAUSSIAN RANDOM WALK MH SAMPLER FOR THE LOG-LINEAR MODEL
# UNDER NON-INFORMATIVE PRIOR

library(mnormt)
p=dim(X)[2]
mod=summary(glm(y~-1+X,family=poisson()))
beta=matrix(0,niter,p)
beta[1,]=as.vector(mod$coeff[,1])
Sigma2=as.matrix(mod$cov.unscaled)
for (i in 2:niter)
{
tildebeta=rmnorm(1,beta[i-1,],scale*Sigma2)
llr=loglinnoinflpost(tildebeta,y,X)-loglinnoinflpost(beta[i-1,],y,X)
if (runif(1)<=exp(llr)) beta[i,]=tildebeta else beta[i,]=beta[i-1,]
}
beta
}

################################################################################

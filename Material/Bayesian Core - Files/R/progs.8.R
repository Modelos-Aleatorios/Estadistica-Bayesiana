# CHAPTER 8: R FUNCTIONS
# 08/01/2007

# IMAGE SEGMENTATION PART

###################################################################################

xvois4=function(x,a,b,col)
{
n=dim(x)[1]
if (a!=n & b!=n) lulu=c(x[a-1,b]==col,x[a+1,b]==col,x[a,b-1]==col,x[a,b+1]==col)
if (a==n & b!=n) lulu=c(x[a-1,b]==col,x[a,b-1]==col,x[a,b+1]==col)
if (a!=n & b==n) lulu=c(x[a-1,b]==col,x[a+1,b]==col,x[a,b-1]==col) 
if (a==n & b==n) lulu=c(x[a-1,b]==col,x[a,b-1]==col) 
sum(lulu)
}

###################################################################################

isingcondloglike=function(beta,x,a,b)
{
n0=xvois4(x,a,b,0)
n1=4-n0
beta*(n0*(x[a,b]==0)+n1*(x[a,b]==1))-log(exp(beta*n0)+exp(beta*n1))
}

###################################################################################

isinglogpseudolike=function(beta,x)
{
lulu=0
n=dim(x)[1]
for (i in 1:n)
{
for (j in 1:n)
{
n0=xvois4(x,i,j,0)
n1=4-n0
lulu=lulu+beta*(n0*(x[i,j]==0)+n1*(x[i,j]==1))-log(exp(beta*n0)+exp(beta*n1))
}
}
lulu
}

###################################################################################

pottscondloglike=function(xtilde,beta,x,a,b)
{
n=rep(0,4)
n[1]=xvois4(x,a,b,1)
n[2]=xvois4(x,a,b,2)
n[3]=xvois4(x,a,b,3)
n[4]=4-n[1]-n[2]-n[3]
beta*n[xtilde]-log(sum(exp(beta*n)))
}

###################################################################################

pottslogpseudolike=function(beta,x)
{
lulu=0
numb=dim(x)[1]
for (i in 1:numb)
{
for (j in 1:numb)
{
n=rep(0,4)
n[1]=xvois4(x,a,b,1)
n[2]=xvois4(x,a,b,2)
n[3]=xvois4(x,a,b,3)
n[4]=4-n[1]-n[2]-n[3]
lulu=lulu+beta*n[x[a,b]]-log(sum(exp(beta*n)))
}
}
lulu
}

###################################################################################

isinggibbs=function(niter,n,beta)
{
x=sample(c(0,1),n^2,prob=c(0.5,0.5),rep=TRUE)
x=matrix(x,n,n)
for (i in 1:niter)
{
echan1=sample(1:n,n,prob=rep(1,n)/n)
echan2=sample(1:n,n,prob=rep(1,n)/n)
for (k in 1:n)
{
for (l in 1:n)
{
n0=xvois4(x,echan1[k],echan2[l],0)
n1=4-n0
x[echan1[k],echan2[l]]=sample(c(0,1),1,prob=c(exp(beta*n0),exp(beta*n1)))
}
}
print(i)
}
x
}

###################################################################################

isinghm=function(niter,n,beta)
{
x=sample(c(0,1),n^2,prob=c(0.5,0.5),rep=TRUE)
x=matrix(x,n,n)
for (i in 1:niter)
{
echan1=sample(1:n,n,prob=rep(1,n)/n)
echan2=sample(1:n,n,prob=rep(1,n)/n)
for (k in 1:n)
{
for (l in 1:n)
{
laccept=1/exp(isingcondloglike(beta,x,echan1[k],echan2[l]))-1
if (runif(1)<=laccept) x[echan1[k],echan2[l]]=1-x[echan1[k],echan2[l]]
}
}
print(i)
}
x
}

###################################################################################

pottsgibbs=function(niter,numb,beta)
{
x=sample(1:4,numb^2,prob=rep(1,4),rep=TRUE)
x=matrix(x,numb,numb)
for (i in 1:niter)
{
echan1=sample(1:numb,numb,prob=rep(1,numb)/numb)
echan2=sample(1:numb,numb,prob=rep(1,numb)/numb)
for (k in 1:numb)
{
for (l in 1:numb)
{
n=rep(0,4)
n[1]=xvois4(x,echan1[k],echan2[l],1)
n[2]=xvois4(x,echan1[k],echan2[l],2)
n[3]=xvois4(x,echan1[k],echan2[l],3)
n[4]=4-n[1]-n[2]-n[3]
x[echan1[k],echan2[l]]=sample(1:4,1,prob=exp(beta*n))
}
}
print(i)
}
x
}

###################################################################################

pottshm=function(niter,numb,beta)
{
x=sample(1:4,numb^2,prob=rep(1,4),rep=TRUE)
x=matrix(x,numb,numb)
for (i in 1:niter)
{
echan1=sample(1:numb,numb,prob=rep(1,numb)/numb)
echan2=sample(1:numb,numb,prob=rep(1,numb)/numb)
for (k in 1:numb)
{
for (l in 1:numb)
{
n=rep(0,4)
n[1]=xvois4(x,echan1[k],echan2[l],1)
n[2]=xvois4(x,echan1[k],echan2[l],2)
n[3]=xvois4(x,echan1[k],echan2[l],3)
n[4]=4-n[1]-n[2]-n[3]
lulu=1:4
lulu=lulu[-x[echan1[k],echan2[l]]]
xtilde=sample(lulu,1,rep(1,3))
laccept=pottscondloglike(xtilde,beta,x,echan1[k],echan2[l])-pottscondloglike(x[echan1[k],echan2[l]],beta,x,echan1[k],echan2[l])
if (runif(1)<=exp(laccept)) x[echan1[k],echan2[l]]=xtilde
}
}
print(i)
}
x
}

###################################################################################

pottshm6=function(niter,numb,beta)
{
x=sample(1:6,numb^2,prob=rep(1,6),rep=TRUE)
x=matrix(x,numb,numb)
s=rep(0,niter)
for (i in 1:niter)
{
echan1=sample(1:numb,numb,prob=rep(1,numb)/numb)
echan2=sample(1:numb,numb,prob=rep(1,numb)/numb)
tutu=0
for (k in 1:numb)
{
for (l in 1:numb)
{
n=rep(0,6)
n[1]=xvois4(x,echan1[k],echan2[l],1)
n[2]=xvois4(x,echan1[k],echan2[l],2)
n[3]=xvois4(x,echan1[k],echan2[l],3)
n[4]=xvois4(x,echan1[k],echan2[l],4)
n[5]=xvois4(x,echan1[k],echan2[l],5)
n[6]=4-n[1]-n[2]-n[3]-n[4]-n[5]
lulu=1:6
lulu=lulu[-x[echan1[k],echan2[l]]]
xtilde=sample(lulu,1,rep(1,5))
laccept=beta*n[xtilde]-beta*n[x[echan1[k],echan2[l]]]
if (runif(1)<=exp(laccept)) x[echan1[k],echan2[l]]=xtilde
tutu=tutu+n[x[echan1[k],echan2[l]]]
}
}
s[i]=tutu
}
s
}

###################################################################################

reconstruct1=function(niter,y,m)
{
numb=dim(y)[1]
x=matrix(0,numb,numb)

mu=matrix(0,niter,6)
sigma2=rep(0,niter)

mu[1,1]=35
mu[1,2]=50
mu[1,3]=65
mu[1,4]=84
mu[1,5]=92
mu[1,6]=120
sigma2[1]=100

beta=rep(1,niter)

xcum=matrix(0,numb^2,6)
n=rep(0,6)

for (i in 2:niter)
{

lvr=0

for (k in 1:numb)
{
for (l in 1:numb)
{
n[1]=xvois4(x,k,l,1)
n[2]=xvois4(x,k,l,2)
n[3]=xvois4(x,k,l,3)
n[4]=xvois4(x,k,l,4)
n[5]=xvois4(x,k,l,5)
n[6]=4-n[1]-n[2]-n[3]-n[4]-n[5]

x[k,l]=sample(1:6,1,prob=exp(beta[i-1]*n)*dnorm(y[k,l],mu[i-1,],sqrt(sigma2[i-1])))
xcum[(k-1)*100+l,x[k,l]]=xcum[(k-1)*100+l,x[k,l]]+1

lvr=lvr+n[x[k,l]]
}
}

mu[i,1]=truncnorm(1,mean(y[x==1]),sqrt(sigma2[i-1]/sum(x==1)),0,mu[i-1,2])
mu[i,2]=truncnorm(1,mean(y[x==2]),sqrt(sigma2[i-1]/sum(x==2)),mu[i,1],mu[i-1,3])
mu[i,3]=truncnorm(1,mean(y[x==3]),sqrt(sigma2[i-1]/sum(x==3)),mu[i,2],mu[i-1,4])
mu[i,4]=truncnorm(1,mean(y[x==4]),sqrt(sigma2[i-1]/sum(x==4)),mu[i,3],mu[i-1,5])
mu[i,5]=truncnorm(1,mean(y[x==5]),sqrt(sigma2[i-1]/sum(x==4)),mu[i,4],mu[i-1,6])
mu[i,6]=truncnorm(1,mean(y[x==6]),sqrt(sigma2[i-1]/sum(x==5)),mu[i,5],255)

sese=sum((y-mu[i,1])^2*(x==1)+(y-mu[i,2])^2*(x==2)+(y-mu[i,3])^2*(x==3)+(y-mu[i,4])^2*(x==4)+(y-mu[i,5])^2*(x==5)+(y-mu[i,6])^2*(x==6))
sigma2[i]=1/rgamma(1,numb^2/2,sese/2)

betatilde=beta[i-1]+runif(1,-0.05,0.05)

# iitilde=1+10*trunc(betatilde)+(betatilde-trunc(betatilde))*10
# ii=1+10*trunc(beta[i-1])+(beta[i-1]-trunc(beta[i-1]))*10

laccept=lvr*(betatilde-beta[i-1])+integrate(lrcst,betatilde,beta[i-1])$value
if (runif(1)<=exp(laccept)) beta[i]=betatilde else beta[i]=beta[i-1]

print(round(c(i,beta[i],mu[i,],sigma2[i]),2))
}

list(beta=beta,mu=mu,sigma2=sigma2,xcum=xcum)
}

###################################################################################

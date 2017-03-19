# Estimacion of a mixture of 3 normal distributions with
# aLL parameters unknown

# A mixture will be represented as a list
# with components k, mu, sig and p

#---------------------------------Main-prgrM-----------------------------------#
# The dataset is called datha
n=length(datha)

meand=mean(datha)
vard=var(datha)
sdd=sqrt(vard)

#------------------------------------------------------------------------------#
# Simulation of n independent random variables
# from a Dirichlet distribution with parameters par
rdirichlet=function(n=1,par)
{

k=length(par)
mat=matrix(0,n,k)
for (i in 1:n) 
{
sim=rgamma(k,shape=par,scale=1)
mat[i,]=sim/sum(sim)
}
mat
}
#------------------------------------------------------------------------------#

#Initialising
k=3
mu=rnorm(k,mean=meand,sd=sdd/3)
sig=1/rgamma(k,shape=10,scale=vard)
p=rdirichlet(par=rep(5.5,k))
mix=list(k=k,p=p,mu=mu,sig=sig)

#------------------------------------------------------------------------------#
# Dirichlet log-prior (1/2,...,1/2)
Dirichlet=function(p){

k=length(p)
ou=-0.5*sum(log(p))+lgamma(k*.5)-k*lgamma(.5)

ou
}

#------------------------------------------------------------------------------#
# Overall empirical type I log-prior
prior=function(mix,meanp=meand,sdp=sdd,varp=vard){

ou=Dirichlet(mix$p)
for (i in 1:mix$k)
   ou=ou+dnorm(mix$mu[i],mean=meanp,sd=sdp,log=T)+
         dgamma(1/mix$sig[i],shape=3,rate=varp,log=T)

ou
}

#------------------------------------------------------------------------------#
# log-Likelihood
like=function(mix){

ou=0

if (n>0){
  for (i in 1:n)
     ou=ou+log(sum(mix$p*dnorm(datha[i],mean=mix$mu,sd=sqrt(mix$sig))))
  }

ou
}

#------------------------------------------------------------------------------#
# Log-posterior distribution 
lpost1=function(mix,meanp=meand,sigp=sdd,varp=vard)
{

like(mix)+prior(mix,meanp,sigp,varp)
}

#------------------------------------------------------------------------------#
# GIBBS sampling
# niter: number of iterations
gibbs=function(niter,max,meanp=meand,sigp=sdd,varp=vard)
{

mix=max

z=rep(0,mix$k)

ssiz=rep(0,mix$k)
nxx=ssiz
ssum=ssiz

mug=matrix(0,nrow=niter,ncol=mix$k)
sigg=mug
prog=mug

for (i in 1:niter)
{
  # Allocs

  for (t in 1:n){

    prob=mix$p*dnorm(datha[t],mean=mix$mu,sd=sqrt(mix$sig))
    z[t]=sample(x=1:mix$k,size=1,prob=prob)	#order(prob)[mix$k]
    }

  # Parameta
  for (j in 1:mix$k){

    ssiz[j]=sum(z==j)
    nxj[j]=sum(as.numeric(z==j)*datha)
    }

  mug[i,]=rnorm(mix$k,mean=((3*meanp/varp)+nxj/mix$sig)/((3/varp)+(ssiz/mix$sig)),sd=1/sqrt((3/varp)+(ssiz/mix$sig)))
  mix$mu=mug[i,]

  for (j in 1:mix$k)
    ssum[j]=sum(as.numeric(z==j)*(datha-mix$mu[j])^2)

  sigg[i,]=1/rgamma(mix$k,shape=10+0.5*ssiz,rate=varp+0.5*(mix$mu-meanp)^2+0.5*ssum)
  mix$sig=sigg[i,]

  prog[i,]=rdirichlet(1,par=ssiz+0.5)
  mix$p=prog[i,]

  }

list(k=k,mu=mug,sig=sigg,p=prog)
}


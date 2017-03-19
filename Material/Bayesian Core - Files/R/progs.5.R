# CHAPTER 5: R FUNCTIONS
# 21/12/2006

################################################################################

pbino=function(nplus)
{

# POSTERIOR PROBABILITIES OF N FOR THE BINOMIAL CAPTURE MODEL
# UNDER UNIFORM PRIOR (N_MAX=400)

prob=c(rep(0,max(nplus,1)-1),1/(max(nplus,1):400+1))
prob/sum(prob)
}

################################################################################

pcapture=function(T,nplus,nc)
{

# POSTERIOR PROBABILITIES OF N FOR THE T-STAGE CAPTURE-RECAPTURE MODEL
# UNDER UNIFORM PRIOR (N_Max=400)

lprob=lchoose(max(nplus,1):400,nplus)+lgamma(T*max(nplus,1):400-nc+1)-lgamma(T*max(nplus,1):400+2)
prob=c(rep(0,max(nplus,1)-1),exp(lprob-max(lprob)))
prob/sum(prob)
}

################################################################################

pdarroch=function(n1,n2,m2)
{

# POSTERIOR PROBABILITIES OF N FOR THE DARROCH MODEL
# UNDER UNIFORM PRIOR (N_Max=400)

prob=c(rep(0,max(n1+n2-m2,1)-1),choose(n1,m2)*choose(max((n1+n2-m2),1):400-n1,n2-m2)/choose(max((n1+n2-m2),1):400,n2))
prob/sum(prob)
}

################################################################################

gibbs1=function(nsimu,T,nplus,nc,lambda)
{

# GIBBS SAMPLING FOR THE T-STAGE CAPTURE-RECAPTURE MODEL
# UNDER POISSON PRIOR

N=rep(0,nsimu)
N[1]=lambda
p=rep(0,nsimu)
p[1]=rbeta(1,nc+1,T*lambda-nc+1)
for (i in 2:nsimu) 
{
N[i]=nplus+rpois(1,lambda*(1-p[i-1])^T)
p[i]=rbeta(1,nc+1,T*N[i]-nc+1)
}
list(N=N,p=p)
}

################################################################################

gibbs2=function(nsimu,n1,c2,c3,N0,r10,r20)
{

# GIBBS SAMPLING FOR THE 2-STAGE OPEN POPULATION MODEL

N=rep(0,nsimu)
p=rep(0,nsimu)
q=rep(0,nsimu)
r1=rep(0,nsimu)
r2=rep(0,nsimu)
N[1]=N0
r1[1]=r10
r2[1]=r20
nplus=n1+c2+c3
for (i in 2:nsimu) 
{
uplus=N[i-1]-r1[i-1]-c2+n1-r1[i-1]-r2[i-1]-c3
p[i]=rbeta(1,nplus+1,uplus+1)
q[i]=rbeta(1,r1[i-1]+r2[i-1]+1,2*n1-2*r1[i-1]-r2[i-1]+1)
N[i]=n1+rnbinom(1,n1,p[i])
mm=min(n1-r2[i-1]-c3,n1-c2)
pq=q[i]/(1+(1-q[i])*(1-p[i]))
pr=lchoose(n1-c2,0:mm)+(0:mm)*log(pq)+lchoose(n1-(0:mm),r2[i-1]+c3)
r1[i]=sample(0:mm,1,prob=exp(pr-max(pr)))
r2[i]=rbinom(1,n1-r1[i]-c3,q[i]/(1+(1-q[i])*(1-p[i])))
}
list(N=N,p=p,q=q,r1=r1,r2=r2)
}

################################################################################

seuil=function(k,n1,c2,c3,r2,q1) 
{

# ACCEPT-REJECT BOUND

choose(n1-c2,k)*0.9^k*choose(n1-k,c3+1)/(choose(n1,k)*choose(n1,c3+1))
} 

################################################################################

ardipper <- function(nsimu,n1,c2,c3,r2,q1)
{

# ACCEPT-REJECT ALGORITHM FOR THE OPEN POPULATION MODEL

echan=1:nsimu
for (i in 1:nsimu)
{
test=TRUE
while (test==TRUE) {
y=rbinom(1,n1,q1)
if (runif(1) <= seuil(y,n1,c2,c3,r2,q1))
{
test=FALSE
echan[i]=y
}
}
}
echan
}

################################################################################

rdirichlet=function(n,par)
{

# PRODUCE n SAMPLES FROM THE SPECIFIED DIRICHLET DISTRIBUTION

k=length(par)
mat=matrix(0,n,k)
for (i in 1:n) 
{
sim=rgamma(k,par,1)
mat[i,]=sim/sum(sim)
}
mat
}

################################################################################

gibbs3=function(nsimu,x,z)
{

# GIBBS SAMPLING FOR THE ARNASON-SCHWARZ CAPTURE-RECAPTURE MODEL

m=max(z)
T=dim(z)[2]
n=dim(z)[1]
p=array(0,c(nsimu,m))
phi=array(0,c(nsimu,m))
psi=array(0,c(m,m,nsimu))
latent=z
for (i in 1:n)
{
for (t in 1:T)
{
if (x[i,t]==0 & sum(x[i,t:T])!=0) latent[i,t]=sample(1:m,1,prob=rep(1,m))
if (x[i,t]==0 & sum(x[i,t:T])==0) latent[i,t]=sample(1:(m+1),1,prob=c(rep(1,m),m))
if (t!=1) if (latent[i,t-1]==m+1) latent[i,t]=m+1
}
}
omega=rep(0,m+1)
for (s in 2:nsimu)
{
for (r1 in 1:m) { for (r2 in 1:(m+1))
{
omega[r2]=sum(latent[,1:(T-1)]==r1 & latent[,2:T]==r2)
}
u=sum(x==1 & latent==r1)
v=sum(x==0 & latent==r1)
p[s,r1]=rbeta(1,1+u,1+v)
phi[s,r1]=rbeta(1,1+sum(omega[1:m]),1+omega[m+1])
psi[r1,,s]=rdirichlet(1,rep(1,m)+omega[1:m])
}
tt=matrix(rep(phi[s,],m),m,byrow=T)
q=rbind(tt*psi[,,s],rep(0,m))
q=cbind(q,1-apply(q,1,sum))
for (i in 1:n)
{
if (x[i,1]==0) latent[i,1]=sample(1:(m+1),1,prob=q[,latent[i,2]]*(1-c(p[s,],0)))
for (t in 2:(T-1))
{
if (x[i,t]==0) latent[i,t]=sample(1:(m+1),1,prob=q[latent[i,t-1],]*q[,latent[i,t+1]]*(1-c(p[s,],0)))
}
if (x[i,T]==0) latent[i,T]=sample(1:(m+1),1,prob=q[latent[i,T-1],]*(1-c(p[s,],0)))
}
}
list(p=p,phi=phi,psi=psi)
}

################################################################################

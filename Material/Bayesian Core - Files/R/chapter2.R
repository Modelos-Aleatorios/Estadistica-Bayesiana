# CHAPTER 2: R FUNCTIONS
# 25/07/2006

################################################################################

inv=function(X)
{

# RETURN THE INVERSE OF THE SYMMETRIC MATRIX X

EV=eigen(X)
EV$vector%*%diag(1/EV$values)%*%t(EV$vector)
}


################################################################################

t1=function(gam)
{

# RETURN THE INDICES OF THE BINARY VECTOR gam WHICH ARE DIFFERENT FROM 0

p=length(gam)
if (sum(gam)!=0) order(gam)[(p-sum(gam)+1):p] else 0
}

################################################################################

invt1=function(ind,p)
{

# CREATE A VECTOR OF DIMENSION p WITH ELEMENTS ind
# EQUAL TO 1 AND THE OTHERS EQUAL TO 0

gam=rep(0,p)
if (sum(ind)!=0) gam[ind]=1
gam
}

################################################################################

newt1=function(gam)
{

# RETURN THE INDICES OF THE BINARY VECTOR gam WHICH ARE DIFFERENT FROM 0
# IN CHARACTER MODE

qgam=sum(gam)
if (qgam!=0) 
{
p=length(gam)
cha=paste(order(gam)[(p-sum(gam)+1):p])
chacha=cha[1]
if (qgam!=1)
{
for (k in 2:length(cha)) chacha=paste(chacha,cha[k],sep=",")
}
}
else chacha="0"
chacha
}

################################################################################

gaga=function(GAM,gam)
{

# RETURN THE ROW INDEX OF THE MATRIX GAM WHICH IS EQUAL TO THE VECTOR gam

bool=T
t=1
while (bool==T)
{
if (sum(GAM[t,]==gam)==length(gam)) bool=F else t=t+1
}
t
}

################################################################################

lpostw=function(gam,y,X,betatilde,c)
{

# RETURN THE POSTERIOR PROBABILITY OF MODEL gam UP TO A NORMALIZING CONSTANT
# UNDER ZELLNER INFORMATIVE G-PRIOR

library(MASS)
n=length(y)
qgam=sum(gam)
t1gam=t1(gam)
Xt1=cbind(rep(1,n),X[,t1gam])
if (qgam!=0) P1=Xt1%*%inv(t(Xt1)%*%Xt1)%*%t(Xt1) else P1=matrix(0,n,n)
-(qgam+1)/2*log(c+1)-n/2*log(t(y)%*%y-c/(c+1)*t(y)%*%P1%*%y-1/(c+1)*t(betatilde)%*%t(cbind(rep(1,n),X))%*%P1%*%cbind(rep(1,n),X)%*%betatilde)
}

#################################################################################

lpostwnoinf=function(gam,y,X)
{

# RETURN THE POSTERIOR PROBABILITY OF MODEL gam UP TO A NORMALIZING CONSTANT
# UNDER ZELLNER NON-INFORMATIVE G-PRIOR

cc=1:100000
library(MASS)
n=length(y)
qgam=sum(gam)
t1gam=t1(gam)
Xt1=cbind(rep(1,n),X[,t1gam])
if (qgam!=0) P1=Xt1%*%inv(t(Xt1)%*%Xt1)%*%t(Xt1) else P1=matrix(0,n,n)
log(sum(cc^(-1)*(cc+1)^(-(qgam+1)/2)*(t(y)%*%y-cc/(cc+1)*t(y)%*%P1%*%y)^(-n/2)))
}

################################################################################

gibbsIP=function(niter,y,X,betatilde,c)
{

# GIBBS SAMPLING UNDER ZELLNER INFORMATIVE G-PRIOR

lga=length(betatilde)-1
gamma=matrix(0,niter,lga)
gamma[1,]=sample(c(0,1),lga,rep=T)
for (i in 1:(niter-1))
{
gamma[i+1,]=gamma[i,]
for (j in 1:lga)
{
if (j==1) {gam0=c(0,gamma[i+1,2:lga]);gam1=c(1,gamma[i+1,2:lga])}
if (j==lga) {gam0=c(gamma[i+1,1:(lga-1)],0);gam1=c(gamma[i+1,1:(lga-1)],1)}
if (j!=1 & j!=lga) {gam0=c(gamma[i+1,1:(j-1)],0,gamma[i+1,(j+1):lga]);gam1=c(gamma[i+1,1:(j-1)],1,gamma[i+1,(j+1):lga])}
pr0=lpostw(gam0,y,X,betatilde,c)
pr1=lpostw(gam1,y,X,betatilde,c)
prob0=exp(pr0-max(pr0,pr1))
prob1=exp(pr1-max(pr0,pr1))
s=prob0+prob1
prob0=prob0/s
prob1=prob1/s
gamma[i+1,j]=sample(c(0,1),1,prob=c(prob0,prob1))
}
print(i)
}
gamma
}

################################################################################

gibbsNIP1=function(niter,y,X)
{

# GIBBS SAMPLING UNDER ZELLNER NON-INFORMATIVE G-PRIOR

lga=dim(X)[2]

lpow=function(gam)
{
lpostwnoinf(gam,y,X)
}

gamma=matrix(0,niter,lga)
gamma[1,]=sample(c(0,1),lga,rep=T)
for (i in 1:(niter-1))
{
gamma[i+1,]=gamma[i,]
for (j in 1:lga)
{
if (j==1) {gam0=c(0,gamma[i+1,2:lga]);gam1=c(1,gamma[i+1,2:lga])}
if (j==lga) {gam0=c(gamma[i+1,1:(lga-1)],0);gam1=c(gamma[i+1,1:(lga-1)],1)}
if (j!=1 & j!=lga) {gam0=c(gamma[i+1,1:(j-1)],0,gamma[i+1,(j+1):lga]);gam1=c(gamma[i+1,1:(j-1)],1,gamma[i+1,(j+1):lga])}
pr0=lpow(gam0)
pr1=lpow(gam1)
prob0=exp(pr0-max(pr0,pr1))
prob1=exp(pr1-max(pr0,pr1))
s=prob0+prob1
prob0=prob0/s
prob1=prob1/s
gamma[i+1,j]=sample(c(0,1),1,prob=c(prob0,prob1))
}
print(i)
}
gamma
}

################################################################################

gibbsNIP2=function(niter,y,X,allmodnoinf)
{

# GIBBS SAMPLING UNDER ZELLNER NON-INFORMATIVE G-PRIOR
# VERSION TO USE WHEN THE MODEL POSTERIOR PROBABILITIES HAVE
# BEEN PREVIOUSLY CALCULATED (allmodnoinf)

library(combinat)

lga=log(length(allmodnoinf))/log(2)
GAM=matrix(0,2^lga,lga)
qq=cumsum(c(choose(lga,0:lga)))
it1=function(ind,p=lga)
{
invt1(ind,p)
}
for (k in 1:lga)
{
tab=t(combn(lga,k))
GAM[(qq[k]+1):qq[k+1],]=t(apply(tab,1,it1))
}

pow=function(gam)
{
tre=order(apply(GAM==t(matrix(gam,lga,2^lga)),1,sum))[2^lga]
allmodnoinf[tre]
}


gamma=matrix(0,niter,lga)
gamma[1,]=sample(c(0,1),lga,rep=T)
for (i in 1:(niter-1))
{
gamma[i+1,]=gamma[i,]
for (j in 1:lga)
{
if (j==1) {gam0=c(0,gamma[i+1,2:lga]);gam1=c(1,gamma[i+1,2:lga])}
if (j==lga) {gam0=c(gamma[i+1,1:(lga-1)],0);gam1=c(gamma[i+1,1:(lga-1)],1)}
if (j!=1 & j!=lga) {gam0=c(gamma[i+1,1:(j-1)],0,gamma[i+1,(j+1):lga]);gam1=c(gamma[i+1,1:(j-1)],1,gamma[i+1,(j+1):lga])}
pr0=log(pow(gam0))
pr1=log(pow(gam1))
prob0=exp(pr0-max(pr0,pr1))
prob1=exp(pr1-max(pr0,pr1))
s=prob0+prob1
prob0=prob0/s
prob1=prob1/s
gamma[i+1,j]=sample(c(0,1),1,prob=c(prob0,prob1))
}
print(i)
}
gamma
}

################################################################################

# R program to estimate the coefficients of an MA(p) model
# stealing the R program for the AR(p) model

# Negative noises are generated, then conditioned upon
# Model parameters are mu, sigma and lambda (vector of roots)
# complex roots are represented by real(z),im(z)

# WArning: the same inversion algorithm as in the AR model
# is used but the sign convention for the parameters of the
# MA model are opposite to those of the AR model, namely
#  	x_t = mu + epsilon_t + sum_i theta_i epsilon_t-i
# hence minuses here and there

# Set-up parameters
p=9		# MA order (minimum value 2)
q=p		# to reduce programming errors
T=550

x=matrix(scan("Dyndata"),ncol=5,byrow=T)[1:(T+50),5] # observations
#hid=x[(T+1):(T+50)] 		   # prediction purposes
x=x[1:T]
#x=scan("MA.3634")
#T=length(x)			   # number of observations

W=10000   		           # number of MCMC steps

mu=mean(x)
varef=var(x)
sig2=varef  #Rough idea

# Initial values = all real, no problem!
lambdareal=2*runif(p)-1
preal=p
pcomp=0
lambdacomp=0
eps=rnorm(p,sd=sqrt(sig2)) #Just about anythin'

#----------------------------------------------------------------------------
llog=function(pr,pc,lr,lc,compsi=T,pepsi=0,eps){
# requires the negative time index epsilons, of dimension q

# Likelihood representation

if (compsi){

# Moving to coefficients from roots
# A generic program that has nothing to
# do with AR nor MA, just polynomials : No Change

Psi=matrix(0,ncol=p,nrow=p+1)
Psi[1,]=1	# Always 1 because of inverse notations

if (pr>0){
  Psi[2,1]=-lr[1]

  if (pr>1){
    for (i in 2:pr)
      Psi[2:(i+1),i]=Psi[2:(i+1),i-1]-lr[i]*Psi[1:i,i-1]
    }
  }

if (pc>0){
  if (pr>0){

    Psi[2,pr+2]=-2*lc[1]+Psi[2,pr]
    Psi[3:(pr+3),pr+2]=(lc[1]^2+lc[2]^2)*Psi[1:(pr+1),pr]-2*lc[1]*Psi[2:(pr+2),pr]+Psi[3:(pr+3),pr]
    }else{

    Psi[2,2]=-2*lc[1];
    Psi[3,2]=(lc[1]^2+lc[2]^2);
    }

  if (pc>2){
  for (i in seq(4,pc,2)){

     pri=pr+i
     prim=pri-2
     Psi[2,pri]=-2*lc[i-1]+Psi[2,prim]
     Psi[3:(pri+1),pri]=(lc[i-1]^2+lc[i]^2)*Psi[1:(pri-1),prim]-2*lc[i-1]*Psi[2:pri,prim]+Psi[3:(pri+1),prim]
     }
     }
  }
    
Psi=Psi[2:(p+1),p] #We do not need Psi[1]=1 

}else{ # no computation of psi

  Psi=pepsi
  }

x=x-mu #getting rid of the mean

# construction of the epsilonhats
# the psi's are the -theta's
# and the heps[1:p] are the negative epsilons

heps=rep(0,T+p)
heps[1:p]=eps	# Simulated ones
for (i in 1:T)
  heps[p+i]=x[i]+sum(rev(Psi)*heps[i:(p+i-1)])

# loglikelihood
# (includes the negative epsilons)
loglike=-((sum(heps^2)/sig2)+(T+p)*log(sig2))/2

x=x+mu #the mean is back

list(ll=loglike,ps=Psi)
}
#----------------------------------------------------------------------------------

# Starting likelihood
llo=llog(pr=preal,pc=pcomp,lr=lambdareal,lc=lambdacomp,eps=eps)

# MCMC reversible  scheme
psis=matrix(0,ncol=p,nrow=W)
mus=rep(0,W)
sigs=rep(0,W)
ncomp=rep(0,W)
llik=rep(0,W)
indacpt=rep(0,10)
epsrec=matrix(0,ncol=p,nrow=W)
preds=matrix(0,ncol=p,nrow=W)

for (m in 1:W){

# (i) Within muv'----------------------------------------------------

# Pick root(s) to modify
# We assume that the p roots are ranked 1) complex 2) real
  ind=sample(1:p,1)

  if (ind<=pcomp){ # modify two conjugate complex roots
     ind=ind-(ind%%2==0) # in case ind is even
			 # since we need to modify both roots

     #theta=2*pi*runif(1)
     #rho=sqrt(runif(1))  # R^2 cdf obtained by Jacobian

     ppropreal=preal
     ppropcomp=pcomp
     lambpropreal=lambdareal
     lambpropcomp=lambdacomp
     lambpropcomp[ind]=lambdacomp[ind]+.05*rnorm(1) #rho*cos(theta)
     lambpropcomp[ind+1]=lambdacomp[ind+1]+.05*rnorm(1) #rho*sin(theta)

     }else{ # modify one real root

      ppropreal=preal
      ppropcomp=pcomp
      lambpropreal=lambdareal
      lambpropcomp=lambdacomp
      lambpropreal[ind-pcomp]=lambdareal[ind-pcomp]+.05*rnorm(1) #2*runif(1)-1
      }

# Acceptance step

  lloprop=llog(pr=ppropreal,pc=ppropcomp,lr=lambpropreal,lc=lambpropcomp,eps=eps)

  if (log(runif(1))<lloprop$ll-llo$ll){

    llo=lloprop
    preal=ppropreal # not necessary but...
    pcomp=ppropcomp
    lambdacomp=lambpropcomp
    lambdareal=lambpropreal
    indacpt[1]=indacpt[1]+1
    }

# (ii) Between muv'-----------------------------------------------------

# Choose to move up or down (unless one boundary is reached)

   if (preal<2){       # no or one real, thus move down

       ppropreal=preal+2
       ppropcomp=pcomp-2

       ind=sample(1:pcomp,1) # complex root to remove
       ind=ind-(ind%%2==0) 

       if (preal==0){
         lambpropreal=2*runif(2)-1
       }else{

	 lambpropreal=c(lambdareal,lambdacomp[ind]+0.01*rnorm(2)) #r2*runif(2)-1)
       }

       lambpropcomp=lambdacomp[((1:pcomp)!=ind)&((1:pcomp)!=(ind+1))]

       coef=(1/2) * (pi/4)       # to compensate for the obligatory move down
       }

   if (pcomp==0){	# no complex, thus move up

       ppropreal=p-2
       ppropcomp=2

       ind=sample(1:p,2)  # real roots to remove
       lambpropreal=lambdareal[((1:p)!=ind[1])&((1:p)!=ind[2])] 

       #theta=2*pi*runif(1)
       #rho=sqrt(runif(1)) 
       lambpropcomp=c(mean(lambdareal[ind])+.01*rnorm(1),.01*rnorm(1)) #rho*c(cos(theta),sin(theta))

       coef=(1/2) * (4/pi)		# same thing
       }

# Uneven weighting of UP and DOWN moves

   if ((preal>1)&&(pcomp>0)){  	# freedom to move up OR down
      
       if (runif(1)<.1){ #DOWN, one less complex root
       
         ppropcomp=pcomp-2
         ppropreal=preal+2

         ind=sample(1:pcomp,1)		#complex root to remove
         ind=ind-(ind%%2==0) 

         if (ppropcomp>0){
           lambpropcomp=lambdacomp[((1:pcomp)!=ind)&((1:pcomp)!=(ind+1))]
         }else{
           lambpropcomp=0 # should be harmless, necessary for llog function
         }

         lambpropreal=c(lambdareal,lambdacomp[ind]+0.05*rnorm(2)) #r2*runif(2)-1)

         coef=9*(1+(preal<2)) * (pi/4)	#in case the new value is a boundary one

       }else{            #UP, one more complex root 

         ppropreal=preal-2
         ppropcomp=pcomp+2

         #theta=2*pi*runif(1)
         #rho=sqrt(runif(1))  

         ind=sample(1:preal,2) 	#real roots to remove
                                                                                                                                         
         if (ppropreal>0){

           lambpropreal=lambdareal[((1:preal)!=ind[1])&((1:preal)!=ind[2])]
	   }else{
          
              lambpropreal=0 #again, should be harmless, necessary for llog function
           }
	
        lambpropcomp=c(lambdacomp,mean(lambdareal[ind])+.05*rnorm(1),.05*rnorm(1)) #rho*cos(theta),rho*sin(theta))

        coef=(4/pi) * (1+(ppropcomp<p-1))/9. 	#in case the new value is a boundary one
      }
    }

  lloprop=llog(pr=ppropreal,pc=ppropcomp,lr=lambpropreal,lc=lambpropcomp,eps=eps)

  if (log(runif(1))<log(coef)+lloprop$ll-llo$ll){

    llo=lloprop
    preal=ppropreal
    pcomp=ppropcomp
    lambdacomp=lambpropcomp
    lambdareal=lambpropreal
    indacpt[2]=indacpt[2]+1
    }

  psis[m,]=llo$ps

# (iii) Hypermove on mu and sigma
# first, change mu (using a flat prior)

  muold=mu
  mu=rnorm(1,mean=muold,sd=sqrt(5*sig2))
    
  lloprop=llog(pr=preal,pc=pcomp,lr=lambdareal,lc=lambdacomp,compsi=F,pepsi=llo$ps,eps=eps)

  if (log(runif(1))>lloprop$ll-llo$ll){
 
    indacpt[3]=1+indacpt[3]
    mu=muold
    }else{

       llo=lloprop
    }

# second, change sig2 (using a noninformative prior)

  if (runif(1)<.3){

  Psi=llo$ps
  heps=rep(0,2*p)
  heps[1:p]=eps       # Simulated ones
  for (j in (p+1):(2*p))
      heps[j]=x[j]-sum(Psi*heps[(j-p):(j-1)])
  varheps=(2*p)*var(heps)

  sig2old=sig2
  #sig2=exp(rnorm(1,mean=log(sig2),sd=sqrt(.2)*sig2))
  sig2=varheps/rgamma(1,2*p)
  difdens=dgamma(sig2/varheps,2*p,log=T)-dgamma(sig2old/varheps,2*p,log=T)+log(sig2)-log(sig2old)

  lloprop=llog(pr=preal,pc=pcomp,lr=lambdareal,lc=lambdacomp,compsi=F,pepsi=llo$ps,eps=eps)

  if (log(runif(1))>lloprop$ll-llo$ll-difdens){		

    sig2=sig2old
    }else{

        indacpt[4]=3+indacpt[4]
        llo=lloprop
     }
  }

# another attempt

  if (runif(1)<.3){

  sig2old=sig2
  sig2=exp(rnorm(1,mean=log(sig2),sd=sqrt(.1*varef)))

  lloprop=llog(pr=preal,pc=pcomp,lr=lambdareal,lc=lambdacomp,compsi=F,pepsi=llo$ps,eps=eps)
  difdens=log(sig2)-log(sig2old) # prior

  if (log(runif(1))>lloprop$ll-llo$ll-difdens){

     sig2=sig2old
     }else{

        llo=lloprop
        indacpt[5]=3+indacpt[5]
     }
  }

# yet another attempt:
# sig20/(1+sum(psi^2)) is the variance of x

  if (runif(1)<.3){

  sig2old=sig2
  thismean=varef/(1+sum(Psi^2))
  sig2=exp(rnorm(1,mean=thismean,sd=sqrt(varef)))

  # Prior and jacobian cancel
  difdens=dnorm(log(sig2),mean=thismean,sd=sqrt(.5*varef),log=T)-dnorm(log(sig2old),mean=thismean,sd=sqrt(.5*varef),log=T)

  lloprop=llog(pr=preal,pc=pcomp,lr=lambdareal,lc=lambdacomp,compsi=F,pepsi=llo$ps,eps=eps)

  if (log(runif(1))>lloprop$ll-llo$ll-difdens){

     sig2=sig2old
     }else{

        llo=lloprop
        indacpt[6]=3+indacpt[6]
     }
  }

# (iv) Updating negative noises
# one epsilon at a time
  
  Psi=llo$ps

  if (runif(1)<.5){

  heps=rep(0,2*p+1)
  keps=rep(0,p)


  for (i in 1:q){ #Actually, 'tis (-q+1):0 ...

# construction of the epsilonhats
    x = x-mu
    heps[1:p]=eps	# Simulated ones
    for (j in (p+1):(2*p+1))
      heps[j]=x[j]+sum(rev(Psi)*heps[(j-p):(j-1)])

    heps[i]=0
    for (j in 1:(q-i+1))
      keps[j]=x[j]+sum(rev(Psi)*heps[j:(j+p-1)])
    x = x+mu

    epsvar = 1/sum(c(1,Psi[i:q]^2))
    epsmean = sum(Psi[i:q]*keps[1:(q-i+1)])*epsvar
    epsmean = epsmean/epsvar
    epsvar = sig2*epsvar

    propeps = rnorm(1,mean=epsmean,sd=sqrt(epsvar))
    epspr=eps
    epspr[i]=propeps
    lloprop=llog(pr=preal,pc=pcomp,lr=lambdareal,lc=lambdacomp,compsi=F,pepsi=Psi,eps=epspr)
    propsal1=dnorm(propeps,mean=epsmean,sd=sqrt(epsvar),log=T)

  # The acceptance probability requires computing the inverse epsilon_hats,
  # sadly enough...

    x = x-mu
    heps[i]=propeps
    for (j in (p+1):(2*p+1))
      heps[j]=x[j]+sum(rev(Psi)*heps[(j-p):(j-1)])

    heps[i]=0
    for (j in 1:(q-i+1))
      keps[j]=x[j]+sum(rev(Psi)*heps[j:(j+p-1)])
    x = x+mu

    epsvar = 1/sum(c(1,Psi[i:q]^2))
    epsmean = sum(Psi[i:q]*keps[1:(q-i+1)])
    epsmean = epsmean*epsvar
    epsvar = sig2*epsvar 
    propsal0=dnorm(eps[i],mean=epsmean,sd=sqrt(epsvar),log=T)

    if (log(runif(1))<lloprop$ll-llo$ll-propsal1+propsal0){ #Big change!

      indacpt[7]=(2/p)+indacpt[7]
      eps[i]=propeps;
      llo=lloprop
      }
  }
  }

# (v) Second attempt
# one epsilon at a time + randowak

  if (runif(1)<.5){

  heps=rep(0,2*p+1)
  keps=rep(0,p)

  for (i in 1:q){ 

  # construction of the epsilonhats
    propeps = rnorm(1,mean=eps[i],sd=0.1*sqrt(sig2))

    heps[1:p]=eps       # Simulated ones
    heps[i]=propeps
    keps=heps[1:p]

    allop=llog(pr=preal,pc=pcomp,lr=lambdareal,lc=lambdacomp,compsi=F,pepsi=Psi,eps=keps)

    if (log(runif(1))<allop$ll-llo$ll){ 

      indacpt[8]=(2/p)+indacpt[8]
      eps[i]=propeps;
      llo=allop
      }
  }
  }

  psis[m,]=llo$ps
  mus[m]=mu
  sigs[m]=sig2
  llik[m]=llo$ll
  ncomp[m]=pcomp
  epsrec[m,]=eps

  # prediction of the next q observations
  # hat x_t = hat mu + sum_i hat theta_i hat epsilon_t-i
  # where hat x_t is the previous predicted value (or zero)

  heps=rep(0,T+2*p)
  heps[1:p]=eps
  hatx=heps
  hatx[1:T]=x-mu

  for (j in (p+1):(p+T))
      heps[j]=x[j-p]+sum(rev(Psi)*heps[(j-p):(j-1)])
  for (j in (p+T+1):(p+T+p)){
      hatx[j-p]=-sum(rev(Psi)*heps[(j-p):(j-1)])
      }

  preds[m,]=hatx[(T+1):(T+p)]+mu
  }

indacpt/W

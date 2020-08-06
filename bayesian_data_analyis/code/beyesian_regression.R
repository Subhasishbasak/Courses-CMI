load(file="data_stock.RData")
library(bayeslm)


# We convert that SBI's savings accound 
# rate at daily level
Rf<-4/252

Nifty50<-data_stock$Nifty50
HDFCBANK<-data_stock$HDFCBANK
RELIANCE<-data_stock$RELIANCE
MARUTI<-data_stock$MARUTI
TCS<-data_stock$TCS

ln_rt_nifty50<-diff(log(Nifty50))*100-Rf
ln_rt_hdfcbnk<-diff(log(HDFCBANK))*100-Rf
ln_rt_rel<-diff(log(RELIANCE))*100-Rf
ln_rt_maruti<-diff(log(MARUTI))*100-Rf
ln_rt_tcs<-diff(log(TCS))*100-Rf

## log-return of the portfolio
ln_r = cbind(ln_rt_hdfcbnk,ln_rt_rel
             ,ln_rt_maruti,ln_rt_tcs)
w = c(0.2,0.3,0.25,0.25)
ln_rt_portf = ln_r%*%w


INR_USD <- data_stock$INR_USD
ln_rt_inr_usd<-diff(log(INR_USD))*100-Rf

Bank_indx <- data_stock$Nifty_Bank_indx
ln_rt_bank_indx <- diff(log(Bank_indx))*100 - Rf


####
## Capital Asset Pricing Model
## OLS Analysis of Regression Model

capm_hdfcbnk<-lm(ln_rt_hdfcbnk~ln_rt_nifty50)
summary(capm_hdfcbnk)


capm_rel<-lm(ln_rt_rel~ln_rt_nifty50)
summary(capm_rel)


capm_maruti<-lm(ln_rt_maruti~ln_rt_nifty50)
summary(capm_maruti)

capm_tcs<-lm(ln_rt_tcs~ln_rt_nifty50)
summary(capm_tcs)

capm_portf<-lm(ln_rt_portf~ln_rt_nifty50)
summary(capm_portf)


plot(ln_rt_nifty50
     ,ln_rt_portf
     ,xlab="Nifty 50"
     ,ylab="Portfolio"
     ,pch=20)
abline(capm_portf
       ,col="blue"
       ,lwd=2)
grid(col="red")

rse<-summary(capm_portf)$sigma
al <-capm_portf$coefficients[1]-2*rse
b <- capm_portf$coefficients[2]
abline(a= al,b= b,col=3,lty=2,lwd=2)
au <-capm_portf$coefficients[1]+2*rse
abline(a=au,b=b,col=3,lty=2,lwd=2)

### Bayesian Regression or Bayes CAPM Lazy Implementation
### with flat prior

fit <- bayeslm(ln_rt_portf ~ ln_rt_nifty50
               , prior = 'horseshoe'
               ,N = 30000, burnin = 1000)

summary(fit)


### Bayesian Regression or Bayes CAPM Detail Implementation

rP<-ln_rt_portf
n<-length(ln_rt_nifty50)
Intercept<-rep(1,n)
f<-cbind(Intercept,ln_rt_nifty50)

## parameters of prior distribution
nms<-c("alpha","beta")
k<-length(nms)
beta_0<-rep(0,k)
Lambda_0<-diag(1,nrow=length(beta_0))
a_0<-0.1
b_0<-0.01


## parameters of posterior distribution
Lambda_pi<-t(f)%*%f+Lambda_0
Lambda_pi_inv<-solve(Lambda_pi)
beta_pi<-Lambda_pi_inv%*%(t(f)%*%rP
                          +Lambda_0%*%beta_0)
rownames(beta_pi)<-nms
a_pi<-a_0+n/2
err <-(rP-f%*%beta_pi)
b_pi<-b_0+0.5*t(err)%*%(err)

## Simulate from Posterior distribution
set.seed(1)
library(mvtnorm)
N.sim<-10000
burnin<-100
sigma<-rep(NA,N.sim)
beta.draws <- matrix(NA
                     , nrow=N.sim
                     , ncol=k)
colnames(beta.draws)<-nms

for(i in 1:(N.sim+burnin)){
  sigma_2_star<-1/rgamma(1,a_pi,b_pi)
  S <- sigma_2_star*Lambda_pi_inv
  beta_star<-rmvnorm(1,mean = beta_pi
                     ,sigma = S)
  if(i >burnin){
    sigma[i-burnin]<-sqrt(sigma_2_star)
    beta.draws[i-burnin,]<-beta_star
  }
}
theta<-cbind(beta.draws,sigma)

head(round(theta,digits = 3))
post_median<-apply(theta,2,median)
post_mean <- apply(theta,2,mean)
post_sd   <- apply(theta,2,sd)
post_band <- apply(theta,2,quantile
                   ,prob=c(0.025,0.975))
post_sumry<-rbind( post_median     
                   ,post_mean
                   ,post_sd
                   ,post_band)
rownames(post_sumry)<-c("median"
                        ,"mean"
                        ,"sd"
                        ,"2.5%"
                        ,"97.5%")
round(post_sumry,digits = 4)


## P(alpha > 0 | Data)
length(theta[theta[,"alpha"]>0
             ,"alpha"])/nrow(theta)


par(mfrow=c(3,2))
plot(ts(theta[,"alpha"])
     ,ylab = "alpha")
plot(density(theta[,"alpha"])
     ,main = ""
     ,xlab = "alpha")
plot(ts(theta[,"beta"])
     ,ylab = "beta")
plot(density(theta[,"beta"])
     ,main = ""
     ,xlab = "beta")
plot(ts(theta[,"sigma"])
     ,ylab = "sigma")
plot(density(theta[,"sigma"])
     ,main = ""
     ,xlab = "sigma")

## Stress Testing

market_drop=-3
rP.sim<-rep(NA,N.sim)
x0<-matrix(c(1,market_drop),nrow=1,ncol = k)
for(i in 1:N.sim){
  ErP <- x0%*%theta[i,1:k]
  rP.sim[i]<-rnorm(1,mean = ErP, sd= theta[i,(k+1)])
}
hist(rP.sim,probability = TRUE,col="orange",main="")
lines(density(rP.sim),lwd=2,col="red")

length(rP.sim[rP.sim<=-3])/N.sim

expected_drop<-mean(rP.sim)
expected_drop

summary(rP.sim)
quantile(rP.sim,probs = c(0.025,0.975))


### Test the assumption
shapiro.test(capm_portf$residuals)

qqnorm(capm_portf$residuals
       ,cex=0.3
       ,pch=20
       ,xlim=c(-4,4)
       ,ylim=c(-4,4)
       ,xlab="idiosyncratic return with fitted Gaussian distribution"
       ,ylab="observed residual or idiosyncratic return"
       ,main="Gaussian Q-Q Plot")
abline(a=0,b=1,col="blue")
grid(col="red")

## Does residual, i.e., idiosyncratic return follow Laplace distribution?

## density of Laplace distribution 
dlaplace<-function(x,mu,lambda){
  exp(-abs(x-mu)/lambda)/(2*lambda)
}
## The negative log-likelihood function
neg_log_likelihood<-function(data,para){
  mu<-para[1]
  lambda<-(para[2])
  f<-dlaplace(x=data
              ,mu=mu
              ,lambda=lambda)
  return(-sum(log(f)))
}

## initial parameters
para.init<-c(0,1)

## residual or idiosyncratic return 
## from capm_portf
e<-summary(capm_portf)$residual

mle<-optim(para.init
           ,neg_log_likelihood
           ,data=e)

mle$par
## The following function simulate from 
## the fitted Laplace distribution
rlaplace<-function(n,mu,lmbda){
  # n : sample size
  # mu: location parameter
  # lambda : scale parameter
  u<-runif(n,-0.5,0.5)
  sgn<-rep(NA,n)
  for(i in 1:n){
    if(u[i]<0) sgn[i]<--1
    if(u[i]>0) sgn[i]<-1
    if(u[i]==0) sgn[i]<-0
  }
  r<-mu-lmbda*sgn*log(1-2*abs(u))
  return(r)
}

## simulate the result
set.seed(4155)
r.sim<-rlaplace(n=1000
                ,mu=mle$par[1]
                ,lmbda=mle$par[2])
qqplot(r.sim,e
       ,xlim=c(-4,4)
       ,ylim=c(-4,4)
       ,cex=0.5
       ,pch=20
       ,xlab="idiosyncratic return with fitted Laplace distribution"
       ,ylab="observed residual or idiosyncratic return"
       ,main="Laplace Q-Q Plot")
abline(a=0,b=1,col="blue")
grid(col="red")

### Bayesian CAPM where Idiosyncratic Return follows Laplace Distribution

log_likelihood <- function(param,y,x){
  a = param[1]
  b = param[2]
  lambda = param[3]
  pred = a + b*x 
  likelihoods = -log(2*lambda)-abs(y-pred)/lambda
  sumll = sum(likelihoods)
  return(sumll)   
}


log_prior <- function(param,x){
  a = param[1]
  b = param[2]
  lambda = param[3]
  a_prior = dcauchy(a,0,1,log = T) 
  b_prior = dcauchy(b,0,1,log = T) 
  scale_prior = dgamma(lambda
                       ,1,1
                       ,log = T)
  
  return(a_prior+b_prior+scale_prior)
}

log_posterior <- function(param,y,x){
  like <- log_likelihood(param=param
                         ,y=y,x=x)
  prior <- log_prior(param=param,x=x)
  post  <- like + prior
  return ( post )
}

### Model Fitting with Metropolis-Hastings

proposalfunction <- function(param,x){
  X=cbind(rep(1,length(x)),x)
  S=param[3]*solve(t(X)%*%X)
  prop<-c(rmvnorm(1
                  ,mean = param[1:2]
                  ,sigma = S)
          ,rgamma(1,param[3]*5,5))
  return(prop)
}

run_metropolis <- function(startvalue, N.sim, burnin){
  iterations = N.sim + burnin
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,],x=x)
    
    probab = exp(log_posterior(param=proposal
                               ,y=y,x=x) 
                 - log_posterior(param=chain[i,]
                                 ,y=y,x=x))
    
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

y=ln_rt_portf
x=ln_rt_nifty50
startvalue = c(0,1,0.1)
N.sim=30000
burnin=1000

set.seed(1)
chain = run_metropolis(startvalue=startvalue
                       ,N.sim=N.sim
                       ,burnin=burnin)
colnames(chain)<-c("alpha","beta","lambda")
chain=chain[(burnin+1):nrow(chain),]

apply(chain,2,mean)
apply(chain,2,sd)
apply(chain,2,quantile,prob=c(0.025
                              ,0.975))

## P(alpha > 0 | Data)
length(chain[chain[,"alpha"]>0
             ,"alpha"])/nrow(chain)

par(mfrow=c(3,2))
plot(ts(chain[,"alpha"])
     ,ylab = "alpha")
plot(density(chain[,"alpha"])
     ,main = ""
     ,xlab = "alpha")
plot(ts(chain[,"beta"])
     ,ylab = "beta")
plot(density(chain[,"beta"])
     ,main = ""
     ,xlab = "beta")
plot(ts(chain[,"lambda"])
     ,ylab = "lambda")
plot(density(chain[,"lambda"])
     ,main = ""
     ,xlab = "lambda")

### Stress Test 

market_drop=-3
rP.sim<-rep(NA,nrow(chain))
x0<-matrix(c(1,market_drop),nrow=1,ncol = k)
for(i in 1:nrow(chain)){
  ErP <- x0%*%chain[i,1:k]
  rP.sim[i]<-rnorm(1,mean = ErP, sd= chain[i,(k+1)])
}
hist(rP.sim,probability = TRUE,col="orange",main="")
lines(density(rP.sim),lwd=2,col="red")
quantile(rP.sim,prob = c(0.025,0.975))






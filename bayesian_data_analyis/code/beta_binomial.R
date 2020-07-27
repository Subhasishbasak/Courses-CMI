## Marketing Effectiveness Model with Beta Binomial

n = 24
y = 17

## Flat prior : pi ~ Unif(0,1)
## Posterior beta(y+1,n-y+1)

pai<-seq(0,1,0.01)
d_post=dbeta(pai,y+1,n-y+1)
plot(pai,d_post,type="l",xlab="pai",ylab="Posterior density")

## Simulate 1000 samples from Beta(18,8)

pai<-rbeta(1000,y+1,n-y+1)
hist(pai,prob=TRUE,col="blue",xlim=c(0,1),main="")
summary(pai)
quantile(pai,prob=c(0.025,0.975))

pai<-seq(0,1,0.01)
lines(pai,dbeta(pai,y+1,n-y+1),col="red",lwd=2)

### Credible Intervals or Highest Posterior Density Region
library(graphics)
LB<-qbeta(0.025,18,8)
UB<-qbeta(0.975,18,8)
cord.pai<-c(LB,seq(LB,UB,0.01),UB)
cord.dpai<-c(0,dbeta(seq(LB,UB,0.01),18,8),0)
curve(dbeta(x,18,8),xlim=c(0,1),xlab="pai",ylab="posterior density")
polygon(cord.pai,cord.dpai,col="blue")


## Four different choices of prior

par(mfrow=c(2,2))
pai<-seq(0,1,0.01)
plot(pai,dbeta(pai,5,5),main="Beta(5,5)",type="l",ylab="")
plot(pai,dbeta(pai,3,10),main="Beta(3,10)",type="l",ylab="")
plot(pai,dbeta(pai,10,3),main="Beta(10,3)",type="l",ylab="")
plot(pai,dbeta(pai,100,30),main="Beta(100,30)",type="l",ylab="")



## Modeling expert's opinion with conjugate beta prior
integrate(dbeta,lower=0.5,upper=0.9,shape1=7,shape2=3)

obj<-function(lower,upper,a,b,p){
  I = integrate(dbeta,lower,upper,shape1=a,shape2=b)$value
  dif = abs(I-p)
  return(dif)
}

obj(lower=0.5,upper=0.9,a=7,b=3,p=0.9)

a<-b<-seq(1,10,0.1)
res<-c(1,1,obj(lower=0.5,upper=0.9,a=1,b=1,p=0.9))
for(i in 1:length(a)){
  for(j in 1:length(b)){
      val<-c(a[i],b[j],obj(lower=0.5,upper=0.9,a=a[i],b=b[j],p=0.9))
      res<-rbind(res,val)
  }
}

round(res[which.min(res[,3]),],2)

integrate(dbeta,0.5,0.9,shape1=9.2,shape2=4.3)


## Posterior Analysis with Expert's Opinion as Beta prior

y<-17
n<-24
a<-9.2
b<-4.3
pai<-seq(0,1,0.01)
plot(pai,dbeta(pai,y+a,n-y+b),type="l",ylab="density",col=2,ylim=c(0,6))
lines(pai,dbeta(pai,a,b),col=4) 
text<-c("Prior","Posterior")
legend(0.01,6,text,col=c("blue","red"),lty=c(1,1))

### Posterior Summary

pai<-rbeta(1000,y+a,n-y+b)
summary(pai)
quantile(pai,prob=c(0.025,0.975))


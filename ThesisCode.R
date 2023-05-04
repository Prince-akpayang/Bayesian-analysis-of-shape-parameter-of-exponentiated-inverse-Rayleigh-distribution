# To compare the MLE and Bayesian Estimates assuming Jeffrey and Gamma Priors. under SELF, PLF and QLF

shapeparameterestimator <- function(n, R, alpha, sigma, a, b){
  
  # To generate random numbers from the Exponentiated Inverse Rayleigh distribution and also get the estimates of alpha
  Estimates <- function(n, alpha, sigma, a, b){
    x <- c()
    for (i in 1:n) {
      u <- runif(1)
      x[i] <- sigma/sqrt(-log(1-(1-u)^(1/alpha)))
    }
    x
    # Computes the various estimates
    
    MLE <- -n/sum(log(1-exp(-(sigma/x)^2)))
    
    # Under Jeffreys' Prior
    
    JSELF <- n/log(prod(1-exp(-(sigma/x)^2))^(-1))
    JPLF <- sqrt((n+1)*n)/log(prod(1-exp(-(sigma/x)^2))^(-1))
    JQLF <- (n-2)/log(prod(1-exp(-(sigma/x)^2))^(-1))
    
    # Under the Gamma Prior
    
    GSELF <- (n+a)/(b+log(prod(1-exp(-(sigma/x)^2))^(-1)))
    GPLF <- sqrt((n+a)*(n+a-1))/(b+log(prod(1-exp(-(sigma/x)^2))^(-1)))
    GQLF <- (n+a-2)/(b+log(prod(1-exp(-(sigma/x)^2))^(-1)))
    
    list(MLE, JSELF, JPLF, JQLF, GSELF, GPLF, GQLF)
    
  }
  
  # To compute the various estimates a large number of times as in monte carlo studies
  
  monte.carloE <- replicate(R, Estimates(n, alpha, sigma, a, b))
  
  # Computes the various posterior risk
  
  posterior.risk <- function(n, alpha, sigma, a, b){
    x <- c()
    for (i in 1:n) {
      u <- runif(1)
      x[i] <- sigma/sqrt(-log(1-(1-u)^(1/alpha)))
    }
    x
    
    # Under Jeffreys' Prior
    
    RJSELF <- n/(log(prod(1-exp(-(sigma/x)^2))^(-1)))^2
    RJPLF <- 2*((sqrt((n+1)*n)- n)/(log(prod(1-exp(-(sigma/x)^2))^(-1))))
    RJQLF <- 1/(n-1)
    
    # Under the Gamma Prior
    
    RGSELF <- (n+a)/((b+log(prod(1-exp(-(sigma/x)^2))^(-1))))^2
    RGPLF <- 2*((sqrt((n+a)*(n+a+1))-(n+a))/(b+log(prod(1-exp(-(sigma/x)^2))^(-1))))
    RGQLF <- 1/(n+a-1)
    
    list(RJSELF, RJPLF, RJQLF, RGSELF, RGPLF,  RGQLF)
    
  }
  
  # To compute the various posterior risk a large number of times as in monte carlo studies
  
  monte.carloR <- replicate(R, posterior.risk(n, alpha, sigma, a, b))
  
  Ests <- matrix(NA, nrow = 4, ncol = 7, 
                 dimnames = list(c("Estimates", "BIAS", "MSE", "POSTERIOR RISK"),
                                 c("MLE", "JSELF", "JPLF", "JQLF","GSELF", "GPLF", "GQLF")))
  
  # Averages out each of the estimates
  Ests[1,1] <- round(mean(unlist(monte.carloE[1,])),4)
  Ests[1,2] <- round(mean(unlist(monte.carloE[2,])),4)
  Ests[1,3] <- round(mean(unlist(monte.carloE[3,])),4)
  Ests[1,4] <- round(mean(unlist(monte.carloE[4,])),4)
  Ests[1,5] <- round(mean(unlist(monte.carloE[5,])),4)
  Ests[1,6] <- round(mean(unlist(monte.carloE[6,])),4)
  Ests[1,7] <- round(mean(unlist(monte.carloE[7,])),4)
  
  # Computes the bias of each of the estimates
  Ests[2,] <- abs(alpha - Ests[1,])
  
  # Computes the MSE of each of the estimates
  Ests[3,1] <- round(sum((alpha - unlist(monte.carloE[1,]))^2)/R, 4)
  Ests[3,2] <- round(sum((alpha - unlist(monte.carloE[2,]))^2)/R, 4)
  Ests[3,3] <- round(sum((alpha - unlist(monte.carloE[3,]))^2)/R, 4)
  Ests[3,4] <- round(sum((alpha - unlist(monte.carloE[4,]))^2)/R, 4)
  Ests[3,5] <- round(sum((alpha - unlist(monte.carloE[5,]))^2)/R, 4)
  Ests[3,6] <- round(sum((alpha - unlist(monte.carloE[6,]))^2)/R, 4)
  Ests[3,7] <- round(sum((alpha - unlist(monte.carloE[7,]))^2)/R, 4)
  
  # compute posterior risk
  Ests[4,1] <- round(NA, 4)
  Ests[4,2] <- round(mean(unlist(monte.carloR[1,])),4)
  Ests[4,3] <- round(mean(unlist(monte.carloR[2,])),4)
  Ests[4,4] <- round(mean(unlist(monte.carloR[3,])),4)
  Ests[4,5] <- round(mean(unlist(monte.carloR[4,])),4)
  Ests[4,6] <- round(mean(unlist(monte.carloR[5,])),4)
  Ests[4,7] <- round(mean(unlist(monte.carloR[6,])),4)
  
  Ests
}

# alpha = 0.5, sigma = 1
set.seed(1452)
shapeparameterestimator(n = 10, R = 10000, alpha = 0.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 20, R = 10000, alpha = 0.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 30, R = 10000, alpha = 0.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 40, R = 10000, alpha = 0.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 60, R = 10000, alpha = 0.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 95, R = 10000, alpha = 0.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 125, R = 10000, alpha = 0.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 150, R = 10000, alpha = 0.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 170, R = 10000, alpha = 0.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 195, R = 10000, alpha = 0.5, sigma = 1, a= 5, b = 3)

# alpha = 1, sigma = 1
set.seed(1452)
shapeparameterestimator(n = 10, R = 10000, alpha = 1, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 20, R = 10000, alpha = 1, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 30, R = 10000, alpha = 1, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 40, R = 10000, alpha = 1, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 60, R = 10000, alpha = 1, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 95, R = 10000, alpha = 1, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 125, R = 10000, alpha = 1, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 150, R = 10000, alpha = 1, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 170, R = 10000, alpha = 1, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 195, R = 10000, alpha = 1, sigma = 1, a= 5, b = 3)

# alpha = 1.5, sigma = 1
set.seed(1452)
shapeparameterestimator(n = 10, R = 10000, alpha = 1.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 20, R = 10000, alpha = 1.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 30, R = 10000, alpha = 1.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 40, R = 10000, alpha = 1.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 60, R = 10000, alpha = 1.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 95, R = 10000, alpha = 1.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 125, R = 10000, alpha = 1.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 150, R = 10000, alpha = 1.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 170, R = 10000, alpha = 1.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 195, R = 10000, alpha = 1.5, sigma = 1, a= 5, b = 3)

# alpha = 2, sigma = 1
set.seed(1452)
shapeparameterestimator(n = 10, R = 10000, alpha = 2, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 20, R = 10000, alpha = 2, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 30, R = 10000, alpha = 2, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 40, R = 10000, alpha = 2, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 60, R = 10000, alpha = 2, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 95, R = 10000, alpha = 2, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 125, R = 10000, alpha = 2, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 150, R = 10000, alpha = 2, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 170, R = 10000, alpha = 2, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 195, R = 10000, alpha = 2, sigma = 1, a= 5, b = 3)

# alpha = 2.5, sigma = 1
set.seed(1452)
shapeparameterestimator(n = 10, R = 10000, alpha = 2.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 20, R = 10000, alpha = 2.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 30, R = 10000, alpha = 2.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 40, R = 10000, alpha = 2.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 60, R = 10000, alpha = 2.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 95, R = 10000, alpha = 2.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 125, R = 10000, alpha = 2.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 150, R = 10000, alpha = 2.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 170, R = 10000, alpha = 2.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 195, R = 10000, alpha = 2.5, sigma = 1, a= 5, b = 3)

# alpha = 3, sigma = 1
set.seed(1452)
shapeparameterestimator(n = 10, R = 10000, alpha = 3, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 20, R = 10000, alpha = 3, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 30, R = 10000, alpha = 3, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 40, R = 10000, alpha = 3, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 60, R = 10000, alpha = 3, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 95, R = 10000, alpha = 3, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 125, R = 10000, alpha = 3, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 150, R = 10000, alpha = 3, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 170, R = 10000, alpha = 3, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 195, R = 10000, alpha = 3, sigma = 1, a= 5, b = 3)

# alpha = 3.5, sigma = 1
set.seed(1452)
shapeparameterestimator(n = 10, R = 10000, alpha = 3.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 20, R = 10000, alpha = 3.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 30, R = 10000, alpha = 3.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 40, R = 10000, alpha = 3.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 60, R = 10000, alpha = 3.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 95, R = 10000, alpha = 3.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 125, R = 10000, alpha = 3.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 150, R = 10000, alpha = 3.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 170, R = 10000, alpha = 3.5, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 195, R = 10000, alpha = 3.5, sigma = 1, a= 5, b = 3)

# alpha = 4, sigma = 1
set.seed(1452)
shapeparameterestimator(n = 10, R = 10000, alpha = 4, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 20, R = 10000, alpha = 4, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 30, R = 10000, alpha = 4, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 40, R = 10000, alpha = 4, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 60, R = 10000, alpha = 4, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 95, R = 10000, alpha = 4, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 125, R = 10000, alpha = 4, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 150, R = 10000, alpha = 4, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 170, R = 10000, alpha = 4, sigma = 1, a= 5, b = 3)
shapeparameterestimator(n = 195, R = 10000, alpha = 4, sigma = 1, a= 5, b = 3)


# Probability density function
deird<-function(x,sigma,alpha){
  p<-(2*alpha*sigma^2)/x^3*exp(-(sigma/x^2))*(1- exp(-(sigma/x)^2))^(alpha-1)
  return(p)
}

# cumulative distribution function
peird<-function(x,sigma,alpha){
  d<-1-(1-exp(-(sigma/x)^2))^alpha
  return(d)
}

# PdfCurves
x <- seq(0, 5, by=.001)
plot(x, deird(x, 1,0.5), type="l", ylim=c(0,5),  ylab="Density",
     main=" ",lwd=3)
lines(x, deird(x, 1,1), col=2,lwd=3)
lines(x, deird(x, 1,1.5), col=3,lwd=3)
lines(x, deird(x, 1,2), col=4,lwd=3)
lines(x, deird(x, 1,2.5), col=5,lwd=3)
lines(x, deird(x, 1,3), col=6,lwd=3)
lines(x, deird(x, 1,3.5), col=7,lwd=3)
lines(x, deird(x, 1,4), col=8,lwd=3)


#legend(3,1.7, legend = c(''exp(1)'', ''exp(2)'', ''exp(3)'', ''exp(4)'',''exp(5)'',''exp(6)'',''exp(7)'',''exp(8)''),

pch=c(20,20,20,20,20,20,20,20)
col=c(1,2,3,4,5,6,7,8) 

legend(par('usr')[2], par('usr')[4], xjust=1,
       c(c(as.expression(substitute(EIRD(sigma==1,alpha==0.5))),
           as.expression(substitute(EIRD(sigma==1,alpha==1))),
           as.expression(substitute(EIRD(sigma==1,alpha==1.5))),
           as.expression(substitute(EIRD(sigma==1,alpha==2))),
           as.expression(substitute(EIRD(sigma==1,alpha==2.5))),
           as.expression(substitute(EIRD(sigma==1,alpha==3))),
           as.expression(substitute(EIRD(sigma==1,alpha==3.5))),
           as.expression(substitute(EIRD(sigma==1,alpha==4))))),
       
       lty=c(2,2,2,2,2,2,2,2), lwd =c(3,3,3,3,3,3,3,3),
       pch=c(20,20,20,20,20,20,20,20), col=c(1,2,3,4,5,6,7,8)) 

####To find the credible interval for alpha assuming the Jeffrey prior###
credintajeff<-function(n,alpha,sigma,startvalue,iterations,BurnIn){
  n=n
  alpha=alpha
  startvalue=startvalue
  iterations=iterations
  BurnIn=BurnIn
  x<-c()
  for (i in 1:n) {
    u<-runif(1)
    x[i]<- sigma/sqrt(-log(1-(1-u)^(1/alpha)))
  }
  x
  
  likelihood<-function(alpha){
    ll<-((alpha^(n))*(prod(1-exp(-(sigma/x)^2))^(alpha-1)))
  }
  summary(likelihood(alpha))
  log(likelihood(alpha))
  
  prior<-function(alpha){
    jefferyPrior<-(1/alpha)
  }
  summary(prior(alpha))
  
  posterior<-function(alpha){
    return(likelihood(alpha)*prior(alpha))
  }
  summary(posterior(alpha))
  
  randnum=c()
  randnum[1] = startvalue    
  for (j in 1:20000) {
    randnum[j+1]<-rnorm(1,alpha,0.5)
    round(randnum,4)
  }
  round(randnum,4)
  
  alphaI = c()
  for (i in 1:iterations){
    
    proposal= randnum[i+10]
    
    probab = posterior(alpha=proposal)/posterior(alpha=randnum[i])
    
    if (probab>1){alphaI[i+1]=randnum[i]
    }else{alphaI[i+1]=alphaI[i]}
  }
  
  alpha_I<-c(alphaI[-(1:BurnIn)])
  return(alpha_I)
}

credintajefffin<-function(alphaI,n=1800,al=0.1){
  S<-sort(alphaI)
  cre.int<-c(S[round((al/2)*n)],S[round((1-(al/2))*n)])
  return(cre.int)
}

{
  SB=
    B=
    set.seed(18)
  n10<-credintajeff(10,B,1,SB,2000,201)
  n20<-credintajeff(20,B,1,SB,2000,201)
  n30<-credintajeff(30,B,1,SB,2000,201)
  n40<-credintajeff(40,B,1,SB,2000,201)
  n60<-credintajeff(60,B,1,SB,2000,201)
  n95<-credintajeff(95,B,1,SB,2000,201)
  n125<-credintajeff(125,B,1,SB,2000,201)
  n150<-credintajeff(150,B,1,SB,2000,201)
  n170<-credintajeff(170,B,1,SB,2000,201)
  n195<-credintajeff(195,B,1,SB,2000,201)
  
  n10CI<-credintajefffin(n10)
  n20CI<-credintajefffin(n20)
  n30CI<-credintajefffin(n30)
  n40CI<-credintajefffin(n40)
  n60CI<-credintajefffin(n60)
  n95CI<-credintajefffin(n95)
  n125CI<-credintajefffin(n125)
  n150CI<-credintajefffin(n150)
  n170CI<-credintajefffin(n170)
  n195CI<-credintajefffin(n195)
  
  ci<-rbind(n10CI,n20CI,n30CI,n40CI,n60CI,n95CI,n125CI,n150CI,n170CI,n195CI)
  credible_interval<-matrix(ci,nrow=10,ncol=2,byrow = F,
                            dimnames=list(c("at n=10","at n=20","at n=30","at n=40","at n=60",
                                            "at n=95","at n=125","at n=150","at n=170","at n=195"),
                                          c("lower","upper")))
  credible_interval
}

{
  b=B
  par(mfrow = c(2,5))
  plot(n10,type = "l", ylab="(assuming a Jeffrey prior)",xlab="number of iteration",
       main = "When sample size=10")
  abline(h=b,col="red")
  abline(h=c(n10CI[1],h=n10CI[2]),col="blue")
  
  plot(n20,type = "l" ,ylab="(assuming a Jeffrey prior)",xlab="number of iteration",
       main = "When sample size=20")
  abline(h=b,col="red")
  abline(h=c(n20CI[1],h=n20CI[2]),col="blue")
  
  plot(n30,type = "l" ,ylab="(assuming a Jeffrey prior)",xlab="number of iteration",
       main = "When sample size=30")
  abline(h=b,col="red")
  abline(h=c(n30CI[1],h=n30CI[2]),col="blue")
  
  plot(n40,type = "l" ,ylab="(assuming a Jeffrey prior)",xlab="number of iteration",
       main = "When sample size=40")
  abline(h=b,col="red")
  abline(h=c(n40CI[1],h=n40CI[2]),col="blue")
  
  plot(n60,type = "l" ,ylab="(assuming a Jeffrey prior)",xlab="number of iteration",
       main = "When sample size=60")
  abline(h=b,col="red")
  abline(h=c(n60CI[1],h=n60CI[2]),col="blue")
  
  plot(n95,type = "l" ,ylab="(assuming a Jeffrey prior)",xlab="number of iteration",
       main = "When sample size=95")
  abline(h=b,col="red")
  abline(h=c(n95CI[1],h=n95CI[2]),col="blue")
  
  plot(n125,type = "l" ,ylab="(assuming a Jeffrey prior)",xlab="number of iteration",
       main = "When sample size=125")
  abline(h=b,col="red")
  abline(h=c(n125CI[1],h=n125CI[2]),col="blue")
  
  plot(n150,type = "l" ,ylab="(assuming a Jeffrey prior)",xlab="number of iteration",
       main = "When sample size=150")
  abline(h=b,col="red")
  abline(h=c(n150CI[1],h=n150CI[2]),col="blue")
  
  plot(n170,type = "l" ,ylab="(assuming a Jeffrey prior)",xlab="number of iteration",
       main = "When sample size=170")
  abline(h=b,col="red")
  abline(h=c(n170CI[1],h=n170CI[2]),col="blue")
  
  plot(n195,type = "l" ,ylab="(assuming a Jeffrey prior)",xlab="number of iteration",
       main = "When sample size=195")
  abline(h=b,col="red")
  abline(h=c(n195CI[1],h=n195CI[2]),col="blue")
}

####To find the credible interval for alpha assuming the gamma prior###
credintagam<-function(n,alpha,sigma,startvalue,iterations,BurnIn){
  n=n
  alpha=alpha
  startvalue=startvalue
  iterations=iterations
  BurnIn=BurnIn
  a=5
  b=3
  
  x<-c()
  for (i in 1:n) {
    u<-runif(1)
    x[i]<- sigma/sqrt(-log(1-(1-u)^(1/alpha)))
  }
  x
  
  likelihood<-function(alpha){
    ll<-((alpha^(n))*(prod(1-exp(-(sigma/x)^2))^(alpha-1)))
  }
  summary(likelihood(alpha))
  log(likelihood(alpha))
  
  prior<-function(alpha){
    GammaPrior<-((b^a)/gamma(a))*(alpha^(a-1))*exp(-alpha*b)
  }
  summary(prior(alpha))
  
  posterior<-function(alpha){
    return(likelihood(alpha)*prior(alpha))
  }
  summary(posterior(alpha))
  
  
  randnum=c()
  randnum[1] = startvalue    
  for (j in 1:10000) {
    randnum[j+1]<-rnorm(1,alpha,0.05)
    round(randnum,4)
  }
  
  alphaI = c()
  for (i in 1:iterations){
    
    proposal= randnum[i+10]
    
    probab = posterior(alpha=proposal)/posterior(alpha=randnum[i])
    
    if (probab>1){alphaI[i+1]=randnum[i]
    }else{alphaI[i+1]=alphaI[i]}
  }
  alpha_I<-c(alphaI[-(1:BurnIn)])
  alpha_I
  
}
credintagamfin<-function(alphaI,n=1800,al=0.1){
  S<-sort(alphaI)
  cre.int<-c(S[round((al/2)*n)],S[round((1-(al/2))*n)])
  return(cre.int)
}

{
  set.seed(11111)
  SA=1.4
  A=1
  n10<-credintagam(10,A,1,SA,2000,201)
  n20<-credintagam(20,A,1,SA,2000,201)
  n30<-credintagam(30,A,1,SA,2000,201)
  n40<-credintagam(40,A,1,SA,2000,201)
  n60<-credintagam(60,A,1,SA,2000,201)
  n95<-credintagam(95,A,1,SA,2000,201)
  n125<-credintagam(125,A,1,SA,2000,201)
  n150<-credintagam(150,A,1,SA,2000,201)
  n170<-credintagam(170,A,1,SA,2000,201)
  n195<-credintagam(195,A,1,SA,2000,201)
  
  n10CI<-credintagamfin(n10)
  n20CI<-credintagamfin(n20)
  n30CI<-credintagamfin(n30)
  n40CI<-credintagamfin(n40)
  n60CI<-credintagamfin(n60)
  n95CI<-credintagamfin(n95)
  n125CI<-credintagamfin(n125)
  n150CI<-credintagamfin(n150)
  n170CI<-credintagamfin(n170)
  n195CI<-credintagamfin(n195)
  
  ci<-rbind(n10CI,n20CI,n30CI,n40CI,n60CI,n95CI,n125CI,n150CI,n170CI,n195CI)
  credible_interval<-matrix(ci,nrow=10,ncol=2,byrow = F,
                            dimnames=list(c("at n=10","at n=20","at n=30","at n=40","at n=60",
                                            "at n=95","at n=125","at n=150","at n=170","at n=195"),
                                          c("lower","upper")))
  credible_interval
}

{
  a=A
  par(mfrow = c(2,5))
  plot(n10,type = "l", ylab="(assuming a gamma prior)",xlab="number of iteration",
       main = "When sample size=10")
  abline(h=a,col="red")
  abline(h=c(n10CI[1],h=n10CI[2]),col="blue")
  
  plot(n20,type = "l" ,ylab="(assuming a gamma prior)",xlab="number of iteration",
       main = "When sample size=20")
  abline(h=a,col="red")
  abline(h=c(n20CI[1],h=n20CI[2]),col="blue")
  
  plot(n30,type = "l" ,ylab="(assuming a gamma prior)",xlab="number of iteration",
       main = "When sample size=30")
  abline(h=a,col="red")
  abline(h=c(n30CI[1],h=n30CI[2]),col="blue")
  
  plot(n40,type = "l" ,ylab="(assuming a gamma prior)",xlab="number of iteration",
       main = "When sample size=40")
  abline(h=a,col="red")
  abline(h=c(n40CI[1],h=n40CI[2]),col="blue")
  
  plot(n60,type = "l" ,ylab="(assuming a gamma prior)",xlab="number of iteration",
       main = "When sample size=60")
  abline(h=a,col="red")
  abline(h=c(n60CI[1],h=n60CI[2]),col="blue")
  
  plot(n95,type = "l" ,ylab="(assuming a gamma prior)",xlab="number of iteration",
       main = "When sample size=95")
  abline(h=a,col="red")
  abline(h=c(n95CI[1],h=n95CI[2]),col="blue")
  
  plot(n125,type = "l" ,ylab="(assuming a gamma prior)",xlab="number of iteration",
       main = "When sample size=125")
  abline(h=a,col="red")
  abline(h=c(n125CI[1],h=n125CI[2]),col="blue")
  
  plot(n150,type = "l" ,ylab="(assuming a gamma prior)",xlab="number of iteration",
       main = "When sample size=150")
  abline(h=a,col="red")
  abline(h=c(n150CI[1],h=n150CI[2]),col="blue")
  
  plot(n170,type = "l" ,ylab="(assuming a gamma prior)",xlab="number of iteration",
       main = "When sample size=170")
  abline(h=a,col="red")
  abline(h=c(n170CI[1],h=n170CI[2]),col="blue")
  
  plot(n195,type = "l" ,ylab="(assuming a gamma prior)",xlab="number of iteration",
       main = "When sample size=195")
  abline(h=a,col="red")
  abline(h=c(n195CI[1],h=n195CI[2]),col="blue")
}
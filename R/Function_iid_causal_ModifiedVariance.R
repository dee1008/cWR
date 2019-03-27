## Di Zhang
## Started: 11/29/2016
## Last Updated: 10/19/2017






####################################################################################################
## ----  Generate Semi-Competing Risk Data ----- ##
##
## ----  Subjects in the same cluster belong to the same membership
## 
## Description:
## (1) use Joint distribution method. Bivariate exponential distribution with Gumbel-Hougaard copula
## (2) generate from marginal and conditional distributions
## (3) it allows one covariate: treatment. One can modify this function to include more covariates.
## (4) assumes the means for MVN are 0.
##
## Input:
##
## n:number of subjects per group 
## n.cluster: number of clusters per group
## dim: number of observations per subject (=2 in our case)
## alpha: association parameter from positive stable copula
## lambdaH: exponential rate of nonterminal events
## lambdaD: exponential rate of terminal events\
## lambdaC: exponential rate of censoring
## etaH: hazard ratio in nonterminal events between groups
## etaD: hazard ratio in terminal events between groups
## etaC: hazard ratio in censoring between groups
## shape: shape parameter for frailty (gamma distribution)
## rate: rate parameter for frailty (gamma distribution)
##
##
#####################################################################################################

gumbel.jeong<-function(n,n.clust,dim,alpha,lambdaH,lambdaD,etaH,etaD,shape,rate)
{
  exprand <- matrix(rexp(dim * n), c(n, dim))
  unifpirand <- runif(n, 0, pi)
  exprand2 <- rexp(n)
  beta <- 1/alpha
  stablerand <- sin((1 - beta) * unifpirand)^((1 - beta)/beta) *
    (sin(beta * unifpirand))/(sin(unifpirand))^(1/beta)
  stablerand <- stablerand/(exprand2^(alpha - 1))
  unifrand <- invphigumbel(exprand/stablerand, alpha) # generating bivariate uniform random variables for marginal survival funtions--(*)

  clust.siz<- n/n.clust
  frail <- rep(rgamma(n.clust,shape=shape,rate=rate),each=clust.siz)
  matrix(c(-log(unifrand[,1])/(frail*lambdaH*exp(-etaH)),-log(unifrand[,2])/(frail*lambdaD*exp(-etaD))),c(n,dim)) # inverting specific forms of survival functions in (*) to create
  # true bivariate event times adjusted for event types and trt groups
}





gen_cluster_data<-function(n.sub=30, n.clust=10, dim=2, alpha=2, lambdaH=0.1, lambdaD=0.08, lambdaC=0.01, etaH=0.2, etaD=0.5, etaC=0.1, shape=5, rate=5){


  group0<-gumbel.jeong(n.sub,n.clust,dim,alpha,lambdaH,lambdaD,0,0,shape,rate)
  group1<-gumbel.jeong(n.sub,n.clust,dim,alpha,lambdaH,lambdaD,etaH,etaD,shape,rate)

  true.t<-rbind(group0,group1)
  temp.data<-cbind(true.t,c(rexp(dim(true.t)[1]/2,lambdaC),rexp(dim(true.t)[1]/2,lambdaC*exp(-etaC))))

  t.obs<-apply(temp.data,1,min)
  delH<-rep(0,dim(true.t)[1])
  delD<-rep(0,dim(true.t)[1])
  delH[temp.data[,1]==t.obs]<-1
  delD[temp.data[,2]==t.obs]<-1

  my.data<-cbind(temp.data,t.obs,delH,delD,rep(0:1,each=dim(true.t)[1]/2))
  y1<-rep(0,dim(true.t)[1])
  y2<-rep(0,dim(true.t)[1])

  my.data.f<-data.frame(cbind(my.data,y1,y2))
  names(my.data.f)<-c("t1","t2","c","t.obs","delH","delD","group","y1","y2")

  indi.1<-(my.data.f$c < my.data.f$t1) & (my.data.f$c < my.data.f$t2)
  my.data.f$y1[indi.1]<-my.data.f$c[indi.1]
  my.data.f$y2[indi.1]<-my.data.f$c[indi.1]

  indi.2<-(my.data.f$t2 < my.data.f$t1) & (my.data.f$t1 < my.data.f$c)
  indi.21<-(my.data.f$t2 < my.data.f$c) & (my.data.f$c < my.data.f$t1)
  my.data.f$y1[indi.2 | indi.21]<-my.data.f$t2[indi.2| indi.21]
  my.data.f$y2[indi.2| indi.21]<-my.data.f$t2[indi.2| indi.21]

  indi.3<-(my.data.f$t1 < my.data.f$c) & (my.data.f$c < my.data.f$t2)
  my.data.f$y1[indi.3]<-my.data.f$t1[indi.3]
  my.data.f$y2[indi.3]<-my.data.f$c[indi.3]

  indi.4<-(my.data.f$t1 < my.data.f$t2) & (my.data.f$t2 < my.data.f$c)
  my.data.f$y1[indi.4]<-my.data.f$t1[indi.4]
  my.data.f$y2[indi.4]<-my.data.f$t2[indi.4]

  my.data.f$delD[indi.4]<-1

  # add cluster information in the data set
  my.data.f$cluster<-rep(1:(2*n.clust),each=n.sub/n.clust)

  names(my.data.f)<-c("time_Non_Fatal","time_Fatal","time_censor","t.obs","delH","delD","treatment","y1","y2","cluster")

  return(my.data.f)

}




####################################################################################################
## ----  Generate Individual Semi-Competing Risk Data ----- ##
##
## 
## Description:
## (1) use Joint distribution method. Bivariate exponential distribution with Gumbel-Hougaard copula
## (2) generate from marginal and conditional distributions
## (3) it allows one covariate: treatment. One can modify this function to include more covariates.
## (4) assumes the means for MVN are 0.
##
## Input:
##
## n:number of subjects per group 
## n.cluster: number of clusters per group
## dim: number of observations per subject (=2 in our case)
## alpha: association parameter from positive stable copula
## lambdaH: exponential rate of nonterminal events
## lambdaD: exponential rate of terminal events\
## lambdaC: exponential rate of censoring
## etaH: hazard ratio in nonterminal events between groups
## etaD: hazard ratio in terminal events between groups
## etaC: hazard ratio in censoring between groups
##
##
#####################################################################################################

gumbel.jeong_individual<-function(n,n.clust,dim,alpha,lambdaH,lambdaD,etaH,etaD)
{   
  exprand <- matrix(rexp(dim * n), c(n, dim))
  unifpirand <- runif(n, 0, pi)
  exprand2 <- rexp(n)
  beta <- 1/alpha
  stablerand <- sin((1 - beta) * unifpirand)^((1 - beta)/beta) * 
    (sin(beta * unifpirand))/(sin(unifpirand))^(1/beta)
  stablerand <- stablerand/(exprand2^(alpha - 1))
  unifrand <- invphigumbel(exprand/stablerand, alpha) # generating bivariate uniform random variables for marginal survival funtions--(*)
  
  clust.siz<- n/n.clust
  frail <- rep(rep(1, n.clust),each=clust.siz)
  matrix(c(-log(unifrand[,1])/(frail*lambdaH*exp(-etaH)),-log(unifrand[,2])/(frail*lambdaD*exp(-etaD))),c(n,dim)) # inverting specific forms of survival functions in (*) to create 
  # true bivariate event times adjusted for event types and trt groups
}





gen_cluster_data_individual<-function(n.sub=30, n.clust=10, dim=2, alpha=2, lambdaH=0.1, lambdaD=0.08, lambdaC=0.01, etaH=0.2, etaD=0.5, etaC=0.1){
  
  
  group0<-gumbel.jeong_individual(n.sub,n.clust,dim,alpha,lambdaH,lambdaD,0,0)
  group1<-gumbel.jeong_individual(n.sub,n.clust,dim,alpha,lambdaH,lambdaD,etaH,etaD)
  
  true.t<-rbind(group0,group1)
  temp.data<-cbind(true.t,c(rexp(dim(true.t)[1]/2,lambdaC),rexp(dim(true.t)[1]/2,lambdaC*exp(-etaC))))
  
  t.obs<-apply(temp.data,1,min)
  delH<-rep(0,dim(true.t)[1])
  delD<-rep(0,dim(true.t)[1])
  delH[temp.data[,1]==t.obs]<-1
  delD[temp.data[,2]==t.obs]<-1
  
  my.data<-cbind(temp.data,t.obs,delH,delD,rep(0:1,each=dim(true.t)[1]/2))
  y1<-rep(0,dim(true.t)[1])
  y2<-rep(0,dim(true.t)[1])
  
  my.data.f<-data.frame(cbind(my.data,y1,y2))
  names(my.data.f)<-c("t1","t2","c","t.obs","delH","delD","group","y1","y2")
  
  indi.1<-(my.data.f$c < my.data.f$t1) & (my.data.f$c < my.data.f$t2)
  my.data.f$y1[indi.1]<-my.data.f$c[indi.1]
  my.data.f$y2[indi.1]<-my.data.f$c[indi.1]
  
  indi.2<-(my.data.f$t2 < my.data.f$t1) & (my.data.f$t1 < my.data.f$c)
  indi.21<-(my.data.f$t2 < my.data.f$c) & (my.data.f$c < my.data.f$t1)
  my.data.f$y1[indi.2 | indi.21]<-my.data.f$t2[indi.2| indi.21]
  my.data.f$y2[indi.2| indi.21]<-my.data.f$t2[indi.2| indi.21]
  
  indi.3<-(my.data.f$t1 < my.data.f$c) & (my.data.f$c < my.data.f$t2)
  my.data.f$y1[indi.3]<-my.data.f$t1[indi.3]
  my.data.f$y2[indi.3]<-my.data.f$c[indi.3]
  
  indi.4<-(my.data.f$t1 < my.data.f$t2) & (my.data.f$t2 < my.data.f$c)
  my.data.f$y1[indi.4]<-my.data.f$t1[indi.4]
  my.data.f$y2[indi.4]<-my.data.f$t2[indi.4]
  
  my.data.f$delD[indi.4]<-1
  
  # add cluster information in the data set
  my.data.f$cluster<-rep(1:(2*n.clust),each=n.sub/n.clust)
  
  names(my.data.f)<-c("time_Non_Fatal","time_Fatal","time_censor","t.obs","delH","delD","treatment","y1","y2","cluster")
  
  return(my.data.f)
  
}













####################################################################################################
## ----  Functions of calculating the true value of log(WR) ----- ##
##
## 
## Description:
## According to the formula from Luo's 2015 paper
##
## Input:
##
##
#####################################################################################################
# double_integral<-function(myfun){
#   
#   out<-integrate(Vectorize(function(y2) { 
#     integrate(Vectorize(function(y1) myfun(y1,y2)), 0, y2)$value
#   }), 0, Inf)
#   
#   return(out)
#   
# }
# 
# 
# 
# # ----------------- Individual setting
# # individual setting: double integration in the numerator
# Ind_doubleNu<-function(y1,y2){
#   loga<-(-((lambdaH*exp(-etaH)*y1)^alpha+(lambdaD*exp(-etaD)*y2)^alpha)^(1/alpha)-lambdaC*exp(-etaC)*y2)
#   logb<-(-((lambdaH*y1)^alpha+(lambdaD*y2)^alpha)^(1/alpha)-lambdaC*y2)
#   logc<- (log(1/alpha)+log(((lambdaH*y1)^alpha+(lambdaD*y2)^alpha)^(1/alpha-1))
#           +log(alpha*(lambdaH*y1)^(alpha-1)*lambdaH)
#   )
#   logd<-log(lambdaC*(1+exp(-etaC)))
# 
#   return(exp(loga + logb + logc + logd))
# }
# 
# 
# # individual setting: double integration in the denominator
# Ind_doubleDe<-function(y1,y2){
#   loga<-(-((lambdaH*exp(-etaH)*y1)^alpha+(lambdaD*exp(-etaD)*y2)^alpha)^(1/alpha)-lambdaC*exp(-etaC)*y2)
#   logb<-(-((lambdaH*y1)^alpha+(lambdaD*y2)^alpha)^(1/alpha)-lambdaC*y2)
#   logc<-(log(1/alpha)+log(((lambdaH*exp(-etaH)*y1)^alpha+(lambdaD*exp(-etaD)*y2)^alpha)^(1/alpha-1))
#          +log(alpha*(lambdaH*exp(-etaH)*y1)^(alpha-1)*lambdaH*exp(-etaH))
#   )
#   logd<-log(lambdaC*(1+exp(-etaC)))
# 
#   return(exp(loga + logb + logc + logd))
# }
# 
# 
# # individual setting: single integration in the numerator
# Ind_singleNu<-function(y2){
#   out<-exp(-lambdaD*exp(-etaD)*y2-lambdaC*exp(-etaC)*y2)*exp(-lambdaD*y2-lambdaC*y2)*lambdaD
#   return(out)
# }
# 
# 
# # individual setting: single integration in the denomenator
# Ind_singleDe<-function(y2){
#   out<-exp(-lambdaD*exp(-etaD)*y2-lambdaC*exp(-etaC)*y2)*exp(-lambdaD*y2-lambdaC*y2)*lambdaD*exp(-etaD)
#   return(out)
# }


# 
# # ----------------- Individual setting with covariates
# # individual setting: double integration in the numerator 
# Ind_doubleNu<-function(y1,y2,z10,z11,z20,z21){
#   loga<-(-((lambdaH*exp(-etaH+t(outcome.H.eta)%*%c(z11,z21))*y1)^alpha+(lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))*y2)^alpha)^(1/alpha)-lambdaC*exp(-etaC)*y2)
#   logb<-(-((lambdaH*exp(t(outcome.H.eta)%*%c(z10,z20))*y1)^alpha+(lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))*y2)^alpha)^(1/alpha)-lambdaC*y2)
#   logc<- (log(1/alpha)+log(((lambdaH*exp(t(outcome.H.eta)%*%c(z10,z20))*y1)^alpha+(lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))*y2)^alpha)^(1/alpha-1))
#           +log(alpha*(lambdaH*exp(t(outcome.H.eta)%*%c(z10,z20))*y1)^(alpha-1)*lambdaH*exp(t(outcome.H.eta)%*%c(z10,z20)))
#   )
#   logd<-log(lambdaC*(1+exp(-etaC)))
#   loge<-log((1/4)*((1/sqrt(2*pi))*exp(-0.5*z10^2))*((1/sqrt(2*pi))*exp(-0.5*z11^2)))
#   
#   return(exp(loga + logb + logc + logd + loge))
# }
# 
# 
# # individual setting: double integration in the denominator
# Ind_doubleDe<-function(y1,y2,z10,z11,z20,z21){
#   loga<-(-((lambdaH*exp(-etaH+t(outcome.H.eta)%*%c(z11,z21))*y1)^alpha+(lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))*y2)^alpha)^(1/alpha)-lambdaC*exp(-etaC)*y2)
#   logb<-(-((lambdaH*exp(t(outcome.H.eta)%*%c(z10,z20))*y1)^alpha+(lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))*y2)^alpha)^(1/alpha)-lambdaC*y2)
#   logc<-(log(1/alpha)+log(((lambdaH*exp(-etaH+t(outcome.H.eta)%*%c(z11,z21))*y1)^alpha+(lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))*y2)^alpha)^(1/alpha-1))
#          +log(alpha*(lambdaH*exp(-etaH+t(outcome.H.eta)%*%c(z11,z21))*y1)^(alpha-1)*lambdaH*exp(-etaH+t(outcome.H.eta)%*%c(z11,z21)))
#   )
#   logd<-log(lambdaC*(1+exp(-etaC)))
#   loge<-log((1/4)*((1/sqrt(2*pi))*exp(-0.5*z10^2))*((1/sqrt(2*pi))*exp(-0.5*z11^2)))
#   
#   return(exp(loga + logb + logc + logd +loge))
# }
# 
# 
# # individual setting: single integration in the numerator
# Ind_singleNu<-function(y2,z10,z11,z20,z21){
#   out<-(
#     exp(-lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))*y2-lambdaC*exp(-etaC)*y2)*exp(-lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))*y2-lambdaC*y2)*lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))
#     *(1/4)*((1/sqrt(2*pi))*exp(-0.5*z10^2))*((1/sqrt(2*pi))*exp(-0.5*z11^2))
#   )
#   return(out)
# }
# 
# 
# # individual setting: single integration in the denomenator
# Ind_singleDe<-function(y2,z10,z11,z20,z21){
#   out<-(
#     exp(-lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))*y2-lambdaC*exp(-etaC)*y2)*exp(-lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))*y2-lambdaC*y2)*lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))
#     *(1/4)*((1/sqrt(2*pi))*exp(-0.5*z10^2))*((1/sqrt(2*pi))*exp(-0.5*z11^2))
#   )
#   
#     return(out)
# }
# 






# # --------------------- cluster setting
# 
# # V2.1
# # cluster setting: single integrationin the numerator
# Clus_singleNu<-function(y2,g1,g2){
#   out<-exp(-g1*numda.D*exp(-beta.D)*y2-numda.C*exp(-beta.C)*y2)*exp(-g2*numda.D*y2-numda.C*y2)*g2*numda.D*((theta.true^theta.true)/gamma(theta.true)*g1^(theta.true-1)*exp(-theta.true*g1))*((theta.true^theta.true)/gamma(theta.true)*g2^(theta.true-1)*exp(-theta.true*g2))
#   return(out)
# }
# 
# # cluster setting: single integration in the denomenator
# Clus_singleDe<-function(y2,g1,g2){
#   out<-exp(-g1*numda.D*exp(-beta.D)*y2-numda.C*exp(-beta.C)*y2)*exp(-g2*numda.D*y2-numda.C*y2)*g2*numda.D*exp(-beta.D)*((theta.true^theta.true)/gamma(theta.true)*g1^(theta.true-1)*exp(-theta.true*g1))*((theta.true^theta.true)/gamma(theta.true)*g2^(theta.true-1)*exp(-theta.true*g2))
#   return(out)
# }
# 
# # cluster setting: double integration in the numerator
# Clus_doubleNu<-function(y1,y2,g1,g2){
#   loga<-(-((g1*numda.H*exp(-beta.H)*y1)^alpha.cor+(g1*numda.D*exp(-beta.D)*y2)^alpha.cor)^(1/alpha.cor)-numda.C*exp(-beta.C)*y2)
#   logb<-(-((g2*numda.H*y1)^alpha.cor+(g2*numda.D*y2)^alpha.cor)^(1/alpha.cor)-numda.C*y2)
#   logc<- (log(((g2*numda.H*y1)^alpha.cor+(g2*numda.D*y2)^alpha.cor)^(1/alpha.cor-1))
#           +log((g2*numda.H*y1)^(alpha.cor-1)*g2*numda.H)
#   )
#   logd<-log(numda.C*(1+exp(-beta.C)))
#   loge<-log((theta.true^theta.true)/gamma(theta.true)*g1^(theta.true-1)*exp(-theta.true*g1))+log(((theta.true^theta.true)/gamma(theta.true)*g2^(theta.true-1)*exp(-theta.true*g2)))
#   
#   return(exp(loga + logb + logc + logd + loge))
# }
# 
# # cluster setting: double integration in the denominator
# Clus_doubleDe<-function(y1,y2,g1,g2){
#   loga<-(-((g1*numda.H*exp(-beta.H)*y1)^alpha.cor+(g1*numda.D*exp(-beta.D)*y2)^alpha.cor)^(1/alpha.cor)-numda.C*exp(-beta.C)*y2)
#   logb<-(-((g2*numda.H*y1)^alpha.cor+(g2*numda.D*y2)^alpha.cor)^(1/alpha.cor)-numda.C*y2)
#   logc<-(log(((g1*numda.H*exp(-beta.H)*y1)^alpha.cor+(g1*numda.D*exp(-beta.D)*y2)^alpha.cor)^(1/alpha.cor-1))
#          +log((g1*numda.H*exp(-beta.H)*y1)^(alpha.cor-1)*g1*numda.H*exp(-beta.H))
#   )
#   logd<-log(numda.C*(1+exp(-beta.C)))
#   loge<-log((theta.true^theta.true)/gamma(theta.true)*g1^(theta.true-1)*exp(-theta.true*g1))+log(((theta.true^theta.true)/gamma(theta.true)*g2^(theta.true-1)*exp(-theta.true*g2)))
#   
#   
#   return(exp(loga + logb + logc + logd + loge))
# }
# 
















####################################################################################################
## ----  Pairwise comparison of comparing x2 < min(y2, y3, x3)  ----- ##
##
## Description:
## (1) x represent treatment group; y represent control group. 
## (2) Row comparison within group.
## (3) Compare each x2, x3 to all the y2, y3.
## (4) It returns a vector of x2 in group A compared to every element of min(y2, y3, x3) in group B.
##
## Input:
##
## xk   : x2 as a scalar.
## yk   : x3 as a list of values.
## miny : min(y2, y3) as a list of values.
##
#####################################################################################################


# com_less<-function(xk,yk,miny){
#   res<-sapply(miny, function(miny) xk<min(yk,miny))
#   return(res)
# }


com_less<-function(xk,yk,miny){
  res<-xk<pmin(yk, miny)
  return(res)
}




####################################################################################################
## ----  Pairwise comparison of comparing x2 >= min(y2, y3, x3), or x2 !< min(y2, y3, x3)  ----- ##
##
## Description:
## (1) x represent treatment group; y represent control group. 
## (2) Row comparison within group.
## (3) Compare each x2, x3 to all the y2, y3.
## (4) It returns a vector of x2 in group A compared to every element of min(y2, y3, x3) in group B.
##
## Input:
##
## xk   : x2 as a scalar.
## yk   : x3 as a list of values.
## miny : min(y2, y3) as a list of values.
##
#####################################################################################################

# 
# com_more<-function(xk,yk,miny){
#   res<-sapply(miny, function(miny) xk>=min(yk,miny))
#   return(res)
# }

com_more<-function(xk,yk,miny){
  res<-xk>=pmin(yk, miny)
  return(res)
}





#############################################################
#
#  Function of comparison of same x(y) with different y(x) from different clusters.
#
# Input:
# i   : vector of sequence from 1 to total number of obs in the fixed group x(y).
# v   : list of one vector of sequence from 1 to number of varied group y(x) - 1.
# win : list of one data frame of all pairwise comparisons. Each column of the data is one obs from the fixed group x(y).
#       The last column of the data is cluster ID. 
#
# Return: a list of comparison matrixes.
#
#############################################################
# x12GTy12<-function(i, v, win){
# 
#   a<-lapply(v, function(v) win[which(win$cluster==v),i])
#   b<-lapply(v, function(v) win[which(win$cluster>v),i])
#   c1<-mapply(function(x,y) sum(outer(x,y,FUN="*")), x=a, y=b)
# 
#   return(c1)
# }

x12GTy12<-function(i, v, nj, win){
  
  w<-win[,i]
  a<-sum(mapply(function(v, nj) sum(outer(w[c(nj:v)], w[-c(nj:v)],FUN="*")), v=v, nj=nj, SIMPLIFY = T))
  
  return(a)
}






#############################################################
#
#  Function of comparison of different x(y) from the same cluster with different y(x) from different clusters.
#
# Input:
# ncluster: number of y(x) clusters. 
# Nj      : vector of number of obs in each x(y) cluster.
# win     : list of one data frame of all pairwise comparisons. Each column of the data is one obs from the fixed group x(y).
#           The last column of the data is cluster ID. 
#
# Return  :  sum of all the pairwise comparison.
#
#############################################################

# xi12GTy12<-function(ncluster, Nj, win){
# 
#   v.row<-c(1:ncluster)
#   v.col<-setdiff(c(1:sum(Nj)), cumsum(Nj))
#   r<-rep(cumsum(Nj), sapply(Nj, function(x) x-1))
#   counter<-unlist(sapply(Nj, function(x) c((x-1):1)))
#   index<-rep(c(1:length(v.col)), counter)
# 
#   a<-sapply(v.col, function(y) lapply(v.row, function(x) win[which(win$cluster==x), y]))
#   cd<-mapply(function(x,y) c((x+1):y), x=v.col, y=r)
# 
#   rm(v.col, r, counter)
# 
#   b<-lapply(cd, function(x)
#     lapply(x, function(z)
#       lapply(v.row, function(h) win[which(win$cluster!=h),z])))
# 
#   rm(cd)
# 
#   out<-sum(
#     unlist(
#       mapply(function(x,y)
#         mapply(function(z,h) sum(outer(unlist(a[z,x]),h,FUN="*")), z=v.row, h=y),
#         x=index, y=unlist(b, recursive=F), SIMPLIFY=F)
#     )
#   )
# 
#   return(out)
# 
# }

xi12GTy12<-function(i, sf, ef, sv, ev, nv, win){
  
  w<-win[,c(sf[i]:ef[i])]
  
  out<-sum(
    sapply(c(1:dim(w)[2]), function(h){
     sum(sapply(c(1:nv), function (j) sum(outer(w[c(sv[j]:ev[j]),h], w[-c(sv[j]:ev[j]),-h], FUN = "*"))))
  }))
 
 return(out) 
}



#############################################################
#
#  Function of comparison of same x(y) with different y(x) from different clusters between two kernels.
#
# Input:
# i      : vector of sequence from 1 to total number of obs in the fixed group x(y).
# v      : list of one vector of sequence from 1 to number of varied group y(x).
# winfix : list of one data frame of all pairwise comparisons for kernel1. Each column of the data is one obs from the fixed group x(y).
#          The last column of the data is cluster ID. 
# winvary: list of one data frame of all pairwise comparisons for kernel2. Each column of the data is one obs from the fixed group x(y).
#          The last column of the data is cluster ID. 
#
# Return: a list of comparison matrixes.
#
#############################################################
# covx12GTy12<-function(i, v, winfix, winvary){
#   
#   a<-lapply(v, function(v) winfix[which(winfix$cluster==v),i])
#   b<-lapply(v, function(v) winvary[which(winvary$cluster!=v),i])
#   c1<-mapply(function(x,y) sum(outer(x,y,FUN="*")), x=a, y=b)
#   
#   return(c1)
# }


covx12GTy12<-function(i, v, nj, winfix, winvary){
  
  wf<-winfix[,i]
  wv<-winvary[,i]
  a<-sum(mapply(function(v, nj) sum(outer(wf[c(nj:v)], wv[-c(nj:v)],FUN="*")), v=v, nj=nj, SIMPLIFY = T))
  
  return(a)
}




#############################################################
#
#  Function of comparison of different x(y) from the same cluster with different y(x) from different clusters between kernels.
#
# Input:
# ncluster: number of y(x) clusters. 
# Nj      : vector of number of obs in each x(y) cluster.
# win     : list of one data frame of all pairwise comparisons. Each column of the data is one obs from the fixed group x(y).
#           The last column of the data is cluster ID. 
#
# Return  :  sum of all the pairwise comparison.
#
#############################################################

# covxi12GTy12<-function(ncluster, Nj, winfix, winvary){
#   
#   v.row<-c(1:ncluster)
#   v.col<-c(1:sum(Nj))
#   r<-rep(cumsum(Nj), Nj)
#   
#   a<-sapply(v.col, function(y) lapply(v.row, function(x) winfix[which(winfix$cluster==x), y]))
#   
#   index1<-c(1:sum(Nj))
#   index2<-rep(mapply(function(x,y) x-y+1, x=cumsum(Nj), y=Nj), Nj)
#   
#   c1<-mapply(function(x,y,z) setdiff(c(z:y),x), x=index1, y=r, z=index2)
#   
#   rm(v.col, r, index1, index2)
#   
#   b<-lapply(c1, function(x) 
#     lapply(x, function(z) 
#       lapply(v.row, function(h) winvary[which(winvary$cluster!=h),z])))
#   
#   rm(c1)
#   
#   
#   out<-sum(
#     unlist(
#       mapply(function(x,y)
#         lapply(y, function(k)
#           mapply(function(z,h) sum(outer(unlist(a[z,x]),h,FUN="*")), z=v.row, h=k)),
#         x=c(1:sum(Nj)), y=b, SIMPLIFY=F)
#     )
#   )
#   
#   return(out)
#   
# }


covxi12GTy12<-function(i, sf, ef, winf, sv, ev, nv, winv){
  
  
  wf<-winf[,c(sf[i]:ef[i])]
  wv<-winv[,c(sf[i]:ef[i])]
  
  out<-sum(
    sapply(c(1:dim(wf)[2]), function(h){
      sum(sapply(c(1:nv), function (j) sum(outer(wf[c(sv[j]:ev[j]),h], wv[-c(sv[j]:ev[j]),-h], FUN = "*"))))
    }))
  
  return(out) 
  
  
  
}










#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  Causal Inference ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#



####################################################################################################
##  ------  Function of finding true causal WR
##
## Input:
## scen: scenario number
## 
## Output:
##  1. etaH, etaD values
####################################################################################################
etaHD<-function(scen){
  etaH<-etaD<-c(0.1, 0.2, 0.3, 0.4, 0.5)
  
  etaH<-rep(etaH, each=5)
  etaD<-rep(etaD, 5)
  
  etaHD<-matrix(cbind(etaH, etaD), ncol=2)
  
  out.H<-etaHD[scen,1]
  out.D<-etaHD[scen,2]
  out<-list(out.H, out.D)
  
  return(out)
  
}






####################################################################################################
##  ------  Simulation scenarios for causal inference 
##
## Input:
## scen: scenario number
## 
## Output:
##  1. etaH
##  2. etaD
##  3. ps.model.phi
####################################################################################################
scenario<-function(scen=1){
  
  #--------------------- etaH and etaD values
  if(scen%in%c(1,5)){
    etaH=0
    etaD=0
    ind_logWR=0
  }else if(scen%in%c(2,6)){
    etaH=0.1
    etaD=0.1
    ind_logWR=0.098061
  }else if(scen%in%c(3,7)){
    etaH=0.3
    etaD=0.2
    ind_logWR=0.223361
  }else if(scen%in%c(4,8)){
    etaH=0.2
    etaD=0.4
    ind_logWR=0.328163
  }
  
  
  #--------------------- PS model: intercept, z1, z2
  if(scen%in%c(1,2,3,4)){
    # correct PS model
    ps.model.phi=c(-0.25, 0.5, 0.5) # this intercept gives equal % trt assignments.
  }else if(scen%in%c(5,6,7,8)){
    # correct PS model
    ps.model.phi=c(-1.15, 0.5, 0.5) # this intercept gives 30% trt assignment
  }
  
  output<-list(etaH, etaD, ind_logWR, ps.model.phi)
  names(output)<-c("etaH", "etaD", "true.logWR", "ps.model.phi")
  
  return(output)
}







####################################################################################################
## ---- Causal Inference:
##      Generate Individual Semi-Competing Risk Data with covariates ----- ##
##
## 
## Description:
## (1) use Joint distribution method. Bivariate exponential distribution with Gumbel-Hougaard copula
## (2) generate from marginal and conditional distributions
## (3) it allows one covariate: treatment. One can modify this function to include more covariates.
## (4) assumes the means for MVN are 0.
##
## Input:
##
## N:number of total subjects 
## n.sub: number of subjects per group (treatment or control)
## dim: number of observations per subject (=2 in our case)
## alpha: association parameter from positive stable copula
## lambdaH: exponential rate of nonterminal events
## lambdaD: exponential rate of terminal events\
## lambdaC: exponential rate of censoring
## etaH: hazard ratio in nonterminal events between groups
## etaD: hazard ratio in terminal events between groups
## etaC: hazard ratio in censoring between groups
##
##
#####################################################################################################

gumbel_causal_individual<-function(n.sub,dim,alpha,lambdaH,lambdaD,etaH,etaD,cov.H, cov.D)
{   
  exprand <- matrix(rexp(dim * n.sub), c(n.sub, dim))
  unifpirand <- runif(n.sub, 0, pi)
  exprand2 <- rexp(n.sub)
  beta <- 1/alpha
  stablerand <- sin((1 - beta) * unifpirand)^((1 - beta)/beta) * 
    (sin(beta * unifpirand))/(sin(unifpirand))^(1/beta)
  stablerand <- stablerand/(exprand2^(alpha - 1))
  unifrand <- invphigumbel(exprand/stablerand, alpha) # generating bivariate uniform random variables for marginal survival funtions--(*)
  
  matrix(c(-log(unifrand[,1])/(lambdaH*exp(-etaH+cov.H)),-log(unifrand[,2])/(lambdaD*exp(-etaD+cov.D))),c(n.sub,dim)) # inverting specific forms of survival functions in (*) to create 
  # true bivariate event times adjusted for event types and trt groups
}





gen_causal_individual<-function(N=100, 
                                dim=2, 
                                alpha=2, 
                                lambdaH=0.1, 
                                lambdaD=0.08, 
                                lambdaC=0.01, 
                                etaH=0.2,
                                outcome.H.eta=c(0.1, 0.3),
                                etaD=0.5,
                                outcome.D.eta=c(0.2, 0.4), 
                                etaC=0.1, 
                                ps.model.phi=c(1, 0.5, 0.5)){
  
  # generate covariates z1 and z2
  z1<-rnorm(N)
  z2<-sample(c(1,0), N, prob=c(0.5, 0.5), replace=TRUE)
  
  # generate treatment status
  # design matrix X
  X<-cbind(rep(1,N), z1, z2)
  # calculate probability of getting treatment for each obs
  p<-exp(X%*%ps.model.phi)/(1+exp(X%*%ps.model.phi))
  trt<-sapply(p, function(x) sample(c(1,0), 1, prob=c(x, (1-x))))
  # number of treatment obs
  m<-sum(trt)
  # number of control obs
  n<-length(trt)-sum(trt)
  # percentage of treatment and control obs
  percent<-c("control%"=n/N, "trt%"=m/N)
  
  # combine data trt, z1, z2
  temp<-cbind(trt, z1, z2)
  # separate treatment temp data
  temp.trt<-temp[(temp[,1]==1),]
  temp.con<-temp[(temp[,1]==0),]
  
  # generate time to events: time to non-fatal and time to fatal events
  group0<-gumbel_causal_individual(n,dim,alpha,lambdaH,lambdaD,0,0,temp.con[,2:3]%*%outcome.H.eta, temp.con[,2:3]%*%outcome.D.eta)
  group1<-gumbel_causal_individual(m,dim,alpha,lambdaH,lambdaD,etaH,etaD, temp.trt[,2:3]%*%outcome.H.eta, temp.trt[,2:3]%*%outcome.D.eta)
  
  # combine time to events and time to censoring
  true.t<-rbind(group0,group1)
  temp.data<-cbind(true.t,c(rexp(n,lambdaC),rexp(m,lambdaC*exp(-etaC))))
  
  t.obs<-apply(temp.data,1,min)
  delH<-rep(0,dim(true.t)[1])
  delD<-rep(0,dim(true.t)[1])
  delH[temp.data[,1]==t.obs]<-1
  delD[temp.data[,2]==t.obs]<-1
  
  my.data<-cbind(temp.data,t.obs,delH,delD,rbind(temp.con, temp.trt))
  y1<-rep(0,n+m)
  y2<-rep(0,n+m)
  
  my.data.f<-data.frame(cbind(my.data,y1,y2))
  names(my.data.f)<-c("t1","t2","c","t.obs","delH","delD","group","z1","z2","y1","y2")
  
  indi.1<-(my.data.f$c < my.data.f$t1) & (my.data.f$c < my.data.f$t2)
  my.data.f$y1[indi.1]<-my.data.f$c[indi.1]
  my.data.f$y2[indi.1]<-my.data.f$c[indi.1]
  
  indi.2<-(my.data.f$t2 < my.data.f$t1) & (my.data.f$t1 < my.data.f$c)
  indi.21<-(my.data.f$t2 < my.data.f$c) & (my.data.f$c < my.data.f$t1)
  my.data.f$y1[indi.2 | indi.21]<-my.data.f$t2[indi.2| indi.21]
  my.data.f$y2[indi.2| indi.21]<-my.data.f$t2[indi.2| indi.21]
  
  indi.3<-(my.data.f$t1 < my.data.f$c) & (my.data.f$c < my.data.f$t2)
  my.data.f$y1[indi.3]<-my.data.f$t1[indi.3]
  my.data.f$y2[indi.3]<-my.data.f$c[indi.3]
  
  indi.4<-(my.data.f$t1 < my.data.f$t2) & (my.data.f$t2 < my.data.f$c)
  my.data.f$y1[indi.4]<-my.data.f$t1[indi.4]
  my.data.f$y2[indi.4]<-my.data.f$t2[indi.4]
  
  my.data.f$delD[indi.4]<-1
  
  
  names(my.data.f)<-c("time_Non_Fatal","time_Fatal","time_censor","t.obs","delH","delD","treatment","z1","z2","y1","y2")
  
  
  
  output<-list(my.data.f, n, m, percent)
  names(output)<-c("data", "#control", "#trt", "assignment%")
  return(output)
  
}





####################################################################################################
## ----  Functions of calculating the true value of log(WR) ----- ##
##
## 
## Description:
## According to the formula from Luo's 2015 paper
##
## Input:
##
##
#####################################################################################################
double_integral<-function(myfun){
  
  out<-integrate(Vectorize(function(y2) { 
    integrate(Vectorize(function(y1) myfun(y1,y2)), 0, y2)$value
  }), 0, Inf)
  
  return(out)
  
}



# # ----------------- Individual setting
# # individual setting: double integration in the numerator
# Ind_doubleNu<-function(y1,y2){
#   loga<-(-((lambdaH*exp(-etaH)*y1)^alpha+(lambdaD*exp(-etaD)*y2)^alpha)^(1/alpha)-lambdaC*exp(-etaC)*y2)
#   logb<-(-((lambdaH*y1)^alpha+(lambdaD*y2)^alpha)^(1/alpha)-lambdaC*y2)
#   logc<- (log(1/alpha)+log(((lambdaH*y1)^alpha+(lambdaD*y2)^alpha)^(1/alpha-1))
#           +log(alpha*(lambdaH*y1)^(alpha-1)*lambdaH)
#   )
#   logd<-log(lambdaC*(1+exp(-etaC)))
#   
#   return(exp(loga + logb + logc + logd))
# }
# 
# 
# # individual setting: double integration in the denominator
# Ind_doubleDe<-function(y1,y2){
#   loga<-(-((lambdaH*exp(-etaH)*y1)^alpha+(lambdaD*exp(-etaD)*y2)^alpha)^(1/alpha)-lambdaC*exp(-etaC)*y2)
#   logb<-(-((lambdaH*y1)^alpha+(lambdaD*y2)^alpha)^(1/alpha)-lambdaC*y2)
#   logc<-(log(1/alpha)+log(((lambdaH*exp(-etaH)*y1)^alpha+(lambdaD*exp(-etaD)*y2)^alpha)^(1/alpha-1))
#          +log(alpha*(lambdaH*exp(-etaH)*y1)^(alpha-1)*lambdaH*exp(-etaH))
#   )
#   logd<-log(lambdaC*(1+exp(-etaC)))
#   
#   return(exp(loga + logb + logc + logd))
# }
# 
# 
# # individual setting: single integration in the numerator
# Ind_singleNu<-function(y2){
#   out<-exp(-lambdaD*exp(-etaD)*y2-lambdaC*exp(-etaC)*y2)*exp(-lambdaD*y2-lambdaC*y2)*lambdaD
#   return(out)
# }
# 
# 
# # individual setting: single integration in the denomenator
# Ind_singleDe<-function(y2){
#   out<-exp(-lambdaD*exp(-etaD)*y2-lambdaC*exp(-etaC)*y2)*exp(-lambdaD*y2-lambdaC*y2)*lambdaD*exp(-etaD)
#   return(out)
# }



# ----------------- Individual setting with covariates
# individual setting: double integration in the numerator
Ind_doubleNu<-function(y1,y2,z10,z11,z20,z21){
  loga<-(-((lambdaH*exp(-etaH+t(outcome.H.eta)%*%c(z11,z21))*y1)^alpha+(lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))*y2)^alpha)^(1/alpha)-lambdaC*exp(-etaC)*y2)
  logb<-(-((lambdaH*exp(t(outcome.H.eta)%*%c(z10,z20))*y1)^alpha+(lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))*y2)^alpha)^(1/alpha)-lambdaC*y2)
  logc<- (log(1/alpha)+log(((lambdaH*exp(t(outcome.H.eta)%*%c(z10,z20))*y1)^alpha+(lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))*y2)^alpha)^(1/alpha-1))
          +log(alpha*(lambdaH*exp(t(outcome.H.eta)%*%c(z10,z20))*y1)^(alpha-1)*lambdaH*exp(t(outcome.H.eta)%*%c(z10,z20)))
  )
  logd<-log(lambdaC*(1+exp(-etaC)))
  loge<-log((1/4)*((1/sqrt(2*pi))*exp(-0.5*z10^2))*((1/sqrt(2*pi))*exp(-0.5*z11^2)))

  return(exp(loga + logb + logc + logd + loge))
}


# individual setting: double integration in the denominator
Ind_doubleDe<-function(y1,y2,z10,z11,z20,z21){
  loga<-(-((lambdaH*exp(-etaH+t(outcome.H.eta)%*%c(z11,z21))*y1)^alpha+(lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))*y2)^alpha)^(1/alpha)-lambdaC*exp(-etaC)*y2)
  logb<-(-((lambdaH*exp(t(outcome.H.eta)%*%c(z10,z20))*y1)^alpha+(lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))*y2)^alpha)^(1/alpha)-lambdaC*y2)
  logc<-(log(1/alpha)+log(((lambdaH*exp(-etaH+t(outcome.H.eta)%*%c(z11,z21))*y1)^alpha+(lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))*y2)^alpha)^(1/alpha-1))
         +log(alpha*(lambdaH*exp(-etaH+t(outcome.H.eta)%*%c(z11,z21))*y1)^(alpha-1)*lambdaH*exp(-etaH+t(outcome.H.eta)%*%c(z11,z21)))
  )
  logd<-log(lambdaC*(1+exp(-etaC)))
  loge<-log((1/4)*((1/sqrt(2*pi))*exp(-0.5*z10^2))*((1/sqrt(2*pi))*exp(-0.5*z11^2)))

  return(exp(loga + logb + logc + logd +loge))
}


# individual setting: single integration in the numerator
Ind_singleNu<-function(y2,z10,z11,z20,z21){
  out<-(
    exp(-lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))*y2-lambdaC*exp(-etaC)*y2)*exp(-lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))*y2-lambdaC*y2)*lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))
    *(1/4)*((1/sqrt(2*pi))*exp(-0.5*z10^2))*((1/sqrt(2*pi))*exp(-0.5*z11^2))
  )
  return(out)
}


# individual setting: single integration in the denomenator
Ind_singleDe<-function(y2,z10,z11,z20,z21){
  out<-(
    exp(-lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))*y2-lambdaC*exp(-etaC)*y2)*exp(-lambdaD*exp(t(outcome.D.eta)%*%c(z10,z20))*y2-lambdaC*y2)*lambdaD*exp(-etaD+t(outcome.D.eta)%*%c(z11,z21))
    *(1/4)*((1/sqrt(2*pi))*exp(-0.5*z10^2))*((1/sqrt(2*pi))*exp(-0.5*z11^2))
  )

    return(out)
}






#############################################################
#
#  Function of Calculating the An
#
# Input:
# phi:  a vector of the estimated parameter values of the PS model
# covMatrix: a matrix of all covariates. The first column is treatment status, which should be coded as 1: treatment, 0: control.
#
# Return  :  the evaluated partial derivative respect to PS model parameters. 
#           attr(,"gradient") stores evaluations (values) for each sample and each               parameter
#
#############################################################
# simple example for demonstration
# fx <- deriv(y ~ b0*x0 + b1*x1 , c("b0", "b1"),
#             function(b0, b1, x = 1:7){} )
# 
# An<-function(phi, x){
#   
#   xdim<-paste0("x",0:1)
#   for(i in 1:length(xdim)){
#     assign(xdim[i], x[,i], inherits=TRUE)
#   }
#   
#   xnam <- paste0("b", 0:1, "*x",0:1)
#   fmla <- as.formula(paste0("y~ ", paste0(xnam, collapse= "+")))
#   argu<-paste0("b", 0:1)
#   fx <- deriv(fmla , c(argu), func=TRUE)
#   
#   out<-parse(text=paste0("fx(",paste0(phi, collapse = ","),")"))
#   return(eval(out))
# }
# 
# 
# An(phi, x)
# 
# 


An<-function(phi, covMatrix, n, m, h){
  
  # separate treatment and control covariates matrix
  conMatrix<-covMatrix[which(covMatrix[,1]==0),-1]
  trtMatrix<-covMatrix[which(covMatrix[,1]==1),-1]
  
  # expand the treatment and control covariate matrix
  trt<-matrix(rep(t(trtMatrix),n), ncol=ncol(trtMatrix), byrow = TRUE)
  con<-matrix(rep(conMatrix, each=m), ncol=ncol(conMatrix))
  
  trt.dim<-paste0("trt",1:ncol(trt))
  for(i in 1:length(trt.dim)){
    assign(trt.dim[i], trt[,i], inherits=TRUE)
  }
  trt.nam<-paste0("phi",1:ncol(trt),"*trt",1:ncol(trt))
  
  con.dim<-paste0("con",1:ncol(con))
  for(i in 1:length(con.dim)){
    assign(con.dim[i], con[,i], inherits=TRUE)
  }   
  con.nam<-paste0("phi",1:ncol(con),"*con",1:ncol(con))
  
  fmla<-as.formula(paste0(
    "y~h/((exp(phi0+", 
    paste0(trt.nam, collapse = "+"),
    ")/(1+exp(phi0+",
    paste0(trt.nam, collapse = "+"),
    ")))*(1-(exp(phi0+",
    paste0(con.nam, collapse = "+"),
    ")/(1+exp(phi0+",
    paste0(con.nam, collapse = "+"),
    ")))))"
  ))
  
  argu<-c("phi0", paste0("phi", 1:ncol(trt)))
  fx <- deriv(fmla , c(argu), func=TRUE)
  
  out<-parse(text=paste0("fx(",paste0(phi, collapse = ","),")"))
  return(eval(out))
}




#############################################################
#
#  Function of calculating the influence function of the PS model: logistic regression: l_tilt
#
# Input:
# phi:  a vector of the estimated parameter values of the PS model
# covMatrix: a matrix of all covariates. All columns are covariates
#
# Return  :  
# - attr(, "gradient") contains all the values of first derivative of log likelihood for each sample.
# - attr(, "hessian") contains all the values of information matrix of log likelihood for each sample. It's 3-dimensional
#
#############################################################
l_tilt<-function(phi, covMatrix){
  
  cov.dim<-paste0("cov",1:ncol(covMatrix))
  for(i in 1:length(cov.dim)){
    assign(cov.dim[i], covMatrix[,i], inherits=TRUE)
  }
  cov.nam<-paste0("phi",1:ncol(covMatrix),"*cov",1:ncol(covMatrix))
  
  #"exp(phi0+X*phi)+yi*(phi0+X*phi)"
  fmla<-as.formula(paste0(
    "y~exp(phi0+",
    paste0(cov.nam, collapse = "+"),
    ")+yi*(phi0+",
    paste0(cov.nam, collapse = "+"),
    ")"
  ))
  
  argu<-c("phi0", paste0("phi", 1:ncol(covMatrix)))
  fx <- deriv(fmla , c(argu), func=TRUE, hessian = TRUE)
  
  out<-parse(text=paste0("fx(",paste0(phi, collapse = ","),")"))
  return(eval(out))
}




################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
######  Causual inference for observation cluster data
################################################################################################################################################################################################################################################################################################


####################################################################################################
##  ------  Simulation scenarios for causal inference for clustered data
##
## Input:
## scen: scenario number
## 
## Output:
##  1. etaH
##  2. etaD
##  3. ps.model.phi
####################################################################################################
scenario_cluster<-function(scen=1){
  
  #--------------------- Sample size
  if(scen%in%c(1:4)){
    # small # clusters, small cluster size
    # total # subjects: N
    N=200
    # total # clusters
    n.cluster=10
  }else if(scen%in%c(5:8)){
    # large # clusters, small cluster size
    # total # subjects: N
    N=1000
    # total # clusters
    n.cluster=50
  }else if(scen%in%c(9:12)){
    # small # clusters, large cluster size
    # total # subjects: N
    N=1000
    # total # clusters
    n.cluster=10
  }else if(scen%in%c(13:16)){
    # total # subjects: N
    N=1000
    # total # clusters
    n.cluster=20
  }
  
  
  
  #--------------------- etaH and etaD values
  if(scen%in%c(1:24)){
    etaH=0
    etaD=0
    ind_logWR=0
  }
  
  
  #--------------------- true PS model: intercept, z1, z2
  if(scen%in%c(1:2,5:6,9:10,13:14,17:18,21:22)){
    # equal % trt assignments
    ps.model.phi=c(-0.2, 0.5, 0.5)
  }else if (scen%in%c(3:4,7:8,11:12,15:16,19:20)){
    # unequal % trt assignments: trt - 30%
    ps.model.phi=c(-1.8, 0.5, 0.5)
  }
  
  
  #--------------------- ICC: cluster effect
  if(scen%in%c(1,3,5,7,9,11,13,15,17,19,21,23)){
    # small cluster effect
    # ICC=0.067, CV=0.258
    shape=15
    rate=15
  }else if(scen%in%c(2,4,6,8,10,12,14,16,18,20,22,24)){
    # large cluster effect
    # ICC=0.2, CV=0.447
    shape=5
    rate=5
  }
  
  
  #--------------------- outcome models
  if(scen%in%c(1:24)){
    # information to simulate time to events
    dim=2 # number of events (fixed)
    alpha=2 # correlation between two events (fixed)
    lambdaH=0.1  # baseline hazard for time to non-fatal event
    lambdaD=0.08  # baseline hazard for time to fatal event
    lambdaC=0.09  # baseline hazard for time to censor
    outcome.H.eta=c(0.1, 0.3) # covariate coefficient on the copula: time to non-fatal event: z1, z2, z1*z2
    outcome.D.eta=c(0.2, 0.4) # covariate coefficient on the copula: time to fatal event: z1, z2, z1*z2
    etaC=0.1  # parameter for treatment difference in censoring
  }
  
  
  
  
  
  output<-list(N, n.cluster, etaH, etaD, ind_logWR, ps.model.phi, shape, rate, dim, alpha, lambdaH, lambdaD, lambdaC, outcome.H.eta, outcome.D.eta, etaC)
  names(output)<-c("Total sample size", "Total # clusters","etaH", "etaD", "true.logWR", "ps.model.phi", "shape", "rate", "# events", "Correlation between events", "lambdaH", "lambdaD", "lambdaC", "outcome.H.eta", "outcome.D.eta", "etaC")
  
  return(output)
}








####################################################################################################
## ---- Causal Inference:
##      Generate cluster Semi-Competing Risk Data with covariates ----- ##
##
## 
## Description:
## (1) use Joint distribution method. Bivariate exponential distribution with Gumbel-Hougaard copula
## (2) generate from marginal and conditional distributions
## (3) it allows one covariate: treatment. One can modify this function to include more covariates.
## (4) assumes the means for MVN are 0.
##
## Input:
##
## N:number of total subjects 
## n.sub: number of subjects per group (treatment or control)
## dim: number of observations per subject (=2 in our case)
## alpha: association parameter from positive stable copula
## lambdaH: exponential rate of nonterminal events
## lambdaD: exponential rate of terminal events\
## lambdaC: exponential rate of censoring
## etaH: hazard ratio in nonterminal events between groups
## etaD: hazard ratio in terminal events between groups
## etaC: hazard ratio in censoring between groups
##
##
#####################################################################################################

gumbel_causal_PScluster<-function(n.sub,dim,alpha,lambdaH,lambdaD,etaH,etaD,cov.H, cov.D, frail)
{   
  exprand <- matrix(rexp(dim * n.sub), c(n.sub, dim))
  unifpirand <- runif(n.sub, 0, pi)
  exprand2 <- rexp(n.sub)
  beta <- 1/alpha
  stablerand <- sin((1 - beta) * unifpirand)^((1 - beta)/beta) * 
    (sin(beta * unifpirand))/(sin(unifpirand))^(1/beta)
  stablerand <- stablerand/(exprand2^(alpha - 1))
  unifrand <- invphigumbel(exprand/stablerand, alpha) # generating bivariate uniform random variables for marginal survival funtions--(*)
  
  matrix(c(-log(unifrand[,1])/(frail*lambdaH*exp(-etaH+cov.H)),-log(unifrand[,2])/(frail*lambdaD*exp(-etaD+cov.D))),c(n.sub,dim)) # inverting specific forms of survival functions in (*) to create 
  # true bivariate event times adjusted for event types and trt groups
}





gen_causal_PScluster<-function(N=100, 
                               n.cluster=50,
                               shape=5,
                               rate=5,
                                dim=2, 
                                alpha=2, 
                                lambdaH=0.1, 
                                lambdaD=0.08, 
                                lambdaC=0.01, 
                                etaH=0.2,
                                outcome.H.eta=c(0.1, 0.3, 0.2),
                                etaD=0.5,
                                outcome.D.eta=c(0.2, 0.4, 0.4), 
                                etaC=0.1, 
                                ps.model.phi=c(1, 0.5, 0.5, 0)){
  
  # generate covariates z1 and z2, z3=z1*z2
  z1<-rnorm(N)
  #z2<-sample(c(1,0), N, prob=c(0.5, 0.5), replace=TRUE)
  z2<-rnorm(N, mean=1, sd=2)
  #z3<-z1*z2
  #z3<-z1^2
  
  # generate clusters
  cluster<-rep(1:n.cluster, each=N/n.cluster)
  # # generate cluster effect : N(0,1)
  # reffect<-rep(rnorm(n.cluster, mean=0, sd=1), each=N/n.cluster)
  # r1<-reffect
  # r2<-reffect
  
  # # generate cluster effect: 
  # reffect<-rep(rgamma(n.cluster,shape=shape,rate=rate), each=N/n.cluster)
  
  # generate cluster effect: use copula to create joint distribution of normal and gamma
  # constructs an elliptical copula
  myCop <- normalCopula(param=c(0.4), dim=2, dispstr="un")
  # creates a multivariate distribution via copula
  myMvd <- mvdc(copula=myCop, margins=c("norm", "gamma"),
                paramMargins=list(list(mean=0,  sd=2),
                                  list(shape=shape, scale=1/rate))
  )
  reffect <- rMvdc(n.cluster, myMvd)
  r1<-rep(reffect[,1], each=N/n.cluster) # cluster effect for treatment
  r2<-rep(reffect[,2], each=N/n.cluster) # cluster effect for outcome

  
  
  # generate treatment status
  # design matrix X
  X<-cbind(rep(1,N), z1, z2)
  # calculate probability of getting treatment for each obs
  p<-exp(X%*%ps.model.phi+r1)/(1+exp(X%*%ps.model.phi+r1))
  trt<-sapply(p, function(x) sample(c(1,0), 1, prob=c(x, (1-x))))
  # number of treatment obs
  m<-sum(trt)
  # number of control obs
  n<-length(trt)-sum(trt)
  # percentage of treatment and control obs
  percent<-c("control%"=n/N, "trt%"=m/N)
  
  
  # combine data trt, z1, z2, z3, cluster ID, cluster effect
  temp<-cbind(trt, z1, z2, cluster, r1, r2)
  # separate treatment temp data
  temp.trt<-temp[(temp[,1]==1),]
  temp.con<-temp[(temp[,1]==0),]
  
  # generate time to events: time to non-fatal and time to fatal events
  group0<-gumbel_causal_PScluster(n,dim,alpha,lambdaH,lambdaD,0,0,
                                  temp.con[,2:3]%*%outcome.H.eta, 
                                  temp.con[,2:3]%*%outcome.D.eta, temp.con[,6])
  group1<-gumbel_causal_PScluster(m,dim,alpha,lambdaH,lambdaD,etaH,etaD, 
                                  temp.trt[,2:3]%*%outcome.H.eta, 
                                  temp.trt[,2:3]%*%outcome.D.eta,temp.trt[,6])
  
  # combine time to events and time to censoring
  true.t<-rbind(group0,group1)
  temp.data<-cbind(true.t,c(rexp(n,lambdaC),rexp(m,lambdaC*exp(-etaC))))
  
  t.obs<-apply(temp.data,1,min)
  delH<-rep(0,dim(true.t)[1])
  delD<-rep(0,dim(true.t)[1])
  delH[temp.data[,1]==t.obs]<-1
  delD[temp.data[,2]==t.obs]<-1
  
  my.data<-cbind(temp.data,t.obs,delH,delD,rbind(temp.con, temp.trt))
  y1<-rep(0,n+m)
  y2<-rep(0,n+m)
  
  my.data.f<-data.frame(cbind(my.data,y1,y2))
  names(my.data.f)<-c("t1","t2","c","t.obs","delta1","delta2","group","z1","z2","cluster","r1", "r2","y1","y2")
  
  indi.1<-(my.data.f$c < my.data.f$t1) & (my.data.f$c < my.data.f$t2)
  my.data.f$y1[indi.1]<-my.data.f$c[indi.1]
  my.data.f$y2[indi.1]<-my.data.f$c[indi.1]
  
  indi.2<-(my.data.f$t2 < my.data.f$t1) & (my.data.f$t1 < my.data.f$c)
  indi.21<-(my.data.f$t2 < my.data.f$c) & (my.data.f$c < my.data.f$t1)
  my.data.f$y1[indi.2 | indi.21]<-my.data.f$t2[indi.2| indi.21]
  my.data.f$y2[indi.2| indi.21]<-my.data.f$t2[indi.2| indi.21]
  
  indi.3<-(my.data.f$t1 < my.data.f$c) & (my.data.f$c < my.data.f$t2)
  my.data.f$y1[indi.3]<-my.data.f$t1[indi.3]
  my.data.f$y2[indi.3]<-my.data.f$c[indi.3]
  
  indi.4<-(my.data.f$t1 < my.data.f$t2) & (my.data.f$t2 < my.data.f$c)
  my.data.f$y1[indi.4]<-my.data.f$t1[indi.4]
  my.data.f$y2[indi.4]<-my.data.f$t2[indi.4]
  
  my.data.f$delD[indi.4]<-1
  
  
  names(my.data.f)<-c("time_Non_Fatal","time_Fatal","time_censor","t.obs","delta1","delta2","treatment","z1","z2","cluster","cluster.trt.effect","cluster.outcome.effect","y1","y2")
  
  
  
  output<-list(my.data.f, n, m, percent)
  names(output)<-c("data", "#control", "#trt", "assignment%")
  return(output)
  
}







###############################################################################################
#
#
#        Pairwise comparisons to compare semi-competing risk data using Win Ratio
#
#  Input:
#  data_c: control data set
#  data_t: treatment data set
#
#  Output:
#  win_t: comparison matrix for kernel 1: trt>control
#  win_c: comparison matrix for kernel 2: trt<control
#
# win_c: 1 if control individual wins within a pair; 0 otherwise.
# win_t: 1 if treatment individual wins within a pair; 0 otherwise.
# If both win_c and win_t are 0, then this pair comparison is indeterminate.
# All the outcome matrixs below are formatted as: First column as first element in control group pairwise comparison to all elements in treatment group. The matrix dimension is #ObsInTrt X #ObsInControl.
###############################################################################################
PairwiseComp<-function(data_c, data_t){
  
  # Treated individual wins:
  
  # First column is the results from comparing the first element in group A to every element in group B, etc.
  
  # @2@
  # For each obs in treatment group, compare time to non-fatal event and time to censor.
  min2_t<-list(pmin(data_t$time_censor, data_t$time_Non_Fatal))
  # All pairwise comparison of H0 < min(H1, C)
  comp2<-matrix(mapply(com_less, xk=data_c$time_Non_Fatal, yk=data_c$time_censor, miny=min2_t), ncol=length(data_c$time_Fatal))
  
  # @3@
  # For each obs in control group, compare time to fatal event and time to censor.
  min3_c<-list(pmin(data_c$time_censor, data_c$time_Fatal))
  # All pairwise comparison of T1 >= min(C, T0)
  comp3<-matrix(t(mapply(com_more, xk=data_t$time_Fatal, yk=data_t$time_censor, miny=min3_c)), ncol=length(data_c$time_Fatal))
  
  comp2<-comp2*comp3 # do this first to save memory
  rm(comp3)
  
  # @1@
  # For each obs in treatment group, compare time to fatal event and time to censor.
  min1_t<-list(pmin(data_t$time_censor, data_t$time_Fatal))
  # All pairwise comparison of T0 < min(C, T1).
  comp1<-matrix(mapply(com_less, xk=data_c$time_Fatal, yk=data_c$time_censor, miny=min1_t), ncol=length(data_c$time_Fatal))
  
  win_t<-comp1 + comp2  # create indicator that treated individual wins
  win_t[win_t==2]<-1  # some comparisons fall into both of the "or" scenarioes. So set 2 to 1.
  rm(comp1, comp2)  # remove objects to save memory
  
  
  
  # Control individual wins:
  
  # @5@
  # For each obs in control group, compare time to non-fatal and time to censor.
  min5_c<-list(pmin(data_c$time_censor, data_c$time_Non_Fatal))
  # All pairwise comparison of H1 < min(H0, C)
  comp5<-matrix(t(mapply(com_less, xk=data_t$time_Non_Fatal, yk=data_t$time_censor, miny=min5_c)), ncol=length(data_c$time_Fatal))
  
  # @6@
  # For each obs in treatment group, compare time to fatal and time to censor.
  min6_t<-min1_t
  # All pairwise comparison of T0 >= min(T1, C)
  comp6<-matrix(mapply(com_more, xk=data_c$time_Fatal, yk=data_c$time_censor, miny=min6_t), ncol=length(data_c$time_Fatal))
  
  comp5<-comp5*comp6
  rm(comp6)
  
  # @4@
  # For each obs in control group, compare time to fatal and time to censor.
  min4_c<-min3_c
  # All pairwise comparison of T1 < min(C, T0)
  comp4<-matrix(t(mapply(com_less, xk=data_t$time_Fatal, yk=data_t$time_censor, miny=min4_c)), ncol=length(data_c$time_Fatal))
  
  win_c<-comp4 + comp5
  win_c[win_c==2]<-1  # some comparisons fall into both of the "or" scenarioes. So set 2 to 1.
  rm(comp4, comp5, min1_t, min2_t, min3_c, min4_c, min5_c, min6_t, data_c, data_t)
  
  output<-list(win_t, win_c)
  return(output)
}







###############################################################################################
#
#                   Calibration method to calculate lambda
#
# Input:
#   lam: initial values of lambda
#   x.cov:   vector of names of covariates
# 
# Output:
#   f: function f is a vector of equations that have the same length as lam
# 
# This function is used for function "BBsolve" - numerical approach to solve lambda. lambda
#
###############################################################################################
calibrationPS_lambda<-function(lam, data, x.cov){
  
  lam1<-lam[1:(length(lam)/2)]
  lam2<-lam[(length(lam)/2+1):length(lam)]
  
  temp<-data
  
  # ni
  nii<-table(temp$cluster)
  ni<-c(rep(nii, nii))
  
  # Aij
  Aij<-c(temp$treatment)
  
  # dij
  dij<-c(temp$raw.ps)
  
  # X: dimension - n*p
  X<-as.matrix(temp[,x.cov])
  
  # alphaij: calculate expression of the new weight
  # calculate denominator sum
  exp1<-exp(X%*%lam1*Aij)
  de.sum1<-aggregate(Aij*dij*exp1, by=list(temp$cluster), FUN=sum)
  de.sum.trt<-rep(de.sum1[,2], nii)
  exp2<-exp(X%*%lam2*(1-Aij))
  de.sum2<-aggregate((1-Aij)*dij*exp2, by=list(temp$cluster), FUN=sum)
  de.sum.con<-rep(de.sum2[,2], nii)
  
  alphaij<-ni*Aij*dij*exp1/de.sum.trt+ni*(1-Aij)*dij*exp2/de.sum.con
  
  
  f<-rep(NA, 2*ncol(X))
  # for treatment group
  f[1:(length(lam)/2)]<-t(X)%*%(Aij*alphaij-rep(1,nrow(temp)))
  f[(length(lam)/2+1):length(lam)]<-t(X)%*%((1-Aij)*alphaij-rep(1,nrow(temp)))
  
  return(f)
  
}





###############################################################################################
#
#                   Calibration method to estimate new PS
#
# Input:
#   lam: initial values of lambda
#   x.cov:   vector of names of covariates
# 
# Output:
#   alphaij: a vector of new estimated weight
#   cal.ps: a vector of new estimated calibrated PS
# 
#
###############################################################################################
calibrationPS<-function(lam, data, x.cov){
  
  lam1<-lam[1:(length(lam)/2)]
  lam2<-lam[(length(lam)/2+1):length(lam)]
  
  temp<-data
  
  # ni
  nii<-table(temp$cluster)
  ni<-c(rep(nii, nii))
  
  # Aij
  Aij<-c(temp$treatment)
  
  # dij
  dij<-c(temp$raw.ps)
  
  # X: dimension - n*p
  X<-as.matrix(temp[,x.cov])
  
  # alphaij: calculate expression of the new weight
  # calculate denominator sum
  exp1<-exp(X%*%lam1*Aij)
  de.sum1<-aggregate(Aij*dij*exp1, by=list(temp$cluster), FUN=sum)
  de.sum.trt<-rep(de.sum1[,2], nii)
  exp2<-exp(X%*%lam2*(1-Aij))
  de.sum2<-aggregate((1-Aij)*dij*exp2, by=list(temp$cluster), FUN=sum)
  de.sum.con<-rep(de.sum2[,2], nii)
  
  alphaij<-ni*Aij*dij*exp1/de.sum.trt+ni*(1-Aij)*dij*exp2/de.sum.con
  
  eij<-((1/alphaij)^Aij)*(1-1/alphaij)^(1-Aij)
  
  out<-list(alphaij, eij)
  
  return(out)
  
}





######################################################################################
#
#  Function of pairwise comparisons using observed data, y1, y2
#
######################################################################################
PairwiseComp_Obs<-function(data_c, data_t){
  # --------- Treatment individual wins:
  # condition (1)
  cond1<-cond1(cat_c=data_c$cat, y2_c=data_c$y2, y2_t=data_t$y2)
  
  # each column is each obs in control group; each row is each obs in treatment group;
  
  # condition (3)
  # (3)_1
  cond3_1<-cond3_1(cat_c=data_c$cat, y2_c=data_c$y2, y2_t=data_t$y2, y1_c=data_c$y1, y1_t=data_t$y1)
  
  # (3)_2
  cond3_2<-cond3_2(cat_c=data_c$cat, cat_t=data_t$cat, y2_c=data_c$y2, y2_t=data_t$y2, y1_c=data_c$y1, y1_t=data_t$y1)
  
  cond3<-cond3_1+cond3_2
  
  win_t<-cond1+cond3
  
  
  # --------- Control individual wins:
  # condition (2)
  cond2<-cond2(cat_t=data_t$cat, y2_c=data_c$y2, y2_t=data_t$y2)
  
  # each column is each obs in control group; each row is each obs in treatment group;
  
  # condition (4)
  # (4)_1
  cond4_1<-cond4_1(cat_t=data_t$cat, y2_c=data_c$y2, y2_t=data_t$y2, y1_c=data_c$y1, y1_t=data_t$y1)
  
  # (4)_2
  cond4_2<-cond4_2(cat_c=data_c$cat, cat_t=data_t$cat, y2_c=data_c$y2, y2_t=data_t$y2, y1_c=data_c$y1, y1_t=data_t$y1)
  
  cond4<-cond4_1+cond4_2
  
  win_c<-cond2+cond4
  
  return(list(win_t, win_c))
}













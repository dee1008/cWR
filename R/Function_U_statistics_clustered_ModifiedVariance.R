## Di Zhang
## Started: 11/29/2016
## Last Updated: 10/19/2017






####################################################################################################
## ----  Generate Semi-Competing Risk Data ----- ##
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
double_integral<-function(myfun){
  
  out<-integrate(Vectorize(function(y2) { 
    integrate(Vectorize(function(y1) myfun(y1,y2)), 0, y2)$value
  }), 0, Inf)
  
  return(out)
  
}



# ----------------- Individual setting
# individual setting: double integration in the numerator 
Ind_doubleNu<-function(y1,y2){
  loga<-(-((numda.H*exp(-beta.H)*y1)^alpha.cor+(numda.D*exp(-beta.D)*y2)^alpha.cor)^(1/alpha.cor)-numda.C*exp(-beta.C)*y2)
  logb<-(-((numda.H*y1)^alpha.cor+(numda.D*y2)^alpha.cor)^(1/alpha.cor)-numda.C*y2)
  logc<- (log(1/alpha.cor)+log(((numda.H*y1)^alpha.cor+(numda.D*y2)^alpha.cor)^(1/alpha.cor-1))
          +log(alpha.cor*(numda.H*y1)^(alpha.cor-1)*numda.H)
  )
  logd<-log(numda.C*(1+exp(-beta.C)))
  
  return(exp(loga + logb + logc + logd))
}


# individual setting: double integration in the denominator
Ind_doubleDe<-function(y1,y2){
  loga<-(-((numda.H*exp(-beta.H)*y1)^alpha.cor+(numda.D*exp(-beta.D)*y2)^alpha.cor)^(1/alpha.cor)-numda.C*exp(-beta.C)*y2)
  logb<-(-((numda.H*y1)^alpha.cor+(numda.D*y2)^alpha.cor)^(1/alpha.cor)-numda.C*y2)
  logc<-(log(1/alpha.cor)+log(((numda.H*exp(-beta.H)*y1)^alpha.cor+(numda.D*exp(-beta.D)*y2)^alpha.cor)^(1/alpha.cor-1))
         +log(alpha.cor*(numda.H*exp(-beta.H)*y1)^(alpha.cor-1)*numda.H*exp(-beta.H))
  )
  logd<-log(numda.C*(1+exp(-beta.C)))
  
  return(exp(loga + logb + logc + logd))
}


# individual setting: single integration in the numerator
Ind_singleNu<-function(y2){
  out<-exp(-numda.D*exp(-beta.D)*y2-numda.C*exp(-beta.C)*y2)*exp(-numda.D*y2-numda.C*y2)*numda.D
  return(out)
}


# individual setting: single integration in the denomenator
Ind_singleDe<-function(y2){
  out<-exp(-numda.D*exp(-beta.D)*y2-numda.C*exp(-beta.C)*y2)*exp(-numda.D*y2-numda.C*y2)*numda.D*exp(-beta.D)
  return(out)
}


# --------------------- cluster setting

# V2.1
# cluster setting: single integrationin the numerator
Clus_singleNu<-function(y2,g1,g2){
  out<-exp(-g1*numda.D*exp(-beta.D)*y2-numda.C*exp(-beta.C)*y2)*exp(-g2*numda.D*y2-numda.C*y2)*g2*numda.D*((theta.true^theta.true)/gamma(theta.true)*g1^(theta.true-1)*exp(-theta.true*g1))*((theta.true^theta.true)/gamma(theta.true)*g2^(theta.true-1)*exp(-theta.true*g2))
  return(out)
}

# cluster setting: single integration in the denomenator
Clus_singleDe<-function(y2,g1,g2){
  out<-exp(-g1*numda.D*exp(-beta.D)*y2-numda.C*exp(-beta.C)*y2)*exp(-g2*numda.D*y2-numda.C*y2)*g2*numda.D*exp(-beta.D)*((theta.true^theta.true)/gamma(theta.true)*g1^(theta.true-1)*exp(-theta.true*g1))*((theta.true^theta.true)/gamma(theta.true)*g2^(theta.true-1)*exp(-theta.true*g2))
  return(out)
}

# cluster setting: double integration in the numerator
Clus_doubleNu<-function(y1,y2,g1,g2){
  loga<-(-((g1*numda.H*exp(-beta.H)*y1)^alpha.cor+(g1*numda.D*exp(-beta.D)*y2)^alpha.cor)^(1/alpha.cor)-numda.C*exp(-beta.C)*y2)
  logb<-(-((g2*numda.H*y1)^alpha.cor+(g2*numda.D*y2)^alpha.cor)^(1/alpha.cor)-numda.C*y2)
  logc<- (log(((g2*numda.H*y1)^alpha.cor+(g2*numda.D*y2)^alpha.cor)^(1/alpha.cor-1))
          +log((g2*numda.H*y1)^(alpha.cor-1)*g2*numda.H)
  )
  logd<-log(numda.C*(1+exp(-beta.C)))
  loge<-log((theta.true^theta.true)/gamma(theta.true)*g1^(theta.true-1)*exp(-theta.true*g1))+log(((theta.true^theta.true)/gamma(theta.true)*g2^(theta.true-1)*exp(-theta.true*g2)))
  
  return(exp(loga + logb + logc + logd + loge))
}

# cluster setting: double integration in the denominator
Clus_doubleDe<-function(y1,y2,g1,g2){
  loga<-(-((g1*numda.H*exp(-beta.H)*y1)^alpha.cor+(g1*numda.D*exp(-beta.D)*y2)^alpha.cor)^(1/alpha.cor)-numda.C*exp(-beta.C)*y2)
  logb<-(-((g2*numda.H*y1)^alpha.cor+(g2*numda.D*y2)^alpha.cor)^(1/alpha.cor)-numda.C*y2)
  logc<-(log(((g1*numda.H*exp(-beta.H)*y1)^alpha.cor+(g1*numda.D*exp(-beta.D)*y2)^alpha.cor)^(1/alpha.cor-1))
         +log((g1*numda.H*exp(-beta.H)*y1)^(alpha.cor-1)*g1*numda.H*exp(-beta.H))
  )
  logd<-log(numda.C*(1+exp(-beta.C)))
  loge<-log((theta.true^theta.true)/gamma(theta.true)*g1^(theta.true-1)*exp(-theta.true*g1))+log(((theta.true^theta.true)/gamma(theta.true)*g2^(theta.true-1)*exp(-theta.true*g2)))
  
  
  return(exp(loga + logb + logc + logd + loge))
}

















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







































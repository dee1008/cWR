
#####################################################################################
#
#   Apply data example for Win Ratio
#
#  Di Zhang
#  2/5/2018
#  Last modified: 3/21/2019
#
#####################################################################################


WR.CRT<-function(treatment, cluster, y1, y2, delta1, delta2, null.WR=1, alpha.sig=0.05){
  # 
# #---------------------------- Data analysis
# # input data set
# data<-read.csv(paste(path.in, "clusterData.csv", sep=""), header = T, stringsAsFactors = F)
# data$treatment[data$treatment==1]<-0
# data$treatment[data$treatment==2]<-1
# 
# 
# # define package input arguements
# treatment=data$treatment
# cluster=data$oinst  # if it's independent analysis without clustering, then set each individual a unique cluster ID. e.g. sequence number from 1 to n. 
# y1=data$Y1
# y2=data$Y2
# delta1=data$delta1
# delta2=data$delta2
# null.WR=1
# alpha.sig=0.05
# 



ll<-list(treatment, cluster, y1, y2, delta1, delta2)
# ---------------- warning message: 
# check if all vectors input have the same length
llen<-sapply(ll, function(x) length(x))
if(sum(llen)!=length(treatment)*length(llen)) stop("Error: Input arguements 'treatment','cluster','y1','y2','delta1', and 'delta2' have to be numeric vectors with the same length")
# check if all vectors input are numeric vectors
lna<-sapply(ll, function(x) is.numeric(x))
if(!(all(lna==T))) stop("Error: All input arguements have to be numeric vectors")
# check if treatment is coded as 0 or 1
if(!(all(treatment%in%c(0,1)))) stop("Error: Treatment status has to be coded as 0 (control) or 1 (treatment)")
# check if cluster ID is a positive interger vector
if(!((all(cluster == floor(cluster)))&(all(cluster>0))))stop("Error: Cluster ID has to be positive integer numbers")
# check if all event times are all non-negative
if(!all(y1>=0)) stop("Error: All time to non-fatal event have to be non-negative")
if(!all(y2>=0)) stop("Error: All time to fatal event have to be non-negative")
# check if delta1 and delta2 are coded as 0 and 1
if(!all(delta1%in%c(0,1))) stop("Error: Vector delta1 can only contain 0 or 1")
if(!all(delta2%in%c(0,1))) stop("Error: Vector delta2 can only contain 0 or 1")





# form a data set
data<-data.frame(treatment, cluster, y1, y2, delta1, delta2)

# rank data by treatment and then by cluster
data<-data[order(data$treatment, data$cluster),]

# obtain number of clusters per comparison group
n<-length(unique(data$cluster[data$treatment==0]))  # number of control clusters
m<-length(unique(data$cluster[data$treatment==1]))   # number of treatment clusters
J<-n+m

# Define some common quantities; 
Nj_c<-as.numeric(table(data$cluster[data$treatment==0]))
Nj_t<-as.numeric(table(data$cluster[data$treatment==1]))
Jm_bar<-mean(Nj_t)
Ln_bar<-mean(Nj_c)
Jm_bar2<-mean(sapply(Nj_t, function(x) x*(x-1)))
Ln_bar2<-mean(sapply(Nj_c, function(x) x*(x-1)))

xy12.den<-m*Jm_bar*((n^2)*(Ln_bar^2)-n*Ln_bar2-n*Ln_bar)
x12y.den<-n*Ln_bar*(m^2*Jm_bar^2-m*Jm_bar2-m*Jm_bar)
xi12GTy12.den<-m*Jm_bar2*(n^2*Ln_bar^2-n*Ln_bar2-n*Ln_bar)
x12GTyi12.den<-n*Ln_bar2*(m^2*Jm_bar^2-m*Jm_bar2-m*Jm_bar)



# Add comparison scenarios for each observation.
data$cat<-ifelse((data$delta1==1)&(data$delta2==1), 1, 
                 ifelse((data$delta1==0)&(data$delta2==1),2,
                        ifelse((data$delta1==0)&(data$delta2==0),3,
                               ifelse((data$delta1==1)&(data$delta2==0),4,NA))))




# (4) Set up criteria for X>Y, X<Y, and X><Y
# Pairwise comparison between all observations from treatment group and control group.


# Seperate patients in treatment group and control group into two data sets.
data_c<-data[which(data$treatment==0),]
data_t<-data[which(data$treatment==1),]


# (5) compare control and treatment observations to determine win or loss.

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


# (3) Calculate U1, U2.
U1<-sum(win_t) / (m*n)
U2<-sum(win_c) / (m*n)

# (4) Calculate theta1_hat, theta2_hat.
theta1<-mean(win_t)
theta2<-mean(win_c)

# (5) Calculate Var(U1), Var(U2);

# Kernel 1----------------------------------------------------------;

### Calculate Var(phi1(x));
# Comparison of same x with different y from different clusters.
xGTy12.out<-x12GTy12_Rcpp_test(t(win_t),m, n, Nj_t, Nj_c)
var.phi1x<-(xGTy12.out/xy12.den) - theta1^2
rm(xGTy12.out)


### Calculate Var(phi1(y));
# Comparison of same y with different x from different clusters.
x12GTy.out<-x12GTy12_Rcpp_test(win_t, n, m, Nj_c, Nj_t)
var.phi1y<-(x12GTy.out/x12y.den) - theta1^2
rm(x12GTy.out)


### Calculate Cov(phi1(x1), phi1(x2));
# Comparison of different x from the same cluster with different y from different clusters.
if (sum(Nj_t)==length(Nj_t)){
  cov.phi1xphix2<-0
}else{
  xi12GTy12.out<-xi12GTy12_Rcpp_test(t(win_t),m, n, Nj_t, Nj_c)
  cov.phi1xphix2<- (xi12GTy12.out/xi12GTy12.den) - theta1^2
  rm(xi12GTy12.out)
}

### Calculate Cov(phi1(y1), phi1(y2));
# Comparison of different y from the same cluster with different x from different clusters.
if (sum(Nj_t)==length(Nj_t)){
  cov.phi1yphiy2<-0
}else{
  x12GTyi12.out<-xi12GTy12_Rcpp_test(win_t, n, m, Nj_c, Nj_t)
  cov.phi1yphiy2<- (x12GTyi12.out/x12GTyi12.den) - theta1^2
  rm(x12GTyi12.out)
}

# Therefore, the Var(U1):
var.U1<-((Ln_bar^2)/m)*(Jm_bar*var.phi1x+Jm_bar2*cov.phi1xphix2) + ((Jm_bar^2)/n)*(Ln_bar*var.phi1y+Ln_bar2*cov.phi1yphiy2)



# For kernel 2----------------------------------------------------------;

### Calculate Var(phi2(x));
# Comparison of same x with different y from different clusters.
xGTy12.out2<-x12GTy12_Rcpp_test(t(win_c),m, n, Nj_t, Nj_c)
var.phi2x<-(xGTy12.out2/xy12.den) - theta2^2
rm(xGTy12.out2)


### Calculate Var(phi2(y));
# Comparison of same y with different x from different clusters.
x12GTy.out2<-x12GTy12_Rcpp_test(win_c, n, m, Nj_c, Nj_t)
var.phi2y<-(x12GTy.out2/x12y.den) - theta2^2
rm(x12GTy.out2)


### Calculate Cov(phi2(x1), phi2(x2));
# Comparison of different x from the same cluster with different y from different clusters.
if (sum(Nj_t)==length(Nj_t)){
  cov.phi1xphix22<-0
}else{
  xi12GTy12.out2<-xi12GTy12_Rcpp_test(t(win_c),m, n, Nj_t, Nj_c)
  cov.phi1xphix22<- (xi12GTy12.out2/xi12GTy12.den) - theta2^2
  rm(xi12GTy12.out2)
}

### Calculate Cov(phi2(y1), phi2(y2));
# Comparison of different y from the same cluster with different x from different clusters.
if (sum(Nj_t)==length(Nj_t)){
  cov.phi1yphiy22<-0
}else{
  x12GTyi12.out2<-xi12GTy12_Rcpp_test(win_c, n, m, Nj_c, Nj_t)
  cov.phi1yphiy22<- (x12GTyi12.out2/x12GTyi12.den) - theta2^2
  rm(x12GTyi12.out2)
}

# Therefore, the Var(U2):
var.U2<-((Ln_bar^2)/m)*(Jm_bar*var.phi2x+Jm_bar2*cov.phi1xphix22) + ((Jm_bar^2)/n)*(Ln_bar*var.phi2y+Ln_bar2*cov.phi1yphiy22)



# For Cov (U1, U2)----------------------------------------------------------;

### Calculate Cov(phi1(x), phi2(x));
# Comparison of same x with different y from different clusters between kernels.
covxGTy12.out<-covx12GTy12_Rcpp_test(t(win_t), t(win_c),m, n, Nj_t, Nj_c)
cov.phi1xphi2x<-(covxGTy12.out/(xy12.den)) - theta1*theta2
rm(covxGTy12.out)

### Calculate Cov(phi1(y), phi2(y));
# Comparison of same y with different x from different clusters between kernels.
covx12GTy.out<-covx12GTy12_Rcpp_test(win_c, win_t, n, m, Nj_c, Nj_t)
cov.phi1yphi2y<-(covx12GTy.out/(x12y.den)) - theta1*theta2
rm(covx12GTy.out)


### Calculate Cov(phi1(x1), phi2(x2));
if (sum(Nj_t)==length(Nj_t)){
  cov.phi1x1phi2x2<-0
}else{
  covxi12GTy12.out<-covxi12GTy12_Rcpp_test(t(win_t), t(win_c),m, n, Nj_t, Nj_c)
  cov.phi1x1phi2x2<-covxi12GTy12.out/(xi12GTy12.den) - theta1*theta2
  rm(covxi12GTy12.out)
}

### Calculate Cov(phi1(y1), phi2(y2));
if (sum(Nj_t)==length(Nj_t)){
  cov.phi1y1phi2y2<-0
}else{
  covx12GTyi12.out<-covxi12GTy12_Rcpp_test(win_t, win_c, n, m, Nj_c, Nj_t)
  cov.phi1y1phi2y2<-covx12GTyi12.out/(x12GTyi12.den) - theta1*theta2
  rm(covx12GTyi12.out)
}

# Therefore, the Cov(U1, U2):
cov.U1U2<-((Ln_bar^2)/m)*(Jm_bar*cov.phi1xphi2x+Jm_bar2*cov.phi1x1phi2x2)+((Jm_bar^2)/n)*(Ln_bar*cov.phi1yphi2y+Ln_bar2*cov.phi1y1phi2y2)


# (6) Use Delta Method to do hypothesis test under null.
test_stat<-(log(U1/U2)-log(null.WR)) / sqrt(var.U1/(U1^2) - 2*cov.U1U2/(U1*U2) + var.U2/(U2^2))


# store estimates
# estimated log(WR)
logWR<-log(U1/U2)
# estimated variance of log(WR) 
EV<-var.U1/(U1^2) - 2*cov.U1U2/(U1*U2) + var.U2/(U2^2)
se<-sqrt(EV)
# 95% CI
ci<-paste("(",round(logWR-1.96*se,3),",",round(logWR+1.96*se,3),")", sep = "")
# calculate P-value
p_value<-(1-pnorm(test_stat))*2
# variance covariance matrix of U1 and U2
var_cov<-matrix(c(var.U1, rep(cov.U1U2,2), var.U2),2,2)


if(sum(Nj_t)==length(Nj_t)){
  title<-"Test of two composite endpoints for independent survival data"
}else{
  title<-"Test of two composite endponits for clustered survival data"
}

output<-list(title, U1, U2, logWR, se, test_stat, ci, p_value, var_cov)
names(output)<-c("name","U1", "U2", "logWR", "se", "z", "ci", "p", "var_cov")

return(output)

}








WR.causal<-function(treatment, cluster, y1, y2, delta1, delta2, x.con=NA, x.char=NA, null.WR=1, alpha.sig=0.05, control=NA, n.boot=200)
{
  
  
  
  #----- Starting the function ... -------------------------------------------------
  
  library(fastDummies)
  library(lme4)
  library(reshape2)
  library(BB)
  library(nleqslv)
  library(MASS)
  
  
  
  ll<-list(treatment, cluster, y1, y2, delta1, delta2)
  
  # ---------------- Warning Message: ---------------------------------
  # check if all vectors input have the same length
  llen1<-sapply(ll, function(x) length(x))
  if(sum(llen1)!=length(treatment)*length(llen1)) stop("Error: Input arguements 'treatment','cluster','y1','y2','delta1', and 'delta2' have to be numeric vectors with the same length")
  # check if all vectors input are numeric vectors
  lna<-sapply(ll, function(x) is.numeric(x))
  if(!(all(lna==T))) stop("Error: Input arguements 'treatment','cluster','y1','y2','delta1', and 'delta2' have to be numeric vectors")
  # check if treatment is coded as 0 or 1
  if(!(all(treatment%in%c(0,1)))) stop("Error: Treatment status has to be coded as 0 (control) or 1 (treatment)")
  # check if cluster ID is a positive interger vector
  if(!((all(cluster == floor(cluster)))&(all(cluster>0))))stop("Error: Cluster ID has to be positive integer numbers")
  # check if all event times are all non-negative
  if(!all(y1>=0)) stop("Error: All time to non-fatal event have to be non-negative")
  if(!all(y2>=0)) stop("Error: All time to fatal event have to be non-negative")
  # check if delta1 and delta2 are coded as 0 and 1
  if(!all(delta1%in%c(0,1))) stop("Error: Vector delta1 can only contain 0 or 1")
  if(!all(delta2%in%c(0,1))) stop("Error: Vector delta2 can only contain 0 or 1")
  
  
  
  # at least one of the covariate matrix should be specified
  if((sum(is.na(x.con))!=0)&(sum(is.na(x.char))!=0)) stop("Error: At least one covariate matrix without any missing should be specified: x.con or x.char.")
  
  if(sum(is.na(x.char))==0){
    
    # check if the covariate matrix is type matrix
    if(!is.matrix(x.char)) stop("Error: Input argument 'x.char' must be a matrix")
    # check if the categorical covariate matrix has characteristic type in each column
    if(!all(apply(x.char, 2, function(x) is.character(x))==TRUE)) stop("Error: Input argument 'x.char' must be a character matrix")
    # check if all covariate matrix have the same number of rows
    llen<-c(llen1, nrow(x.char))
    if(sum(llen)!=length(treatment)*length(llen)) stop("Error: The row number of the input arguements 'x.char' has to be the same as the length of 'treatment','cluster','y1','y2','delta1' or 'delta2' ")
    
    
    if(sum(is.na(x.con))==0){
      
      # check if the covariate matrix is type matrix
      if(!is.matrix(x.con)) stop("Error: Input argument 'x.con' must be a matrix")
      # check if the continuous covariate matrix has numeric type in each column
      if(!all(apply(x.con, 2, function(x) is.numeric(x))==TRUE)) stop("Error: Input argument 'x.con' must be a numeric matrix")
      # check if all covariate matrix have the same number of rows
      llen<-c(llen1, nrow(x.con))
      if(sum(llen)!=length(treatment)*length(llen)) stop("Error: The row number of the input arguements 'x.con' has to be the same as the length of 'treatment','cluster','y1','y2','delta1' or 'delta2' ")
      
    }
  }else if(sum(is.na(x.con))==0){
    # check if the covariate matrix is type matrix
    if(!is.matrix(x.con)) stop("Error: Input argument 'x.con' must be a matrix")
    # check if the continuous covariate matrix has numeric type in each column
    if(!all(apply(x.con, 2, function(x) is.numeric(x))==TRUE)) stop("Error: Input argument 'x.con' must be a numeric matrix")
    # check if all covariate matrix have the same number of rows
    llen<-c(llen1, nrow(x.con))
    if(sum(llen)!=length(treatment)*length(llen)) stop("Error: The row number of the input arguements 'x.con' has to be the same as the length of 'treatment','cluster','y1','y2','delta1' or 'delta2' ")
  }
  
  
  
  
  #--------------------------------------------------------------------------
  
  
  
  #---------- Analysis ----------------
  if(sum(is.na(x.char))==0){
    # define covariate matrix column names
    x.char.name<-paste0("char", 1:ncol(x.char))
    colnames(x.char)<-x.char.name
    
    
    # convert character variables to factor variable
    for(yy in 1:length(x.char.name)){
      x.char[,x.char.name[yy]]<-as.factor(x.char[,x.char.name[yy]])
    }
    
    # create dummy variables for character variables
    dummy<-dummy_cols(x.char[,], remove_first_dummy = TRUE)
    dummy<-dummy[,-which(colnames(dummy)%in%x.char.name)]
    colnames(x.char)<-x.char.name<-colnames(dummy)
    x.char<-dummy
    
    if(sum(is.na(x.con))==0){
      data.cov<-data.frame(x.con, x.char)
      x.con.name<-paste0("con", 1:ncol(x.con))
      colnames(data.cov)<-c(x.con.name, x.char.name)
    }else{
      data.cov<-data.frame(x.char)
      colnames(data.cov)<-c(x.char.name)
    }
    
  }else if(sum(is.na(x.con))==0){
    data.cov<-data.frame(x.con)
    x.con.name<-paste0("con", 1:ncol(x.con))
    colnames(data.cov)<-c(x.con.name)
  }
  
  
  
  
  data<-data.frame(treatment, cluster, y1, y2, delta1, delta2, data.cov)
  
  # check missingness
  # delete any missing observations
  if(!all(complete.cases(data)==TRUE)){
    data<-data[complete.cases(data),]
    print(paste0("Warning: Missing observations are deleted, and resulting ", nrow(data), " observations in the analysis"))
  }
  
  
  
  # define basic quantities
  N<-nrow(data) # total # subjects
  list.cluster<-unique(data$cluster) # list of cluster ID
  n.cluster<-length(unique(data$cluster)) # total # clusters
  ind_logWR<-log(null.WR) # null log WR.
  n<-table(data$treatment)[1]  # # control patients
  m<-table(data$treatment)[2]  # # trt patients
  assign.percent<-m/N   # % trt assignment
  x=colnames(data.cov)
  
  
  
  # ---- Pairwise comparisons --------------
  # Add comparison scenarios for each observation.
  data$cat<-ifelse((data$delta1==1)&(data$delta2==1), 1, 
                   ifelse((data$delta1==0)&(data$delta2==1),2,
                          ifelse((data$delta1==0)&(data$delta2==0),3,
                                 ifelse((data$delta1==1)&(data$delta2==0),4,NA))))
  
  
  
  
  # (4) Set up criteria for X>Y, X<Y, and X><Y
  # Pairwise comparison between all observations from treatment group and control group.
  
  
  # Seperate patients in treatment group and control group into two data sets.
  data_c<-data[which(data$treatment==0),]
  data_t<-data[which(data$treatment==1),]
  
  
  # (5) compare control and treatment observations to determine win or loss.
  comp<-PairwiseComp_Obs(data_c=data_c, data_t=data_t)
  win_t<-comp[[1]]
  win_c<-comp[[2]]
  
  
  #-----------------------------------------------------------------------------------
  
  
  
  # analysis
  if(length(unique(cluster))==nrow(data)){
    
    title<-"Causal inference of two composite endpoints for independent survival data"
    ########## For independent subjects without clustering ######################
    
    
    # fit PS model
    # create formula for later use
    fmla<-formula(paste0("treatment~", paste0(x, collapse = "+")))
    # calculate the propensity score (correctly specified)
    fit<-glm(fmla, family=binomial(link='logit'), data=data)
    phi<-summary(fit)$coeff[,1]  # PS model coefficient: intercept, z1, z2
    dm<-as.matrix(cbind(rep(1, N),data[,x]))
    est.p<-exp(dm%*%phi)/(1+exp(dm%*%phi))
    data$ps<-est.p
    
    # separate treatment and control ps
    ps_c<-data[which(data$treatment==0),"ps"]
    ps_t<-data[which(data$treatment==1),"ps"]
    
    #----------------------- Point Estimation -------------------------------------
    # create matrix of control and PS
    ps_cMatrix<-matrix(rep(ps_c, m), nrow=m, byrow = T)
    ps_tMatrix<-matrix(rep(ps_t, n), nrow=m)
    
    # calculate IPW estimator U1 and U2
    # First kernel: phi_1(X>Y)
    U1_IPW1<-0.5*(1/choose(m+n,2))*sum(win_t/((1-ps_cMatrix)*ps_tMatrix))
    # Second kernel: phi_1(X<Y)
    U2_IPW1<-0.5*(1/choose(m+n,2))*sum(win_c/((1-ps_cMatrix)*ps_tMatrix))
    # log(WR)
    logWR<-log(U1_IPW1/U2_IPW1)
    
    
    #----------------------- Variance Estimation -------------------------------------
    
    # ----- Influence function of U1_IPW
    h0_Y<-apply(win_t, 2, mean)
    h1_Y<-apply(win_t, 1, mean)
    g1_IP<-c(h1_Y/ps_t,h0_Y/(1-ps_c))*0.5
    first<-2*(g1_IP-U1_IPW1)
    
    
    # pass in covariate matrix: first column - treatment indicator (1=trt, 0=control)
    covMatrix<-as.matrix(data[, c("treatment",x)])
    h<<-c(win_t)
    Ani<-An(phi, covMatrix, n, m, h)
    An.f<-apply(attr(Ani,"gradient"), 2, sum)/choose(m+n,2)
    
    # calculate the influence function of the PS model: logistic regression: l_tilt
    covMatrix<-as.matrix(data[, x])
    yi<<-data[, "treatment"]
    # ---- first derivative of log likelihood function
    li<-attr(l_tilt(phi, covMatrix), "gradient")
    # attr(, "gradient") contains all the values of first derivative of log likelihood for each sample
    # ---- information matrix
    ii<-attr(l_tilt(phi, covMatrix), "hessian")
    # attr(, "hessian") contains all the values of information matrix of log likelihood for each sample. It's 3-dimensional
    
    l_tilt.f<-sapply(1:(n+m), function(x) ginv(ii[x,,])%*%li[x,])
    
    second<-c(t(An.f)%*%l_tilt.f)
    
    # influence function of U1_IPW
    inf.U1_IPW<-first+second
    var.U1_IPW<-mean(inf.U1_IPW^2)/(n+m)
    
    
    
    
    # ----- Influence function of U2_IPW
    h0_Y<-apply(win_c, 2, mean)
    h1_Y<-apply(win_c, 1, mean)
    g1_IP<-c(h1_Y/ps_t,h0_Y/(1-ps_c))*0.5
    first<-2*(g1_IP-U1_IPW1)
    
    
    # pass in covariate matrix: first column - treatment indicator (1=trt, 0=control)
    covMatrix<-as.matrix(data[, c("treatment",x)])
    h<<-c(win_c)
    Ani<-An(phi, covMatrix, n, m, h)
    An.f<-apply(attr(Ani,"gradient"), 2, sum)/choose(m+n,2)
    
    # calculate the influence function of the PS model: logistic regression: l_tilt
    covMatrix<-as.matrix(data[,x])
    yi<<-data[, "treatment"]
    # ---- first derivative of log likelihood function
    li<-attr(l_tilt(phi, covMatrix), "gradient")
    # attr(, "gradient") contains all the values of first derivative of log likelihood for each sample
    # ---- information matrix
    ii<-attr(l_tilt(phi, covMatrix), "hessian")
    # attr(, "hessian") contains all the values of information matrix of log likelihood for each sample. It's 3-dimensional
    
    l_tilt.f<-sapply(1:(n+m), function(x) ginv(ii[x,,])%*%li[x,])
    
    second<-c(t(An.f)%*%l_tilt.f)
    
    # influence function of U2_IPW
    inf.U2_IPW<-first+second
    var.U2_IPW<-mean(inf.U2_IPW^2)/(n+m)
    
    ### Therefore, the estimated asymtotic variance of estimated log(WR) is
    se<-sqrt(mean((inf.U1_IPW/U1_IPW1-inf.U2_IPW/U2_IPW1)^2)/(n+m))
    
    test_stat<-(logWR-ind_logWR)/se
    
    U1<-U1_IPW1
    U2<-U2_IPW1
    
    lam.out<-NA
    # data info
    info<-paste0("Total ", N, " subjects. Control patients = ", n, " treatment patients = ", m)
    list.cluster<-NA
    
    
    
    
  }else{
    title<-"Causal inference of two composite endpoints for cluster-dependent survival data"
    ########## For cluster-dependent subjects with clustering ######################
    
    
    #--------------- exclude clusters that have only one membership
    ti<-table(data$treatment, data$cluster)
    tt<-t(ti)
    tt[tt==0]<-NA
    tt<-as.numeric(colnames(ti)[complete.cases(tt)])
    
    data3<-data[which(data$cluster%in%tt),]  # data3 is the data set with all clusters have both memberships
    
    N<-nrow(data3) # total # subjects
    list.cluster<-unique(data3$cluster) # list of selected cluster ID
    n.cluster<-length(unique(data3$cluster)) # total # clusters
    n<-table(data3$treatment)[1]  # # control patients
    m<-table(data3$treatment)[2]  # # trt patients
    
    # rank the data3 by cluster ID and by treatment
    data3<-data3[order(data3$cluster, data3$treatment),]
    
    # re-name cluster ID to make it have consecutive cluster ID
    for(j in 1:length(unique(data3$cluster))){
      data3[which(data3$cluster==unique(data3$cluster)[j]), "cluster"]<-j
    }
    
    
    #------------------- Estimate PS 
    # create formula for later use
    fmla<-formula(paste0("treatment~", paste0(x, collapse = "+")))
    # fit PS model, without cluster effect: standard logistic regression with iid r.v.
    fit<-glm(fmla, family=binomial(link='logit'), data=data3)
    phi<-summary(fit)$coeff[,1]  # PS model coefficient: intercept, z1, z2
    dm<-as.matrix(cbind(rep(1, nrow(data3)),data3[,x]))
    est.p<-exp(dm%*%phi)/(1+exp(dm%*%phi))
    data3$raw.ps<-est.p
    
    
    # estimate PS using calibration technique.
    fmla.boot<-formula(paste0("treatment~0+", paste0(x, collapse = "+")))
    fit.boot<-glm(fmla.boot, family=binomial(link='logit'), data=data3)
    lam1<-summary(fit.boot)$coeff[,1]
    lam<-rep(lam1, 2)
    #lam<-BBsolve(par=lam, data=data3, x=x, fn=calibrationPS_lambda, control = list(maxit=maxit, tol=tol))$par
    if(sum(is.na(control))){
      lam.out<-nleqslv(lam, fn=calibrationPS_lambda,data=data3, x.cov=x)
    }else{
      lam.out<-nleqslv(lam, control=control, fn=calibrationPS_lambda,data=data3, x.cov=x)
    }
    lam<-lam.out$x
    # print convergence status
    if(lam.out$termcd==1){
      print("Successful Convergence.")
    }else if(lam.out$termcd==2){
      print("Lambda values within tolerance.")
    }else{
      print("Unsuccessful Convergence. Check '$convergence' and R package 'nleqslv' for details.")
    }
    
    
    
    # calculate the new estimated PS
    cal<-calibrationPS(lam=lam, data=data3, x.cov=x)
    data3$cal.weight<-cal[[1]]
    data3$cal.ps<-cal[[2]]
    
    
    
    
    #--------------- Outcome model: stratification
    # subset data by clusters
    aa<-split.data.frame(data3, data3$cluster)
    
    # calculate cluster specific estimates
    U1<-rep(NA, length(aa))
    U2<-rep(NA, length(aa))
    tw<-rep(NA, length(aa))  # total weight in each cluster
    for(i in 1:length(aa)){
      # data set for cluster i
      temp<-aa[[i]]
      
      # re-calculate pairwise comparison of win_t and win_c
      temp_c<-temp[which(temp$treatment==0),]
      temp_t<-temp[which(temp$treatment==1),]
      
      comp<-PairwiseComp_Obs(data_c=temp_c, data_t=temp_t)
      temp.win_t<-comp[[1]]
      temp.win_c<-comp[[2]]
      
      # separate treatment and control weight
      w_c<-temp[which(temp$treatment==0),"cal.weight"]
      w_t<-temp[which(temp$treatment==1),"cal.weight"]
      
      temp.m<-length(w_t)
      temp.n<-length(w_c)
      
      #----------------------- Point Estimation -------------------------------------
      # create matrix of control and treatment weight
      w_cMatrix<-matrix(rep(w_c, temp.m), nrow=temp.m, byrow = T)
      w_tMatrix<-matrix(rep(w_t, temp.n), nrow=temp.m)
      
      # calculate IPW estimator U1 and U2
      # First kernel: phi_1(X>Y)
      U1[i]<-0.5*(1/choose(nrow(temp),2))*sum(temp.win_t*(w_cMatrix*w_tMatrix))
      # Second kernel: phi_1(X<Y)
      U2[i]<-0.5*(1/choose(nrow(temp),2))*sum(temp.win_c*(w_cMatrix*w_tMatrix))
      
      # calculate weight for each observation
      tw[i]<-sum(temp$cal.weight)
    }
    
    # calculate the clustered estimator
    U1_IPW5<-sum(tw*U1)/sum(tw)
    U2_IPW5<-sum(tw*U2)/sum(tw)
    logWR<-log(U1_IPW5/U2_IPW5)
    
    
    #--- Variance estimation using Bootstrap -----
    print("Estimating variance using bootstrap on clusters ...")
    boot.logWR5<-rep(NA, n.boot)
    for(b in 1:n.boot){
      ind.cluster<-sample(unique(data$cluster), n.cluster, replace=T)
      
      for(i in 1:length(ind.cluster)){
        temp<-data[data$cluster==ind.cluster[i],]
        temp$cluster<-i
        if(i==1){
          bootdata<-temp
        }else{
          bootdata<-rbind(bootdata, temp)
        }
      }
      
      dc<-bootdata[which(bootdata$treatment==0),]
      dt<-bootdata[which(bootdata$treatment==1),]
      boot.N<-nrow(bootdata)
      boot.m<-nrow(dt)
      boot.n<-nrow(dc)
      
      #--------------- exclude clusters that have only one membership
      ti<-table(bootdata$treatment, bootdata$cluster)
      tt<-t(ti)
      tt[tt==0]<-NA
      tt<-as.numeric(colnames(ti)[complete.cases(tt)])
      
      data3<-bootdata[which(bootdata$cluster%in%tt),]  # data3 is the data set with all clusters have both memberships
      
      # rank the data3 by cluster ID and by treatment
      data3<-data3[order(data3$cluster, data3$treatment),]
      
      # re-name cluster ID to make it have consecutive cluster ID
      for(j in 1:length(unique(data3$cluster))){
        data3[which(data3$cluster==unique(data3$cluster)[j]), "cluster"]<-j
      }
      
      
      #------------------- Estimate PS 
      # fit PS model, without cluster effect: standard logistic regression with iid r.v.
      fit<-glm(fmla, family=binomial(link='logit'), data=data3)
      phi<-summary(fit)$coeff[,1]  # PS model coefficient: intercept, z1, z2
      dm<-as.matrix(cbind(rep(1, nrow(data3)),data3[,x]))
      est.p<-exp(dm%*%phi)/(1+exp(dm%*%phi))
      data3$raw.ps<-est.p
      
      # calculate initial values separately in treatment and control group 
      fmla.boot<-formula(paste0("treatment~0+", paste0(x, collapse = "+")))
      fit.boot<-glm(fmla.boot, family=binomial(link='logit'), data=data3)
      lam1<-summary(fit.boot)$coeff[,1]
      lam<-rep(lam1, 2)
      #lam<-BBsolve(par=lam, data=data3, x=x, fn=calibrationPS_lambda, control = list(maxit=maxit, tol=tol))$par
      if(sum(is.na(control))){
        lam.out.boot<-nleqslv(lam, fn=calibrationPS_lambda,data=data3, x.cov=x)
      }else{
        lam.out.boot<-nleqslv(lam, control=control, fn=calibrationPS_lambda,data=data3, x.cov=x)
      }
      lam<-lam.out.boot$x
      # # print convergence status
      # if(lam.out$termcd==1){
      #   print("Successful Convergence.")
      # }else if(lam.out$termcd==2){
      #   print("Lambda values within tolerance.")
      # }else{
      #   print("Unsuccessful Convergence. Check '$convergence' and R package 'nleqslv' for details.")
      # }
      
      # calculate the new estimated PS
      cal<-calibrationPS(lam=lam, data=data3, x=x)
      data3$cal.weight<-cal[[1]]
      data3$cal.ps<-cal[[2]]
      
      
      
      
      #--------------- Outcome model: stratification
      # subset data by clusters
      aa<-split.data.frame(data3, data3$cluster)
      
      # calculate cluster specific estimates
      U1<-rep(NA, length(aa))
      U2<-rep(NA, length(aa))
      tw<-rep(NA, length(aa))  # total weight in each cluster
      for(i in 1:length(aa)){
        # data set for cluster i
        temp<-aa[[i]]
        
        # re-calculate pairwise comparison of win_t and win_c
        temp_c<-temp[which(temp$treatment==0),]
        temp_t<-temp[which(temp$treatment==1),]
        
        comp<-PairwiseComp_Obs(data_c=temp_c, data_t=temp_t)
        temp.win_t<-comp[[1]]
        temp.win_c<-comp[[2]]
        
        # separate treatment and control weight
        w_c<-temp[which(temp$treatment==0),"cal.weight"]
        w_t<-temp[which(temp$treatment==1),"cal.weight"]
        
        temp.m<-length(w_t)
        temp.n<-length(w_c)
        
        #----------------------- Point Estimation -------------------------------------
        # create matrix of control and treatment weight
        w_cMatrix<-matrix(rep(w_c, temp.m), nrow=temp.m, byrow = T)
        w_tMatrix<-matrix(rep(w_t, temp.n), nrow=temp.m)
        
        # calculate IPW estimator U1 and U2
        # First kernel: phi_1(X>Y)
        U1[i]<-0.5*(1/choose(nrow(temp),2))*sum(temp.win_t*(w_cMatrix*w_tMatrix))
        # Second kernel: phi_1(X<Y)
        U2[i]<-0.5*(1/choose(nrow(temp),2))*sum(temp.win_c*(w_cMatrix*w_tMatrix))
        
        # calculate weight for each observation
        tw[i]<-sum(temp$cal.weight)
      }
      
      # calculate the clustered estimator
      boot.U1<-sum(tw*U1)/sum(tw)
      boot.U2<-sum(tw*U2)/sum(tw)
      boot.logWR5[b]<-log(boot.U1/boot.U2)
      
      
    }
    
    # bootstrap SE
    se<-sqrt(var(boot.logWR5))
    
    test_stat<-(logWR-ind_logWR)/se
    
    indicator<-ifelse(abs(test_stat)>qnorm(1-alpha.sig/2), 1, 0)
    
    U1<-U1_IPW5
    U2<-U2_IPW5
    
    # data info
    info<-paste0("Total ", N, " subjects with ", n.cluster, " clusters. Total control patients = ", n, " total treatment patients = ", m)
    
    
    
    
  }
  
  
  # 95% CI
  ci<-paste("(",round(logWR-1.96*se,3),",",round(logWR+1.96*se,3),")", sep = "")
  # calculate P-value
  p_value<-(1-pnorm(abs(test_stat)))*2
  
  
  
  
  output<-list(title, info, list.cluster, U1, U2, logWR, se, test_stat, ci, p_value, lam.out)
  names(output)<-c("name","info", "clusters","U1", "U2", "logWR", "se", "z", "ci", "p-value", "convergence")
  
  return(output)
  
}











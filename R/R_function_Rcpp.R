
#####################################################################################
#
#   Apply data example for Win Ratio
#
#  Di Zhang
#  2/5/2018
#
#####################################################################################


cWR<-function(treatment, cluster, y1, y2, delta1, delta2, null.WR=1, alpha.sig=0.05){
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
if(sum(llen)!=length(treatment)*length(llen)) stop("Input arguements 'treatment','cluster','y1','y2','delta1', and 'delta2' have to be numeric vectors with the same length")
# check if all vectors input are numeric vectors
lna<-sapply(ll, function(x) is.numeric(x))
if(!(all(lna==T))) stop("All input arguements have to be numeric vectors")
# check if treatment is coded as 0 or 1
if(!(all(treatment%in%c(0,1)))) stop("Treatment status has to be coded as 0 (control) or 1 (treatment)")
# check if cluster ID is a positive interger vector
if(!((all(cluster == floor(cluster)))&(all(cluster>0))))stop("Cluster ID has to be positive integer numbers")
# check if all event times are all non-negative
if(!all(y1>=0)) stop("All time to non-fatal event have to be non-negative")
if(!all(y2>=0)) stop("All time to fatal event have to be non-negative")
# check if delta1 and delta2 are coded as 0 and 1
if(!all(delta1%in%c(0,1))) stop("Vector delta1 can only contain 0 or 1")
if(!all(delta2%in%c(0,1))) stop("Vector delta2 can only contain 0 or 1")





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





























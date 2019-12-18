library(moments)
library(ggplot2)
library(tidyr)
library(dplyr)
library(matrixStats)
library(abc)
library(np)
library(ks)
library(statip)
set.seed(1234)
rm(list=ls())

N <- 25000 
# Number of Monte Carlo for ABC
B <- 100
# Number of replications for the diagnostic (used to calculate empirical quantile)
T <- 100                                                                   
# Sample size
J<-5
# Number of different values to use for sigma

theta1 <- c(2,1,1,1,1) 
# Vector of parameters for true model. 
# The first value of theta1 is used to get the value of t_n "under correct specifiation". 
# That is, the first entry in the loop does not have any misspecifiaction (sigmaSeq[1]=1).
# We tried 0,1,2 for the first entry of theta1 and the results were not sensative to this value.

sigmaSeq <- seq(1,J,1)
# After this, we go into the loop all the way to sigmaSeq[5]=5. 


# Degree of Misspecification Experiment %
z <- matrix(0,T,N) 
Rep_rej_preds1 = array(0,c(B,3,J))
Rep_rejLN_preds3 = array(0,c(B,3,J))
Rep_rejNN_preds4 = array(0,c(B,3,J))
Rep_rej_preds1_justMean = array(0,c(B,3,J))
Rep_rejLN_preds3_justMean =  array(0,c(B,3,J))

for (j in 1:J){
  
for (b in 1:B){
# PRIOR specificaiton
mu_p <- rnorm(N)*5 


#%%%%%%%%%Generate real "observed data" with misspecifcation %%%%%%%%%%%%%%%
eps <- rnorm(T)
y <- vector()


y<-theta1[j]+sigmaSeq[j]*eps 
y1 <- mean(y) 
y2 <- mean((y-y1)^2) 

eta_y <- c(y1, y2) 

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  % First loop to get the ABC posterior that we will ater draw from 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#  %%%%%%%%% Generate fake "observed data" without misspecifcation %%%%%%%%%%%%

#% Generate random number between 1:N %
center_val <- sample(1:N,1)
sim_vals <- sample(1:N,B)

s0=function(x) mean(x)
s1=function(x) mean((x-s0(x))^2)

y1_m1 <- s0(y) 
y2_m1 <- s1(y)
eta_y <- c(y1_m1, y2_m1) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta_z <- matrix(0,nrow=N, ncol = 2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Simulate the errors: do once of dimension T,N
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu <- matrix( rnorm(T*N,mean=0,sd=1), T, N) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Simulate the simulated data: do once of dimension T,N
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_mat <-matrix(mu_p,nrow=T,ncol=N,byrow=TRUE)
z <- mu_mat+nu 

eta_z[,1] <- apply(z,2,s0)
eta_z[,2] <- apply(z,2,s1)

d_rho1 <- eta_z[,1]-eta_y[1] 
d_rho2 <- eta_z[,2]-eta_y[2] 
d_rho1_c<- as.matrix(cbind(d_rho1, d_rho2))
d_rho_c<-apply(abs(d_rho1_c),1,sum)

tol1 <- quantile(d_rho_c,.01) 
d_keep1 <- d_rho_c<=tol1
mu11 <- mu_p[d_keep1] 
mu_class<-density(mu11) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STEP 5   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Generating observed distances
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Thetas<- cbind.data.frame(mu_p)


rej=abc(target=eta_y, param=Thetas, sumstat=eta_z, tol=.01, method="rejection")
lin_abc=abc(target=eta_y, param=Thetas, sumstat=eta_z, tol=.01, method="loclinear", hcorr=FALSE,transf=c("none"))



# Recording the results for the different functions in the experiment: h(\theta)=\theta and h(\theta)=(\theta^2,\theta^3)'
# For ABC-AR
Rep_rej_preds1[b,,j] = cbind(mean(mu11),mean(mu11^2),mean(mu11^3) )
# For ABC-Reg
Rep_rejLN_preds3[b,,j] = cbind(mean((lin_abc$adj.values)),mean((lin_abc$adj.values)^2),mean((lin_abc$adj.values)^3) )


}

}
# Now we take the difference of the statistics
# This first vals1 contains the results for the correct specifiation case and is used to build the value of t_n
# [,k,j] refers to the k-th array elements across the replications under the j-th experiment value for sigma^2
# For instance [,2,1] is the value of the second element of the function h(\theta)=(\theta,\theta^2,\theta^3)' 
# for the first value of sigmaSeq, i.e., for sigmaSeq=1 (correct specification)

vals1 = cbind(Rep_rej_preds1[,2,1]-Rep_rejLN_preds3[,2,1],Rep_rej_preds1[,3,1]-Rep_rejLN_preds3[,3,1])
vals2 = cbind(Rep_rej_preds1[,2,2]-Rep_rejLN_preds3[,2,2],Rep_rej_preds1[,3,2]-Rep_rejLN_preds3[,3,2])
vals3 = cbind(Rep_rej_preds1[,2,3]-Rep_rejLN_preds3[,2,3],Rep_rej_preds1[,3,3]-Rep_rejLN_preds3[,3,3])
vals4 = cbind(Rep_rej_preds1[,2,4]-Rep_rejLN_preds3[,2,4],Rep_rej_preds1[,3,4]-Rep_rejLN_preds3[,3,4])
vals5 = cbind(Rep_rej_preds1[,2,5]-Rep_rejLN_preds3[,2,5],Rep_rej_preds1[,3,5]-Rep_rejLN_preds3[,3,5])

# Calculating the norm we will compare across correc and incorrect experiments
valss1 = sqrt(T)*sqrt(rowSums((vals1^2)))
valss2 = sqrt(T)*sqrt(rowSums((vals2^2)))
valss3 = sqrt(T)*sqrt(rowSums((vals3^2)))
valss4 = sqrt(T)*sqrt(rowSums((vals4^2)))
valss5 = sqrt(T)*sqrt(rowSums((vals5^2)))

# Getting the diagnostic value to compare against as the quantile of the simulated distirbution of the statistic
# in the case of h(\theta)=(\theta^2,\theta^3)'
t_n <- as.numeric(quantile(valss1,.95))
# Calculating rejection frequencies 
Reject_frequ <- c(sum(valss2>t_n),sum(valss3>t_n),sum(valss4>t_n),sum(valss5>t_n))/B


# Plotting the outputs of the statistic across the different values of sigma
nameS <-  c("sigma^2=1","sigma^2=2","sigma^2=3","sigma^2=4","sigma^2=5")
par(mfrow=c(1,2))
nameS <-  c("sigma2=1","sigma2=2")
boxplot(cbind(valss1, valss2),names=nameS,outline=FALSE,col = c("red","royalblue2"),main="Panel A")
nameS <-  c("sigma2=2","sigma2=3")
boxplot(cbind(valss2,valss3),names=nameS,outline=FALSE,col = c("royalblue2","green"),main="Panel B")

### Repeating the experiment for just the mean (i.e., h(\theta)=\theta) ###

vals1_mean = cbind(Rep_rej_preds1[,1,1]-Rep_rejLN_preds3[,1,1])
vals2_mean = cbind(Rep_rej_preds1[,1,2]-Rep_rejLN_preds3[,1,2])
vals3_mean = cbind(Rep_rej_preds1[,1,3]-Rep_rejLN_preds3[,1,3])
vals4_mean = cbind(Rep_rej_preds1[,1,4]-Rep_rejLN_preds3[,1,4])
vals5_mean = cbind(Rep_rej_preds1[,1,5]-Rep_rejLN_preds3[,1,5])
valss1_mean = sqrt(T)*sqrt(rowSums((vals1_mean^2)))
valss2_mean = sqrt(T)*sqrt(rowSums((vals2_mean^2)))
valss3_mean = sqrt(T)*sqrt(rowSums((vals3_mean^2)))
valss4_mean = sqrt(T)*sqrt(rowSums((vals4_mean^2)))
valss5_mean = sqrt(T)*sqrt(rowSums((vals5_mean^2)))


par(mfrow=c(1,1))
#tikz( 'Box_Stat_dist_compare_mean.tex', width=10, height=4 )

par(mfrow=c(1,2))
nameS <-  c("sigma2=1","sigma2=2")
boxplot(cbind(valss1_mean, valss2_mean),names=nameS,outline=FALSE,col = c("red","royalblue2"),main="Panel A")
nameS <-  c("sigma2=2","sigma2=3")
boxplot(cbind(valss2_mean,valss3_mean),names=nameS,outline=FALSE,col = c("royalblue2","green"),main="Panel B")

# Getting the diagnostic value to compare against as the quantile of the simulated distirbution of the statistic
# in the case of h(\theta)=(\theta)'

t_n <- as.numeric(quantile(valss1_mean,.95))
# Calculating rejection frequency 
Reject_frequ <- c(sum(valss2_mean>t_n),sum(valss3_mean>t_n),sum(valss4>t_n),sum(valss5>t_n))/B



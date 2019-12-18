library(moments)
library(lazyeval)
library(ggplot2)
library(tidyr)
library(dplyr)
library(matrixStats)
library(abc)
library(np)
library(ks)
library(tikzDevice)
rm(list=ls())
set.seed(5)

N <- 25000 
# Number of Monte Carlo for ABC
B <- 1000
# Number of Replications
T <- 100                       
# Sample size
m <- 10
# Reps for ABC-RegN
theta1 <- 1 
# True value of the parameters


sigma <- seq(1,5,1)
# Degree of Misspecification Experiment %
z <- matrix(0,T,N) 

# Allocating
Rep_rej = array(0,c(B,4,length(sigma)))
Rep_cov_preds1 = matrix(0,B,4)
Rep_rejLN = array(0,c(B,4,length(sigma)))
Rep_rejLN_Corr = array(0,c(B,4,length(sigma)))
Rep_rejLN_New = array(0,c(B,4,length(sigma)))
Rep_rejLN_New_Corr = array(0,c(B,4,length(sigma)))
Rep_rejNN = array(0,c(B,4,length(sigma)))
Rep_rejNN_Corr = array(0,c(B,4,length(sigma)))

for (k in 1:length(sigma))
{  


for(j in 1:B){
    # PRIOR specificaiton
    mu_p <- rnorm(N)*5 
    Sigma<-sqrt(sigma[k])
    
    #%%%%%%%%%Generate real "observed data" with misspecifcation %%%%%%%%%%%%%%%
    eps <- rnorm(T)
    y <- vector()
    
    
    y<-theta1+Sigma*eps 
    y1 <- mean(y) 
    y2 <- mean((y-y1)^2) 
    
    eta_y <- c(y1, y2) 
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #  % First loop to get the ABC posterior that we will ater draw from 
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    #  %%%%%%%%% Generate fake "observed data" without misspecifcation %%%%%%%%%%%%
    

    s0=function(x) mean(x)
    s1=function(x) mean((x-s0(x))^2)
    
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
 
     lin_abc <- abc(target=eta_y, param=Thetas, sumstat=eta_z, tol=.01, method ="loclinear", hcorr=FALSE, transf=c("none"))
     lin_abc_corr <- abc(target=eta_y, param=Thetas, sumstat=eta_z, tol=.01, method ="loclinear", hcorr=TRUE, transf=c("none"))
     nn_abc <- abc(target=eta_y, param=Thetas, sumstat=eta_z, tol=.01, method ="neuralnet", hcorr=FALSE, transf=c("none"))
     nn_abc_corr <- abc(target=eta_y, param=Thetas, sumstat=eta_z, tol=.01, method ="neuralnet", hcorr=TRUE, transf=c("none"))
 

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #% ABC-RegN
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ## Simulated New summary stat at post mean ##
     nu_new <- matrix( rnorm(T*m,mean=0,sd=1), T, m) 
     z_new <- mean(mu11)+nu_new
     eta_new = cbind(mean(apply(z_new,2,s0)),mean(apply(z_new,2,s1)))
     
    lin_abc_new=abc(target=eta_new, param=Thetas, sumstat=eta_z, tol=.01, method="loclinear", hcorr=FALSE,transf=c("none"))
    lin_abc_new_corr=abc(target=eta_new, param=Thetas, sumstat=eta_z, tol=.01, method="loclinear", hcorr=TRUE,transf=c("none"))
    sum_lin <- summary(lin_abc, intvl = .95)
    
    Rep_rej[j, ,k] = cbind(mean(mu11),mean((mu11-mean(mu11))^2 ),quantile(mu11,.025),quantile(mu11,.975))
    Rep_rejLN[j, ,k] = cbind(mean((lin_abc$adj.values)),var((lin_abc$adj.values)),quantile((lin_abc$adj.values),.025),quantile((lin_abc$adj.values),.975))
    Rep_rejLN_Corr[j, ,k] = cbind(mean((lin_abc_corr$adj.values)),var((lin_abc_corr$adj.values)),quantile((lin_abc_corr$adj.values),.025),quantile((lin_abc_corr$adj.values),.975))
    Rep_rejLN_New[j, ,k] = cbind(mean((lin_abc_new$adj.values)),var((lin_abc_new$adj.values)),quantile((lin_abc_new$adj.values),.025),quantile((lin_abc_new$adj.values),.975))
    Rep_rejLN_New_Corr[j, ,k] = cbind(mean((lin_abc_new_corr$adj.values)),var((lin_abc_new_corr$adj.values)),quantile((lin_abc_new_corr$adj.values),.025),quantile((lin_abc_new_corr$adj.values),.975))
    Rep_rejNN[j, ,k] = cbind(mean((nn_abc$adj.values)),var(nn_abc$adj.values),quantile((nn_abc$adj.values),.025),quantile((nn_abc$adj.values),.975))
    Rep_rejNN_Corr[j, ,k] = cbind(mean((nn_abc_corr$adj.values)),var(nn_abc_corr$adj.values),quantile((nn_abc_corr$adj.values),.025),quantile((nn_abc_corr$adj.values),.975))

}


}


### Calculating Coverages ###

# sigma2=1
Mat1 = (cbind(Rep_rej[,3,1]<1,Rep_rej[,4,1]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej1 = length(Mat2[Mat2==TRUE])/B

# sigma2=2
Mat1 = (cbind(Rep_rej[,3,2]<1,Rep_rej[,4,2]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej2 = length(Mat2[Mat2==TRUE])/B

# sigma2=3
Mat1 = (cbind(Rep_rej[,3,3]<1,Rep_rej[,4,3]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej3 = length(Mat2[Mat2==TRUE])/B

# sigma2=4
Mat1 = (cbind(Rep_rej[,3,4]<1,Rep_rej[,4,4]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej4 = length(Mat2[Mat2==TRUE])/B

# sigma2=5
Mat1 = (cbind(Rep_rej[,3,5]<1,Rep_rej[,4,5]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej5 = length(Mat2[Mat2==TRUE])/B

COV_REJ <- cbind(Cov_rej1,Cov_rej2,Cov_rej3,Cov_rej4,Cov_rej5)

# sigma2=1
Mat1 = (cbind(Rep_rejLN[,3,1]<1,Rep_rejLN[,4,1]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej1 = length(Mat2[Mat2==TRUE])/B

# sigma2=2
Mat1 = (cbind(Rep_rejLN[,3,2]<1,Rep_rejLN[,4,2]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej2 = length(Mat2[Mat2==TRUE])/B

# sigma2=3
Mat1 = (cbind(Rep_rejLN[,3,3]<1,Rep_rejLN[,4,3]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej3 = length(Mat2[Mat2==TRUE])/B

# sigma2=4
Mat1 = (cbind(Rep_rejLN[,3,4]<1,Rep_rejLN[,4,4]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej4 = length(Mat2[Mat2==TRUE])/B

# sigma2=5
Mat1 = (cbind(Rep_rejLN[,3,5]<1,Rep_rejLN[,4,5]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej5 = length(Mat2[Mat2==TRUE])/B

COV_REJLN <- cbind(Cov_rej1,Cov_rej2,Cov_rej3,Cov_rej4,Cov_rej5)


# sigma2=1
Mat1 = (cbind(Rep_rejLN_New[,3,1]<1,Rep_rejLN_New[,4,1]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej1 = length(Mat2[Mat2==TRUE])/B

# sigma2=2
Mat1 = (cbind(Rep_rejLN_New[,3,2]<1,Rep_rejLN_New[,4,2]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej2 = length(Mat2[Mat2==TRUE])/B

# sigma2=3
Mat1 = (cbind(Rep_rejLN_New[,3,3]<1,Rep_rejLN_New[,4,3]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej3 = length(Mat2[Mat2==TRUE])/B

# sigma2=4
Mat1 = (cbind(Rep_rejLN_New[,3,4]<1,Rep_rejLN_New[,4,4]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej4 = length(Mat2[Mat2==TRUE])/B

# sigma2=5
Mat1 = (cbind(Rep_rejLN_New[,3,5]<1,Rep_rejLN_New[,4,5]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej5 = length(Mat2[Mat2==TRUE])/B

COV_REJLN_New <- cbind(Cov_rej1,Cov_rej2,Cov_rej3,Cov_rej4,Cov_rej5)


# sigma2=1
Mat1 = (cbind(Rep_rejNN[,3,1]<1,Rep_rejNN[,4,1]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej1 = length(Mat2[Mat2==TRUE])/B

# sigma2=2
Mat1 = (cbind(Rep_rejNN[,3,2]<1,Rep_rejNN[,4,2]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej2 = length(Mat2[Mat2==TRUE])/B

# sigma2=3
Mat1 = (cbind(Rep_rejNN[,3,3]<1,Rep_rejNN[,4,3]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej3 = length(Mat2[Mat2==TRUE])/B

# sigma2=4
Mat1 = (cbind(Rep_rejNN[,3,4]<1,Rep_rejNN[,4,4]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej4 = length(Mat2[Mat2==TRUE])/B

# sigma2=5
Mat1 = (cbind(Rep_rejNN[,3,5]<1,Rep_rejNN[,4,5]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej5 = length(Mat2[Mat2==TRUE])/B

COV_REJNN <- cbind(Cov_rej1,Cov_rej2,Cov_rej3,Cov_rej4,Cov_rej5)

# Hetero Versions

# sigma2=1
Mat1 = (cbind(Rep_rejLN_Corr[,3,1]<1,Rep_rejLN_Corr[,4,1]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej1 = length(Mat2[Mat2==TRUE])/B

# sigma2=2
Mat1 = (cbind(Rep_rejLN_Corr[,3,2]<1,Rep_rejLN_Corr[,4,2]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej2 = length(Mat2[Mat2==TRUE])/B

# sigma2=3
Mat1 = (cbind(Rep_rejLN_Corr[,3,3]<1,Rep_rejLN_Corr[,4,3]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej3 = length(Mat2[Mat2==TRUE])/B

# sigma2=4
Mat1 = (cbind(Rep_rejLN_Corr[,3,4]<1,Rep_rejLN_Corr[,4,4]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej4 = length(Mat2[Mat2==TRUE])/B

# sigma2=5
Mat1 = (cbind(Rep_rejLN_Corr[,3,5]<1,Rep_rejLN_Corr[,4,5]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej5 = length(Mat2[Mat2==TRUE])/B

COV_REJLN_Corr <- cbind(Cov_rej1,Cov_rej2,Cov_rej3,Cov_rej4,Cov_rej5)


# sigma2=1
Mat1 = (cbind(Rep_rejLN_New_Corr[,3,1]<1,Rep_rejLN_New_Corr[,4,1]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej1 = length(Mat2[Mat2==TRUE])/B

# sigma2=2
Mat1 = (cbind(Rep_rejLN_New_Corr[,3,2]<1,Rep_rejLN_New_Corr[,4,2]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej2 = length(Mat2[Mat2==TRUE])/B

# sigma2=3
Mat1 = (cbind(Rep_rejLN_New_Corr[,3,3]<1,Rep_rejLN_New_Corr[,4,3]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej3 = length(Mat2[Mat2==TRUE])/B

# sigma2=4
Mat1 = (cbind(Rep_rejLN_New_Corr[,3,4]<1,Rep_rejLN_New_Corr[,4,4]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej4 = length(Mat2[Mat2==TRUE])/B

# sigma2=5
Mat1 = (cbind(Rep_rejLN_New_Corr[,3,5]<1,Rep_rejLN_New_Corr[,4,5]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej5 = length(Mat2[Mat2==TRUE])/B

COV_REJLN_New_Corr <- cbind(Cov_rej1,Cov_rej2,Cov_rej3,Cov_rej4,Cov_rej5)


# sigma2=1
Mat1 = (cbind(Rep_rejNN_Corr[,3,1]<1,Rep_rejNN_Corr[,4,1]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej1 = length(Mat2[Mat2==TRUE])/B

# sigma2=2
Mat1 = (cbind(Rep_rejNN_Corr[,3,2]<1,Rep_rejNN_Corr[,4,2]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej2 = length(Mat2[Mat2==TRUE])/B

# sigma2=3
Mat1 = (cbind(Rep_rejNN_Corr[,3,3]<1,Rep_rejNN_Corr[,4,3]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej3 = length(Mat2[Mat2==TRUE])/B

# sigma2=4
Mat1 = (cbind(Rep_rejNN_Corr[,3,4]<1,Rep_rejNN_Corr[,4,4]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej4 = length(Mat2[Mat2==TRUE])/B

# sigma2=5
Mat1 = (cbind(Rep_rejNN_Corr[,3,5]<1,Rep_rejNN_Corr[,4,5]>1))
Mat2 <- Mat1[,1]==Mat1[,2]
Cov_rej5 = length(Mat2[Mat2==TRUE])/B

COV_REJNN_Corr <- cbind(Cov_rej1,Cov_rej2,Cov_rej3,Cov_rej4,Cov_rej5)



COVS <- rbind(COV_REJ,COV_REJLN_New,COV_REJLN,COV_REJNN,COV_REJLN_New_Corr,COV_REJLN_Corr,COV_REJNN_Corr)
### Plotts ##

#par(mfrow = c(1, 2))
#boxplot(cbind(Rep_rej_preds1[,1],Rep_rejLN_preds3[,1],Rep_rejNN_preds4[,1]), names=Names, col=c("red","green","royalblue2"), main="sigma^2=1")
#tikz( 'Box_reg_new_all', width=7, height=3.5 )
#boxplot(cbind(Rep_rej_preds1[,1],Rep_rejLN_preds3[,1],Rep_rejNN_preds4[,1]), names=Names, col=c("red","green","royalblue2"))

# Plot in main paper
Names=cbind("AR", "RegN", "Reg")
#
#tikz( 'Box_all_MainPaper.tex', width=8, height=3.5 )
par(mfrow = c(1, 5))
boxplot(cbind(Rep_rej[,1,1],Rep_rejLN_New[,1,1],Rep_rejLN[,1,1]), names=Names, col=c("red","green","royalblue2"), main="sigma2=1")
boxplot(cbind(Rep_rej[,1,2],Rep_rejLN_New[,1,2],Rep_rejLN[,1,2]), names=Names, col=c("red","green","royalblue2"), main="sigma2=2")
boxplot(cbind(Rep_rej[,1,3],Rep_rejLN_New[,1,3],Rep_rejLN[,1,3]), names=Names, col=c("red","green","royalblue2"), main="sigma2=3")
boxplot(cbind(Rep_rej[,1,4],Rep_rejLN_New[,1,4],Rep_rejLN[,1,4]), names=Names, col=c("red","green","royalblue2"), main="sigma2=4")
boxplot(cbind(Rep_rej[,1,5],Rep_rejLN_New[,1,5],Rep_rejLN[,1,5]), names=Names, col=c("red","green","royalblue2"), main="sigma2=5")

#dev.off()

# Plots in appendix: No Hetero correction
Names=cbind("AR", "RegN", "Reg","NN")
#
#tikz( 'Box_all_SuppAppNoHetero.tex', width=10, height=4.0 )
par(mfrow = c(1, 5))
boxplot(cbind(Rep_rej[,1,1],Rep_rejLN_New[,1,1],Rep_rejLN[,1,1],Rep_rejNN[,1,1]), names=Names, col=c("red","green","royalblue2"), main="sigma2=1")
boxplot(cbind(Rep_rej[,1,2],Rep_rejLN_New[,1,2],Rep_rejLN[,1,2],Rep_rejNN[,1,2]), names=Names, col=c("red","green","royalblue2"), main="sigma2=2")
boxplot(cbind(Rep_rej[,1,3],Rep_rejLN_New[,1,3],Rep_rejLN[,1,3],Rep_rejNN[,1,3]), names=Names, col=c("red","green","royalblue2"), main="sigma2=3")
boxplot(cbind(Rep_rej[,1,4],Rep_rejLN_New[,1,4],Rep_rejLN[,1,4],Rep_rejNN[,1,4]), names=Names, col=c("red","green","royalblue2"), main="sigma2=4",outline=FALSE)
boxplot(cbind(Rep_rej[,1,5],Rep_rejLN_New[,1,5],Rep_rejLN[,1,5],Rep_rejNN[,1,5]), names=Names, col=c("red","green","royalblue2"), main="sigma2=5",outline=FALSE)
#dev.off()

# Plots in appendix: Hetero correction
Names=cbind("AR", "RegNC", "RegC","NNC")
#
#tikz( 'Box_all_SuppAppHetero.tex', width=11, height=4.0 )
par(mfrow = c(1, 5))
boxplot(cbind(Rep_rej[,1,1],Rep_rejLN_New_Corr[,1,1],Rep_rejLN_Corr[,1,1],Rep_rejNN_Corr[,1,1]), names=Names, col=c("red","green","royalblue2"), main="sigma2=1")
boxplot(cbind(Rep_rej[,1,2],Rep_rejLN_New_Corr[,1,2],Rep_rejLN_Corr[,1,2],Rep_rejNN_Corr[,1,2]), names=Names, col=c("red","green","royalblue2"), main="sigma2=2")
boxplot(cbind(Rep_rej[,1,3],Rep_rejLN_New_Corr[,1,3],Rep_rejLN_Corr[,1,3],Rep_rejNN_Corr[,1,3]), names=Names, col=c("red","green","royalblue2"), main="sigma2=3")
boxplot(cbind(Rep_rej[,1,4],Rep_rejLN_New_Corr[,1,4],Rep_rejLN_Corr[,1,4],Rep_rejNN_Corr[,1,4]), names=Names, col=c("red","green","royalblue2"), main="sigma2=4",outline=FALSE)
boxplot(cbind(Rep_rej[,1,5],Rep_rejLN_New_Corr[,1,5],Rep_rejLN_Corr[,1,5],Rep_rejNN_Corr[,1,5]), names=Names, col=c("red","green","royalblue2"), main="sigma2=5",outline=FALSE)
#dev.off

# Average Posterior Variances
var_1<-cbind(median((Rep_rej[,2,1])),median((Rep_rej[,2,2])),median((Rep_rej[,2,3])),median((Rep_rej[,2,4])),median((Rep_rej[,2,5])))
var_2<-cbind(median((Rep_rejLN_New[,2,1])),median((Rep_rejLN_New[,2,2])),median((Rep_rejLN_New[,2,3])),median((Rep_rejLN_New[,2,4])),median((Rep_rejLN_New[,2,5])))
var_3<-cbind(median((Rep_rejLN[,2,1])),median((Rep_rejLN[,2,2])),median((Rep_rejLN[,2,3])),median((Rep_rejLN[,2,4])),median((Rep_rejLN[,2,5])))
var_4<-cbind(median((Rep_rejNN[,2,1])),median((Rep_rejNN[,2,2])),median((Rep_rejNN[,2,3])),median((Rep_rejNN[,2,4])),median((Rep_rejNN[,2,5])))
var_5<-cbind(median((Rep_rejLN_New_Corr[,2,1])),median((Rep_rejLN_New_Corr[,2,2])),median((Rep_rejLN_New_Corr[,2,3])),median((Rep_rejLN_New_Corr[,2,4])),median((Rep_rejLN_New_Corr[,2,5])))
var_6<-cbind(median((Rep_rejLN_Corr[,2,1])),median((Rep_rejLN_Corr[,2,2])),median((Rep_rejLN_Corr[,2,3])),median((Rep_rejLN_Corr[,2,4])),median((Rep_rejLN_Corr[,2,5])))
var_7<-cbind(median((Rep_rejNN_Corr[,2,1])),median((Rep_rejNN_Corr[,2,2])),median((Rep_rejNN_Corr[,2,3])),median((Rep_rejNN_Corr[,2,4])),median((Rep_rejNN_Corr[,2,5])))

Var<- rbind(var_1,var_2,var_3,var_4,var_5,var_6,var_7)

# Average Credible set Length
l1<-cbind(median((Rep_rej[,4,1])-(Rep_rej[,3,1])),median((Rep_rej[,4,2])-(Rep_rej[,3,2])),median((Rep_rej[,4,3])-(Rep_rej[,3,3])),median((Rep_rej[,4,4])-(Rep_rej[,3,4])),median((Rep_rej[,4,5])-(Rep_rej[,3,5])))
l2<-cbind(median((Rep_rejLN_New[,4,1])-(Rep_rejLN_New[,3,1])),median((Rep_rejLN_New[,4,2])-(Rep_rejLN_New[,3,2])),median((Rep_rejLN_New[,4,3])-(Rep_rejLN_New[,3,3])),median((Rep_rejLN_New[,4,4])-(Rep_rejLN_New[,3,4])),median((Rep_rejLN_New[,4,5])-(Rep_rejLN_New[,3,5])))
l3<-cbind(median((Rep_rejLN[,4,1])-(Rep_rejLN[,3,1])),median((Rep_rejLN[,4,2])-(Rep_rejLN[,3,2])),median((Rep_rejLN[,4,3])-(Rep_rejLN[,3,3])),median((Rep_rejLN[,4,4])-(Rep_rejLN[,3,4])),median((Rep_rejLN[,4,5])-(Rep_rejLN[,3,5])))
l4<-cbind(median((Rep_rejNN[,4,1])-(Rep_rejNN[,3,1])),median((Rep_rejNN[,4,2])-(Rep_rejNN[,3,2])),median((Rep_rejNN[,4,3])-(Rep_rejNN[,3,3])),median((Rep_rejNN[,4,4])-(Rep_rejNN[,3,4])),median((Rep_rejNN[,4,5])-(Rep_rejNN[,3,5])))
l5<-cbind(median((Rep_rejLN_New_Corr[,4,1])-(Rep_rejLN_New_Corr[,3,1])),median((Rep_rejLN_New_Corr[,4,2])-(Rep_rejLN_New_Corr[,3,2])),median((Rep_rejLN_New_Corr[,4,3])-(Rep_rejLN_New_Corr[,3,3])),median((Rep_rejLN_New_Corr[,4,4])-(Rep_rejLN_New_Corr[,3,4])),median((Rep_rejLN_New_Corr[,4,5])-(Rep_rejLN_New_Corr[,3,5])))
l6<-cbind(median((Rep_rejLN_Corr[,4,1])-(Rep_rejLN_Corr[,3,1])),median((Rep_rejLN_Corr[,4,2])-(Rep_rejLN_Corr[,3,2])),median((Rep_rejLN_Corr[,4,3])-(Rep_rejLN_Corr[,3,3])),median((Rep_rejLN_Corr[,4,4])-(Rep_rejLN_Corr[,3,4])),median((Rep_rejLN_Corr[,4,5])-(Rep_rejLN_Corr[,3,5])))
l7<-cbind(median((Rep_rejNN_Corr[,4,1])-(Rep_rejNN_Corr[,3,1])),median((Rep_rejNN_Corr[,4,2])-(Rep_rejNN_Corr[,3,2])),median((Rep_rejNN_Corr[,4,3])-(Rep_rejNN_Corr[,3,3])),median((Rep_rejNN_Corr[,4,4])-(Rep_rejNN_Corr[,3,4])),median((Rep_rejNN_Corr[,4,5])-(Rep_rejNN_Corr[,3,5])))

Length<- rbind(l1,l2,l3,l4,l5,l6,l7)


# Median Quantiles: 975
qU1<-cbind(median((Rep_rej[,4,1])),median((Rep_rej[,4,2])),median((Rep_rej[,4,3])),median((Rep_rej[,4,4])),median((Rep_rej[,4,5])))
qU2<-cbind(median((Rep_rejLN_New[,4,1])),median((Rep_rejLN_New[,4,2])),median((Rep_rejLN_New[,4,3])),median((Rep_rejLN_New[,4,4])),median((Rep_rejLN_New[,4,5])))
qU3<-cbind(median((Rep_rejLN[,4,1])),median((Rep_rejLN[,4,2])),median((Rep_rejLN[,4,3])),median((Rep_rejLN[,4,4])),median((Rep_rejLN[,4,5])))
qU4<-cbind(median((Rep_rejNN[,4,1])),median((Rep_rejNN[,4,2])),median((Rep_rejNN[,4,3])),median((Rep_rejNN[,4,4])),median((Rep_rejNN[,4,5])))
qU5<-cbind(median((Rep_rejLN_New_Corr[,4,1])),median((Rep_rejLN_New_Corr[,4,2])),median((Rep_rejLN_New_Corr[,4,3])),median((Rep_rejLN_New_Corr[,4,4])),median((Rep_rejLN_New_Corr[,4,5])))
qU6<-cbind(median((Rep_rejLN_Corr[,4,1])),median((Rep_rejLN_Corr[,4,2])),median((Rep_rejLN_Corr[,4,3])),median((Rep_rejLN_Corr[,4,4])),median((Rep_rejLN_Corr[,4,5])))
qU7<-cbind(median((Rep_rejNN_Corr[,4,1])),median((Rep_rejNN_Corr[,4,2])),median((Rep_rejNN_Corr[,4,3])),median((Rep_rejNN_Corr[,4,4])),median((Rep_rejNN_Corr[,4,5])))

Q975<-rbind(qU1,qU2,qU3,qU4,qU5,qU6,qU7)

# Median Quantiles: 025
qL1<-cbind(median((Rep_rej[,3,1])),median((Rep_rej[,3,2])),median((Rep_rej[,3,3])),median((Rep_rej[,3,4])),median((Rep_rej[,3,5])))
qL2<-cbind(median((Rep_rejLN_New[,3,1])),median((Rep_rejLN_New[,3,2])),median((Rep_rejLN_New[,3,3])),median((Rep_rejLN_New[,3,4])),median((Rep_rejLN_New[,3,5])))
qL3<-cbind(median((Rep_rejLN[,3,1])),median((Rep_rejLN[,3,2])),median((Rep_rejLN[,3,3])),median((Rep_rejLN[,3,4])),median((Rep_rejLN[,3,5])))
qL4<-cbind(median((Rep_rejNN[,3,1])),median((Rep_rejNN[,3,2])),median((Rep_rejNN[,3,3])),median((Rep_rejNN[,3,4])),median((Rep_rejNN[,3,5])))
qL5<-cbind(median((Rep_rejLN_New_Corr[,3,1])),median((Rep_rejLN_New_Corr[,3,2])),median((Rep_rejLN_New_Corr[,3,3])),median((Rep_rejLN_New_Corr[,3,4])),median((Rep_rejLN_New_Corr[,3,5])))
qL6<-cbind(median((Rep_rejLN_Corr[,3,1])),median((Rep_rejLN_Corr[,3,2])),median((Rep_rejLN_Corr[,3,3])),median((Rep_rejLN_Corr[,3,4])),median((Rep_rejLN_Corr[,3,5])))
qL7<-cbind(median((Rep_rejNN_Corr[,3,1])),median((Rep_rejNN_Corr[,3,2])),median((Rep_rejNN_Corr[,3,3])),median((Rep_rejNN_Corr[,3,4])),median((Rep_rejNN_Corr[,3,5])))

Q025<-rbind(qL1,qL2,qL3,qL4,qL5,qL6,qL7)



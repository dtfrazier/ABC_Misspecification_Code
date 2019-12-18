library(moments)
library(ggplot2)
library(tidyr)
library(dplyr)
library(matrixStats)
library(reshape)
library(ggplot2)
library(tikzDevice)
rm(list=ls())

N <- 25000 
# Number of MC trials
T <- 100             
# Sample size
M = 50
# Number of grid points for the diagnostic
J = 10
# Number of different values to use for the experiment (will determine number of sigma values to try)

theta1 <- 1 
# True value of parameter
sigma <- seq(1,J,by=(3-1)/(J-1)) 
# grid of sigma values for the experiment

# Allocating matrices
Vals = matrix(0,M,J)
Val_Line = matrix(0,M,J)
Eps_T2 = matrix(0,M,J)

eps <- rnorm(T)
# Generating errors for observed data just once

mu_p <- rnorm(N)*5 
# PRIOR specificaiton

nu <- matrix( rnorm(T*N,mean=0,sd=1), T, N) 
# Generating errors for simulated data just once

for (q in 1:J)
{
  z <- matrix(0,T,N) 
  #%%%%%%%%% Generate real "observed data" with misspecifcation %%%%%%%%%%%%%%%

  y <- vector()
  
  
  y<-theta1+sqrt(sigma[q])*eps 
  y1 <- mean(y) 
  y2 <- mean((y-y1)^2)-1 
  
  eta_y <- c(y1, y2) 
  

  mu_mat <-matrix(mu_p,nrow=T,ncol=N,byrow=TRUE)
  z <- mu_mat+nu 
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #% Generating observed distances
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  d_rho1 <- matrix(0,nrow=N,ncol=1) 
  d_rho2 <- matrix(0,nrow=N,ncol=1) 
  d_rho1_c<- matrix(0,nrow=N,ncol=1)
  
  for(k in 1:N){
    zz <- z[,k] 
    zz1 <- mean(zz) 
    zz2 <- (mean((zz-zz1)^2))-1 
    
    eta_zz <- c(zz1, zz2) 
    
    
    d_rho1[k]<-sqrt(sum(((eta_y-eta_zz)^2))) 
  }
  

  Min_obs <- min(d_rho1)
  # Starting tolerance level:
  tol1 <- quantile(d_rho1,.1) 
  # End tolerance level:
  tol_end <- quantile(d_rho1,.0001)

  # Generating a sequence of tolerances
    eps_T2 = seq(from =   tol_end, to = tol1, by = ((tol1 -   tol_end)/(M - 1)))
  # Recording them across the experiments
    Eps_T2[,q] = eps_T2
    vals1 = matrix(0,M,1)
    
# Generating the realized acceptance rates. 
    for(i in 1:M)
    { new = d_rho1<eps_T2[M+1-i]
    vals1[i] = sum(new, na.rm="TRUE")/N
    }
    
# Recording realized values across the experiments
    Vals[,q]= (vals1)
# Theoretical relationship that should be in evidence if the model is correctly specified
    Val_Line[,q] = seq(from = min(vals1), to = max(vals1), by = ((max(vals1) - min(vals1))/(M - 1)))
  }
 
# Ploting the relationships across the experiments: Forward axes
 # tikz( 'Normal_graph_diagnostic.tex', width=11, height=4.0 )
  par(mfrow=c(3,3))
  matplot(cbind(Val_Line[,1],rev(Vals[,1])),Eps_T2[,1],xlab="alpha", ylab="eps",type = "l",lwd=2)
  matplot(cbind(Val_Line[,2],rev(Vals[,2])),Eps_T2[,2],xlab="alpha", ylab="eps",type = "l",lwd=2)
  matplot(cbind(Val_Line[,3],rev(Vals[,3])),Eps_T2[,3],xlab="alpha", ylab="eps",type = "l",lwd=2)
  matplot(cbind(Val_Line[,4],rev(Vals[,4])),Eps_T2[,4],xlab="alpha", ylab="eps",type = "l",lwd=2)
  matplot(cbind(Val_Line[,5],rev(Vals[,5])),Eps_T2[,5],xlab="alpha", ylab="eps",type = "l",lwd=2)
  matplot(cbind(Val_Line[,6],rev(Vals[,6])),Eps_T2[,6],xlab="alpha", ylab="eps",type = "l",lwd=2)
  matplot(cbind(Val_Line[,7],rev(Vals[,7])),Eps_T2[,7],xlab="alpha", ylab="eps",type = "l",lwd=2)
  matplot(cbind(Val_Line[,8],rev(Vals[,8])),Eps_T2[,8],xlab="alpha", ylab="eps",type = "l",lwd=2)
  matplot(cbind(Val_Line[,9],rev(Vals[,9])),Eps_T2[,9],xlab="alpha", ylab="eps",type = "l",lwd=2)
  dev.off() 
  # PLotting reversiing the axes. 
  # Ploting the relationships across the experiments: Reverse axes
  
  par(mfrow=c(3,3))
  matplot(Eps_T2[,1],cbind(Val_Line[,1],rev(Vals[,1])),ylab="alpha", xlab="eps",type = "l",lwd=2)
  matplot(Eps_T2[,2],cbind(Val_Line[,2],rev(Vals[,2])),ylab="alpha", xlab="eps",type = "l",lwd=2)
  matplot(Eps_T2[,3],cbind(Val_Line[,3],rev(Vals[,3])),ylab="alpha", xlab="eps",type = "l",lwd=2)
  matplot(Eps_T2[,4],cbind(Val_Line[,4],rev(Vals[,4])),ylab="alpha", xlab="eps",type = "l",lwd=2)
  matplot(Eps_T2[,5],cbind(Val_Line[,5],rev(Vals[,5])),ylab="alpha", xlab="eps",type = "l",lwd=2)
  matplot(Eps_T2[,6],cbind(Val_Line[,6],rev(Vals[,6])),ylab="alpha", xlab="eps",type = "l",lwd=2)
  matplot(Eps_T2[,7],cbind(Val_Line[,7],rev(Vals[,7])),ylab="alpha", xlab="eps",type = "l",lwd=2)
  matplot(Eps_T2[,8],cbind(Val_Line[,8],rev(Vals[,8])),ylab="alpha", xlab="eps",type = "l",lwd=2)
  matplot(Eps_T2[,9],cbind(Val_Line[,9],rev(Vals[,9])),ylab="alpha", xlab="eps",type = "l",lwd=2)
  
  
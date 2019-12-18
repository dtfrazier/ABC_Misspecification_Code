library(abc)
library(readr)
library(abctools)

TPrior <- read_csv("TPrior", 
                   col_names = FALSE)

SSY1 <- read_csv("SSY1", col_names = FALSE)
SSY2 <- read_csv("SSY2", col_names = FALSE)
SSX1 <- read_csv("SSX1", col_names = FALSE)
SSX2 <- read_csv("SSX2", col_names = FALSE)

#Eps <- read_csv("Eps", col_names = FALSE)

sigmas = seq(0.5,5,.05)
B<-length(sigmas)

ABC_LocLinRect = matrix(0,B,4)
ABC_Rej = matrix(0,B,4)

for (i in 1:B)
{
  
  eta_y = cbind(SSY1[i,],SSY2[i,])
  eta_z = cbind(SSX1[,i],SSX2[,i])
  Thetas = TPrior
  
  d_rho1 <- eta_z[,1]-eta_y[1,1] 
  d_rho2 <- eta_z[,2]-eta_y[1,2] 
  d_rho1_c<- as.matrix(cbind(d_rho1, d_rho2))
  d_rho_c<-apply(abs(d_rho1_c),1,sum)
  
  tol1 <- quantile(d_rho_c,.01) 
  d_keep1 <- d_rho_c<=tol1
  mu11 <- Thetas[d_keep1,] 
  
   
  
lin_abcR <- abc(target=eta_y, param=Thetas, sumstat=eta_z, tol=.01, method ="loclinear")

ABC_Rej[i,] = cbind(mean(mu11$X1),mean((mu11$X1-mean(mu11$X1))^2 ),quantile(mu11$X1,.025),quantile(mu11$X1,.975))
ABC_LocLinRect[i,] = cbind(mean((lin_abcR$adj.values)),var((lin_abcR$adj.values)),quantile((lin_abcR$adj.values),.025),quantile((lin_abcR$adj.values),.975))
 
   
}



df <- data.frame(x=sigmas,val=cbind(ABC_Rej[,1],ABC_LocLinRect[,1]))
ggplot(data = df, aes(x)) + geom_line(aes(y=val.1, colour="ABC-AR"))+geom_line(aes(y=val.2, colour="ABC-Reg"))




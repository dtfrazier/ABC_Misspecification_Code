% 
rng(12)
clear
clc
load('Normal_data_exper_Fixed_eps.mat')
% values for sigma:  sigmas = 0.5:0.05:5;
% How the errors are darwn: r= randn(n,1);
% Number of ABC reps: M=25000;
% 
% 
% errors for ABC simulation:   eps = randn(length(r),M);
% prior for theta:   theta_prior = randn(M,1).*5;
% ABC tolerance:  tol=.01;
for i = 1:length(sigmas)
    
     sigma_true = sigmas(i)
    y = 1 + sqrt(sigma_true)*r;

  z = ones(length(y),1)*theta_prior'+eps;
   ssx=[mean(z)',var(z)'];
   ssy = [mean(y),var(y)];
   
SSX(:,:,i) = ssx;
SSY(:,:,i) = ssy;
THETA_PRIOR(:,:,i) = theta_prior;

   temp = ssx-ones(M,1)*ssy;
   
   dist = sqrt(sum(temp.^2,2));
   
tempN=[dist, (1:1:length(dist))'];
B = sortrows(tempN,1);
    
   
level=floor(tol*size(B,1));

BB(:,i) = B(level,1);


theta=theta_prior(B(1:level,2),:);
acc_ssx = ssx(B(1:level,2),:);

ACC_theta(:,i) = theta;
ACC_SSX1(:,i) = acc_ssx(:,1);
ACC_SSX2(:,i) = acc_ssx(:,2);


X=[ones(length(acc_ssx),1), (acc_ssx - ones(length(acc_ssx),1)*ssy)];
beta = (X'*X)\X'*theta;

ACC_theta_reg(:,i) = theta_reg;
Mean_theta(i,1) = mean(theta);

%% Local-Linear regression adjustment: Epanechnikov Kernel (accepted draws) %%

% EP Kernel
 w_i = .75.*(1-((B(1:level,1)./B(level,1))).^2)./(B(level,1));
% Triangle Kernel
%  w_i = 1.*(1-abs((B(1:level,1)./B(level,1))))./(B(level,1));

temp_X = acc_ssx-ones(length(acc_ssx),1)*ssy;
temp_X = [ones(size(acc_ssx,1),1) temp_X];
X = bsxfun(@times,temp_X,sqrt(w_i));
Y = bsxfun(@times,theta,sqrt(w_i));


fit = lscov(temp_X,theta,w_i);


beta_hat = (X'*X)\X'*Y;
 theta_NL_reg = theta - (acc_ssx - ones(length(acc_ssx),1)*ssy)*beta_hat(2:3);

 Mean_theta_WReg(i,1) = beta_hat(1);
ACC_theta_Wreg(:,i) = theta_NL_reg;




end

% Generate Figre 1
% plot(sigmas,Mean_theta,'blue',sigmas,Mean_theta_WReg,'black')

for i=1:length(sigmas)
    
    qL_ar(:,i) = quantile(ACC_theta(:,i),.025);
    qU_ar(:,i) = quantile(ACC_theta(:,i),.975);
    
    qL_reg(:,i) = quantile(ACC_theta_reg(:,i),.025);
    qU_reg(:,i) = quantile(ACC_theta_reg(:,i),.975);
    
    
end

SSX1 = squeeze(SSX(:,1,:));
SSX2 = squeeze(SSX(:,2,:));


for i=1:length(sigmas)
    
    [ft1,xt1] = ksdensity(ACC_theta(:,i));
    
    Ft1(:,i) = ft1;
    Xt1(:,i) = xt1;
    
    [ftr1,xtr1] = ksdensity(ACC_theta_reg(:,i));
    
    Ftr1(:,i) = ftr1;
    Xtr1(:,i) = xtr1;


end



clear
clc

load('Accepted_data.mat')
seq_length = length(sigmas);
seq = 1:2:seq_length;
cmap = colormap(summer(length(seq)));


for i= 1:length(seq)
    ss = seq(i);
    [fm1,xm1] = ksdensity(ACC_theta(:,ss));
    Fm1(:,i) = fm1;
    Xm1(:,i) = xm1;
    [fm2,xm2] = ksdensity(ACC_theta_Wreg(:,ss));
     Fm2(:,i) = fm2;
    Xm2(:,i) = xm2;

    
    hold on
     subplot(1,2,1)
    plot(Xm1(:,i),Fm1(:,i),'Color',cmap(i,:))
    
    hold on
    subplot(1,2,2)
    plot(Xm2(:,i),Fm2(:,i),'Color',cmap(i,:))
% 
%     hold on
%     subplot(2,2,3)
%     plot(Xv1(:,i),Fv1(:,i),'Color',app_map(i,:))
% 
%     hold on
%     subplot(2,2,4)
%     plot(Xv2(:,i),Fv2(:,i),'Color',app_map(i,:))

end

% cmap = colormap(hot(length(seq))) ; %Create Colormap
 cbh = colorbar ; %Create Colorbar
 vals=linspace(1,5,11);%  cbh.Ticks = sigmas; %Create 8 ticks from zero to 1
 cbh.TickLabels = num2cell(vals) ;   
 
 
 for i = 1:length(sigmas)
     qL_ar(:,i) = quantile(ACC_theta(:,i),.025);
     qU_ar(:,i) = quantile(ACC_theta(:,i),.975);
     
       qL_Reg(:,i) = quantile(ACC_theta_Wreg(:,i),.025);
     qU_Reg(:,i) = quantile(ACC_theta_Wreg(:,i),.975);
 
 end
 
 % plot(sigmas,Mean_theta,'blue',sigmas,Mean_theta_Reg,'red',sigmas,Mean_theta_NLReg,'black',sigmas,Mean_theta_WReg,'green')

 
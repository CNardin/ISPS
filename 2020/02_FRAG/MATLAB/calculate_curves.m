%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%________________________________ ETH ZURICH __________________________________
%_____Probabilistic seismic risk analysis and management for civil system______
%______________________________________________________________________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================================
%               marco broccardo (bromarco@ethz.ch)
%  Department of Civil and Environmental Engineering, IBK, ETH Zurich
%  Copyright(c) ETH Zurich. All Rights Reserved.
%===============================================================================
%===============================================================================
%========================== Script =============================================
clear all 
clc 
g=9.816;
%%
%%%%%%%%%%%%% load your hazard curve  %%%%%
% data dowloaded from 
% http://geohazards.usgs.gov/hazardtool/curves.php?format=2&lat=37.64903&lon= -122.34375&site=760&period=0p00
% data are given in mean annual rate of exceedence 
% IM= PGA

l_IM=load('results.csv');      
pga_d  = l_IM(1,:);
lm_pga_d = l_IM(2,:);

pga = 0:0.001:pga_d(numel(pga_d));
log_P_pga = spline(pga_d,log(lm_pga_d),pga);
lm_pga = exp(log_P_pga);      % rate of PGA
d_lm_pga = abs(diff(lm_pga)); % differential rate of PGA; 
%
%

close all
figure('OuterPosition',[100 100 650 460]);
loglog(pga,lm_pga,'linewidth',2.5);
set(gca,'fontname','Times')
ylabel('$\lambda(pga)$','Interpreter', 'latex')
xlabel('$pga [g]$','Interpreter', 'latex')
set(gca,'FontSize',20)
ylim([10^-5 0.5])
grid on

pga = pga(1:(numel(pga)-1)); %redefine pga for computational purposes 



%%   Loss values for damage states ductile system %%% 
 

mean_cov_loss =[0.05 0.2 0.55 0.9 0.95;
                 0.5 0.5 0.20 0.10 0.05]; % first row mean second row loss 
             
             
 
alpha_loss = (1-mean_cov_loss(1,:))./mean_cov_loss(2,:).^2 - mean_cov_loss(1,:);
beta_loss = alpha_loss.*(1-mean_cov_loss(1,:))./mean_cov_loss(1,:);
%
%
loss= 0:0.005:1;
close all
figure('OuterPosition',[100 100 560 420]);
plot(loss,1-betacdf(loss,alpha_loss(1),beta_loss(1)))
hold on 
plot(loss,1-betacdf(loss,alpha_loss(2),beta_loss(2)))
hold on 
plot(loss,1-betacdf(loss,alpha_loss(3),beta_loss(3)))
hold on 
plot(loss,1-betacdf(loss,alpha_loss(4),beta_loss(4)))
hold on 
plot(loss,1-betacdf(loss,alpha_loss(5),beta_loss(5)))
%
set(gca,'FontSize',20)
%title('ii)','Interpreter','Latex')
set(gca,'fontname','Times')
ylabel('$P(L>l|ds)$','Interpreter', 'latex')
xlabel('$l$[\% building value loss]','Interpreter', 'latex')
hhh=legend('$DS1$','$DS2$','$DS3$','$DS4$','$DS_c$');
set(hhh,'Interpreter','Latex');
set(hhh,'Location','SouthEast')
set(hhh,'FontSize',18)
legend boxoff




%% %%%%%%%%%%%%%%%%%%%%%%   Ductile system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the parameter of the fragility function for the ductile system
b = load('Ductile_par'); 
ds_par = b.Par;


%% Compute the P(DS|IM)


% Hint
p_ds0 = 1-normcdf((log(pga/ds_par(1,1)))/ds_par(1,2));
p_ds1 = abs(normcdf((log(pga/ds_par(1,1)))/ds_par(1,2)) - normcdf((log(pga/ds_par(2,1)))/ds_par(2,2)));
p_ds2 = abs(normcdf((log(pga/ds_par(2,1)))/ds_par(2,2)) - normcdf((log(pga/ds_par(3,1)))/ds_par(3,2)));
p_ds3 = abs(normcdf((log(pga/ds_par(3,1)))/ds_par(3,2)) - normcdf((log(pga/ds_par(4,1)))/ds_par(4,2)));
p_ds4 = abs(normcdf((log(pga/ds_par(4,1)))/ds_par(4,2)) - normcdf((log(pga/ds_par(5,1)))/ds_par(5,2)));
p_ds5 = abs(normcdf((log(pga/ds_par(5,1)))/ds_par(5,2)));


%% Rate of damage states
lm_ds0 = p_ds0*d_lm_pga';
lm_ds1 = p_ds1*d_lm_pga';
lm_ds2 = p_ds2*d_lm_pga';
lm_ds3 = p_ds3*d_lm_pga';
lm_ds4 = p_ds4*d_lm_pga';
lm_ds5 = p_ds5*d_lm_pga';


%% Loss curve

 lm_loss= (1-betacdf(loss,alpha_loss(1),beta_loss(1)))*lm_ds1 + (1-betacdf(loss,alpha_loss(2),beta_loss(2)))*lm_ds2+...
          (1-betacdf(loss,alpha_loss(3),beta_loss(3)))*lm_ds3 + (1-betacdf(loss,alpha_loss(4),beta_loss(4)))*lm_ds4+...
          (1-betacdf(loss,alpha_loss(4),beta_loss(5)))*lm_ds5;

figure


close all
figure('OuterPosition',[100 100 650 460]);
semilogy(10*loss,lm_loss)
set(gca,'fontname','Times')
ylabel('$\lambda(l)$','Interpreter', 'latex')
xlabel('$l$ [\$]','Interpreter', 'latex')
set(gca,'FontSize',20)
ylim([2*10^-5 0.1])
xlim([0.0 10])
grid on



%% %%%%%%%%%%%%%%%%%%%%%%   Fragile system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%% Compute the P(DS|IM)
p_ds0_f = 1-normcdf((log(pga/ds_par(1,1)))/ds_par(1,2));
p_ds1_f = abs(normcdf((log(pga/ds_par(1,1)))/ds_par(1,2)) - normcdf((log(pga/ds_par(2,1)))/ds_par(2,2)));
p_dsc_f = abs(normcdf((log(pga/ds_par(2,1)))/ds_par(2,2)));

lm_ds0_f = p_ds0_f*d_lm_pga';
lm_ds1_f = p_ds1_f*d_lm_pga';
lm_dsc_f = p_dsc_f*d_lm_pga';


mean_cov_loss_f =[0.05 0.95;
                   0.5 0.05]; % first row mean second row loss 

alpha_loss_f = (1-mean_cov_loss_f(1,:))./mean_cov_loss_f(2,:).^2 - mean_cov_loss_f(1,:);
beta_loss_f = alpha_loss_f.*(1-mean_cov_loss_f(1,:))./mean_cov_loss_f(1,:);
loss= 0:0.005:1;
close all
figure('OuterPosition',[100 100 560 420]);
plot(loss,1-betacdf(loss,alpha_loss_f(1),beta_loss_f(1)))
hold on 
plot(loss,1-betacdf(loss,alpha_loss_f(2),beta_loss_f(2)))
%
set(gca,'FontSize',20)
%title('ii)','Interpreter','Latex')
set(gca,'fontname','Times')
ylabel('$P(L>l|ds)$','Interpreter', 'latex')
xlabel('$l$[\% building value]','Interpreter', 'latex')
hhh=legend('$DS_1$','$DS_c$');
set(hhh,'Interpreter','Latex');
set(hhh,'Location','SouthEast')
set(hhh,'FontSize',18)
legend boxoff

 lm_loss_f= (1-betacdf(loss,alpha_loss_f(1),beta_loss_f(1)))*lm_ds1_f + (1-betacdf(loss,alpha_loss_f(2),beta_loss_f(2)))*lm_dsc_f;

 close all
figure('OuterPosition',[100 100 650 460]);
semilogy(10*loss,lm_loss,10*loss,lm_loss_f);
set(gca,'fontname','Times')
ylabel('$\lambda(l)$','Interpreter', 'latex')
xlabel('$l$ [\$]','Interpreter', 'latex')
set(gca,'FontSize',20)
ylim([2*10^-5 0.1])
xlim([0.0 10])
legend('ductile','fragile')
grid on
 
 
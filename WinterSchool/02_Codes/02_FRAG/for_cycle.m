%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%______________________________ UNITN _____________________________________
%____Department of Civil, Environmental and Mechanical Engineering ________
%__________________________________________________________________________
%___________________________ ETH ZURICH ___________________________________
%___Probabilistic seismic risk analysis and management for civil system____
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
%               Chiara Nardin (chiara.nardin@unitn.it)
%  Department of Civil, Environmental and Mechanical Engineering, UNITN
%  Trento - Italy
%               Marco Broccardo (marco.broccardo@unitn.it)
%  Department of Civil and Environmental Engineering, IBK, ETH Zurich
%  Copyright(c) ETH Zurich. All Rights Reserved.
%==========================================================================
%% ------------------------------------------------------------------------
%_________________________: FRAGILITY ANALYSIS :___________________________
%% ------------------------------------------------------------------------

%% General Setup
close all
clear all
clc
set(0,'DefaultFigureColor',[1 1 1])
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
g=9.816;
workDir = (cd)
imagesDir = fullfile(workDir,'images');
%% General flag/debug
 flagGM = 1; % 1 = true; 0 = false


%%  Ground motions 
Ground_motions = load('accelrot_cellarray.mat'); 
%
NN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19, 21, 23, 24,...
      25, 26, 28, 29, 30, 31, 32, 33, 34, 36, 37]; % number ID of the SS ground motions
%  
DT = [0.01 , 0.01 , 0.005, 0.005, 0.005, 0.005, 0.02 , 0.02, 0.02, 0.01,...
      0.01 , 0.01 , 0.02 , 0.02 , 0.02 , 0.02 , 0.02 , 0.02, 0.01, 0.01,...
      0.005, 0.005, 0.005, 0.005, 0.02 , 0.005, 0.005, 0.01, 0.01, 0.01,...
      0.01 , 0.005, 0.005, 0.02 , 0.005, 0.01 , 0.01]; % integration time step 

  
%% MDOF Properties 
% MDOF_properties_BW_MRF 
MDOF_properties_BW_BF

%% Structural behaviour
prompt = 'Type  -le - or - bw - for linear elastic behaviour or Bouc Wen hysteretic model: ';
system_type = input(prompt,'s')
        if isempty(system_type)
            disp('Fragile system: linear elastic behaviour assumed.')
            str = 'le';
        end
switch system_type
    case 'le'
        %% Bouc Wen Model parameters        
        n=1.5;
        Mat.ALPHA(1)=alphaaa;
        Mat.ALPHA(2)=alphaaa;
        Mat.ALPHA(3)=alphaaa;
        Mat.ALPHA(4)=alphaaa;
        Mat.ALPHA(5)=alphaaa;        
        Mat.KO(1)=ks(1,1);
        Mat.KO(2)=ks(2,2);
        Mat.KO(3)=ks(3,3);
        Mat.KO(4)=ks(4,4);
        Mat.KO(5)=ks(5,5);
        Mat.N(1:Mat.nelem)=n;
        Mat.gamma(1:Mat.nelem)=1/(2*(uy)^n);
        Mat.beta(1:Mat.nelem)=Mat.gamma(1);
        Mat.AO(1:Mat.nelem)=1;
        Mat.DELTAA(1:Mat.nelem)=0;
        Mat.DELTAv(1:Mat.nelem)=0;
        Mat.DELTAeta(1:Mat.nelem)=0;
        %% Limit state for collapse
        edpc = 0.04; % 1% interstorey drift fixed to 4[m] hight
    case 'bw'
        %% Bouc Wen Model parameters        
        n=1.5;
        Mat.ALPHA(1)=alphaaa;
        Mat.ALPHA(2)=alphaaa;
        Mat.ALPHA(3)=alphaaa;
        Mat.ALPHA(4)=alphaaa;
        Mat.ALPHA(5)=alphaaa;        
        Mat.KO(1)=ks(1,1);
        Mat.KO(2)=ks(2,2);
        Mat.KO(3)=ks(3,3);
        Mat.KO(4)=ks(4,4);
        Mat.KO(5)=ks(5,5);
        Mat.N(1:Mat.nelem)=n;
        Mat.gamma(1:Mat.nelem)=1/(2*(uy)^n);
        Mat.beta(1:Mat.nelem)=Mat.gamma(1);
        Mat.N(1:Mat.nelem)=n;
        Mat.gamma(1)=1/(2*(uy)^n);
        Mat.beta(1:Mat.nelem)=Mat.gamma(1);
        Mat.AO(1:Mat.nelem)=1;
        Mat.DELTAA(1:Mat.nelem)=01*2*10^-12;
        Mat.DELTAv(1:Mat.nelem)=01*10^-8;
        Mat.DELTAeta(1:Mat.nelem)=01*10^-8;
        %% Limit state for collapse
        edpc = 0.12; % 3% interstorey drift fixed to 4[m] hight
                edpc = 0.04; % 1% interstorey drift fixed to 4[m] hight
    otherwise
        disp('No structural type identified.')
end


%% Limit States
% LS = [0.50 0.75 1 2 3]*4/100; %MRF
LS = [0.30 0.50 1 2]*4/100; %BF
for ls_i = 1:numel(LS)
ls_val = LS(ls_i)

structuralAnalysis

ls_i = ls_i + 1;
end

return


%% Plot fragility curves for every identified limit state
for ls_i=1:numel(LS)

    matName = strcat('LS',num2str(ls_i),'.mat');
    load(matName)
    thetaName = strcat('theta_mle_DS',num2str(ls_i));
    betaName = strcat('beta_mle_DS',num2str(ls_i));
    eval(strcat(thetaName,'=mu_IDA_t;'));
    eval(strcat(betaName,'=sigma_IDA_t;'));

end

% Par=[theta_mle_DS1, beta_mle_DS1;
%   theta_mle_DS2, beta_mle_DS2;
%   theta_mle_DS3, beta_mle_DS3;
%   theta_mle_DS4, beta_mle_DS4;
%   theta_mle_DS5, beta_mle_DS5];
%save('MRF_par','Par')
Par=[theta_mle_DS1, beta_mle_DS1;
  theta_mle_DS2, beta_mle_DS2;
  theta_mle_DS3, beta_mle_DS3;
  theta_mle_DS4, beta_mle_DS4];
%save('BF_par','Par')

%% compute fragility functions using estimated parameters
xlim_sup=4;
x_vals = 0:0.01:xlim_sup; % IM levels to plot fragility function 

p_mle_DS1 = normcdf((log(x_vals/theta_mle_DS1))/beta_mle_DS1);
p_mle_DS2 = normcdf((log(x_vals/theta_mle_DS2))/beta_mle_DS2); 
p_mle_DS3 = normcdf((log(x_vals/theta_mle_DS3))/beta_mle_DS3); 
p_mle_DS4 = normcdf((log(x_vals/theta_mle_DS4))/beta_mle_DS4); 
% p_mle_DS5 = normcdf((log(x_vals/theta_mle_DS5))/beta_mle_DS5); 


figure('OuterPosition',[100 100 1100 450]);
plot(x_vals,p_mle_DS1, '-', 'linewidth', 2, 'color', [0    0.4470    0.7410])
hold on 
plot(x_vals,p_mle_DS2, '-', 'linewidth', 2)
hold on 
plot(x_vals,p_mle_DS3, '.-', 'linewidth', 2)
hold on 
plot(x_vals,p_mle_DS4, '-', 'linewidth', 2)
hold on 
% plot(x_vals,p_mle_DS5, '-', 'linewidth', 2)
xlim([0 xlim_sup])
title('Fragility Curves','Interpreter','Latex')
set(gca,'fontname','Times')
ylabel('$P(DS>ds|pga)$','Interpreter', 'latex')
xlabel('$pga [g]$','Interpreter', 'latex')
hhh=legend('$DS1$','$DS2$','$DS3$','$DS4$','$DS_c$');
set(hhh,'Interpreter','Latex');
set(hhh,'Location','SouthEast')
set(hhh,'FontSize',18)
title('ii)','Interpreter','Latex')
grid on
legend boxoff
grid on
tightfig;




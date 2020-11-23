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
MDOF_properties_BW_MRF

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
%% Beginning of possible function ----Z    
clear EDP SCALE
ii=0; % counter 
neq = numel(NN); % number of earthquakes
PGA(1:neq) = 0;

for nn =1: numel(NN)
    n = NN(nn);
    disp('nn'); disp(nn);
%     a_g_norm = Ground_motions.SS_Frot_1{n}/max(abs(Ground_motions.SS_Frot_1{n})); % normalized ground motion (read only SS= strike slip)
    a_g_norm = Ground_motions.SS_Frot_1{n};
    PGA(nn) = max(abs(Ground_motions.SS_Frot_1{n}));  % PGA of the earthquake                                           
    edp = 0;
    EDP = [];
    SCALE = [];
    scale = 0.40;                                      % starting scaling point 
    Mat.dt = DT(n);                                   % read the time delt t for the given earthquake
    Mat.t = ((1:numel(a_g_norm ))-1)*Mat.dt;          % time for the givent earthquake
    %
    if flagGM == 1
        figure
        plot(Mat.t,a_g_norm)
    end
    %
    i=0;                                              % Number of time history analysis 
    tic
    while edp<edpc %&& scale <= 2.2
          i = i + 1;
          scale = scale+0.10;                         % increment the scale factor
          %% Initial condition         
          Mat.dFe=zeros(Mat.NDOF,numel(a_g_norm));    % Preallocation for the load for the time series (you do not need to touch this)
          a_g = a_g_norm*scale;                       % Scaled ground motion
          Mat.Fe=Mat.M*Mat.r'*a_g'*g;
          %% Computation response 
          [HistVarBw]=ResponceMDF_Bw(Mat); 
          edp = max(abs(HistVarBw.eps(1,:)));
          EDP(i) = edp; %#ok<SAGROW> store the EDP fpr each time histroy analysis 
          SCALE(i) = scale; %#ok<SAGROW> strore the scale factor for each time history analysis
          disp('EDP')
          disp(edp)
          disp(' ')

%% plot the sequence
          if nn == 2
            figure('OuterPosition',[100 100 1200 500]);
            subplot(2,5,[1 2 3])
            plot(Mat.t,a_g_norm*scale,'Color',[105,105,105]/255)
            hold on 
            xlabel('$t$ [s]')
            ylabel('$\ddot x_g(t)$[g]')
            ylim([-1.5 1.5])
            subplot(2,5,[6 7 8])
            plot(Mat.t,HistVarBw.eps(1,:),'Color',[105,105,105]/255)
            ylim([-0.13 0.13])
            xlabel('$t$ [s]')
            ylabel('$x(t)$[m]')
            subplot(2,5,[4 5 9 10])
            plot(EDP,SCALE,'o-','Color',[169,169,169]/255)
            hold on 
            plot(edp,scale,'o','Color',[255, 0, 0]/255,'MarkerFaceColor','r','MarkerSize',7)
            hold on 
            grid on 
            xlim([0 0.13])
            ylim([0 1.5])
            xlabel('$edp$')
            tightfig
            ptex = [imagesDir 'IM2_' num2str(i)];
            saveas(gcf,ptex,'jpg')
          end
    end
    time_t(nn) = toc;
    disp('time')
    disp(time_t(nn));
    Numb(nn) = i
    disp('scale'); disp(scale);
    EDP_gm{nn} = EDP; %#ok<SAGROW> store the all seqeunce of edp
    SCALE_gm{nn} = SCALE; %#ok<SAGROW> store the all scale of EDP
    PGA_gm{nn} = max(abs(a_g)); %#ok<SAGROW> store the all scale of EDP     
    figure(999)
    plot(EDP,SCALE);
    hold on
    clear EDP SCALE 
    s_plotQvsID
    
end  



%% Plot IDA

close all 
figure('OuterPosition',[100 100 800 500]);
for n=1:numel(NN)
    IM_t(n) = max(SCALE_gm{n});
    disp('')
    plot(EDP_gm{n},SCALE_gm{n},'Color',[105,105,105]/255)
    hold on
    plot(ls_val,max(SCALE_gm{n}),'o','Color',[255, 0, 0]/255,'MarkerFaceColor','r','MarkerSize',7)
    hold on
end
xlim([0 0.15])
ylim([0 3])
xlim([0 ls_val])
xlabel('$edp$')
ylabel('$im=pga [g]$')
grid on 
%tightfig
isd_I(i)=max(abs(HistVarBw.eps(1,:)));

IM_t_c = sort(IM_t);

%% Untruncated IDA
[parmhat,parmci] = lognfit(IM_t_c ,0.01);
mu_IDA = parmhat(1);
sigma_IDA = parmhat(2);

%% 
IM_max = 2.2;
IM_trunc = IM_t_c(IM_t_c < IM_max); % take only the results with IM < IM_max
eq_over = sum(IM_t_c >= IM_max);    % number of analyses reached IM_max without collapsing

% Maximum likelihood fit, using equation 17-18 Lecture notes 8
[mu_IDA_t, sigma_IDA_t ] = truncated_ida(IM_trunc, IM_max, eq_over);

string2save = strcat('LS',num2str(ls_i),'.mat');
save(string2save,'mu_IDA_t','sigma_IDA_t','mu_IDA','sigma_IDA','IM_t_c','IM_trunc','IM_max','eq_over',...
                 'IM_t','EDP_gm','SCALE_gm','PGA_gm')
             
return             
%% Plot fragility curves for every identified limit state
Par=[theta_mle_DS1, beta_mle_DS1;
  theta_mle_DS2, beta_mle_DS2;
  theta_mle_DS3, beta_mle_DS3;
  theta_mle_DS4, beta_mle_DS4;
  theta_mle_DS5, beta_mle_DS5];
%save('Ductile_par','Par')

%% compute fragility functions using estimated parameters
x_vals = 0:0.01:4; % IM levels to plot fragility function 

p_mle_DS1 = normcdf((log(x_vals/theta_mle_DS1))/beta_mle_DS1);
p_mle_DS2 = normcdf((log(x_vals/theta_mle_DS2))/beta_mle_DS2); 
p_mle_DS3 = normcdf((log(x_vals/theta_mle_DS3))/beta_mle_DS3); 
p_mle_DS4 = normcdf((log(x_vals/theta_mle_DS4))/beta_mle_DS4); 
p_mle_DS5 = normcdf((log(x_vals/theta_mle_DS5))/beta_mle_DS5); 


figure('OuterPosition',[100 100 1100 450]);
plot(x_vals,p_mle_DS1, '-', 'linewidth', 2, 'color', [0    0.4470    0.7410])
hold on 
plot(x_vals,p_mle_DS2, '-', 'linewidth', 2)
hold on 
plot(x_vals,p_mle_DS3, '.-', 'linewidth', 2)
hold on 
plot(x_vals,p_mle_DS4, '-', 'linewidth', 2)
hold on 
plot(x_vals,p_mle_DS5, '-', 'linewidth', 2)
xlim([0 4])
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







%______________PROBABILISTIC SEISMIC HAZARD ANALYSIS _____________________
close all
clear all
clc
%Set up visualization
set(0,'DefaultFigureColor',[1 1 1])
set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',20)
%Set up general purpouse
global maindir zs9dir imagesdir
maindir = fullfile('D:\Dati\GoogleDrive_Uni_DD\GIT\ISPS\2020\01_PSHA\00_MATLAB');
zs9dir = fullfile(maindir,'\Italia');
imagesdir = fullfile(maindir,'\images');
cd(maindir)

%% ----------- I. Earthquake source characterization ------------------- %%
prompt = 'Has source characterization already been evaluted? Type Y/N [Y]: ';
flag_choose = input(prompt,'s');
if isempty(flag_choose)
    flag_choose = 'Y';
    load('seismic_source.mat') % già compilata la sorgente
    
elseif flag_choose == 'N'
    %Set up the source site
    IDS=            xls2mtl('Dati','Data','L7','Q7')';
    nsorg=          length(IDS);                            %Source number
    lambdaS=        xls2mtl('Dati','Data','F8','F9');       %rate of exceedance
    mminS=          xls2mtl('Dati','Data','E8','E9');       %minimum magnitude
    mmaxS=          xls2mtl('Dati','Data','D8','D9');       %maximum magnitude
    bS=             xls2mtl('Dati','Data','C8','C9');       %coefficient b
    hS=             xls2mtl('Dati','Data','G8','G9');       %Location: shallow intermediate or deep
    SS=             xls2mtl('Dati','Data','H8','H9');       %Site: 1 = rigid
    fS=             xls2mtl('Dati','Data','I8','H9');       %Fault mechanism
    
    [xsor,ysor]=        coordinate_source('Dati','Coordinate','B3','S4');  %Draw polygons of zonations
    
    dati=[IDS,lambdaS,mminS,mmaxS,bS,hS,SS,fS];
    SOR=cell(nsorg,1);
    save('seismic_source.mat','dati','xsor','ysor','SOR','IDS', ...
        'nsorg','lambdaS','mminS','mmaxS','bS','hS','SS','fS')
else
    disp('Loading seismic_source.mat')
    load('seismic_source.mat')                          % già compilata la sorgente
    
end

% Visualize seismic area selected
prompt = 'Select to change visualization of seismic area? Type Y/N [N]: ';
flag_choose = input(prompt,'s');
if isempty(flag_choose)
    flag_choose = 'N';
end
cd(zs9dir)
mappa

%% ----------- (N.B.): Global definitions            ------------------- %%
%% a. Nfault - defining faults characteristics
% Select the number of faults involved in calculations
flag_choose = 1;
prompt = 'How many sources? Type 1/2 [1]: ';
flag_choose = input(prompt);

N_faults = flag_choose;
% Magnitude distribution parameters
for nof = 1:N_faults
    %Fault 1
    M_min = dati(nof,3);
    b  = -dati(nof,5);
    lambda_M_min = dati(nof,2);
    M_max = dati(nof,4);
    % collect varibles
    M_PAR{nof} = [M_min M_max b];
    LAMBDA_MIN(nof)=lambda_M_min;
    M_step = (M_max-M_min)/1000;   %0.00161 - 0.0023
    mag = M_min:M_step:(M_max-M_step);
end
%% b. IM - definition
INT_step = 0.01; INT_start = 0.01; INT_end = 1.51;
INT=(INT_start:INT_step:(INT_end-INT_step));                       % Discretization of intensity measure PGA in [g]
imeasure = INT(1):INT_step:INT(end);

%% c. R - definition
flag_plot = 0;  % 1 =y , 0 = n
% Set up the grid coordinate and discretization
Rmax = 0.51;    % 50km
Rmin = 0.;      % site
R_step = 0.01;% step 1 km
rho = Rmin:R_step:(Rmax-R_step);
theta_step = pi/25;
theta = 0:theta_step:(2*pi);
% initialize coordinate in polar system
grid_cp = zeros(size(theta,2)*size(rho,2),2);
for i =1:size(rho,2)
    if i == 1
        for j = 1:size(theta,2)
            grid_cp(j,1) = rho(i);
            grid_cp(j,2) = theta(j);
        end
    else
        for j = 1:size(theta,2)
            grid_cp(j+(i-1)*(size(rho,2)),1) = rho(i);
            grid_cp(j+(i-1)*(size(rho,2)),2) = theta(j);
            if flag_plot ==1
                figure(1002);
                hold on
                plot(grid_cc(j+(i-1)*(size(rho,2)),1),grid_cc(j+(i-1)*(size(rho,2)),2),'b*')
            end
        end
    end
end
% initialize coordinate in cartesian system
grid_cc = zeros(size(theta,2)*size(rho,2),2);
[grid_cc(:,1),grid_cc(:,2)] = pol2cart(grid_cp(:,2),grid_cp(:,1));
r_distance = rho;
% Mapping
if flag_plot ==1
    figure(1003)
    subplot(1,2,1)
    plot(grid_cp(:,1),grid_cp(:,2),'r*')
    xlabel('Cartesian coordinates')
    subplot(1,2,2)
    plot(grid_cc(:,1),grid_cc(:,2),'r*')
    xlabel('Polar coordinates')
    sgtitle('Mapping','fontsize',20)
end

%% ----------- Integral of Hazard                     ------------------ %%
% initialization of cell array
cont_plot = 1; flag_plot = 1; % 1 =y , 0 = n
P_i = 0;
Prim_mag = cell(numel(N_faults),numel(r_distance),numel(INT));
lambda_mag = cell(numel(N_faults),numel(r_distance));

for k = 1:N_faults % foreach source
    % Magnitude's PDF
    M_par = M_PAR{k};    %read Fault parameters
    lambda_M_min = LAMBDA_MIN(k);
    
    m0 = M_par(1);
    mu = M_par(2);
    b = M_par(3);
    beta = b*log(10);
    
    % PDF of M
    Nm = 1000; % No. of discretized points between m0 and mu for numerical integration
    M = linspace(m0,mu,Nm);
    fMc = beta*exp(-beta*(M-m0))/(1-exp(-beta*(mu-m0)));
    dM=M(2)-M(1);
    Mdiscr=fMc*dM;
    FMc = (1-10.^(-b*(M-m0)))/(1-10.^(-b*(mu-m0)));
    
    PMm = zeros(1,Nm);
    FM = zeros(1,Nm);
    for iifm = 2:(Nm)
        FM(iifm)= FM(iifm-1)+Mdiscr(iifm-1);
        PMm(iifm) = (FM(iifm)-FM(iifm-1));
    end
    
    %     if flag_plot ==1
    %         figure;
    %         subplot(1,2,1)
    %         plot(M,Mdiscr)
    % %         xlim([m0 mu+1]); xlabel('Magnitude'); ylabel('$f_M (m)$')
    %         subplot(1,2,2)
    %         histogram(PMm,'normalization','probability')  % GR bounded
    % %         xlim([m0 mu+1]); xlabel('Magnitude'); ylabel('$P(M = m)$')
    %         sgtitle('Continuous vs discrete PDF of magnitude distribution','fontsize',20)
    %     end
    
    if flag_plot ==1
        figure;
        subplot(1,2,1)
        plot(M,fMc)
        xlabel('magnitude'); ylabel('$PDF: f_M (m)$')
        
        subplot(1,2,2)
        hold on
        plot(M,FMc); plot(M,FM)
        xlabel('magnitude'); ylabel('$CDF: F_M (m)$')
        legend({'cont','discr'},'Location','southeast','interpreter','latex')
        legend('boxoff')
        box on
    end
    
    if abs(sum(Mdiscr)-1)>0.02
        warning('Numerical integration is off the safty boundary of 2%')
        disp('Numerical integration:')
        disp(sum(Mdiscr))
    end
    Mdiscr=Mdiscr/sum(Mdiscr);
    
    if flag_plot ==1
        figure;
        semilogy(M,Mdiscr)  % GR bounded
        xlabel('magnitude'); ylabel('GR recurrence law')
        xlim([m0 mu])
    end
    
    lambda_m = lambda_M_min*(exp(-beta*(M-m0))-exp(-beta*(mu-m0)))/(1-exp(-beta*(mu-m0)));
    
    if flag_plot ==1
        figure;
        semilogy(M,lambda_m)  % mean annual rate of exceedance
        xlabel('magnitude'); ylabel('$\lambda_m$, rate of magn. $> m $')
        xlim([m0 mu]); grid minor
    end
    
    % PDF of R
    rmax = Rmax*100;    % 50km
    rmin = Rmin.*100;   % site
    
    r_step = R_step;
    r = rmin:r_step:(rmax-r_step);
    Fr = r.^2./rmax.^2;
    
    fr = 2*r./rmax.^2;
    DFr = r_step;
    frdiscr = zeros(1,numel(Fr));
    Frdiscr = zeros(1,numel(Fr));
    for idr = 2:numel(Fr)
    frdiscr(idr) = -(Fr(idr-1)-Fr(idr))./(DFr);
    Frdiscr(idr) = Frdiscr(idr-1)+Fr(idr);
    end
    
    if flag_plot ==1
        figure;
        subplot(1,2,1) 
        hold on
        plot(r,fr); plot(r,frdiscr);
        xlabel('r'); ylabel('$PDF: f_R (r)$')
        legend({'cont','discr'},'Location','southeast','interpreter','latex')
        legend('boxoff')
        box on

        subplot(1,2,2)
        hold on
        plot(r,Fr);plot(r,Frdiscr./max(Frdiscr))
        xlabel('r'); ylabel('$CDF: F_R (r)$')
        legend({'cont','discr'},'Location','southeast','interpreter','latex')
        legend('boxoff')
        box on
    end
    
    
    % Beginning LOOP
    Ptemp = cell(numel(r_distance),numel(INT),numel(mag));
    Ptemp_IMRM = zeros(1,numel(imeasure));
    for iir = 2:(numel(r_distance)) %Rmin:R_step:Rmax % per ogni R
        r_i = r_distance(iir);
        disp('r = ')
        disp(r_i*100)
        
        for im = 1:(numel(INT)) %INT(1):INT_step:INT(end) % per ogni IM
            im_i = INT(im);
            %         Ptemp_RM = zeros(1,numel(r_distance));
            %             disp('IM = ')
            %             disp(im_i);
            
            Ptemp_im = zeros(1,numel(mag));
            Ptemp_m = zeros(1,numel(mag));
            for iimag = 1:numel(mag) % per ogni magnitude
                mag_i = mag(iimag);
                [mean_im, sigma_im] = GMPE(mag_i,r_i);
                Ptemp_im(iimag) =  1 - cdf('Normal',log(im_i),log(mean_im),sigma_im); % P(IM>im | M=mj, R=ri)
                Ptemp_m(iimag) = Mdiscr(iimag);  % P(M=mj)
                Ptemp_r(iimag) = frdiscr(iir);  % P(R=ri)
                P_i(iimag) = Ptemp_im(iimag)*Ptemp_m(iimag)*Ptemp_r(iimag); % P(IM>im | M=mj, R=ri)*P(M=mj)*P(R=ri)
                %P_i = Mdiscr*Pint_M_R*frdiscr';
                %Pint_M_R=1-normcdf(log(int),log(mean_im),sigma_im);
                %P_i=P_i+sum(Pint_M_R.*frdiscr)*Mdiscr(i);
                Ptemp{iir,im,iimag} = P_i(iimag);
            end
            Prim_mag{k,iir,im} = sum(cell2mat(squeeze(Ptemp(iir,im,:))));
            
            
        end % per ogni im
        lambda_mag{k,iir-1} = [cell2mat(squeeze(Prim_mag(k,iir,:)))]*lambda_M_min;
        P_ex_time1{k,iir-1} = 1-exp(-lambda_mag{k,iir-1}*1);                         % Annual Probability of exeedance
        P_ex_time50{k,iir-1} = 1-exp(-lambda_mag{k,iir-1}*50);                       % 50 years Probability of exceedance
        P_ex_time475{k,iir-1} = 1-exp(-lambda_mag{k,iir-1}*475);                     % 475 years Probability of exceedance
        
        
        if rem(iir-1,10) == 0
            figure(101);
            hold on
            set(gca, 'XScale', 'log', 'YScale', 'log');
            h1(cont_plot) = plot(INT,lambda_mag{k,iir-1}');  % mean annual rate of exceedance
            title('Annual frequency of exceedance');
            xlabel('im'); ylabel('$\lambda (IM>im,1 year)$');
            grid on
            Vr = [120 200 1898 3900 200];
            legendInfo{cont_plot} = ['R = ' num2str(r_i*100) 'km'];
            legend(legendInfo);
%             if iir == numel(r_distance) %last run
%                 plot(INT,ones(1,numel(INT))/Vr(1),'--r','linewidth',1.1,'Color',[128 128 0]/255)
%                 txt2 = 'frequent';
%                 text(INT(2),1/Vr(1),txt2,'FontSize',16)
%                 plot(INT,ones(1,numel(INT))/Vr(2),'--r','linewidth',1.1,'Color',[255,140,0]/255)
%                 txt2 = 'occasional';
%                 text(INT(2),1/Vr(2),txt2,'FontSize',16)
%                 plot(INT,ones(1,numel(INT))/Vr(3),'--r','linewidth',1.1,'Color',[255,69,0]/255)
%                 txt2 = 'rare';
%                 text(INT(2),1/Vr(3),txt2,'FontSize',16)
%                 plot(INT,ones(1,numel(INT))/Vr(4),'--r','linewidth',1.1,'Color',[255,69,0]/255)
%                 txt2 = 'very rare';
%                 text(INT(2),1/Vr(4),txt2,'FontSize',16)
%             end
            
            figure(102);
            hold on
            set(gca, 'XScale', 'log', 'YScale', 'log');
            h2(cont_plot) = plot(INT,P_ex_time1{k,iir-1});  % annual hazard curve
            title('Annual hazard curve');
            xlabel('im'); ylabel('$P(IM>im, 1 year)$');
            grid on
            legend(legendInfo);
            
            figure(103);
            hold on
            set(gca, 'XScale', 'log', 'YScale', 'log');
            h3(cont_plot) = plot(INT,P_ex_time50{k,iir-1}); % 50 years hazard curve
            title('50 years hazard curve');
            xlabel('im'); ylabel('$P(IM>im, 50 year)$');
            grid on
            legend(legendInfo);
            
            figure(104);
            hold on
            set(gca, 'XScale', 'log', 'YScale', 'log');
            h3(cont_plot) = plot(INT,P_ex_time475{k,iir-1}); % 50 years hazard curve
            title('475 years hazard curve');
            xlabel('im'); ylabel('$P(IM>im, 475 year)$');
            grid on
            legend off
            if iir == numel(r_distance) %last run
                plot(INT,ones(1,numel(INT))/Vr(1),'--r','linewidth',1.1,'Color',[255 69 0]/255)
                txt2 = 'frequent';
                text(INT(2),1/Vr(1),txt2,'FontSize',16)
                plot(INT,ones(1,numel(INT))/Vr(2),'--r','linewidth',1.1,'Color',[255,69,0]/255)
                txt2 = 'occasional';
                text(INT(2),1/Vr(2),txt2,'FontSize',16)
                plot(INT,ones(1,numel(INT))/Vr(3),'--r','linewidth',1.1,'Color',[255,69,0]/255)
                txt2 = 'rare';
                text(INT(2),1/Vr(3),txt2,'FontSize',16)
                plot(INT,ones(1,numel(INT))/Vr(4),'--r','linewidth',1.1,'Color',[255,69,0]/255)
                txt2 = 'very rare';
                text(INT(2),1/Vr(4),txt2,'FontSize',16)
            end
            legend(legendInfo);

            
            
            cont_plot = cont_plot +1;
        end
        
        %         Ptemp_IMRM(1,im) = Ptemp_RM;
    end % per ogni R
end % per ogni sorgente
disp('Hazard analysis completed')
disp('* ----------------------------- *')


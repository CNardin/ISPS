edpc = ls_val;

clear EDP SCALE
ii=0; % counter 
neq = numel(NN); % number of earthquakes
PGA(1:neq) = 0;

for nn =1: numel(NN)
    n = NN(nn);
    disp('nn'); disp(nn);
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

          disp('EDP')
          disp(edp)
          disp(' ')

    end
    time_t(nn) = toc;
    disp('time')
    disp(time_t(nn));
    Numb(nn) = i
    disp('scale'); disp(scale);

    figure(999)
    plot(EDP,SCALE);
    hold on
    clear EDP SCALE 
%     s_plotQvsID
    
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


[parmhat,parmci] = lognfit(IM_t_c ,0.01);
mu_IDA = parmhat(1);
sigma_IDA = parmhat(2);


IM_max = 2.2;
IM_trunc = IM_t_c(IM_t_c < IM_max); % take only the results with IM < IM_max
eq_over = sum(IM_t_c >= IM_max);    % number of analyses reached IM_max without collapsing

% Maximum likelihood fit, using equation 17-18 Lecture notes 8
[mu_IDA_t, sigma_IDA_t ] = truncated_ida(IM_trunc, IM_max, eq_over);

string2save = strcat('LS',num2str(ls_i),'.mat');
save(string2save,'mu_IDA_t','sigma_IDA_t','mu_IDA','sigma_IDA','IM_t_c','IM_trunc','IM_max','eq_over',...
                 'IM_t','EDP_gm','SCALE_gm','PGA_gm')
             
          




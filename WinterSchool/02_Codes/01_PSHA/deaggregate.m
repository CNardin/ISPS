temp = [];
icol = 1;

for icol = 1:numel(P_ex_time50)
temp = [temp  P_ex_time50{1,icol}(:)];
end

PSHA_P50 = sum(temp,2);
lambdaINGV=xlsread('ntc.csv','ntc','A1:D9');
plotPSHA_P50 = interp1(imeasure,PSHA_P50,lambdaINGV(:,3))
S={lambdaNTC(:,3),lambdaNTC(:,1),lambdaNTC(:,3),plotPSHA_P50};
log_plot_dati(S,"Confronto",["INGV","PSHA"])


% dea2 = [];
% magn ={}; dist ={};
% save('deaggregate.mat','dea2','magn','dist')
% clear all; clc
load('deaggregate.mat')
magn = {'3.5-4.0', '4.0-4.5', '4.5-5.0','5.0-5.5', '5.5-6.0','6.0-6.5', '6.5-7.0','7.0-7.5', '7.5-8.0','8.0-8.5','8.5-9.0' };
dist = {'0-10', '10-20','20-30', '30-40','40-50', '50-60','60-70', '70-80','80-90', ...
        '90-100','100-110','110-120','120-130','130-140','140-150','150-160','160-170','170-180','180-190','190-200'};

    
figure;
set(gcf,'Color',[1 1 1])
b=bar3(dea2);
xlabel('Magnitude [-]')
ylabel('Distance [km]')
set(gca,'xTickLabel',magn,'XtickLabelRotation',-30)
set(gca,'YTick',[1:1:20],'yticklabel',dist,'YTickLabelRotation',30);
set(gca,'Fontsize',10)
colorbar
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
colormap(hot)
hold on


N_pcolor = dea2;
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(0,11,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(0,20,size(N_pcolor,1)); % Rows of N_pcolor

h = pcolor(xl,yl,N_pcolor);
colormap('jet') % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];





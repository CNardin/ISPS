% dea2 = [];
% magn ={}; dist ={};
save('deaggregate.mat','dea2','magn','dist')
clear all; clc
load('deaggregate.mat')
magn = {'3.5-4.0', '4.0-4.5', '4.5-5.0','5.0-5.5', '5.5-6.0','6.0-6.5', '6.5-7.0','7.0-7.5', '7.5-8.0','8.0-8.5','8.5-9.0' };
dist = {'0-10', '10-20','20-30', '30-40','40-50', '50-60','60-70', '70-80','80-90', ...
        '90-100','100-110','110-120','120-130','130-140','140-150','150-160','160-170','170-180','180-190','190-200'};

fontsize = 10;
fontname = 'Helvetica';
myscale = 1;
    
    
figure;
set(gcf,'Color',[1 1 1],'Units','centimeters','Position',[0,0,19,19]/myscale)
set(gca,'FontSize',fontsize,'FontName',fontname);
bar3(dea2);
xlabel('Magnitude [-]')
ylabel('Distance [km]')
set(gca,'xTickLabel',magn)
set(gca,'YTick',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],'yticklabel',dist);

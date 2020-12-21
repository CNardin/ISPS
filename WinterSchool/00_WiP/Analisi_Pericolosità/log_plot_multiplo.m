function log_plot_multiplo(x,y,titolo)
%stampa in un piano bi-logaritmico

fontsize = 10;
fontname = 'Helvetica';
myscale = 1;
nplot=size(y);

figure()
set(gcf,'Color',[1 1 1],'Units','centimeters','Position',[0,0,19,19]/myscale)
set(gca,'FontSize',fontsize,'FontName',fontname);
for j=1:nplot(2)
    if j==1
        loglog(x,y(:,j),'k-s','Markersize',6,'MarkerFaceColor','auto','LineWidth',2);
        hold on
    else
%     loglog(x,y(:,j),'k-s','Markersize',6,'MarkerFaceColor','auto','LineWidth',0.5);
    loglog(x,y(:,j),'Markersize',6,'MarkerFaceColor','auto','LineWidth',0.5);
    hold on
    end
end
legend('Pericolosità totale','Pericolosità sorgenti')
xlabel('PGA [g]','FontSize',fontsize,'FontName',fontname);
ylabel('Tasso medio di eccedenza annuo \lambda ','FontSize',fontsize,'FontName',fontname);
title(titolo)
grid on
box on
xlim([10^(-2) 1])
ylim([10^(-5) 10^(-1)])
end
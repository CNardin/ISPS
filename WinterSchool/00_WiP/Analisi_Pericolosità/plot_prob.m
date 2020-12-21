function plot_prob(x,y,titolo,xlab,ylab)
%stampa in 2D

fontsize = 10;
fontname = 'Helvetica';
myscale = 0.8;
nplot=size(y);

figure()
set(gcf,'Color',[1 1 1],'Units','centimeters','Position',[0,0,19,19]/myscale)
set(gca,'FontSize',fontsize,'FontName',fontname);
for j=1:nplot(2)
    if j==1
        plot(x,y(:,j),'k-s','Markersize',6,'MarkerFaceColor','auto','LineWidth',2);
        hold on
    else
    plot(x,y(:,j),'k-s','Markersize',6,'MarkerFaceColor','auto','LineWidth',0.5);
    hold on
    end
end
xlabel(xlab,'FontSize',fontsize,'FontName',fontname);
ylabel(ylab,'FontSize',fontsize,'FontName',fontname);
title(titolo)
grid on
box on

end
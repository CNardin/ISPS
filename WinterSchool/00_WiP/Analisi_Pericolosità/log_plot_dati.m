function log_plot_dati(S,titolo,legenda)
%stampa in un piano bi-logaritmico

fontsize = 10;
fontname = 'Helvetica';
myscale = 1.;
[~,m]=size(S);

figure()
set(gcf,'Color',[1 1 1],'Units','centimeters','Position',[0,0,19,19]/myscale)
set(gca,'FontSize',fontsize,'FontName',fontname);

for i=1:2:m-1
    A=S(i);
    B=S(i+1);
    loglog(A{:}, B{:},'LineWidth',2);
%         loglog(A{:}, B{1,1}(:,6)+ B{1,1}(:,3)+ B{1,1}(:,5),'LineWidth',2);
    hold on
end
xlabel('PGA [g]','FontSize',fontsize,'FontName',fontname);
ylabel('Tasso medio di eccedenza annuo \lambda ','FontSize',fontsize,'FontName',fontname);
title(titolo)
legend(legenda)
grid on
box on

end
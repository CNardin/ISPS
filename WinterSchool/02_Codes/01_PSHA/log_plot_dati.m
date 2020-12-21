function log_plot_dati(S,titolo,legenda)
%stampa in un piano bi-logaritmico

[~,m]=size(S);

figure()
set(gcf,'Color',[1 1 1])

for i=1:2:m-1
    A=S(i);
    if i == m-1
    B=S(i+1);
    loglog(A{:}, B{:}/2,'LineWidth',2);
    else
    B=S(i+1);       
    loglog(A{:}, B{:},'LineWidth',2);
    end
    hold on        
end
xlabel('PGA [g]');
ylabel('Tasso medio di eccedenza annuo $\lambda$ ');
title(titolo)
legend(legenda)
grid on
box on
% xlim([0.05 0.45])

end
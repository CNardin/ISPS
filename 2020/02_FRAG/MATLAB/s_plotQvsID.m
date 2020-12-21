figure()
Mat.h = 4; %[m]
plot(HistVarBw.eps(1,:),HistVarBw.Q(1,:)./10^3)
xlabel('ID [m]', 'Interpreter','Latex')
ylabel('Q [kN]', 'Interpreter','Latex')

% switch system_type
%     case 'le'
%         xline(edpc,'--r',{'Collapse','LS'},'LineWidth',1.5); xline(-edpc,'--r',{'Collapse','LS'},'LineWidth',1.5);
%     case 'bw'
%         xline(edpc,'--r',{'Collapse','LS'},'LineWidth',1.5); xline(-edpc,'--r',{'Collapse','LS'},'LineWidth',1.5);
%     otherwise
%         disp('No system type defined.')
%         warning('You should stop the analyses: no sense output.')
%         return
% end

% figure
% plot(HistVarBw.u(1,:))
% hold on
% plot(HistVarBw.eps(1,:))

% plot([LS( 0.12],[-2500 2500],'--','linewidth', 0.8,'color',[0.64 0.08 0.18])
% plot([-0.12 -0.12],[-2500 2500],'--','linewidth', 0.8,'color',[0.64 0.08 0.18])
% % if ls_i == numel(LS)
% % for ls_i = 1:numel(LS)
% % ls_val = LS(ls_i);
% % string2plot = strcat('D# ',num2str(ls_i),' Limit state');
% % xline(ls_val,'--r',string2plot,'LineWidth',1.5); xline(-ls_val,'--r',string2plot,'LineWidth',1.5);
% % ls_i = ls_i + 1;
% % end
% % end
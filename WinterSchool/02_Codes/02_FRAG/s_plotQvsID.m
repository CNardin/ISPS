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

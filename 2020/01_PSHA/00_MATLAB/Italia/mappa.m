switch flag_choose
    case 'N'
        %% Load previously saved plot
        cd(imagesdir)
        openfig('fig1000.fig');
        openfig('fig1001.fig');
        cd(maindir)
    otherwise
        %% Plot in real time
        % Plot in real time info about ZS9 all
        filename_map = 'PRO.shp';
        close all
        
        S = shaperead(filename_map);
        figure(1000)
        hold on
        mapshow(S,'LineStyle',':','FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
        
        filename = 'zs9.shp';
        ZS9 = shaperead(filename);
        selectSource= [ZS9(20,:)  ...
            ZS9(24,:)  ...
            ];
        mapshow(ZS9,'FaceColor', [1 1 1], 'EdgeColor', 'black');
        mapshow(selectSource,'FaceColor', [0.4 1 0.6], 'EdgeColor', 'black');
        
        Palmoli.x = 14.253;
        Palmoli.y = 41.927;
        h1= plot(Palmoli.x,Palmoli.y,'LineStyle','none','Marker','^','MarkerSize',5,'Color','r');
        set(h1, 'markerfacecolor', get(h1, 'color')); % Use same color to fill in markers
        h2 = viscircles([Palmoli.x,Palmoli.y], 3.0,'LineStyle','--','color','r');
        hold off
        xlabel('Longitude');ylabel('Latitude');title('ZS9')
        
        % Plot info only for the simple source case
        figure(1001)
        hold on
        mapshow(S,'LineStyle',':','FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
        Palmoli.x = 14.483;
        Palmoli.y = 41.627;
        h1= plot(Palmoli.x,Palmoli.y,'LineStyle','none','Marker','^','MarkerSize',5,'Color','r');
        set(h1, 'markerfacecolor', get(h1, 'color')); % Use same color to fill in markers
        
        h3 = viscircles([Palmoli.x,Palmoli.y], 0.5,'LineStyle','--','color','b');
        xlabel('Longitude');ylabel('Latitude');title('Palmoli, Lon. 14.483 - Lat. 41.627');
        grid on
        grid minor
        % Create textbox
        annotation(figure(1001),'textbox',[0.626071428571427 0.585190475282218 0.15377012304371 0.0619047628130232],...
            'String',{'R = 50 km'},'LineStyle','none','Interpreter','latex');
        cd(maindir)
end

%*************************ANALISI DI PERICOLOSITA'*************************

%**********************************NOTE************************************
%Inserire punti delle coordinate traslati in modo che sito sia origine
%    Lat=[37.0819722222222; 41.939023]; %Siracusa ; Palmoli
%    Lon=[15.2852777777778; 14.58157];
%    [x,y,utmzone] = deg2utm(Lat,Lon);
%    fprintf('%7.0f ',[x,y]); fprintf(utmzone)


%*********************************IPOTESI***********************************
%Calcolo delle aree con metodo di Montecarlo (lento a convergere)
%Legge di attenuazione di Sabetta e Pugliese

%***********************************UDM*************************************
%Acc    [g]
%L      [km]

%*************************STRUTTURA VARIABILI*******************************
%Vettore sorgente è composto così
%SORGENTE=[ID(1x1); lambda{1x1};    Pmodel{nm x nr x 3};   P(M=m) {nm x 2};    P(R=r) {nr x 2}]

                    %({r1,m1,p}    ...     {r1,mn,p})
%Pmodel=P(IM>x|m,r)=%(...                       ... )
                    %({rn,m1,p}    ...     {rn,mn,p})

                    %({m1,p})
%Pmagn=P(M=m)=      %( ...  )
                    %({mn,p})

                    
%*************************RESET MACCHINA***********************************
clear all
close all
clc
%Set up visualization
set(0,'DefaultFigureColor',[1 1 1])
set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',20)
global flagDebug
flagDebug = 1;                   
%****************************DATI******************************************
nomeFoglio = 'convalida'; % 'convalida' - 'Foglio1'
IDS=            xls2mtl('sorgenti','A2','A30',nomeFoglio);
nsorg=          length(IDS);                     %Numero delle sorgenti
lambdaS=        xls2mtl('sorgenti','B2','B30',nomeFoglio);
mminS=          xls2mtl('sorgenti','C2','C30',nomeFoglio);
mmaxS=          xls2mtl('sorgenti','D2','D30',nomeFoglio);
bS=             xls2mtl('sorgenti','E2','E30',nomeFoglio);
hS=             xls2mtl('sorgenti','F2','F30',nomeFoglio);
SS=             xls2mtl('sorgenti','G2','G30',nomeFoglio);
fS=             xls2mtl('sorgenti','H2','H30',nomeFoglio);

[xsor,ysor]=        coordinate_sorgenti('sorgenti','L2','AY30',nomeFoglio);

dati=[IDS,lambdaS,mminS,mmaxS,bS,hS,SS,fS];
SOR=cell(nsorg,1);


%*************************DEFINIZIONI GRIGLIE*******************************
%DISTANZE
nr=40;
rmin=0.001;
rmax=350.001;
pr=(rmax-rmin)/(nr-1);

%MAGNITUDO CON CONTROLLO PASSI
nm=40;
mmin=min(mminS);
mmax=max(mmaxS);
pm=(mmax-mmin)/(nm-1);
if pm>0.05                                                                  %Limite consigliato 0.1 ma risultati con almeno 25% di errore
    disp('Passo magnitudo elevato')
end

%PGA
pgamin=0.03;
pgamax=0.5;
npga=11;
ppga=(pgamax-pgamin)/(npga-1);


%*************************DEFINIZIONE SORGENTI*****************************
for l=1:nsorg
    SOR{l}=sorgente(dati(l,:),rmin,nr,pr,mmin,nm,pm,xsor(l,:),ysor(l,:));  %Definizione delle quantità fisse per ogni sorgente
end

%*************************TASSI ECCEDENZA*********************************
lambdasor=NaN(npga,2,nsorg);                                               %Contiene la pericolosità sismica di tutte le sorgenti

for k=1:npga
    pgalim=pgamin+ppga*(k-1);
    for l=1:nsorg
        SOR{l}=aggiorna_sorgente(rmin,nr,pr,mmin,nm,pm,pgalim,SOR{l},dati(l,:));%Aggiorna le probabilità di eccedenza di ogni sorgente
        lambdasor(k,:,l)=[pgalim,tasso_eccedenza(SOR{l})];                 %Tasso di eccedenza di ogni sorgente
    end
end

lambdaplot=lambdasor(:,2,:);
lambdaplot=lambdaplot(:,:);
lambdaplot=[sum(lambdaplot,2),lambdaplot];


%******************************STAMPA**************************************
log_plot_multiplo(lambdasor(:,1,1),lambdaplot,'Pericolosità sismica totale')

lambdaNTC=xlsread('ntc.csv','ntc','A1:D9');
S={lambdaNTC(:,3),lambdaNTC(:,1),lambdasor(:,1,1),lambdaplot(:,1)};
log_plot_dati(S,"Confronto",["INGV","PSHA"])
S={lambdaNTC(:,4),lambdaNTC(:,1),lambdasor(:,1,1),lambdaplot};
log_plot_dati(S,"Confronto",["INGV","PSHA"])


%%
xx = SOR{1,1}{1,3}(1,:,2)'; % m
yy = SOR{1,1}{1,3}(:,1,1);  % r
zz =  SOR{1,1}{1,3}(:,:,3);  % p



figure
bar3(zz)
xlabel('Magnitude [-]')
ylabel('Distance [km]')
set(gca,'XTickLabel',{'3-4','4-5','5-6','6-7'})
set(gca,'YTickLabel',[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350])


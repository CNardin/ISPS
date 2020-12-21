function fdp=FDP_discreta_distanza(rmin,nr,pr,xv,yv)
global flagDebug

%Calcola la FDP discreta dellla sitanza con il metodo di Montecarlo.


%**************************COORDINATE POLIGONO*****************************
xv=[xv,xv(1)];                  %chiudo il poligono
yv=[yv,yv(1)];


%**************************DEFINIZIONE GRIGLIE*****************************
%Campionamento Montecarlo (trovo il rettangolo minimo che inscrive il poligono)
lx=abs(max(xv)-min(xv));    
ly=abs(max(yv)-min(yv));
xmin=min(xv);
ymin=min(yv);

area=polyarea(xv,yv);           %calcolo area del poligono


% %Stampa del sistema come controllo
if flagDebug == 1
xsample=[xmin,  xmin+lx,  xmin+lx,      xmin,       xmin];
ysample=[ymin,  ymin,     ymin+ly,      ymin+ly,    ymin];
figure(1000001)
plot(xv,yv,'b',0,0,'O')
hold on
plot(xsample,ysample,'r')
hold off
end

%****************************MONTECARLO************************************
n=1*10^5;                       %Numero di punti
cdf=NaN(nr,2);                  %Inizializzo il vettore della CDF
fdp=NaN(nr,2);                  %Inizializzo il vettore della FDP
d=NaN(n,2);                     %Inizializzo il vettore delle distanze
flag2=false(n,1);               %Inizializzo il vettore dei conteggi


%****************************GENERAZIONE PUNTI*****************************
xq=xmin+(lx*rand(n,1));                 %Creo una serie di punti nel quadrato di interesse con coordinate casuali
yq=ymin+(ly*rand(n,1));                 %Creo una serie di punti nel quadrato di interesse con coordinate casuali
flag1=inpolygon(xq,yq,xv,yv);           %Trovo i punti nel e sul poligono
  
for i=1:n
    d(i)=(xq(i)^2+yq(i)^2)^0.5;
end

for l=1:nr
    r=rmin+pr*(l-1);                        %Itero su cerchi sempre maggiori
    
    for i=1:n
        if d(i)<=r
            flag2(i)=true(1);
        end
    end
    
    flag3=and(flag1,flag2);                 %Trovo i punti che soddisfano sia il cerchio che il poligono
    misura=sum(flag3);                      %Conto i punti
    freq=misura/n;                          %Frequenza relativa dei bunti che soddisfano le relazioni
    areaqua=lx*ly;                          %Area del rettangolo di campionamento
    areain=freq*areaqua;                    %Area della di interesse
    p=areain/area;                          %Probabilità che d<r
    cdf(l,:)=[r,p];

    if flagDebug == 1
        if l==1 || rem(l,5)==0                        %Stampa dei punti dentro e fuori solo come controllo
            stampa_montecarlo(xv,yv,flag3,xq,yq)
        end
    end
end


%**************************CALCOLO PDF CUMULATA*****************************
for i=1:nr   
    if i>1
        p=(cdf(i,2)-cdf(i-1,2));        %Dovrei divididere e moltiplicare per il passo
    else
        p=cdf(1,2);
    end
    fdp(i,:)=[p,cdf(i,1)];
end


%***************************STAMPE E MESSAGGI******************************
err=((sum(flag1)/n*lx*ly)-area)/area*100;
disp(strcat('Errore stimato metodo di Montecarlo= ',num2str(err),'%'))

if flagDebug ==1
    plot_prob(cdf(:,1),cdf(:,2),'Funzione cumulata di probabilita delle distanze','r [km]','Pr(R<r)')
end

end

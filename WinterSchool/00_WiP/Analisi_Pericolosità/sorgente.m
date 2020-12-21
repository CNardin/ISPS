function [SORGENTE]=sorgente(dati,rmin,nr,pr,mmin,nm,pm,xv,yv)
%Genera la matrice che caratterizza ogni sorgente

%******************************INIZIALIZZA*********************************
ID=         dati(1);
lambda=     dati(2);
Pmodel=     NaN(nr,nm,3);
Pmagn=      NaN(nm,2);

%*************************VALUTAZIONI SORGENTI*****************************
%Definisco una finestra di valori di distanza e magnitudo e per ciascuna
%valuto le probabilità che il sisma abbia quella distanza dalla sorgente
%e quella magnitudo i valori di targa sono quelli del limite inferiore
%della finestra

rinf=rmin;
minf=mmin;

Pdist=FDP_discreta_distanza(rmin,nr,pr,xv,yv);                                 %Aggiorno il valore della FDP della distanza

for i=1:nr
    rsup=rmin+pr*(i);
    for j=1:nm
        msup=mmin+pm*(j);
        Pmagn(j,:)=[FDP_discreta_magnitudo(minf,msup,dati),minf];               %Aggiorno il valore della FDP della magnitudo
        Pmodel(i,j,:)=[rinf,minf,NaN];                                          %Definisco la matrice della FDP ma non aggiorno le probabilità
        minf=msup;
    end
    minf=mmin;      %Ripristino il valore iniziale della magnitudo minima per rifare la procedura
    rinf=rsup;
end

SORGENTE={ID,lambda,Pmodel,Pmagn,Pdist};

end
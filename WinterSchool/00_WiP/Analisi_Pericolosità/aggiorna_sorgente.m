function [SORGENTE]=aggiorna_sorgente(rmin,nr,pr,mmin,nm,pm,pgalim,SORGENTE,dati)
%Aggiorna la matrice che caratterizza ogni sorgente

%******************************INIZIALIZZA*********************************
Pmodel=     NaN(nr,nm,3);

%*************************VALUTAZIONI SORGENTI*****************************
%Aggiorno solo le probabilità di superamento

rinf=rmin;
minf=mmin;

for i=1:nr
    rsup=rmin+pr*(i);
    for j=1:nm
        msup=mmin+pm*(j);
        Pmodel(i,j,:)=[rinf,minf,FDP_leggediattenuazione(pgalim,rinf,minf,dati)];    %Aggiorno la matrice della FDP
        minf=msup;
    end
    minf=mmin;                                                                  %Ripristino il valore iniziale della magnitudo minima per rifare la procedura
    rinf=rsup;
end

SORGENTE(1,3)={Pmodel};

end
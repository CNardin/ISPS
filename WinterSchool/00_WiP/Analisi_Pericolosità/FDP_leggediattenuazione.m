function p=FDP_leggediattenuazione(pgalim,r,m,dati)
%Trova la probabilità di eccedenza del valore di pga limite

%Il vetore dati lo tengo nel caso mi serva per cambiere le leggi di
%attenuazione


%[media,sigma]=leggeattenuazione_cornel(r,m);    %Trovo i parametri della distribuzione con la legge di Cornell
[media,sigma]=leggeattenuazione_SP(r,m,dati);    %Trovo i parametri della distribuzione con la legge di Sabetta Pugliese
p=1-normcdf((log10(pgalim)-media)/sigma);         %Trovo la probabilità di eccedenza della pga limite 

end
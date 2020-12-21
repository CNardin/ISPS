function [lambda]=tasso_eccedenza(S)
%Trova il taso di eccedenza medio annuo del valore soglia

lambdas=S{2};                       %Estrae il tasso di superamento del valore minimo della sorgente
M=S{1,3};                           %Matrice della probabilità di superamento assegnati distanza e magnitudo
V1=S{1,4};                          %Vettore di probabilità della magnitudo
V2=S{1,5};                          %Vettore di probabilità della distanza

Pcomb=M(:,:,3);                     %Estrae solo le probabilità (da ottimizzare)
Pm=V1(:,1);                         %Estrae solo le probabilità (da ottimizzare)
Pr=V2(:,1);                         %Estrae solo le probabilità (da ottimizzare)

lambda=lambdas*Pr'*Pcomb*Pm;        %Prodotto dei termini vedi formula (1.23) del Backer

end
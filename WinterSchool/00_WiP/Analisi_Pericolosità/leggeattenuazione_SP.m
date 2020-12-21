function [log10PGAm,sigmapga]=leggeattenuazione_SP(R,M,dati)
%Valuta la legge di attenuazione di Sabetta Pugliese (1996)


Fr_eq=1.15;
Fn_eq=0.89;
Fss_eq=0.94;

%Coeff. della legge di attenuazione
a=-1.562;               %Termine costante
b=0.306;                %Termine magnitudo
c=-1;                   %Coeff. distanza
e1=0.169;               %Coeff. di sito
e2=0;                   %Coeff. di sito
h=5.8;
sigmapga=0.173;
S1=dati(7);
S2=S1;

switch dati(8)
    case 1
        cfaglia=Fn_eq;          %1 faglia normale (N)
    case 2
        cfaglia=Fr_eq;          %2 faglia inversa (S)
    case 3
        cfaglia=Fss_eq;         %3 trascorrrenti (SS)
    case 4
        cfaglia=1;              %4 indeterminato
end
    
    log10PGAm=cfaglia*(a+b*M+c*log10((R^2+h^2)^(1/2))+e1*S1+e2*S2);

end
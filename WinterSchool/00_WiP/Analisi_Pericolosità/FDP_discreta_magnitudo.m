function [pmag]=FDP_discreta_magnitudo(minf,msup,dati)
%Genera il vettore con la densità di probabilità discreta della magnitudo.
%Il valore di targa è il limite inferiore

mmin=   dati(3);    
mmax=   dati(4);
b=      dati(5);

if minf>=mmin && msup<=mmax
    num1=1-10^(-b*(minf-mmin));
    num2=1-10^(-b*(msup-mmin));
    den= 1-10^(-b*(mmax-mmin));
    Fsup=num2/den;
    Finf=num1/den;
    pmag=Fsup-Finf;
else
    pmag=0;

end


end
function [lnPGAm,sigmapga]=leggeattenuazione_cornel(R,M)
%Valuta la legge di attenuazione

sigmapga=0.57;
lnPGAm=-0.152+0.859*M-1.803*log(R+25);

end
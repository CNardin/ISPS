function [input]=xls2mtl(nomefile,cella1,cellan,nomeFoglio)

    if nargin>3
        nameFoglio = nomeFoglio;
    else
        nameFoglio = 'Foglio1';
    end

filename=strcat(nomefile,'.xlsm');
strcat(cella1,':',cellan);
input= xlsread(filename,nameFoglio,strcat(cella1,':',cellan),'basic');

end
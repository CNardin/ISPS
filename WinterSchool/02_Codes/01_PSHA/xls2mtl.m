function [input]=xls2mtl(nomefile,foglio,cella1,cellan)

filename=strcat(nomefile,'.xlsm');
sheet = foglio;
input= xlsread(filename,sheet,strcat(cella1,':',cellan));

end
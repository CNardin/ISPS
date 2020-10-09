function [input]=xls2mtl(nomefile,foglio,cella1,cellan)

filename=strcat(nomefile,'.xlsx');
sheet = foglio;
input= xlsread(filename,sheet,strcat(cella1,':',cellan));

end
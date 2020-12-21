function [input]=xls2mtl(nomefile,cella1,cellan)

filename=strcat(nomefile,'.xlsx');
strcat(cella1,':',cellan);
input= xlsread(filename,'coordinate',strcat(cella1,':',cellan));

end
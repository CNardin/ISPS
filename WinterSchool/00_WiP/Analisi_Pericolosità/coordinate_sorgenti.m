function [xsor,ysor]=coordinate_sorgenti(nomefile,cella1,cellan,nomeFoglio)

mat=xls2mtl(nomefile,cella1,cellan,nomeFoglio);
[n,m]=size(mat);

xsor=NaN(n,m/2);
ysor=NaN(n,m/2);
k=1;
for i=1:n
    for j=1:m/2
        xsor(i,k)=mat(i,2*j-1);
        ysor(i,k)=mat(i,2*j);
        k=k+1;
    end
    k=1;
end
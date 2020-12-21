clear all
close all
clc


c1=[4,19,29,40,47,60,73,82,91,108,120,130,143,153,161,173,185,194,204,212,232,245,253,264,273,280,287,307,314,333,347,354,366,376,383,391];
n=length(c1);
cn=[c1(2:n)-3,399];

for i=1:n
    cella1=c1(i);
    cellan=cn(i);

    rangein=strcat('A',num2str(cella1),':','B',num2str(cellan));
    rangeout=strcat('C',num2str(cella1),':','E',num2str(cellan));
    coordinate=xlsread('coordinate','coordinate',rangein);
    [x,y,utmzone] = deg2utm(coordinate(:,2),coordinate(:,1));
    out=table(x,y,utmzone);
    writetable(out,'coordinate.xlsx','Sheet','coordinate','Range',rangeout,'WriteVariableNames',false)
end

for i=1:n
    disp(['Check for ZS9 # ',num2str(i),' out of #', num2str(n)])
    cella1=c1(i);
    cellan=cn(i);
    
    rangein=strcat('F',num2str(cella1),':','G',num2str(cellan));
    rangeout=strcat('A',num2str(i));
    in=xlsread('coordinate','coordinate',rangein);
    
    m=length(in+1);
    out=zeros(1,2*m);
    
    if m==cellan-cella1+1    
        out(1)=900+i;
    else
        out(1)=-1;
    end
   
    if m~=2
    for j=1:m   
        out(2*j)=in(j,1)/1000;
        out(2*j+1)=in(j,2)/1000;
    end
    end
    
    output=table(out);
    writetable(output,'coordinate.xlsx','Sheet','orizzontale','Range',rangeout,'WriteVariableNames',false)
end


disp('Finito')

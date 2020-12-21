function stampa_montecarlo(xv,yv,flag3,xq,yq)    

figure(999) %plot evolution of target the proper area
plot(xv,yv) 
axis equal
plot(0,0,'O')
hold on
plot(xq(flag3),yq(flag3),'r.') 
plot(xq(~flag3),yq(~flag3),'c.')
hold off
    
end
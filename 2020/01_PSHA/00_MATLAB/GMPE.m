function [mean_im, sigma_im] = GMPE(magnitude,R_distance)
%Cornell attenuation law
mean_lnPGA = -0.152+0.859*magnitude -1.803*log(R_distance+25);
mean_PGA = exp(mean_lnPGA); 

mean_im = mean_PGA;
sigma_im = 0.57;
end


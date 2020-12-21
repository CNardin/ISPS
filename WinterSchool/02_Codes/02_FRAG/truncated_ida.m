function [mu_IDA, sigma_IDA] = truncated_ida(IM_i, IM_max, num_exceed)
%
% This function has been adapted form Baker 2013. It fits a lognormal CDF to observed probability of collapse 
% data by maximaze the loglikelihood based on equation (17) of Lecture 7 of
% your lecture notes
%
%
%
% INPUTS:
% IM_i          vector [1xn]             IM levels at which a ground motion
% wasexceeded
% IM_max        scalar                   Maximum IM level for which analysis was performed
% num_exceed 	scalar                   Number of analyses for which the IM_i
%                                           exceeds IM_max)
% 
% OUTPUTS:
% mu_IDA        scalar                   lognormal mean of the fragility function
% beta_IDA      scalar                   lognormal standard deviation of fragility function


% Initial guess for the fragility function parameters theta and beta. 
% These initial choices should not need revision in most cases, but they 
% could be altered if needed.
x0 = [0.5 0.4]; 

options = optimset('MaxFunEvals',1000, 'GradObj', 'off'); %maximum 1000 iterations, gradient of the function not provided
x = fminsearch(@mle, x0, options, IM_i, IM_max, num_exceed) ;
mu_IDA = x(1);
sigma_IDA = x(2);


% objective function to be optimized
function [ neg_loglik ] = mle( x0, IM_i, IM_max, num_exceed )

if (x0(1)<=0) || (x0(2)<=0)
    neg_loglik = 10^4; % penalize negative parameters 
else
    loglik_untr = sum(log(lognpdf(IM_i, x0(1), x0(2))));
    loglik_tr = num_exceed * (log(1-logncdf(IM_max, x0(1), x0(2))));
    neg_loglik = -(loglik_untr + loglik_tr); % Equation 18
end



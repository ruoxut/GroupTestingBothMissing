function [ res ] = f_r_mis_1D(x_1,d,rx,rd,t_r)
% Misspecified parametric missing data model for models (i) and (ii) with a 1-d covariate.
% Input:
% x_1: value of X_1; 
% d: value of d;
% rx: value of missing indicator R^X;
% rd: value of missing indicator R^D;
% t_r: parameters for f_R^X,R^D|X,D;
% Output:
% res: probability of f_R^X,R^D|X,D(rx,rd|x_1,d;t_r);
% Author: Ruoxu Tan; date: 1/Sep/2022; Matlab version: R2020a.

if (rx==0 && rd==0)
    res = exp(t_r(1))./(1+exp(t_r(1))+exp(t_r(2).*x_1.^2+t_r(3))+exp(t_r(4).*d+t_r(5)));
elseif (rx==1 && rd==0)
    res = exp(t_r(2).*x_1.^2+t_r(3))./(1+exp(t_r(1))+exp(t_r(2).*x_1.^2+t_r(3))+exp(t_r(4).*d+t_r(5)));
elseif (rx==0 && rd==1)
    res = exp(t_r(4).*d+t_r(5))./(1+exp(t_r(1))+exp(t_r(2).*x_1.^2+t_r(3))+exp(t_r(4).*d+t_r(5)));
else
    res = 1./(1+exp(t_r(1))+exp(t_r(2).*x_1.^2+t_r(3))+exp(t_r(4).*d+t_r(5)));
end

end


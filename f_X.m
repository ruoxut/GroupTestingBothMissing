function [ res ] = f_X(x_1,x_2,t_x) 
% Bivariate normal density function evaluated at (x_1,x_2).
% Input:
% x_1: value of X_1;
% x_2: value of X_2;  
% t_r: parameters for f_X;
% Output:
% res: probability of f_X(x_1,x_2;t_x);
% Author: Ruoxu Tan; date: 1/Sep/2022; Matlab version: R2020a.

mu_X_1 = t_x(1);
mu_X_2 = t_x(2);
var_X_1 = t_x(3);
cov_X_12 = t_x(4);
var_X_2 = t_x(5);

rho = cov_X_12/(sqrt(var_X_1*var_X_2));
res = exp( ((x_1-mu_X_1).^2./var_X_1 - 2.*rho.*(x_1-mu_X_1).*(x_2-mu_X_2)./(sqrt(var_X_1*var_X_2)) + (x_2-mu_X_2).^2./var_X_2)./(-2*(1-rho^2)) ) ./ (2*pi*sqrt(var_X_1*var_X_2*(1-rho^2)));

end


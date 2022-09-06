function [ res ] = f_r_v(x_1,x_2,d,rx,rd,t_r)
% Correctly specified parametric missing data model for model (v) with a 2-d covariate.
% Input:
% x_1: value of X_1;
% x_2: value of X_2;
% d: value of d;
% rx: value of missing indicator R^X;
% rd: value of missing indicator R^D;
% t_r: parameters for f_R^X,R^D|X,D;
% Output:
% res: probability of f_R^X,R^D|X,D(rx,rd|x_1,x_2,d;t_r);
% Author: Ruoxu Tan; date: 1/Sep/2022; Matlab version: R2020a.

if (rx==0 && rd==0)
    res = exp(t_r(1).*x_2+t_r(2).*d+t_r(3))./(1+exp(t_r(1).*x_2+t_r(2).*d+t_r(3))+exp(t_r(4).*x_1+t_r(5).*d+t_r(6))+exp(t_r(7).*x_1+t_r(8).*x_2+t_r(9)));
elseif (rx==1 && rd==0)
    res = exp(t_r(4).*x_1+t_r(5).*d+t_r(6))./(1+exp(t_r(1).*x_2+t_r(2).*d+t_r(3))+exp(t_r(4).*x_1+t_r(5).*d+t_r(6))+exp(t_r(7).*x_1+t_r(8).*x_2+t_r(9)));
elseif (rx==0 && rd==1)
    res = exp(t_r(7).*x_1+t_r(8).*x_2+t_r(9))./(1+exp(t_r(1).*x_2+t_r(2).*d+t_r(3))+exp(t_r(4).*x_1+t_r(5).*d+t_r(6))+exp(t_r(7).*x_1+t_r(8).*x_2+t_r(9)));
else
    res = 1./(1+exp(t_r(1).*x_2+t_r(2).*d+t_r(3))+exp(t_r(4).*x_1+t_r(5).*d+t_r(6))+exp(t_r(7).*x_1+t_r(8).*x_2+t_r(9)));
end


end


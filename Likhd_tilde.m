function [ nelogL ] = Likhd_tilde( X,Y_star_tilde,RX,RD,n,s_1,s_2,t_fix,t_est,opt )
% Negative log-likelihood function with the grouped outcome $\tilde{Y}^star$, i.e.
% the setting at (2.2) in the paper.
% Input:
% X: N*2 2-d covariate;
% Y_star_tilde: J*1 grouped outcome;
% RX: N*1 missing indicator for X_2;
% RD: N*1 missing indicator for D;
% n: J*1 group sizes;
% s_1: 1 - specificity;
% s_2: 1 - sensitivity;
% t_fix: 1*2 fixed parameters for f_X;
% t_est: 1*14 parameters to be estimated, (1:3) for f_X, (4:6) for p, (7:14) for f_RX,RD|X,D;
% opt: choices of the model for p, 1: model (iii) in the paper; 2: model (iv) in the paper.
% Output:
% nelogL: negative log-likelihood.
% Author: Ruoxu Tan; date: 1/Sep/2022; Matlab version: R2020a.

n1 = [0;n];

mu_X_1 = t_fix(1);
var_X_1 = t_fix(2);

mu_X_2 = t_est(1);
cov_X_12 = t_est(3)*sqrt(var_X_1);
var_X_2 = t_est(2)^2+t_est(3)^2;
 
% Model of p
if opt == 1
    p = @(x_1,x_2) 1./(8+exp(t_est(4).*x_1+t_est(5).*x_2+t_est(6)));
elseif opt == 2
    p = @(x_1,x_2) 1./(1+exp(t_est(4).*x_1.^2+t_est(5).*x_2+t_est(6))); 
end

% Parameters for f_X and f_r
t_x = [mu_X_1,mu_X_2,var_X_1,cov_X_12,var_X_2];
t_r = t_est(7:end);

res = 0; 

parfor j = 1:length(Y_star_tilde)
    n1_par = n1;
    RX_par = RX;
    X_1_par = X(:,1);
    X_2_par = X(:,2);
    
    f_z = 0;
    f_o = 0; 
    % We compute log(f_W) first and then take exp(log(f_W_D))
    for i = sum(n1_par(1:j))+1:sum(n1_par(1:j+1))
        if RX_par(i) == 1 
            f_z = f_z + log(f_X(X_1_par(i),X_2_par(i),t_x)) + log(f_r(X_1_par(i),X_2_par(i),0,1,1,t_r)) + log(1-p(X_1_par(i),X_2_par(i)));% log(f_W_D_0)
            f_o = f_o + log(f_X(X_1_par(i),X_2_par(i),t_x)) + log(f_r(X_1_par(i),X_2_par(i),1,1,1,t_r).*p(X_1_par(i),X_2_par(i)) + f_r(X_1_par(i),X_2_par(i),0,1,1,t_r).*(1-p(X_1_par(i),X_2_par(i))) );% log(f_W_D_0-f_W_D_1)
        else
            inti_z = @(x) f_X(X_1_par(i),x,t_x).*(1-p(X_1_par(i),x)).*f_r(X_1_par(i),x,0,0,1,t_r);
            f_z = f_z + log(integral(inti_z,-Inf,Inf));% log(\int f_W_D_0 dx)
            inti_o = @(x) f_X(X_1_par(i),x,t_x).*(f_r(X_1_par(i),x,1,0,1,t_r).*p(X_1_par(i),x) + f_r(X_1_par(i),x,0,0,1,t_r).*(1-p(X_1_par(i),x)) );
            f_o = f_o + log(integral(inti_o,-Inf,Inf));% log(\int (f_W_D_0-f_W_D_1) dx)
        end 
    end
        
    if Y_star_tilde(j) == 0 
        res = res + log( (1-s_1).*exp(f_z) + s_2.*(exp(f_o)-exp(f_z)) );% Take exp back
    else
        res = res + log( s_1.*exp(f_z) + (1-s_2).*(exp(f_o)-exp(f_z)) );
    end
      
end

N_obs = sum(RD);

for i = (N_obs+1):length(RD)
    % Compute the likelihood of the groups where R^D=0
    if RX(i) == 1
        res = res + log(f_X(X(i,1),X(i,2),t_x).* ( p(X(i,1),X(i,2)).*f_r(X(i,1),X(i,2),1,1,0,t_r)+(1-p(X(i,1),X(i,2))).*f_r(X(i,1),X(i,2),0,1,0,t_r) )); 
    else
        inti_one = @(x) f_X(X(i,1),x,t_x).* ( p(X(i,1),x).*f_r(X(i,1),x,1,0,0,t_r)+(1-p(X(i,1),x)).*f_r(X(i,1),x,0,0,0,t_r) );
        res = res + log(integral(inti_one,-Inf,Inf)); 
    end
end

nelogL = -res;

end


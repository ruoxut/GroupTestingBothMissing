function [ nelogL ] = Likhd_ign( X,Y_star,RX,RD,n,s_1,s_2,t_fix,t_est,opt )
% Negative log-likelihood function ignoring missing data mechanism
% with the grouped outcome $Y^star$, i.e. the setting at (2.3) in the paper.
% Input:
% X: N*2 2-d covariate;
% Y_star: J*1 grouped outcome;
% RX: N*1 missing indicator for X_2;
% RD: N*1 missing indicator for D;
% n: J*1 group sizes;
% s_1: 1 - specificity;
% s_2: 1 - sensitivity;
% t_fix: 1*2 fixed parameters for f_X;
% t_est: 1*6 parameters to be estimated, (1:3) for f_X, (4:6) for p;
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

% Parameters for f_X
t_x = [mu_X_1,mu_X_2,var_X_1,cov_X_12,var_X_2];

res = 0;
parfor j = 1:length(Y_star)
   if Y_star(j) == -1
       for i = sum(n1(1:j))+1:sum(n1(1:j+1))
           if RX(i) == 1
               res = res + log(f_X(X(i,1),X(i,2),t_x));
           else 
               res = res + log(normpdf(X(i,1),t_x(1),sqrt(t_x(3))));
           end
       end
   else
       f_0 = 0;
       f_1 = 0;
       % We compute log(f_W) first and then take exp(log(f_W_D))
       for i = sum(n1(1:j))+1:sum(n1(1:j+1))
           if RX(i) == 1
               if RD(i) == 1
                   f_0 = f_0 + log(1-p(X(i,1),X(i,2))) + log(f_X(X(i,1),X(i,2),t_x));% log(f_W_D_0)
               else
                   f_0 = f_0 + log(f_X(X(i,1),X(i,2),t_x));
               end
               f_1 = f_1 + log(f_X(X(i,1),X(i,2),t_x));% log(f_W_D_1)
           else
               if RD(i) == 1
                   inti_0 = @(x) (1-p(X(i,1),x)).*f_X(X(i,1),x,t_x);
                   f_0 = f_0 + log(integral(inti_0,-Inf,Inf));% log(\int f_W_D_0 dx)
               else
                   f_0 = f_0 + log(normpdf(X(i,1),t_x(1),sqrt(t_x(3))));
               end
               f_1 = f_1 + log(normpdf(X(i,1),t_x(1),sqrt(t_x(3))));% log(\int f_W_D_1 dx)
           end
       end
       
       if Y_star(j) == 0
           res = res + log( (1-s_1).*exp(f_0) + s_2.*(exp(f_1)-exp(f_0)) );% Take exp back
       else
           res = res + log( s_1.*exp(f_0) + (1-s_2).*(exp(f_1)-exp(f_0)) );
       end
   end
end
    
nelogL = -res;


end


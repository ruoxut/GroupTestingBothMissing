function [ nelogL ] = Likhd_v_ign( X,Y_star,RX,RD,n,s_1,s_2,t_est )
% Negative log-likelihood function ignoring missing data mechanism
% with the grouped outcome $Y^star$, i.e. the setting at (2.3) for model (v) in the paper.
% Input:
% X: N*2 2-d covariate;
% Y_star: J*1 grouped outcome;
% (RX,RD): N*2 missing indicators; 
% n: J*1 group sizes;
% s_1: 1 - specificity;
% s_2: 1 - sensitivity; 
% t_est: 1*8 parameters to be estimated, (1:5) for f_X, (6:8) for p; 
% Output:
% nelogL: negative log-likelihood.
% Author: Ruoxu Tan; date: 1/Sep/2022; Matlab version: R2020a.

n1 = [0;n];

mu_X_1 = t_est(1);
var_X_1 = t_est(2);

mu_X_2 = t_est(3);
cov_X_12 = t_est(5)*sqrt(var_X_1);
var_X_2 = t_est(4)^2+t_est(5)^2;
 
% Model of p 
p = @(x_1,x_2) 1./(1+exp(t_est(6).*x_1.^2+t_est(7).*x_2+t_est(8))); 

% Parameters for f_X
t_x = [mu_X_1,mu_X_2,var_X_1,cov_X_12,var_X_2];

res = 0;
parfor j = 1:length(Y_star)
   if Y_star(j) == -1
       for i = sum(n1(1:j))+1:sum(n1(1:j+1)) 
           res = res + log(f_X(X(i,1),X(i,2),t_x));
       end
   else
       f_0 = 0;
       f_1 = 0;
       % We compute log(f_W) first and then take exp(log(f_W_D))
       for i = sum(n1(1:j))+1:sum(n1(1:j+1))
           if RX(i) == 0 && RD(i) == 0
               inti_0 = @(x) (1-p(x,X(i,2))).*f_X(x,X(i,2),t_x);
               f_0 = f_0 + log(integral(inti_0,-Inf,Inf)); 
               f_1 = f_1 + log(normpdf(X(i,2),t_x(2),sqrt(t_x(5))));
           elseif RX(i) == 1 && RD(i) == 0
               inti_0 = @(x) (1-p(X(i,1),x)).*f_X(X(i,1),x,t_x);
               f_0 = f_0 + log(integral(inti_0,-Inf,Inf)); 
               f_1 = f_1 + log(normpdf(X(i,1),t_x(1),sqrt(t_x(3))));
           elseif RX(i) == 0 && RD(i) == 1
               f_0 = f_0 + log(f_X(X(i,1),X(i,2),t_x));
               f_1 = f_1 + log(f_X(X(i,1),X(i,2),t_x));
           else
               f_0 = f_0 + log(1-p(X(i,1),X(i,2))) + log(f_X(X(i,1),X(i,2),t_x));% log(f_W_D_0)
               f_1 = f_1 + log(f_X(X(i,1),X(i,2),t_x));% log(f_W_D_1)
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


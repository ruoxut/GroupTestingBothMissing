function [ nelogL ] = Likhd_v( X,Y_star,RX,RD,n,s_1,s_2,t_est )
% Negative log-likelihood function with the grouped outcome $Y^star$, i.e.
% the setting at (2.3) for model (v) in the paper.
% Input:
% X: N*2 2-d covariate;
% Y_star: J*1 grouped outcome;
% (RX,RD): N*2 missing indicators; 
% n: J*1 group sizes;
% s_1: 1 - specificity;
% s_2: 1 - sensitivity;
% t_est: 1*17 parameters to be estimated, (1:5) for f_X, (6:8) for p, (9:17) for f_R^X,R^D|X,D; 
% Output:
% nelogL: Negative log-likelihood.
% Author: Ruoxu Tan; date: 1/Sep/2022; Matlab version: R2020a.

n1 = [0;n];

mu_X_1 = t_est(1);
var_X_1 = t_est(2);

mu_X_2 = t_est(3);
cov_X_12 = t_est(5)*sqrt(var_X_1);
var_X_2 = t_est(4)^2+t_est(5)^2;
 
% Model of p 
p = @(x_1,x_2) 1./(1+exp(t_est(6).*x_1.^2+t_est(7).*x_2+t_est(8)));
 
% Parameters for f_X and f_r
t_x = [mu_X_1,mu_X_2,var_X_1,cov_X_12,var_X_2];
t_r = t_est(9:end);

res = 0;
parfor j = 1:length(Y_star)
    f_z = 0;
    f_o = 0; 
    % We compute log(f_W) first and then take exp(log(f_W_D))
    if Y_star(j) == -1
        for i = sum(n1(1:j))+1:sum(n1(1:j+1)) 
            res = res + log(f_X(X(i,1),X(i,2),t_x)) + ...
                  log( f_r_v(X(i,1),X(i,2),1,0,1,t_r).*p(X(i,1),X(i,2))...
                        + f_r_v(X(i,1),X(i,2),0,0,1,t_r).*(1-p(X(i,1),X(i,2))) );% log(f_W_D_-1)
        end
    else
        for i = sum(n1(1:j))+1:sum(n1(1:j+1))
            if RX(i) == 0 && RD(i) == 0
                inti_z = @(x) f_X(x,X(i,2),t_x).*(1-p(x,X(i,2))).*f_r_v(x,X(i,2),0,0,0,t_r);
                f_z = f_z + log(integral(inti_z,-Inf,Inf));
                inti_o = @(x) f_X(x,X(i,2),t_x).*(f_r_v(x,X(i,2),1,0,0,t_r).*p(x,X(i,2)) + f_r_v(x,X(i,2),0,0,0,t_r).*(1-p(x,X(i,2))) );
                f_o = f_o + log(integral(inti_o,-Inf,Inf));
            elseif RX(i) == 1 && RD(i) == 0
                inti_z = @(x) f_X(X(i,1),x,t_x).*(1-p(X(i,1),x)).*f_r_v(X(i,1),x,0,1,0,t_r);
                f_z = f_z + log(integral(inti_z,-Inf,Inf));
                inti_o = @(x) f_X(X(i,1),x,t_x).*(f_r_v(X(i,1),x,1,1,0,t_r).*p(X(i,1),x) + f_r_v(X(i,1),x,0,1,0,t_r).*(1-p(X(i,1),x)) );
                f_o = f_o + log(integral(inti_o,-Inf,Inf));
            else
                f_z = f_z + log(f_X(X(i,1),X(i,2),t_x)) + log(1-p(X(i,1),X(i,2))) + log(f_r_v(X(i,1),X(i,2),0,RX(i),RD(i),t_r));
                f_o = f_o + log(f_X(X(i,1),X(i,2),t_x)) + log(f_r_v(X(i,1),X(i,2),1,RX(i),RD(i),t_r).*p(X(i,1),X(i,2)) + f_r_v(X(i,1),X(i,2),0,RX(i),RD(i),t_r).*(1-p(X(i,1),X(i,2))) );
            end 
        end
        
        if Y_star(j) == 0 
            res = res + log( (1-s_1).*exp(f_z) + s_2.*(exp(f_o)-exp(f_z)) );% Take exp back
        else
            res = res + log( s_1.*exp(f_z) + (1-s_2).*(exp(f_o)-exp(f_z)) );
        end
    end
      
end

nelogL = -res;

end


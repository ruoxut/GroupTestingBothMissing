function [ nelogL ] = Likhd_1D( X,Y_star,RX,RD,n,s_1,s_2,t_est,opt )
% Negative log-likelihood function with the grouped outcome $Y^star$, i.e.
% the setting at (2.3) in the paper.
% Input:
% X: N*1 1-d covariate;
% Y_star: J*1 grouped outcome;
% RX: N*1 missing indicator for X;
% RD: N*1 missing indicator for D;
% n: J*1 group sizes;
% s_1: 1 - specificity;
% s_2: 1 - sensitivity;
% t_est: 1*9 (resp 10) parameters to be estimated, (1:2) for f_X, (3:4) (resp 5) for p, (5:9) (resp (6:10)) for f_RX,RD|X,D
% with model (i) (resp (ii)) in the paper;
% opt: choices of the model for p, 1: model (i) in the paper; 2: model (ii) in the paper.
% Output:
% nelogL: negative log-likelihood.
% Author: Ruoxu Tan; date: 1/Sep/2022; Matlab version: R2020a.

n1 = [0;n];

% Model of p
if opt == 1
    p = @(x_1) 1./(1+exp(t_est(3).*x_1+t_est(4)));
elseif opt == 2
    p = @(x_1) 1./(1+exp(t_est(3).*x_1.^2+t_est(4).*x_1+t_est(5)));
end

% Parameters for f_X and f_r_1D
t_x = [t_est(1),exp(t_est(2))];
if opt == 1
    t_r = t_est(5:end);
elseif opt == 2
    t_r = t_est(6:end);
end
    
res = 0;
parfor j = 1:length(Y_star)
    f_z = 0;
    f_o = 0; 
    % We compute log(f_W) first and then take exp(log(f_W_D))
    if Y_star(j) == -1
        for i = sum(n1(1:j))+1:sum(n1(1:j+1))
            if RX(i) == 1
                res = res + log(normpdf(X(i),t_x(1),t_x(2))) + ...
                      log(f_r_1D(X(i),1,1,0,t_r).*p(X(i)) + f_r_1D(X(i),0,1,0,t_r).*(1-p(X(i))) );% log(f_W_D_-1)
            else
                inti_m = @(x) normpdf(x,t_x(1),t_x(2)).* (f_r_1D(x,1,0,0,t_r).*p(x) + f_r_1D(x,0,0,0,t_r).*(1-p(x)));
                res = res + log(integral(inti_m,-Inf,Inf));% log(\int f_W_D_-1 dx)
            end
        end
    else
        for i = sum(n1(1:j))+1:sum(n1(1:j+1))
            if RX(i) == 1
                if RD(i) == 1
                    f_z = f_z + log(normpdf(X(i),t_x(1),t_x(2))) + log(1-p(X(i))) + log(f_r_1D(X(i),0,1,1,t_r));% log(f_W_D_0)
                else
                    f_z = f_z + log(normpdf(X(i),t_x(1),t_x(2))) + log(f_r_1D(X(i),1,1,0,t_r).*p(X(i)) + f_r_1D(X(i),0,1,0,t_r).*(1-p(X(i))) );% log(f_W_D_0)
                end
                f_o = f_o + log(normpdf(X(i),t_x(1),t_x(2))) + log(f_r_1D(X(i),1,1,RD(i),t_r).*p(X(i)) + f_r_1D(X(i),0,1,RD(i),t_r).*(1-p(X(i))) );% log(f_W_D_0-f_W_D_1)
            else
                if RD(i) == 1
                    inti_z1 = @(x) normpdf(x,t_x(1),t_x(2)).*(1-p(x)).*f_r_1D(x,0,0,1,t_r);
                    f_z = f_z + log(integral(inti_z1,-Inf,Inf));% log(\int f_W_D_0 dx)
                else
                    inti_z2 = @(x) normpdf(x,t_x(1),t_x(2)).*(f_r_1D(x,1,0,0,t_r).*p(x) + f_r_1D(x,0,0,0,t_r).*(1-p(x)) );
                    f_z = f_z + log(integral(inti_z2,-Inf,Inf));% log(\int f_W_D_0 dx)
                end
                inti_o = @(x) normpdf(x,t_x(1),t_x(2)).*(f_r_1D(x,1,0,RD(i),t_r).*p(x) + f_r_1D(x,0,0,RD(i),t_r).*(1-p(x)) );
                f_o = f_o + log(integral(inti_o,-Inf,Inf));% log(\int (f_W_D_0-f_W_D_1) dx)
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


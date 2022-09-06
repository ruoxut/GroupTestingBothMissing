function [ nelogL ] = Likhd_tilde_ign_1D( X,Y_star,RX,RD,n,s_1,s_2,t_est,opt )
% Negative log-likelihood function ignoring missing data mechanism 
% with the grouped outcome $\tilde{Y}^star$, i.e. the setting at (2.2) in the paper.
% Input:
% X: N*1 1-d covariate;
% Y_star: J*1 grouped outcome;
% RX: N*1 missing indicator for X;
% RD: N*1 missing indicator for D;
% n: J*1 group sizes;
% s_1: 1 - specificity;
% s_2: 1 - sensitivity;
% t_est: 1*4 (resp 5) parameters to be estimated, (1:2) for f_X, (3:4) (resp 5) for p with model (i) (resp (ii)) in the paper;
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

% Parameters for f_X
t_x = [t_est(1),exp(t_est(2))];

res = 0;
parfor j = 1:length(Y_star)  
   f_0 = 0;
   f_1 = 0;
   % We compute log(f_W) first and then take exp(log(f_W_D))
   for i = sum(n1(1:j))+1:sum(n1(1:j+1))
       if RX(i) == 1
           f_0 = f_0 + log(1-p(X(i))) + log(normpdf(X(i),t_x(1),t_x(2)));% log(f_W_D_0)
           f_1 = f_1 + log(normpdf(X(i),t_x(1),t_x(2)));% log(f_W_D_0-f_W_D_1)
       else
           inti_0 = @(x) (1-p(x)).*normpdf(x,t_x(1),t_x(2));
           f_0 = f_0 + log(integral(inti_0,-Inf,Inf));% log(\int f_W_D_0 dx)
       end
   end
   if Y_star(j) == 0
       res = res + log( (1-s_1).*exp(f_0) + s_2.*(exp(f_1)-exp(f_0)) );% Take exp back
   else
       res = res + log( s_1.*exp(f_0) + (1-s_2).*(exp(f_1)-exp(f_0)) );
   end
end

N_obs = sum(RD);

for i = (N_obs+1):length(RD)
    % Compute the likelihood of the groups where R^D=0
    if RX(i) == 1
        res = res + log(normpdf(X(i),t_x(1),t_x(2)));
    end
end

nelogL = -res;

end


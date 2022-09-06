% Simulations for MNAR X and D with the non-grouped outcome $Y$ in the paper
% A. Delaigle and R. Tan, Group testing regression analysis with covariates
% and specimens subject to missingness.
% This is to reproduce the simulation results of models (iii) and (iv) for the estimator
% $\hat{p}_{UMLE}$ in the paper.
% Author: Ruoxu Tan; date: 1/Sep/2022; Matlab version: R2020a.

n_rep = 200; % Number of repeated simulations
J_sample = [500 1000 2000]; % Number of groups

% Interested range for the 2-d covariate 
x_1_range = linspace(-1.5,1.5,200)';
x_2_range = x_1_range;
[x_1_eva,x_2_eva] = meshgrid(x_1_range,x_2_range); 

% To store the results of ISE and parameter estimates
ISE = zeros(n_rep,1);
Opt_res = cell(n_rep,1); 

s_1 = 0.01; % 1 - specificity
s_2 = 0.15; % 1 - sensitivity

opt = 1;% the options for the model
if opt == 1
    % Model (iii)
    p_eva = 1./(8+exp(8.*x_1_eva-8.*x_2_eva+8));
    p_fun = @(x_1,x_2) 1./(8+exp(8.*x_1-8.*x_2+8));
elseif opt == 2
    % Model (iv)
    p_eva = 1./(1+exp(x_1_eva.^2+0.5.*x_2_eva+1));
    p_fun = @(x_1,x_2) 1./(1+exp(x_1.^2+0.5.*x_2+1));
end

for s = 1:3

    J = J_sample(s); 

    % Uncomment if using grouping (A) ------------------
    N = 6*J;
    J = N;
    % --------------------------------------------------

    % Uncomment if using grouping (B) ------------------
    %N = 5*J;
    %J = N;
    % --------------------------------------------------

    n = ones(J,1); 

    for n_rep_ind = 1:n_rep
       rng(33*n_rep_ind)

       % Generate data
       X = mvnrnd([0,0],[0.75^2 0.5^2;0.5^2 0.75^2],N);
       D = binornd(ones(N,1),p_fun(X(:,1),X(:,2)));
       [ RX,RD ] = MisDM( X,D );

       X(RX==0,2) = NaN;
       D(RD==0) = NaN;
       n1 = [0;n];

       % Generate grouped status
       D_star = zeros(J,1);
       for j = 1:J
           RD_j = RD(sum(n1(1:j))+1:sum(n1(1:j+1)));
           if sum(RD_j==0) == length(RD_j)
               D_star(j) = -1;
           else
               D_star(j) = max(D(sum(n1(1:j))+1:sum(n1(1:j+1))),[],'omitnan');
           end
       end
       
       % Generate grouped test results
       Y_star = zeros(J,1);
       Y_star(D_star==-1) = -1;
       Y_star(D_star==0) = binornd(ones(sum(D_star==0),1),s_1);
       Y_star(D_star==1) = binornd(ones(sum(D_star==1),1),1-s_2);

       % MLE
       t_fix = [mean(X(:,1)) var(X(:,1))];% Fixed parameters
       
       % Initial choices for the parameters to be estimated
       t_ini3 = cov(X(RX==1,1),X(RX==1,2));
       t_ini3 = t_ini3(1,2)/sqrt(var(X(:,1)));
       t_ini2 = sqrt(var(X(RX==1,2))-t_ini3^2);
       t_ini = [mean(X(RX==1,2)) t_ini2 t_ini3 zeros(1,11)];
       t_ini_naive = [mean(X(RX==1,2)) t_ini2 t_ini3 zeros(1,3)];

       options = optimoptions(@fminunc,'Display','off');

       % $\hat{p}_{UMLE}$
       nelogL = @(t_est) Likhd( X,Y_star,RX,RD,n,s_1,s_2,t_fix,t_est,opt );
       t_est = fminunc(nelogL,t_ini,options);
       Opt_res{n_rep_ind} = t_est;% Estimated parameters
       if opt == 1
           p_hat = 1./(8+exp(t_est(4).*x_1_eva+t_est(5).*x_2_eva+t_est(6)));
       elseif opt == 2
           p_hat = 1./(1+exp(t_est(4).*x_1_eva.^2+t_est(5).*x_2_eva+t_est(6)));
       end
       sqer = (p_hat-p_eva).^2;
       ISE(n_rep_ind) = trapz(x_2_range,trapz(x_1_range,sqer));

    end
    
    % Store results
    fname = sprintf('211104_XY_(%d)_ungr_J%d',opt,J_sample(s));
    save(fname,'ISE','Opt_res');

    % Print quartile ISEs
    ISE_qt = quantile(ISE,[0.25,0.5,0.75]).*1000;
 
    fprintf([fname ': ']);
    fprintf('%0.2f (%0.2f) \n',ISE_qt(2),ISE_qt(3)-ISE_qt(1))

end
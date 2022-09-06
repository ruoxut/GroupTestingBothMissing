% Simulations for MNAR X and D with the grouped outcome $Y^star$, i.e. the setting at (2.3), in the paper
% A. Delaigle and R. Tan, Group testing regression analysis with covariates
% and specimens subject to missingness.
% This is to reproduce the simulation results of model (v) for the estimators
% $\hat{p}_{MLE,2}^{LOG}$, $\hat{p}_{MLE,2}^{MISP}$ and $\hat{p}_{MLE,2}^{IGN}$ in the paper.
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
ISE_missp = zeros(n_rep,1);
Opt_res_missp = cell(n_rep,1);
ISE_ign = zeros(n_rep,1);
Opt_res_ign = cell(n_rep,1);

s_1 = 0.01; % 1 - specificity
s_2 = 0.15; % 1 - sensitivity

% Model (v)
p_eva = 1./(1+exp(x_1_eva.^2+0.5.*x_2_eva+1));
p_fun = @(x_1,x_2) 1./(1+exp(x_1.^2+0.5.*x_2+1));


for s = 1:3

    J = J_sample(s);
    
    % Uncomment if using grouping (A) ------------------
    N = J*6;
    n = zeros(J,1);
    n(1:J/2) = 4;
    n(J/2+1:end) =8;
    % --------------------------------------------------
    
    % Uncomment if using grouping (B) ------------------
    %N = J*5;
    %n = ones(J,1).*5;
    % --------------------------------------------------
    
    for n_rep_ind = 1:n_rep
        rng(33*n_rep_ind)

        % Generate data
        X = mvnrnd([0,0],[0.75^2 0.5^2;0.5^2 0.75^2],N);
        D = binornd(ones(N,1),p_fun(X(:,1),X(:,2)));
        [ RX,RD ] = MisDM_v( X,D );

        X(RX==0&RD==0,1) = NaN;
        X(RX==1&RD==0,2) = NaN;
        D(RX==0&RD==1) = NaN;
        n1 = [0;n];

        % Generate grouped status
        D_star = zeros(J,1);
        for j = 1:J
            RD_j = RX==0&RD==1;
            RD_j = RD_j(sum(n1(1:j))+1:sum(n1(1:j+1)));
            if sum(RD_j==1) == length(RD_j)
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
        
        % Initial choices for the parameters to be estimated
        t_ini3 = cov(X(:,1),X(:,2),'omitrows');
        t_ini3 = t_ini3(1,2)/sqrt(var(X(:,1),'omitnan'));
        t_ini2 = sqrt(var(X(:,2),'omitnan')-t_ini3^2);
        t_ini = [mean(X(:,1),'omitnan') var(X(:,1),'omitnan') mean(X(:,2),'omitnan') t_ini2 t_ini3 zeros(1,12)];
        t_ini_naive = t_ini(1:8);
 
        options = optimoptions(@fminunc,'Display','off','MaxFunctionEvaluations',2500);
        
        % $\hat{p}_{MLE,2}^{LOG}$  
        nelogL = @(t_est) Likhd_v( X,Y_star,RX,RD,n,s_1,s_2,t_est);
        t_est = fminunc(nelogL,t_ini,options);
        Opt_res{n_rep_ind} = t_est;% Estimated parameters 
        p_hat = 1./(1+exp(t_est(6).*x_1_eva.^2+t_est(7).*x_2_eva+t_est(8)));
        sqer = (p_hat-p_eva).^2;
        ISE(n_rep_ind) = trapz(x_2_range,trapz(x_1_range,sqer));

        % $\hat{p}_{MLE,2}^{MISP}$
        nelogL = @(t_est) Likhd_v_missp( X,Y_star,RX,RD,n,s_1,s_2,t_est);
        t_est_missp = fminunc(nelogL,t_ini,options);
        Opt_res_missp{n_rep_ind} = t_est_missp;% Estimated parameters
        p_hat_missp = 1./(1+exp(t_est_missp(6).*x_1_eva.^2+t_est_missp(7).*x_2_eva+t_est_missp(8)));
        sqer_missp = (p_hat_missp-p_eva).^2;
        ISE_missp(n_rep_ind) = trapz(x_2_range,trapz(x_1_range,sqer_missp));

        % $\hat{p}_{MLE,2}^{IGN}$
        nelogL = @(t_est) Likhd_v_ign( X,Y_star,RX,RD,n,s_1,s_2,t_est);
        t_est_ign = fminunc(nelogL,t_ini_naive,options);
        Opt_res_ign{n_rep_ind} = t_est_ign;% Estimated parameters 
        p_hat_ign = 1./(1+exp(t_est_ign(6).*x_1_eva.^2+t_est_ign(7).*x_2_eva+t_est_ign(8)));
        sqer_ign = (p_hat_ign-p_eva).^2;
        ISE_ign(n_rep_ind) = trapz(x_2_range,trapz(x_1_range,sqer_ign));

    end
    
    % Store results
    fname = sprintf('220802_XY_(v)_A_J%d',J_sample(s));
    save(fname,'ISE','Opt_res','ISE_missp','Opt_res_missp','ISE_ign','Opt_res_ign');

    % Print quartile ISEs
    ISE_qt = quantile(ISE,[0.25,0.5,0.75]).*1000;
    ISE_qt_missp = quantile(ISE_missp,[0.25,0.5,0.75]).*1000;
    ISE_qt_ign = quantile(ISE_ign,[0.25,0.5,0.75]).*1000;

    fprintf([fname ': ']);
    fprintf('%0.2f (%0.2f) & %0.2f (%0.2f) & %0.2f (%0.2f)\n',...
        ISE_qt(2),ISE_qt(3)-ISE_qt(1),ISE_qt_missp(2),ISE_qt_missp(3)-ISE_qt_missp(1),...
        ISE_qt_ign(2),ISE_qt_ign(3)-ISE_qt_ign(1))

end
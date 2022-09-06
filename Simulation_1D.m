% Simulations for MNAR X and D with the grouped outcome $Y^star$, i.e. the setting at (2.3), in the paper
% A. Delaigle and R. Tan, Group testing regression analysis with covariates
% and specimens subject to missingness.
% This is to reproduce the simulation results of models (i) and (ii) for the estimators
% $\hat{p}_{MLE,2}^{LOG}$, $\hat{p}_{MLE,2}^{MISP}$ and $\hat{p}_{MLE,2}^{IGN}$ in the paper.
% Author: Ruoxu Tan; date: 1/Sep/2022; Matlab version: R2020a.

n_rep = 200; % Number of repeated simulations
J_sample = [500 1000 2000]; % Number of groups

% Interested range for the 1-d covariate 
x_1_range = linspace(-1.5,1.5,200)';  
x_1_eva = x_1_range;

% To store the results of ISE and parameter estimates
ISE = zeros(n_rep,1);
Opt_res = cell(n_rep,1);
ISE_missp = zeros(n_rep,1);
Opt_res_missp = cell(n_rep,1);
ISE_ign = zeros(n_rep,1);
Opt_res_ign = cell(n_rep,1);

s_1 = 0.01; % 1 - specificity
s_2 = 0.15; % 1 - sensitivity

opt = 2; % the options for the model
if opt == 1
    % Model (i)
    p_eva = 1./(1+exp(2.*x_1_range+3));
    p_fun = @(x_1) 1./(1+exp(2.*x_1+3));
elseif opt == 2
    % Model (ii)
    p_eva = 1./(1+exp(x_1_range.^2-3.*x_1_range+3));
    p_fun = @(x_1) 1./(1+exp(x_1.^2-3.*x_1+3));
end

for s = 1:3

    J = J_sample(s);

    % Uncomment if using grouping (A) ------------------
    %N = J*6;
    %n = zeros(J,1);
    %n(1:J/2) = 4;
    %n(J/2+1:end) =8;
    % --------------------------------------------------

    % Uncomment if using grouping (B) ------------------
    N = J*5;
    n = ones(J,1).*5;
    % --------------------------------------------------

    for n_rep_ind = 1:n_rep
        rng(33*n_rep_ind)

        % Generate data
        X = normrnd(0,0.75,[N,1]);
        D = binornd(ones(N,1),p_fun(X));
        [ RX,RD ] = MisDM_1D( X,D );

        X(RX==0) = NaN;
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
        if opt == 1
            t_ini = [mean(X(RX==1)) log(std(X(RX==1))) zeros(1,7)];
            t_ini_naive = [mean(X(RX==1)) log(std(X(RX==1))) zeros(1,2)];
        elseif opt ==2
            t_ini = [mean(X(RX==1)) log(std(X(RX==1))) zeros(1,8)];
            t_ini_naive = [mean(X(RX==1)) log(std(X(RX==1))) zeros(1,3)];
        end
        
        options = optimoptions(@fminunc,'Display','off');
        
        % $\hat{p}_{MLE,2}^{LOG}$ 
        nelogL = @(t_est) Likhd_1D( X,Y_star,RX,RD,n,s_1,s_2,t_est,opt );
        t_est = fminunc(nelogL,t_ini,options); 
        Opt_res{n_rep_ind} = t_est;% Estimated parameters
        if opt == 1
            p_hat = 1./(1+exp(t_est(3).*x_1_eva+t_est(4))); 
        elseif opt == 2
            p_hat = 1./(1+exp(t_est(3).*x_1_eva.^2+t_est(4).*x_1_eva+t_est(5))); 
        end
        
        sqer = (p_hat-p_eva).^2;
        ISE(n_rep_ind) = trapz(x_1_range,sqer);

        % $\hat{p}_{MLE,2}^{MISP}$ 
        nelogL = @(t_est) Likhd_missp_1D( X,Y_star,RX,RD,n,s_1,s_2,t_est,opt );
        t_est_missp = fminunc(nelogL,t_ini,options);
        Opt_res_missp{n_rep_ind} = t_est_missp;% Estimated parameters
        if opt == 1
            p_hat_missp = 1./(1+exp(t_est_missp(3).*x_1_eva+t_est_missp(4))); 
        elseif opt == 2
            p_hat_missp = 1./(1+exp(t_est_missp(3).*x_1_eva.^2+t_est_missp(4).*x_1_eva+t_est_missp(5))); 
        end
        
        sqer_missp = (p_hat_missp-p_eva).^2;
        ISE_missp(n_rep_ind) = trapz(x_1_range,sqer_missp);

        % $\hat{p}_{MLE,2}^{IGN}$
        nelogL = @(t_est) Likhd_ign_1D( X,Y_star,RX,RD,n,s_1,s_2,t_est,opt );
        t_est_ign = fminunc(nelogL,t_ini_naive,options);
        Opt_res_ign{n_rep_ind} = t_est_ign;% Estimated parameters
        if opt == 1
            p_hat_ign = 1./(1+exp(t_est_ign(3).*x_1_eva+t_est_ign(4))); 
        elseif opt == 2
            p_hat_ign = 1./(1+exp(t_est_ign(3).*x_1_eva.^2+t_est_ign(4).*x_1_eva+t_est_ign(5))); 
        end
        
        sqer_ign = (p_hat_ign-p_eva).^2;
        ISE_ign(n_rep_ind) = trapz(x_1_range,sqer_ign);

    end
    
    % Store results
    fname = sprintf('220728_XY_1D_(%d)_B_J%d',opt,J_sample(s));
    save(fname,'ISE','Opt_res','ISE_missp','Opt_res_missp','ISE_ign','Opt_res_ign');

    % Print quartile ISEs
    ISE_qt = quantile(ISE,[0.25,0.5,0.75]).*10000;
    ISE_qt_missp = quantile(ISE_missp,[0.25,0.5,0.75]).*10000;
    ISE_qt_ign = quantile(ISE_ign,[0.25,0.5,0.75]).*10000;

    fprintf([fname ': ']);
    fprintf('%0.2f (%0.2f) & %0.2f (%0.2f) & %0.2f (%0.2f)\n',...
        ISE_qt(2),ISE_qt(3)-ISE_qt(1),ISE_qt_missp(2),ISE_qt_missp(3)-ISE_qt_missp(1),...
        ISE_qt_ign(2),ISE_qt_ign(3)-ISE_qt_ign(1))

end
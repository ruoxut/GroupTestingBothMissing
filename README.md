# GroupTestingBothMissing
Simulation code for Delaigle, A. and Tan, R. (2023). Group testing regression analysis with covariates  and specimens subject to missingness. Statistics in Medicine. 42(6), 731-744.

To reproduce the simulation, run

1, Simulation_tilde_1D.m for models (i) and (ii) with \hat{p}_{MLE,1}^{LOG}, \hat{p}_{MLE,1}^{MISP} and \hat{p}_{MLE,1}^{IGN};

2, Simulation_1D.m for models (i) and (ii) with \hat{p}_{MLE,2}^{LOG}, \hat{p}_{MLE,2}^{MISP} and \hat{p}_{MLE,2}^{IGN};

3, Simulation_ungr_1D.m for models (i) and (ii) with \hat{p}_{UMLE};

4, Simulation_tilde.m for models (iii) and (iv) with \hat{p}_{MLE,1}^{LOG}, \hat{p}_{MLE,1}^{MISP} and \hat{p}_{MLE,1}^{IGN};

5, Simulation.m for models (iii) and (iv) with \hat{p}_{MLE,2}^{LOG}, \hat{p}_{MLE,2}^{MISP} and \hat{p}_{MLE,2}^{IGN};

6, Simulation_ungr.m for models (iii) and (iv) with \hat{p}_{UMLE};

7, Simulation_tilde_v.m for model (v) with \hat{p}_{MLE,1}^{LOG}, \hat{p}_{MLE,1}^{MISP} and \hat{p}_{MLE,1}^{IGN};

8, Simulation_v.m for model (v) with \hat{p}_{MLE,2}^{LOG}, \hat{p}_{MLE,2}^{MISP} and \hat{p}_{MLE,2}^{IGN};

9, Simulation_ungr_v.m for model (v) with \hat{p}_{UMLE}.

Remarks:

1, opt is the option for the model, e.g., opt = 1 and 2 in Simulation.m means models (iii) and (iv), respectively.

2, Un/comment the marked parts to select grouping (A) or (B).

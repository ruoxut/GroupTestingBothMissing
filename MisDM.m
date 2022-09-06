function [ RX,RD ] = MisDM( X,D )
% True missing data mechanism model for models (iii) and (iv) with the 2-d covariate.
% Input:
% X: N*2 1-d covariate;
% D: N*1 true individual status;
% Output:
% RX: N*1 missing indicator for X_2;
% RD: N*1 missing indicator for D.
% Author: Ruoxu Tan; date: 1/Sep/2022; Matlab version: R2020a.

n = length(X);
RX = zeros(n,1);
RD = RX;
for i = 1:n
    p00 = exp(0.5*X(i,1)-3)/(1+exp(0.5*X(i,1)-3)+exp(X(i,1)+X(i,2)-1.8)+exp(-0.5*X(i,1)+D(i)-1.5));% P(R^X=0,R^D=0|X,D)
    p10 = exp(X(i,1)+X(i,2)-1.8)/(1+exp(0.5*X(i,1)-3)+exp(X(i,1)+X(i,2)-1.8)+exp(-0.5*X(i,1)+D(i)-1.5));% P(R^X=1,R^D=0|X,D)
    p01 = exp(-0.5*X(i,1)+D(i)-1.5)/(1+exp(0.5*X(i,1)-3)+exp(X(i,1)+X(i,2)-1.8)+exp(-0.5*X(i,1)+D(i)-1.5));% P(R^X=0,R^D=1|X,D)
    p11 = 1/(1+exp(0.5*X(i,1)-3)+exp(X(i,1)+X(i,2)-1.8)+exp(-0.5*X(i,1)+D(i)-1.5));% P(R^X=1,R^D=1|X,D)
    
    p = [p00,p10,p01,p11];
    r = mnrnd(1,p);
    if r(2)==1
        RX(i)=1;
    elseif r(3)==1
        RD(i)=1;
    elseif r(4)==1
        RX(i)=1;
        RD(i)=1;
    end    
end

end


function [param,funH,funA]=tolparam

% param
param=[];

param.e=double(0.0005); %0.01; 
param.a=double(0.01); % entry of H2O2 into cell

param.b=0.05;  %degradation of H2O2 
%param.n1=3; % Hill function for H2O2 degradation
%param.k1=0.15*1/90; % half  mu for H2O2 degradation rate

param.g=0.005; %Antiox synthesis rate

param.mu0=1/90 ; % growth rate in absence of H2O2
param.d=0.0025; %rate of change of growth rate 

param.n=2; % hill function for growth inhibition
param.k=0.3; % half concentration for H2O2 growth inhibition


param.k2=0.1; % PKA inhibition 

% functions for analytical solving 

%funH= @(H,A,I) param.e + param.a * I - param.b * A; %+ param.a * (I-H)+  * param.k.^param.n / (param.k.^param.n + H.^param.n);

syms funH(H,A,I)
syms funA(H,A)

funH(H,A,I) = param.e + param.a * I - param.b * A;

funA(H,A) = param.g * H * param.k.^param.n / (param.k.^param.n + H.^param.n)  - param.d*A; 

% * param.mu0* param.k.^param.n / (param.k.^param.n + H.^param.n); 
%funM= @(H,A,M) param.d * (param.k2.^param.n2*param.mu0./ (param.k2.^param.n2 + H.^param.n2)  - M);

% init conditions
param.IC =double( [ param.d*param.e/(param.b*param.g) param.e/param.b]); 
%param.IC = [0 0.27 0 0.9 0.1];

param.IT= [];
param.I=[];
param.Tspan = []; % Solve from t=1 to t=5






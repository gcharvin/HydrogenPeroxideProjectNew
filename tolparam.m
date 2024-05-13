function [param,funH,funA,funM]=tolparam

% param
param=[];

param.e=0.01; 
param.a=0.1; % entry of H2O2 into cell

param.b=0.1;  %degradation of H2O2 
param.n1=3; % Hill function for H2O2 degradation
param.k1=0.15*1/90; % half  mu for H2O2 degradation rate

param.g=10; %Antiox synthesis rate

param.mu0=1/90 ; % growth rate in absence of H2O2
param.d=0.1; %rate of change of growth rate 

param.n2=2; % hill function for growth inhibition
param.k2=0.2; % half concentration for H2O2 growth inhibition

% functions for analytical solving 
funH= @(H,A,M,I) param.e + param.a * (I-H) - param.b * A * M.^param.n1 ./ (param.k1.^param.n1 + M.^param.n1); 
funA= @(H,A,M) param.g * M * H - M * A; 
funM= @(H,A,M) param.d * (param.k2.^param.n2*param.mu0./ (param.k2.^param.n2 + H.^param.n2)  - M);

% init conditions
param.IC = [0 0 param.mu0]; 
%param.IC = [0 0.27 0 0.9 0.1];

param.IT= [];
param.I=[];
param.Tspan = []; % Solve from t=1 to t=5






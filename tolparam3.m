function [param,funH,funA,funQ,funR]=tolparam3

% param
param=[];
param.tscreen=300; % masks the first x frames to allow for the equilibrium to settle. 

% H2O2 detoxification 
param.e=0.02; %0.01; 
param.a=0.1; % double(0.01); % entry of H2O2 into cell
param.b=1;  %degradation of H2O2 
param.g=1; %Antiox synthesis rate
param.kh=1; % half activation of Yap1 regulon 
param.nh=1; % non linearity in activiation of Yap1 regulon
param.a0=0; % offset in Yap1 regulon activation 

% growth signalling
param.mu0=1/90 ; % growth rate in absence of H2O2
param.k=0.6; %0.5 half concentration for H2O2 growth inhibition

% protection and toxicity parameters 
param.k2=0.0004; % reverse reaction from Q to R
param.kq=0.03; % formaton of portected state
param.e2=0.1; % formation of R
param.kr=1; % transiton to cell death

% ODE model definition

syms funH(H,A,I)
syms funA(H,A)
syms funQ(H,Q,R)
syms funR(H,Q,R)

funH(H,A,I) = param.e + param.a * (I-H) - param.b * A * H * max(0,(1 - param.b *H / param.k));

funA(H,A) = param.mu0 *max(0, (1 - param.b *H / param.k)) * (param.a0+param.g * (H.^param.nh/(H.^param.nh+(param.kh)^param.nh) - A)); 

funQ(H,Q,R)= (param.kq * H * R  / param.k  - (param.mu0   + param.k2) * Q ) * max(0,1 - param.b *H / param.k);

funR(H,Q,R)=  (param.e2 + param.k2 * Q - ( param.mu0 + param.kq * H / param.k ) * R ) * max(0,1 - param.b *H / param.k) ;

% initial conditions

%param.IC =double( [ param.d*param.e/(param.b*param.g) param.e/param.b]); 
param.IC =double( [ 0 0 0 0]); 

param.IT= [];
param.I=[];
param.Tspan = [];

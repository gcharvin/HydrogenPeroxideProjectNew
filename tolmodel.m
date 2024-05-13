function [T Y]=tolmodel(param,funH,funA,funM)
% experimental version to test features of the model
% vartiables
% I, H , A, X, Xagg, X*
% growth rate mu is proportional to X: mu= k* X



options=odeset('NonNegative',1:3,'RelTol',1e-5,'AbsTol',[1e-5 1e-5 1e-5],'MaxOrder',2);

[T Y] = ode15s(@(t,y) myode(t,y,param,funH,funA,funM),param.Tspan,param.IC,options);

% ode15s much faster than ode45

function dydt = myode(t,y,param, funH, funA, funM)

%int,in
in = interp1(param.IT,param.I,t);

% Interpolate the data set (ft,f) at time t

dydt(1) = funH(y(1),y(2),y(3),in);
dydt(2) = funA(y(1),y(2),y(3));  
dydt(3) = funM(y(1),y(2),y(3));

dydt=dydt';



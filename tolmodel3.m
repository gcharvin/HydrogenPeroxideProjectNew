function [T Y]=tolmodel3(param,funH,funA,funQ,funR)
% experimental version to test features of the model
% vartiables
% I, H , A, Q, R, R*

options=odeset('NonNegative',1:4,'RelTol',1e-3,'AbsTol',[1e-3 1e-3 1e-3 1e-3],'MaxOrder',2);

%[T Y] = ode45(@(t,y) myode(t,y,param,funH,funA),param.Tspan,param.IC,options);
[T, Y] = ode15s(@(t,y) myode(t,y,param,funH,funA,funQ,funR), double(param.Tspan), double(param.IC), options);
% ode15s much faster than ode45

function dydt = myode(t,y,param, funH, funA,funQ,funR)

%int,i
in = interp1(param.IT,param.I,t);

% Interpolate the data set (ft,f) at time t

dydt(1) = funH(y(1),y(2),in);

if y(1)>param.k % growth rate is zero, cell is frozen, only chemical reaction can occur
 dydt(2)=0;
 dydt(3)=0;
 dydt(4)=0;
else
dydt(2) = funA(y(1),y(2));  
dydt(3) = funQ(y(1),y(3),y(4));
dydt(4) = funR(y(1),y(3),y(4));  
end


dydt=double(dydt');




function nullclines=tolanalytics(display,param)



if nargin<2
[param,funH,funA,funM]=tolparam;
param.I=0.8;
end

if nargin
    display=0;
end

% nullclines analysis 

M1=param.mu0*(0.01:0.01:1.5);
H1= (param.e + param.a * param.I) ./ (param.a + (param.b * param.g) .* M1.^param.n1 ./ (param.k1.^param.n1 + M1.^param.n1));

H2= 0.005:0.01:1.5*param.I;
M2= param.k2.^param.n2*param.mu0./ (param.k2.^param.n2 + H2.^param.n2);

nullclines.M1=M1;
nullclines.H1=H1;
nullclines.M2=M2;
nullclines.H2=H2;


if display==1
figure, plot(H1,M1,'Color','r','lineWidth',2);
hold on; 
plot(H2,M2,'Color','b','lineWidth',2);

% figure, loglog(H1,M1,'Color','r');
% hold on; 
% loglog(H2,M2,'Color','b');


x=0.02:0.005:1.5*param.I;
y=param.mu0*(0.005:0.005:1.5);
% 
 [x,y] = meshgrid(x,y);
 

% 
 v = param.d*(param.k2.^param.n2*param.mu0./ (param.k2.^param.n2 + x.^param.n2)-y);
 
 u = param.e + param.a*(param.I-x) - param.b * param.g .* x .* y.^param.n1 ./(y.^param.n1 + param.k1.^param.n1) ;
% 

 quiver(x,y,u,v,'MaxHeadSize',0.);
 
 set(gca,'XScale','log');
 set(gca,'YScale','log');



xlabel('[H2O2] (mM)');
ylabel('Growth rate (min^{-1})');
% funH= @(H,A,M,I) param.a * I - param.b * A * M.^param.n1 ./ (param.k1.^param.n1 + M.^param.n1); 
% funA= @(H,A,M) param.g * H - M * A; 
% funM= @(H,A,M) param.d * (param.mu0./ (param.k2.^param.n2 + H.^param.n2)  - M);
end


function nullclines=tolanalytics2(display,param)



if nargin<2
[param,funH,funA]=tolparam2;
param.I=0.5;
end

if nargin==0
    display=0;
end

% nullclines analysis 

H= 0:0.01:10; %1.5*param.I;

A1= ones(1,length(H)).*(param.e + param.a * (param.I) ) ./param.b; % ( param.b *  param.k.^param.n ./ (param.k.^param.n + H.^param.n));

A2= (H* param.g * param.k.^param.n ./ (param.k.^param.n + H.^param.n)) ./param.d;

nullclines.H=H;
nullclines.A1=A1;
nullclines.A2=A2;
%nullclines.H2=H2;


if display==1
figure, plot(H,A1,'Color','r','lineWidth',2);
hold on; 
plot(H,A2,'Color','b','lineWidth',2);

% figure, loglog(H1,M1,'Color','r');
% hold on; 
% loglog(H2,M2,'Color','b');


x=0.0:0.1:1.2*param.I;
y=param.I*(param.a/param.b)*(0:0.1:5);

% 
 [x,y] = meshgrid(x,y);
 
% 

 v = (param.g* x - param.mu0*y) * param.k.^param.n ./(x.^param.n + param.k.^param.n);
 
 u = param.e + param.a*(param.I-x) - param.b *  y .* param.k.^param.n ./(x.^param.n + param.k.^param.n) ;
% 

 quiver(x,y,u,v,'MaxHeadSize',0.1);
 
 %set(gca,'XScale','log');
 %set(gca,'YScale','log');

 xlim([0 max(x(:))])
 ylim([0 max(y(:))])


xlabel('[H2O2] (mM)');
ylabel('Antioxydants (A.U.)');
% funH= @(H,A,M,I) param.a * I - param.b * A * M.^param.n1 ./ (param.k1.^param.n1 + M.^param.n1); 
% funA= @(H,A,M) param.g * H - M * A; 
% funM= @(H,A,M) param.d * (param.mu0./ (param.k2.^param.n2 + H.^param.n2)  - M);
end


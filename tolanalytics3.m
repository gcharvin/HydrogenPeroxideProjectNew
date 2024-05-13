function nullclines=tolanalytics3(display,param)
% computes and drwas nucclines 


if nargin<2
[param,funH,funA]=tolparam3;
param.I=0.2;
end

if nargin==0
    display=0;
end

% nullclines analysis 

H= 0.01:0.005:2; %1.5*param.I;
A1= (param.e + param.a * (param.I-H))./(param.b * H .* (1 - param.b *H / param.k)); % ( param.b *  param.k.^param.n ./ (param.k.^param.n + H.^param.n));
A2= param.a0+(H.^param.nh./(H.^param.nh+param.kh.^param.nh))* param.g; %* param.k.^param.n ./ (param.k.^param.n + H.^param.n)) ./param.d;

nullclines.H=H;
nullclines.A1=A1;
nullclines.A2=A2;
%nullclines.H2=H2;


if display==1
figure, plot(H,A1,'Color','r','lineWidth',2);
hold on; 
plot(H,A2,'Color','b','lineWidth',2);

xlabel('[H2O2] (mM)');
ylabel('Antioxydants (A.U.)');

set(gca,'XScale','log','YScale','log')
xlim([0.01 10]);
ylim([0.01 10]);
end


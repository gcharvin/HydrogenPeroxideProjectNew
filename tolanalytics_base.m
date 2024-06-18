function nullclines=tolanalytics_base(display,param)
% computes and drwas nucclines 

if nargin<2
[param,funH,funA]=tolparam3;
param.I=0.1;
end

if nargin==0
    display=0;
end

% nullclines analysis 

H= 0.01:0.005:2; %1.5*param.I;
A1= (param.e + param.a * (param.I-H)) .* (1 + param.b *H / param.k) ./(param.b * H); % ( param.b *  param.k.^param.n ./ (param.k.^param.n + H.^param.n));
%A2= param.a0+(H.^param.nh./(H.^param.nh+param.kh.^param.nh))* param.g; %* param.k.^param.n ./ (param.k.^param.n + H.^param.n)) ./param.d;
A2=  H* param.g; %* param.k.^param.n ./ (param.k.^param.n + H.^param.n)) ./param.d;

nullclines.H=H;
nullclines.A1=A1;
nullclines.A2=A2;
%nullclines.H2=H2;


if display==1
figure('Color','w'), plot(H,A1,'Color','r','lineWidth',2);
hold on; 
plot(H,A2,'Color','b','lineWidth',2);

xx=param.k/param.b;

line([xx xx],[0.01 10],'LineWidth',1,'LineStyle','--','Color','k');
xlabel('[H2O2] (mM)');
ylabel('Antioxodants (A.U.)');

set(gca,'XScale','log','YScale','log','FontSize',16)
xlim([0.02 2]);
ylim([0.02 5]);
end


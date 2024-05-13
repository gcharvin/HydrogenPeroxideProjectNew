function [outT,outY,options]=tolstep(display)


if nargin==0
    display=0;
end

display=1;

% define the time-sequence for input (external H2O2)
maxeT=400;
np=maxeT;

[param,funH,funA,funM]=tolparam; % initialize paremeters for WT;
param.nullclines=[];

outT={};
outY={};
options=[];

irange=0.1*[1  3 5 10]; %step
duration=240;
tstart=200;
tend=maxeT;

options.tstart=tstart;
options.tend=tend;
options.param=param;
options.irange=irange;

%[I, IT]=makePulse(irange,duration,tstart,maxeT);
[I IT]=makeStep(irange,tstart,maxeT);

%IT= linspace(0,maxeT,np);

% I=zeros(length(irange),np);
% rampe
% IT=I;
%
% for i=1:length(irange);
%     IT(i,:)= linspace(0,maxeT,maxeT);
%     I(i,1:np)=irange(i)*(1:np)/np;
% end


%param.I=I;
%param.IT=IT;
param.Tspan=[0 maxeT];

%irange=[0.0005 0.001 0.002]; % ramp experiment

%mutants={'WT','zwf1','trx1/2','zwf1 trx1/2' 'pde2'};

mutants={'WT'};

options.mutants=mutants;

%mutants={'WT'};
%param2=param;

  
for j=1:numel(mutants)
    disp(['Simulating strain: ' mutants{j}]);
     
      
    for i=1:numel(irange)
        
        disp(['Simulating condition: ' num2str(i) '/' num2str(length(irange))]);
        
        param2(i,j)=param;
        param2(i,j).IT=IT(i,:);
        param2(i,j).I=I(i,:);
        
        parama=param;
        parama.I= irange(i);
        
        param2(i,j).nullclines= tolanalytics(0,parama);
        %         if strcmp(mutants{j},'zwf1')
        %             param2(i,j).b=0.005;
        %             %param2(i,j).beta_agg=0.0005;
        %         end
        %         if strcmp(mutants{j},'trx1/2')
        %             param2(i,j).alpha_agg=0.005;
        %             param2(i,j).b=0.005;
        %             param2(i,j).beta_agg=0.0005;
        %         end
        %         if strcmp(mutants{j},'zwf1 trx1/2')
        %             param2(i,j).alpha_agg=0.005;
        %             param2(i,j).b=0.005;
        %             param2(i,j).beta_agg=0.0005;
        %         end
        %
        %         if strcmp(mutants{j},'pde2')
        %             param2(i,j).K=10;
        %            % param2(i,j).b=0.005;
        %            % param2(i,j).beta_agg=0.0005;
        %         end
        
        [T,Y]=tolmodel(param2(i,j),funH,funA,funM);
        
        outT{i,j}=T;
        outY{i,j}=Y;
        
        
    end
    
    
   
    if display
         hplot=[];
        %col=shallowColormap(length(irange));
        
        
        
        
        hf=figure('Color','w','Position',[100+200*(j-1) 100 600 1200]);
        col=colormap(lines(length(irange)));
        
        p=panel;
        p.pack('v',{1/4 1/4 1/4 1/4});
        
       
        for i=1:length(irange)
            displaySim(hf,p,outT{i,j},outY{i,j},param2(i,j),col(i,:));
        end
        
        p.de.margin=7;
        p.fontsize=20;
        
        p(1).select();
        title(mutants{j});
        xlim([0 maxeT]);
        %ylim([0 2]);
        set(gca,'FontSize',14,'XTickLabel',{});
        ylabel('I');
        
        p(2).select();
        xlim([0 maxeT]);
        set(gca,'FontSize',14,'XTickLabel',{});
        ylabel('H');
        xlim([0 maxeT]);
        ylim([0 1]);
        
        p(3).select();
        set(gca,'FontSize',14);
        ylabel('A');
        xlim([0 maxeT]);
        ylim([0 2]);
        
        set(gca,'FontSize',14,'XTickLabel',{});
        
        p(4).select();
        set(gca,'FontSize',14);
        ylabel('mu');
        xlim([0 maxeT]);
        %ylim([0 1]);
        set(gca,'FontSize',14,'XTickLabel',{});
  
        xlim([0 maxeT]);
        
        xlabel('Time (min.)');

        p.fontsize=20;
        p.marginbottom=20;
        p.marginleft=25;
        %ylim([0 2]);
        
        xlabel('Time (min.)');
        
        
        % now plot trajectories in the H / mu hyperplan
   
        ht=figure('Color','w','Position',[100+200*(j-1) 300 600 600]);
        
        nullclines= param2(end,j).nullclines;
        plot(nullclines.H2,nullclines.M2,'Color','k','LineWidth',2,'LineStyle','--'); hold on
        
        for i=length(irange):-1:1
            nullclines= param2(i,j).nullclines;
            
            plot(nullclines.H1,nullclines.M1,'Color',col(i,:),'LineStyle','--'); hold on
            
            T=outT{i,j};
            pix=find(T>np/2-10,1,'first');
            Y=outY{i,j};
            
            hplot(i)=plot(Y(pix:end,1),Y(pix:end,3),'LineWidth',3,'Color',col(i,:)); hold on;
        end
        
        nullclines= param2(end,j).nullclines;
        xlim([0.9*min(nullclines.H2) 1.1*max(nullclines.H2)]);
        ylim([0.9*min(nullclines.M2) 1.1*max(nullclines.M2)]);
        
        
        set(gca,'XScale','log','FontSize',20);
        %   set(gca,'YScale','log','FontSize',20);
        
        xlabel('[H2O2] (mM)');
        ylabel('Growth rate (min^{-1})');
        
        % now plot the survical fraction = f(time) & f(dose) and gowth fold change
        % during stress , and fitness
        
        
        %         hf=figure('Color','w','Position',[600 100 600 900]);
        %         p=panel;
        %         p.pack('v',{1/3 1/3 1/3});
        %
        %         p.de.margin=7;
        %         p.fontsize=20;
        %
        %         p(1).select();
        %         xlim([0 maxeT]);
        %
        %         %pix=find(param(1).I>0,1,'first');
        %         %t0=param2(1,1).IT(pix);
        %
        %         %pix2=find(param(1).I>0,1,'last');
        %         %t1=param2(1,1).IT(pix2);
        %
        %         t0=tstart;
        %         t1=tend;
        %
        %         for i=1:length(irange)
        %
        %             color=col(i,:);
        %             time=outT{i,j};
        %             pt0=find(time>=t0,1,'first');
        %             pt1=find(time<t1,1,'last');
        %
        %             mortality=  outY{i,j}(:,5);
        %             growth=     outY{i,j}(:,3);
        %
        %             fgrowth=exp(cumtrapz(time(pt0:pt1),param2(i,j).mu.*growth(pt0:pt1)));
        %             fgrowth=fgrowth./fgrowth(1);
        %             %aa=time(pt0:pt1)
        %
        %             %size(pt0:pt1), size(fgrowth)
        %             p(1).select();
        %             semilogy(time(pt0:pt1),fgrowth,'Color',color,'LineWidth',2); hold on;
        %
        %             fmortality=exp(cumtrapz(time(pt0:pt1),-param2(i,j).mu.*mortality(pt0:pt1)));
        %             fmortality=fmortality./fmortality(1);
        %             %aa=time(pt0:pt1)
        %
        %             %size(pt0:pt1), size(fgrowth)
        %             p(2).select();
        %             semilogy(time(pt0:pt1),fmortality,'Color',color,'LineWidth',2); hold on;
        %
        %             fitness=exp(cumtrapz(time(pt0:pt1),param2(i,j).mu.*(growth(pt0:pt1)-mortality(pt0:pt1))));
        %             fitness=fitness./fitness(1);
        %             %aa=time(pt0:pt1)
        %
        %             %size(pt0:pt1), size(fgrowth)
        %             p(3).select();
        %             semilogy(time(pt0:pt1),fitness,'Color',color,'LineWidth',2); hold on;
        %
        %
        %         end
        %
        %         p(1).select();
        %         title(mutants{j});
        %         ylabel('Growth fold change');
        %         xlim([0 maxeT]);
        %         ylim([1 50]);
        %
        %         p(2).select();
        %         ylabel('Survival');
        %         xlim([0 maxeT]);
        %         ylim([0.05 1]);
        %
        %         p(3).select();
        %         ylabel('Fitness');
        %         xlim([0 maxeT]);
        %         ylim([0.1 100]);
        %
        %         p.fontsize=20;
        %         p.marginbottom=20;
        %         p.marginleft=25;
        %
        %         xlabel('Time (min.)');
        
        
        % plot tolerance and resistance
        
        figure('Color','w','Position',[100 100 1600 400]);
        
        
        %param
        x=0.9*min(nullclines.H2):0.005:1.1*max(nullclines.H2);
        y=param.mu0*(0.005:0.005:1.5);
        %
        [x,y] = meshgrid(x,y);
        
        ax=subplot(1,3,1);
        %z=1;
        c=y;
        
        
        h=pcolor(x,y,c);
        
               cl=caxis;
        set(h, 'EdgeColor', 'none');
        hc=colorbar;
        
        for i=length(irange):-1:1
            copyobj(hplot(i),ax);
        end
        
        ylabel(hc,'Growth rate (A.U.)');
        
        xlim([0.9*min(nullclines.H2) 1.1*max(nullclines.H2)]);
        ylim([0.9*min(nullclines.M2) 1.1*max(nullclines.M2)]);
        
        set(gca,'XScale','log','FontSize',20);
        
        xlabel('[H2O2] (mM)');
        ylabel('Growth rate (min^{-1})');
        title('Resistance');
        
        ax=subplot(1,3,2);
        
        %a=3*cl(2);
        %tol=(y.*x).^3./(1*0.01)^3;
        %tol=a*y.^4.*x./(0.01.^4*1); 
        
        tol=0.1*x.*(y./0.01).^6;
        %tol=800*(x.*y).^2;
        c=tol;
        
     
        
        h=pcolor(x,y,c);
        set(h, 'EdgeColor', 'none');
        hc=colorbar;
        
        for i=length(irange):-1:1
            copyobj(hplot(i),ax);
        end
        
           caxis(5*cl);
           
        ylabel(hc,'Death rate (A.U.)');
        
        xlim([0.9*min(nullclines.H2) 1.1*max(nullclines.H2)]);
        ylim([0.9*min(nullclines.M2) 1.1*max(nullclines.M2)]);
        
        set(gca,'XScale','log','FontSize',20);
        
        xlabel('[H2O2] (mM)');
        ylabel('Growth rate (min^{-1})');
        title('Tolerance');
        
        ax=subplot(1,3,3);
        
        c=y-tol;
        
        h=pcolor(x,y,c);
        
        for i=length(irange):-1:1
            copyobj(hplot(i),ax);
        end
        
           caxis([-2*cl(2) cl(2)]);
           
        set(h, 'EdgeColor', 'none');
        hc=colorbar;
        ylabel(hc,'Fitness (A.U.)');
        
        xlim([0.9*min(nullclines.H2) 1.1*max(nullclines.H2)]);
        ylim([0.9*min(nullclines.M2) 1.1*max(nullclines.M2)]);
        
        set(gca,'XScale','log','FontSize',20);
        
        xlabel('[H2O2] (mM)');
        ylabel('Growth rate (min^{-1})');
        title('Fitness');
        
    end
    
    
end


function displaySim(hf,p,T,Y,param,col)

figure(hf);

p(1).select();
plot(param.IT,param.I,'LineWidth',2,'Color',col); hold on;

p(2).select();
plot(T,Y(:,1),'LineWidth',2,'Color',col); hold on;

p(3).select();
plot(T,Y(:,2),'LineWidth',2,'Color',col); hold on;

p(4).select();
plot(T,Y(:,3),'LineWidth',2,'Color',col); hold on;


function [I, IT]=makePulse(irange,duration,tstart,maxeT)

%irange=[0.0005 0.001 0.002]; % ramp experiment

%col=colormap(lines(length(irange)));
I=zeros(length(irange),maxeT);
IT=I;

for i=1:length(irange)
    IT(i,:)= linspace(0,maxeT,maxeT);
    I(i,tstart:tstart+duration-1)=irange(i); %step
end

function [I, IT]=makeStep(irange,tstart,maxeT)

%irange=[0.0005 0.001 0.002]; % ramp experiment

%col=colormap(lines(length(irange)));
I=zeros(length(irange),maxeT);
IT=I;

for i=1:length(irange)
    IT(i,:)= linspace(0,maxeT,maxeT);
    I(i,tstart:end)=irange(i); %step
end


function y=repress(x,K,n)

y=K^n./(K^n + x.^n);

function cmap=shallowColormap(n)

% n is the number of colors

switch n
    case 1
        cmap=[1 0 0];
        
    case 2
        cmap=[1 0 0 ; 0 1 0];
        
    case 3
        cmap=[1 0 0 ; 0 1 0; 0 0 1];
        
    otherwise
        
        if n<8
            cmap=jet(10);
        else
            
            cmap=jet(n);
        end
        
        
        tmp=cmap;
        tmp(2:2:end,:)=cmap(1:size(cmap,1)/2,:);
        tmp(1:2:end-1,:)=cmap(size(cmap,1)/2+1:end,:);
        %cmap=tmp(end:-1:1,:);
        cmap=tmp;
        
        if n<8
            cmap=cmap(1:n,:);
        end
        
        if n==9
            cmap(9,:)=[1 0.5 0.5];
        end
end

cmap=[0 0 0; cmap];





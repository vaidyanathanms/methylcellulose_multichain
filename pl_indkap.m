function flag_success = pl_indkap(nchains,trialnum,initdist,config,tplot,kappain)
%% Color scheme

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];
brown = [0.6 0.2 0];violet = [0.5,0,0.5];gray = [0.75 0.75 0.75];
p4clr = {orange,'c',green,gold};
p8clr = {green,'r',gray,'b',orange,'m','c',gold};
lsty = {'-','--',':'};

%% Start plotting data
% Plot kappa
hz = figure; %kappa
hold on
box on
set(gca,'FontSize',20)
xlabel('Time ($\tau$)','FontSize',20,'Interpreter','Latex')
ylabel('$\kappa^2_{i}$','FontSize',20,'Interpreter','Latex')
for i = 1:nchains
    if nchains == 4
        plot(tplot(:,1),kappain(:,i),'Color',p4clr{i},'LineWidth',2,'LineStyle',lsty{1});
    elseif nchains == 8
        plot(tplot(:,1),kappain(:,i),'Color',p8clr{i},'LineWidth',2,'LineStyle',lsty{1});
    end
    %legendinfo{i} = ['Chain Number: ' num2str(i)];
end
xlim([0 1.2*max(tplot)]);
xline=1:max(1.2*max(tplot));
yline=0.25*ones(length(xline),1);
plot(xline,yline,'Color','k','LineWidth',2,'LineStyle',lsty{2});
%legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','best')
%legend boxoff
saveas(hz,sprintf('../allfigures/kappa_%d_%d_%s_%s',nchains,trialnum,initdist,config),'png')
clear legendinfo;


% Plot kappa as a function of frame number
hz = figure; %kappa
hold on
box on
set(gca,'FontSize',20)
xlabel('Frame Number','FontSize',20,'Interpreter','Latex')
ylabel('$\kappa^2_{i}$','FontSize',20,'Interpreter','Latex')
for i = 1:nchains
    if nchains == 4
        plot(1:length(tplot(:,1)),kappain(:,i),'Color',p4clr{i},'LineWidth',2,'LineStyle',lsty{1});
    elseif nchains == 8
        plot(1:length(tplot(:,1)),kappain(:,i),'Color',p8clr{i},'LineWidth',2,'LineStyle',lsty{1});
    end
    %legendinfo{i} = ['Chain Number: ' num2str(i)];
end
xlim([0 1.2*max(length(tplot(:,1)))]);
xline=1:max(1.2*max(tplot));
yline=0.25*ones(length(xline),1);
plot(xline,yline,'Color','k','LineWidth',2,'LineStyle',lsty{2});
%legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','best')
%legend boxoff
saveas(hz,sprintf('../allfigures/kappavsframenum_%d_%d_%s_%s',nchains,trialnum,initdist,config),'png')
clear legendinfo;


% Plot Delta lambdaz in normal plot
hz = figure; %lambdaz(t)-lambdaz(t-1)
hold on
box on
set(gca,'FontSize',20)
xlabel('$\Delta t$ ($\tau$)','FontSize',20,'Interpreter','Latex')
ylabel('$\Delta \kappa$','FontSize',20,'Interpreter','Latex')
tlen = length(tplot(:,1));
tnew = zeros(tlen-1,1);
eigznew = zeros(tlen-1,1);
for j = 1:tlen-1
    tnew(j,1) = tplot(j+1,1)-tplot(1,1);
    for i = 1:nchains
        eigznew(j,i) = abs(kappain(j+1,i)-kappain(j,i));
    end
end

for i = 1:nchains
    if nchains == 4
        plot(tnew(:,1),eigznew(:,i),'Color',p4clr{i},'LineWidth',2,'LineStyle',lsty{1});
    elseif nchains == 8
        plot(tnew(:,1),eigznew(:,i),'Color',p8clr{i},'LineWidth',2,'LineStyle',lsty{1});
    end
    legendinfo{i} = ['Chain Number: ' num2str(i)];
end
xlim([0 1.2*max(tplot)]);
legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','best')
legend boxoff
saveas(hz,sprintf('../allfigures/kappadiff_%d_%d_%s_%s',nchains,trialnum,initdist,config),'png')
clear legendinfo;

% Plot Delta kappa in subplot

hz = figure; %delta_kappa = kappa_z-avg(kappa_z) in subplot
hold on
box on
set(gca,'FontSize',16)
h = zeros(4,1);
for i = 1:4
    h(i) = subplot(4,1,i);
    if nchains == 4
        plot(tplot(:,1)*10^-4,abs(kappain(:,i)-kappain(1,i)),'Color',p4clr{i},'LineWidth',2,'LineStyle',lsty{1});
        ylim([0 1.1*max(kappain(:,1)-kappain(1,i))]);
    elseif nchains == 8
        plot(tplot(:,1)*10^-4,kappain(:,i)-kappain(1,i),'Color',p8clr{i},'LineWidth',2,'LineStyle',lsty{1});
        ylim([0 1.1*max(kappain(:,1)-kappain(1,i))]);
    end
end

set(h(1),'xticklabel',[]);
set(h(2),'xticklabel',[]);
set(h(3),'xticklabel',[]);

pos=get(h,'position');
bottom=pos{4}(2);
top=pos{1}(2)+pos{1}(4);
plotspace=top-bottom;

pos{4}(4)=plotspace/4;
pos{3}(4)=plotspace/4;
pos{2}(4)=plotspace/4;
pos{1}(4)=plotspace/4;

pos{1}(2)=bottom + 3*plotspace/4;
pos{2}(2)=bottom + 2*plotspace/4;
pos{3}(2)=bottom + plotspace/4;


set(h(1),'position',pos{1},'FontSize',16);
set(h(2),'position',pos{2},'FontSize',16);
set(h(3),'position',pos{3},'FontSize',16);
set(h(4),'position',pos{4},'FontSize',16);
box on


xlabel('$\Delta t \times 10^{-4}$ ($\tau$)','FontSize',20,'Interpreter','Latex')
ylabel('$\Delta \kappa$','FontSize',20,'Interpreter','Latex')
box(h(1),'on');
box(h(2),'on');
box(h(3),'on');
saveas(hz,sprintf('../allfigures/kappadiffsubplot_%d_%d_%s_%s',...
    nchains,trialnum,initdist,config),'png')
flag_success = 1;
end
function flag_success = pl_indeig(nchains,trialnum,initdist,config,tplot,eigzinplot)
%% Color scheme

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];
brown = [0.6 0.2 0];violet = [0.5,0,0.5]; gray = [0.75 0.75 0.75];
p4clr = {orange,'c',green,gold};
p8clr = {green,'r',gray,'b',orange,'m','c',gold};
lsty = {'-','--',':'};

%% Start plotting data
% Plot lambda z
hz = figure; %lambdaz
hold on
box on
set(gca,'FontSize',20)
xlabel('Time ($\tau$)','FontSize',20,'Interpreter','Latex')
ylabel('$\lambda_{z}$','FontSize',20,'Interpreter','Latex')
for i = 1:nchains
    if nchains == 4
        plot(tplot(:,1),eigzinplot(:,i),'Color',p4clr{i},'LineWidth',2,'LineStyle',lsty{1});
    elseif nchains == 8
        plot(tplot(:,1),eigzinplot(:,i),'Color',p8clr{i},'LineWidth',2,'LineStyle',lsty{1});
    end
    legendinfo{i} = ['Chain Number: ' num2str(i)];
end
xlim([0 1.2*max(tplot)]);
legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','best')
legend boxoff
saveas(hz,sprintf('./Figure_Results/lambdaZ_%d_%d_%s_%s',nchains,trialnum,initdist,config),'png')
clear legendinfo;

% Plot Delta lambdaz in normal plot
hz = figure; %lambdaz(t)-lambdaz(t-1)
hold on
box on
set(gca,'FontSize',20)
xlabel('$\Delta t$ ($\tau$)','FontSize',20,'Interpreter','Latex')
ylabel('$\Delta \lambda_{z}$','FontSize',20,'Interpreter','Latex')
tlen = length(tplot(:,1));
tnew = zeros(tlen-1,1);
eigznew = zeros(tlen-1,1);
for j = 1:tlen-1
    tnew(j,1) = tplot(j+1,1)-tplot(1,1);
    for i = 1:nchains
        eigznew(j,i) = abs(eigzinplot(j+1,i)-eigzinplot(j,i));
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
saveas(hz,sprintf('./Figure_Results/lambdaZdiff_%d_%d_%s_%s',nchains,trialnum,initdist,config),'png')
clear legendinfo;

% Plot Delta lambdaz in subplot

hz = figure; %delta_lambda in subplot
hold on
box on
set(gca,'FontSize',16)
h = zeros(4,1);
for i = 1:4
    h(i) = subplot(4,1,i);
    if nchains == 4
        plot(tplot(:,1)*10^-4,abs(eigzinplot(:,i)-mean(eigzinplot(:,i))),'Color',p4clr{i},'LineWidth',2,'LineStyle',lsty{1});
        ylim([0 1.1*max(abs(eigzinplot(:,i)-mean(eigzinplot(:,i))))]);
    elseif nchains == 8
        plot(tplot(:,1)*10^-4,abs(eigzinplot(:,i)-mean(eigzinplot(:,i))),'Color',p8clr{i},'LineWidth',2,'LineStyle',lsty{1});
        ylim([0 1.1*max(abs(eigzinplot(:,i)-mean(eigzinplot(:,i))))]);
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
ylabel('$\Delta \lambda_{z}$','FontSize',20,'Interpreter','Latex')
box(h(1),'on');
box(h(2),'on');
box(h(3),'on');
saveas(hz,sprintf('./Figure_Results/lambdaZdiffsubplot_%d_%d_%s_%s',...
    nchains,trialnum,initdist,config),'png')
flag_success = 1;
end
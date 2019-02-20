%% To analyze all the systems in one curve

clc;
clear;
close all;
format long;

%% Color scheme

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];
brown = [0.6 0.2 0];violet = [0.5,0,0.5];gray = [0.75 0.75 0.75];
p4clr = {orange,'c',green,gold};
p8clr = {green,'r',gray,'b',orange,'m','c',gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Flags

indkappa = 1;
globfac  = 0;

%% Input Details
% Make sure the trial_arr order matches the order of the arch_arr and
% initdist_arr

nchains = 4;

if nchains == 4
    trial_arr    = [4] %[61,63,6,4,5,63,63]; %[5] %[6]  %
    arch_arr     = {'block'} %{'block','block','alternate','block','block','alternate','alternate'}; %{'block'} %{'alternate'} %
    initdist_arr = {'0.05'} %{'0.05','0.03','0.07','0.05','0.05','0.05','0.07'}; %{'0.05'} %{'0.07'}%
    labelchar = 'B4';
    delxax = 0.001;
elseif nchains == 8
    trial_arr    = [6,7,6,7];%,6,63,63];
    arch_arr     = {'block','block','alternate','alternate'};
    initdist_arr = {'0.07','0.03','0.07','0.07'};%,'0.07','0.05','0.07'};
    delxax = 0.001;
end


%% Basic Checks

len_alldata = length(trial_arr);
len_archarr = length(arch_arr);
len_distarr = length(initdist_arr);

if len_alldata ~= len_archarr || len_alldata ~= len_distarr
    errorMessage=sprintf('Error: Different length of arrays');
    uiwait(warndlg(errorMessage));
    return;
else
    fprintf('Analyzing Data... \n');
end



%% Read Data and Plot

colorbase = hsv(20);
lenbase   = length(colorbase);
colorarr  = colorbase;
lencolor  = length(colorarr);

if indkappa == 1
    
    yax_val  = 0.2;
    
    hz = figure; %kappa
    hold on
    box on
    set(gca,'FontSize',20)
    ylabel('Time ($\tau$)','FontSize',20,'Interpreter','Latex')
    %xlabel('Chain ID','FontSize',20,'Interpreter','Latex')
    
    
    for tval = 1:len_alldata
        
        trialnum = trial_arr(tval);
        initdist = initdist_arr{tval};
        config   = arch_arr{tval};
        
        dirname = sprintf('../trial_alldata/n%d_t%d_%s_%s',nchains,trialnum,...
            config,initdist);
        
        if ~isdir(dirname)
            errorMessage=sprintf('Error: Folder\t%s does not exist',dirname);
            uiwait(warndlg(errorMessage));
            return;
        else
            fprintf('Analyzing\t%s\n',dirname);
        end
        
        for chid = 1:nchains
            
            fylename = strcat(dirname,sprintf('/indshapefac_chID_%d.dat',chid));
            alldata  = importdata(fylename);
            tplot    = alldata.data(:,1);
            kappasq  = alldata.data(:,2);
            yax_arr  = yax_val*ones(length(tplot),1);
            
            for pllen = 1:length(tplot)
                clrrow = 1 + floor(kappasq(pllen,1)*lencolor);
                clrstr = [colorarr(clrrow,1) colorarr(clrrow,2) colorarr(clrrow,3)];
                if nchains == 4 && len_alldata > 1
                    if(tplot(pllen,1) < 10000)
                        plot(yax_arr(pllen,1),tplot(pllen,1),'Color',clrstr,...
                            'Marker','square','MarkerFaceColor',clrstr,'MarkerSize',16);
                    end
                elseif nchains == 4 && len_alldata == 1
                    if tplot(pllen,1) > 10500 && tplot(pllen,1) < 12000
                        plot(yax_arr(pllen,1),tplot(pllen,1),'Color',clrstr,...
                            'Marker','square','MarkerFaceColor',clrstr,'MarkerSize',75);
                    end
                elseif nchains == 8
                    if(tplot(pllen,1) < 10000)
                        plot(yax_arr(pllen,1),tplot(pllen,1),'Color',clrstr,...
                            'Marker','square','MarkerFaceColor',clrstr,'MarkerSize',12);
                    end
                end
            end
            
            if len_alldata ~= 1
                yax_val = yax_val + 2*delxax;
            else
                yax_val = yax_val + 0.0001*delxax;
            end
            clear yax_arr
            
        end
        
        if len_alldata ~= 1
            yax_val = yax_val + 4*delxax;
        else
            yax_val = yax_val + 0.4*delxax;
        end
        clear yax_arr
        

        
    end
    
    
    
    set(gca,'XTickLabel',[]);
    set(gca,'XTick',[]);
    
    if nchains == 4 && len_alldata ~=1
        xticks(0.202:12*delxax:0.202+16*delxax*(len_alldata-1));
        xticklabels({'B1','B2','B3','B4','A1','A2','A3'});
        xlim([0.201-6*delxax 0.202+12*delxax*(len_alldata)])
        ylim([-300 10                                                                                                                                                                                                                                   300]);
    elseif nchains == 4 && len_alldata ==1
        set(gca,'XTickLabel',[]);
        set(gca,'XTick',[]);
        %         xlim([0.200 0.20005])
    elseif nchains == 8
        xticks(0.206:22*delxax:0.202+32*delxax*(len_alldata-1));
        xticklabels({'B1','B2','A1','A2'});
        xlim([0.201-3*delxax 0.202+20*delxax*(len_alldata)])
    end
    
    if len_alldata == 1
        colormap(colorarr);
        c=colorbar;
        c.Limits = [0 1];
        c.FontSize = 20;
        c.Location = 'southoutside';
        c.Label.String = '$\kappa_{i}^{2}$';
        c.Label.FontSize = 20;
        c.Label.Interpreter = 'Latex';
    end
    
    
    if len_alldata ~= 1
        saveas(hz,sprintf('../allfigures/nocolorbar_allchaindata_%d',nchains),'png')
    else
        saveas(hz,sprintf('../allfigures/colorbar_allchaindata_%s',labelchar),'png')
    end
    %     saveas(hz,sprintf('../allfigures/allchaindata_%d',nchains),'pdf')
    
end

%% Global shape anisotropy factor for each system

% if globfac == 1
%
%     yax_val  = 0.2;
%
%     hz = figure; %kappa
%     hold on
%     box on
%     set(gca,'FontSize',20)
%     xlabel('Time ($\tau$)','FontSize',20,'Interpreter','Latex')
%     ylabel('System ID','FontSize',20,'Interpreter','Latex')
%
%
%     for tval = 1:len_alldata
%
%         trialnum = trial_arr(tval);
%         initdist = initdist_arr{tval};
%         config   = arch_arr{tval};
%
%         dirname = sprintf('../trial_alldata/n%d_t%d_%s_%s',nchains,trialnum,...
%             config,initdist);
%
%         if ~isdir(dirname)
%             errorMessage=sprintf('Error: Folder\t%s does not exist',dirname);
%             uiwait(warndlg(errorMessage));
%             return;
%         else
%             fprintf('Analyzing\t%s\n',dirname);
%         end
%
%         fylename = strcat(dirname,'/globshapefac.dat');
%         alldata  = importdata(fylename);
%         tplot    = alldata.data(:,1);
%         kappasq  = alldata.data(:,2);
%         yax_arr  = yax_val*ones(length(tplot),1);
%
%         for pllen = 1:length(tplot)
%
%             clrrow = 1 + floor(kappasq(pllen,1)*lencolor);
%             clrstr = [colorarr(clrrow,1) colorarr(clrrow,2) colorarr(clrrow,3)];
%             plot(tplot(pllen,1),yax_arr(pllen,1),'Color',clrstr,...
%                 'Marker','square','MarkerFaceColor',clrstr,'MarkerSize',5);
%
%         end
%
%         yax_val = yax_val + 0.5;
%         clear yax_arr
%
%     end
%
%     set(gca,'YTickLabel',[]);
%     set(gca,'YTick',[]);
%     colormap(colorarr);
%     c=colorbar;
%     c.Limits = [0 1];
%     c.FontSize = 20;
%     c.Location = 'northoutside';
%     c.Label.String = '$\kappa_{i}^{2}$';
%     c.Label.FontSize = 16;
%     c.Label.Interpreter = 'Latex';
%     saveas(hz,sprintf('../allfigures/globchaindata_%d',nchains),'png')
% %     saveas(hz,sprintf('../allfigures/globchaindata_%d',nchains),'pdf')
%
% end

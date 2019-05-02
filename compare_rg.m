%% To compare the rg and temperature between Langevin and NVT systems

clear;
clc;
close all;
format long;

%% Color Styles

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];
pclr = {'r','b',green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Input Data

n = 200;
nchains = 1;
dt_nvt  = 0.0005;
dt_lang = 0.001;

%% Load files

nvt_fylename = sprintf('./n_%d/nvt/jobsingle.sh_nvt',n);
fTemp_nvt = fopen(nvt_fylename,'r');
if exist(nvt_fylename, 'file') ~= 2
    fprintf('File does not exist \t%s\n',nvt_fylename);
end


lang_fylename = sprintf('./n_%d/lang/jobsingle.sh_lang',n);
fTemp_lang = fopen(lang_fylename,'r');
if exist(lang_fylename, 'file') ~= 2
    fprintf('File does not exist \t%s\n',lang_fylename);
end

frg_nvt  = fopen(sprintf('./n_%d/nvt/rg_nvt.dat',n),'r');
if frg_nvt <= 0
    fprintf('Rg file does not exist for NVT \n');
end

frg_lang = fopen(sprintf('./n_%d/lang/rg_lang.dat',n),'r');
if frg_lang <= 0
    fprintf('Rg file does not exist for Langevin \n');
end

%% Find Rg values and load into arrays

data_nvt_rg  = textscan(frg_nvt,'%f%f','Headerlines',3);
data_lang_rg = textscan(frg_lang,'%f%f','Headerlines',3);

rg_nvt  = cell2mat(data_nvt_rg);
rg_lang = cell2mat(data_lang_rg);

i = 1;k = 1;
rg_time_nvt = zeros(floor(length(rg_nvt(:,1))/(nchains+1)),2);
while i <= length(rg_nvt(:,1))
    tval = rg_nvt(i,1); nch_val = rg_nvt(i,2);
    rgavg = 0;
    for j = 1:nch_val
        rgavg = rgavg + rg_nvt(i+j,2)^2;
    end
    rg_time_nvt(k,1) = tval;
    rg_time_nvt(k,2) = mean(sqrt(rgavg));
    k = k + 1;
    i = i + nch_val + 1;
end

i = 1;k = 1;
rg_time_lang = zeros(floor(length(rg_lang(:,1))/(nchains+1)),2);
while i <= length(rg_lang(:,1))
    tval = rg_lang(i,1); nch_val = rg_lang(i,2);
    rgavg = 0;
    for j = 1:nch_val
        rgavg = rgavg + rg_lang(i+j,2)^2;
    end
    rg_time_lang(k,1) = tval;
    rg_time_lang(k,2) = mean(sqrt(rgavg));
    k = k + 1;
    i = i + nch_val + 1;
end

%% Find temperature column and load into arrays

fprintf('Analyzing Temperature of NVT data \n');
cntvals = 0;
tempfind = -1;

temp_time_nvt = zeros(floor(0.9*length(rg_time_nvt(:,1))),2);
tline = strtrim(fgetl(fTemp_nvt));

equilsteps = min(rg_time_nvt(:,1));
maxsteps   = max(rg_time_nvt(:,1));

while ~feof(fTemp_nvt)
    
    strarr = strsplit(tline);
    
    %Find box dimensions if flag is not 1 - for new cases and when restart occurs
    if tempfind ~= 1 && strcmp(strarr{1},'Step')
        
        str2arr = strsplit(tline);
        lenstr = length(str2arr);
        
        tempcol = -1;
        
        for sval = 1:lenstr
            if strcmp(str2arr{sval},'Temp')
                tempcol = sval;
            end
        end
        
        tempfind = 1;
        
    elseif tempfind == 1 && all(ismember(strarr{1}, '0123456789+-.eEdD'))
        
        stepval = str2double(strarr{1});
        
        if stepval < equilsteps %% Skip until equilibration
            
            tline = strtrim(fgetl(fTemp_nvt));
            continue;
            
        end
        
        if (tempcol == -1)
            fprintf('Temp Column Flags = %d at %s\n',tempcol,nvt_fylename)
            break;
        end
        
        % Find values until end of that given run cycle or until end of
        % file
        
        while ~feof(fTemp_nvt)
            
            if isnan(str2double(strarr{1}))
                
                tline = strtrim(fgetl(fTemp_nvt));
                strarr = strsplit(tline);
                continue;
                
            elseif all(ismember(strarr{2}, '0123456789+-.eEdD'))
                
                
                if(str2double(strarr{1})>=equilsteps && str2double(strarr{1}) <= maxsteps)
                    
                    cntvals = cntvals + 1;
                    temp_time_nvt(cntvals,1) = str2double(strarr{1});
                    temp_time_nvt(cntvals,2) = str2double(strarr{tempcol});
                    
                end
                
                
                tline = strtrim(fgetl(fTemp_nvt));
                strarr = strsplit(tline);
                
               
            else
                
                flagold = -1;
                tempfind = -1; %Refind temp for new cycles
                break;
                
            end
            
        end
        
    elseif ~feof(fTemp_nvt)
        
        tline = strtrim(fgetl(fTemp_nvt));
        
    end
    
    
end


fprintf('Analyzing Temperature of Langevin data \n')
cntvals = 0;
tempfind = -1;

temp_time_lang = zeros(floor(0.9*length(rg_time_lang(:,1))),2);
tline = strtrim(fgetl(fTemp_lang));

equilsteps = min(rg_time_lang(:,1));
maxsteps   = max(rg_time_lang(:,1));

while ~feof(fTemp_lang)
    
    strarr = strsplit(tline);
    
    %Find box dimensions if flag is not 1 - for new cases and when restart occurs
    if tempfind ~= 1 && strcmp(strarr{1},'Step')
        
        str2arr = strsplit(tline);
        lenstr = length(str2arr);
        
        tempcol = -1;
        
        for sval = 1:lenstr
            if strcmp(str2arr{sval},'Temp')
                tempcol = sval;
            end
        end
        
        tempfind = 1;
        
    elseif tempfind == 1 && all(ismember(strarr{1}, '0123456789+-.eEdD'))
        
        stepval = str2double(strarr{1});
        
        if stepval < equilsteps %% Skip until equilibration
            
            tline = strtrim(fgetl(fTemp_lang));
            continue;
            
        end
        
        if (tempcol == -1)
            fprintf('Temp Column Flags = %d at %s\n',tempcol,lang_fylename)
            break;
        end
        
        % Find values until end of that given run cycle or until end of
        % file
        
        while ~feof(fTemp_lang)
            
            if isnan(str2double(strarr{1}))
                
                tline = strtrim(fgetl(fTemp_lang));
                strarr = strsplit(tline);
                continue;
                
            elseif all(ismember(strarr{2}, '0123456789+-.eEdD'))
                
                
                if(str2double(strarr{1})>equilsteps && str2double(strarr{1}) <= maxsteps)
                    
                    cntvals = cntvals + 1;
                    temp_time_lang(cntvals,1) = str2double(strarr{1});
                    temp_time_lang(cntvals,2) = str2double(strarr{tempcol});
                    
                end
                
                
                tline = strtrim(fgetl(fTemp_lang));
                strarr = strsplit(tline);
                
               
            else
                
                flagold = -1;
                tempfind = -1; %Refind temp for new cycles
                break;
                
            end
            
        end
        
    elseif ~feof(fTemp_lang)
        
        tline = strtrim(fgetl(fTemp_lang));
        
    end
    
    
end




%% Plot Temperature-time and Rg-time

h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('Time ($\tau$)','FontSize',20,'Interpreter','Latex')
ylabel('$R_g$ ($\sigma$)','FontSize',20,'Interpreter','Latex')
plot(dt_nvt*(rg_time_nvt(:,1)-rg_time_nvt(1,1)),rg_time_nvt(:,2),'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3})
plot(dt_lang*(rg_time_lang(:,1)-rg_time_lang(1,1)),rg_time_lang(:,2),'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3})
legendinfo{1} = 'NVT Thermostat';
legendinfo{2} = 'Langevin Thermostat';
legend(legendinfo,'FontSize',20,'Location','best')
legend boxoff


h2 = figure;
hold on
box on
set(gca,'FontSize',16)
set(gca,'yscale','log')
xlabel('Time ($\tau$)','FontSize',20,'Interpreter','Latex')
ylabel('Temperature ($K$)','FontSize',20,'Interpreter','Latex')
plot(dt_nvt*(temp_time_nvt(:,1)-temp_time_nvt(1,1)),temp_time_nvt(:,2),'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3})
plot(dt_lang*(temp_time_lang(:,1)-temp_time_lang(1,1)),temp_time_lang(:,2),'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3})
legendinfo{1} = 'NVT Thermostat';
legendinfo{2} = 'Langevin Thermostat';
legend(legendinfo,'FontSize',20,'Location','best')
legend boxoff


%% Prepare for phase-plot

minl = min(length(rg_time_lang(:,1)),length(temp_time_lang(:,1)));
temp_rg_lang(:,1) = temp_time_lang(1:minl,2);
temp_rg_lang(:,2) = rg_time_lang(1:minl,2);

minl = min(length(rg_time_nvt(:,1)),length(temp_time_nvt(:,1)));
temp_rg_nvt(:,1) = temp_time_nvt(1:minl,2);
temp_rg_nvt(:,2) = rg_time_nvt(1:minl,2);

    

%% Plot Temperature-Rg

h3 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('Temperature ($K$)','FontSize',20,'Interpreter','Latex')
ylabel('$R_g$ ($\sigma$)','FontSize',20,'Interpreter','Latex')
plot(temp_rg_nvt(:,1),temp_rg_nvt(:,2),'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3})
plot(temp_rg_lang(:,1),temp_rg_lang(:,2),'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3})
legendinfo{1} = 'NVT Thermostat';
legendinfo{2} = 'Langevin Thermostat';
legend(legendinfo,'FontSize',20,'Location','best')
legend boxoff

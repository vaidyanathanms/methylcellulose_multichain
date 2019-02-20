%% To analyze the shape anisotropy factor and the eigenvalues
% Many details are imported from anaeig.m
% Refer to timedata_sheet.csv/timesheet.dat for frames and time details

clear;
clc;
close all;
format long;

%% Color scheme

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];
brown = [0.6 0.2 0];violet = [0.5,0,0.5];gray = [0.75 0.75 0.75];
p4clr = {orange,'c',green,gold};
p8clr = {green,'r',gray,'b',orange,'m','c',gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Input Details

config   = 'alternate';  %block/alternate/jammed
trialarr = [6];
nmons    = 1000;
nchains  = 8;
initdist = '0.07';

%% Initialize all flags

flagindeig =1;
flagindkap = 1;
globkap    = 0;

%% Initialize all arrays

tval = -1*ones(20000,nchains);
eigx = -1*ones(20000,nchains);
eigy = -1*ones(20000,nchains);
eigz = -1*ones(20000,nchains);
globeigx = -1*ones(20000,1);
globeigy = -1*ones(20000,1);
globeigz = -1*ones(20000,1);

%% Main Calculation Begins Here

%First check whether path exists. If yes, find how many dump files exist

for trcnt = 1:length(trialarr)
    
    trialnum = trialarr(trcnt);
    
    dirname = sprintf('../trial%d/n%d/%s/dist_%s',trialnum,nchains,config...
        ,initdist);
    
    if ~isdir(dirname)
        errorMessage=sprintf('Error: Folder\t%s does not exist',dirname);
        uiwait(warndlg(errorMessage));
        return;
    else
        fprintf('Analyzing\t%s\n',dirname);
    end
    
    fylepattern = fullfile(dirname,'dump*.lammpstrj'); %Change to reqd pattern.
    dumplist    = dir(fylepattern);
    
    nfiles = length(dumplist);
    fprintf('Total number of files: \t%g\n',nfiles);
    
    %Cross refer with the doc from excel sheet to obtain deltaT
    timefyle = strcat(dirname,'/timesheet.dat');
    if exist(timefyle,'file') ~= 2
        fprintf('Timefile does not exist \t%s\n', timefyle);
        fprintf('Check Global Excel Sheet for details \n');
        fprintf('If not, create timesheet.dat from the jobdetails \n');
        continue;
    else
        ftime   = fopen(timefyle,'r');
        tline   = strtrim(fgetl(ftime)); %Find all columns
        strarr  = strsplit(tline);
        for colcnt = 1:length(strarr)
            if strcmp(strarr{colcnt},'trial')
                trialcol = colcnt;
            elseif strcmp(strarr{colcnt},'nchains')
                chaincol = colcnt;
            elseif strcmp(strarr{colcnt},'dist')
                distcol = colcnt;
            elseif strcmp(strarr{colcnt},'dumpname')
                dumpcol = colcnt;
            elseif strcmp(strarr{colcnt},'delta_t')
                deltatcol = colcnt;
            elseif strcmp(strarr{colcnt},'delta_nfr')
                delnframecol = colcnt;
            else
                fprintf('%s\t%s\n','Unknown Keyword',strarr{colcnt});
            end
        end
    end

    deltaT  = zeros(nfiles,1);
    deltafr = zeros(nfiles,1);

    
    %Store all deltaT and Number of Frames
    while ~feof(ftime)
        tline   = strtrim(fgetl(ftime)); %Read rows line by line
        strarr  = strsplit(tline);
        if str2double(strarr{chaincol}) ~= nchains || str2double(strarr{trialcol}) ...
                ~= trialnum
            fprintf('chaincol: \t %s\n',strarr{chaincol},'nchains: \t%g',nchains);
            fprintf('trialcol: \t %s\n',strarr{trialcol},'trialnum: \t%g',trialnum);
            errorMessage=sprintf('Error: Data do not match with input in %s\n', timefyle);
            uiwait(warndlg(errorMessage));
            return;
        else
            flagdumpfyle = -1;rowcnt = 1;
            while flagdumpfyle == -1 && rowcnt <= nfiles %Find row of dumpfyle
                if strcmp(dumplist(rowcnt).name,strarr{dumpcol})
                    flagdumpfyle = 1;
                    deltaT(rowcnt,1)  = str2double(strarr{deltatcol});
                    deltafr(rowcnt,1) = str2double(strarr{delnframecol});
                else
                    rowcnt = rowcnt+1;
                end
            end
            if flagdumpfyle == -1
                errorMessage=sprintf('Error: Did not find in timesheet.dat%s\n',...
                    strarr{dumpcol});
                uiwait(warndlg(errorMessage));
                return;
            end
        end
    end
    
    %Now Process files one by one
    timecntrold = 0;timecntrnew=0;acttimeold = 0;
    for fylecnt = 1:nfiles
        
        baseFileName = dumplist(fylecnt).name;
        fullFileName = fullfile(dirname, baseFileName);
        
        %Analyzing individual eigenvalues
        if flagindeig == 1 || flagindkap == 1
            resultdir = strcat(dirname,sprintf('/outresults_%s_%d_%s',initdist,...
                nchains,config));
            anafylename = strcat(resultdir,sprintf('/indeig_%s',dumplist(fylecnt).name));
            if exist(anafylename,'file') ~= 2
                fprintf('Eigenvalue file does not exist \t%s\n', anafylename);
                continue;
            else
                fprintf('Analyzing dumpfile\t%s\n',anafylename);
                [outdata, nframes] = read_analyze_indeig(anafylename,nchains);
            end
            if nframes == 0
                fprintf('No/Corrupted data in \t %s\n',anafylename);
                continue;
            end
            fprintf('Number of frames in %s:\t%g\n',anafylename,nframes);
            for framecnt = 1:nframes
                timecntrnew = timecntrold + framecnt;
                tval(timecntrnew,1) = acttimeold + framecnt*deltaT(fylecnt,1)*deltafr(fylecnt,1);
                eigx(timecntrnew,:) = outdata(framecnt,:,1);
                eigy(timecntrnew,:) = outdata(framecnt,:,2);
                eigz(timecntrnew,:) = outdata(framecnt,:,3);
            end
            clear outdata
        end
        
        %Analyzing global anisotropy values
        if globkap == 1
            resultdir = strcat(dirname,sprintf('/outresults_%s_%d_%s',initdist,...
                nchains,config));
            anafylename = strcat(resultdir,sprintf('/globeig_%s',dumplist(fylecnt).name));
            if exist(anafylename,'file') ~= 2
                fprintf('Eigenvalue file does not exist \t%s\n', anafylename);
                continue;
            else
                fprintf('Analyzing dumpfile\t%s\n',anafylename);
                [outdata, nframes] = analyze_kappa(anafylename,nchains);
            end
            if nframes == 0
                fprintf('No/Corrupted data in \t %s\n',anafylename);
                continue;
            end
            for framecnt = 1:nframes
                timecntrnew = timecntrold + framecnt;
                tval(timecntrnew,1) = acttimeold + framecnt*deltaT(fylecnt,1)*deltafr(fylecnt,1);
                globeigx(timecntrnew,1) = outdata(framecnt,1);
                globeigy(timecntrnew,1) = outdata(framecnt,2);
                globeigz(timecntrnew,1) = outdata(framecnt,3);
            end
            clear outdata
        end
        acttimeold  = acttimeold + nframes*deltafr(fylecnt,1)*deltaT(fylecnt,1);
        timecntrold = timecntrold + nframes;
    end
end

%% Analyze Global Anisotropy and Local Anisotropy

if globkap == 1
    glob_asphere  = globeigz(1:timecntrold) - 0.5*(globeigx(1:timecntrold)+globeigy(1:timecntrold));
    glob_acylind  = globeigy(1:timecntrold) - globeigx(1:timecntrold);
    glob_radgyr   = globeigx(1:timecntrold) + globeigy(1:timecntrold) + globeigz(1:timecntrold);
    glob_shapefac = (glob_asphere.^2 + 0.75*glob_acylind.^2)./glob_radgyr.^2;
end

%% Save to plot arrays

tplot    = tval(1:timecntrold,1);
eigxplot = eigx(1:timecntrold,:);
eigyplot = eigy(1:timecntrold,:);
eigzplot = eigz(1:timecntrold,:);

%clear tval eigx eigy eigz

%% Write log file
flog = fopen(sprintf('../allfile_data/log_%s_%d_%s.dat',initdist,nchains,config),'w');
fprintf(flog,'Trialnumber/Nchains/initdist/Config:%d\t%d\t%s\t%s\n',...
    trialnum,nchains,initdist,config);
fprintf(flog,'Number of files  analyzed:\t%d\n',nfiles);
fprintf(flog,'Name of files analyzed:\n');
for i = 1:nfiles
    fprintf(flog,'%d\t%s\n',i,dumplist(i).name);
end
fprintf(flog,'Total number of frames analyzed:\t%d\n',timecntrold);
fprintf(flog,'Total time analyzed:\t%g\n',acttimeold);
fclose(flog);

%% Plot Data

if globkap
    hz = figure; %kappa
    hold on
    box on
    set(gca,'FontSize',20)
    xlabel('Time ($\tau$)','FontSize',20,'Interpreter','Latex')
    ylabel('$\kappa^2_{g}$','FontSize',20,'Interpreter','Latex')
    if nchains == 4
        plot(tplot,glob_shapefac,'Color','k','LineWidth',2,'LineStyle',lsty{1});
    elseif nchains == 8
        plot(tplot,glob_shapefac,'Color','k','LineWidth',2,'LineStyle',lsty{1});
    end
    xline=1:max(1.2*max(tplot));
    yline=0.25*ones(length(xline),1);
    plot(xline,yline,'Color','k','LineWidth',2,'LineStyle',lsty{2});
    xlim([0 1.2*max(tplot)]);
    saveas(hz,sprintf('../allfigures/globshapefac_%d_%d_%s_%s',nchains,trialnum,initdist,config),'png')
end

if flagindeig
    sval = pl_indeig(nchains,trialnum,initdist,config,tplot,eigzplot);
end

if flagindkap
    asphere  = zeros(timecntrold,nchains);
    acylind  = zeros(timecntrold,nchains);
    radgyr   = zeros(timecntrold,nchains);
    shapefac = zeros(timecntrold,nchains);
    
    for i = 1:nchains
        asphere(:,i)  = eigzplot(:,i) - 0.5*(eigxplot(:,i)+eigyplot(:,i));
        acylind(:,i)  = eigyplot(:,i) - eigxplot(:,i);
        radgyr(:,i)   = eigxplot(:,i) + eigyplot(:,i) + eigzplot(:,i);
        shapefac(:,i) = (asphere(:,i).^2 + 0.75*acylind(:,i).^2)./radgyr(:,i).^2;
    end
    sval = pl_indkap(nchains,trialnum,initdist,config,tplot,shapefac); 
end



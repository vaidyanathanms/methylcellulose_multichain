function [outdata, nframes] = read_analyze_indeig(anafylename,nchains)
nframes = 0; dumoutdata = -1*ones(10000,nchains,3);
fana = fopen(anafylename,'r');flagstep = 0;
tline   = strtrim(fgetl(fana));
strarr  = strsplit(tline);
firstline = -1;
while ~feof(fana)
    if length(strarr) == 2
        flagstep = 1;
        if firstline == -1
            tstep  = 1;
            firstline = 1;
            fstep = str2double(strarr{1});
        else
            tstep = str2double(strarr{1}) - fstep + 1;
        end
        nchdum = str2double(strarr{2});
        nframes = nframes + 1;
        if nchdum ~= nchains
            fprintf('Warning: Number of chains do not match: %d on %d\t for %s\n',...
                tstep,anafylename)
        end
        tline   = strtrim(fgetl(fana));
        strarr  = strsplit(tline);
    elseif length(strarr) == 4
        if flagstep == 0
            fprintf('Warning: Did not find timestep %d for %s\n',...
                tstep,anafylename)
        end
        for chcnt = 1:nchdum
            chnum = str2double(strarr{1});
            dumoutdata(tstep,chnum,1) = str2double(strarr{2});
            dumoutdata(tstep,chnum,2) = str2double(strarr{3});
            dumoutdata(tstep,chnum,3) = str2double(strarr{4});
            if ~feof(fana)
                tline   = strtrim(fgetl(fana));
                strarr  = strsplit(tline);
            else
                fprintf('Finished analyzing \t%s\n with %d frames\n',...
                    anafylename,nframes);
                break;
            end
        end
        flagstep = 0;
    else
        errorMessage=sprintf('Error: Unknown line in data%s\n',tline);
        uiwait(warndlg(errorMessage));
        return;
    end
end
outdata = zeros(tstep,nchains,3);
cntr = 1;
while min(dumoutdata(cntr,:,:)) ~= -1
    outdata(cntr,:,:) = dumoutdata(cntr,:,:);
    cntr = cntr+1;
end
end
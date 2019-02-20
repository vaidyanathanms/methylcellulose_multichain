function [outdata, nframes] = analyze_kappa(anafylename,nchains)
    dumoutdata = -1*ones(10000,3);
    fana = fopen(anafylename,'r');
    tline   = strtrim(fgetl(fana));
    strarr  = strsplit(tline);
    nframes   = 0;
    while ~feof(fana)
        if length(strarr) ~= 4
            fprintf('ERROR: Unknown line %s\n',tline);
            nframes = 0;
            break;
        else
            if str2double(strarr{1}) ~= nchains
                fprint('ERROR: Mismatch in number of chains%d\t%d\n',nchains,...
                    str2double(strarr{1}));
            end
            nframes = nframes + 1;
            dumoutdata(nframes,1) = str2double(strarr{2});
            dumoutdata(nframes,2) = str2double(strarr{3});
            dumoutdata(nframes,3) = str2double(strarr{4});
            if ~feof(fana)
                tline   = strtrim(fgetl(fana));
                strarr  = strsplit(tline);
            else
                fprintf('Finished analyzing \t%s\n with %d frames\n',...
                    anafylename,nframes);
                break;
            end
        end
    end
    outdata = zeros(nframes,3);
    cntr = 1;
    while min(dumoutdata(cntr,:)) ~= -1
        outdata(cntr,:) = dumoutdata(cntr,:);
        cntr = cntr+1;
    end
end
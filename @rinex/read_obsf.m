% DDG: obtained and modified from PPPH by Berkay Bahadur
% function to read multiconstellation RINEX files
%

function [ obs,inf ] = read_obsf(files,options)

[ver] = rinex.r_rnxvers(files.rinex);

if ver>=3
    
    [inf] = rinex.r_rnxheadv3(files.rinex);
    
    if options.dcb == 1
        [dcb] = rinex.r_dcb(files.dcb);
        
        [obs] = rinex.r_rnxobsv3(files.rinex,inf,options,dcb);
    else
        [obs] = rinex.r_rnxobsv3(files.rinex,inf,options);
    end
elseif ver>=2
    
    [inf] = rinex.r_rnxheadv2(files.rinex);
    
    if options.dcb == 1
        [dcb] = rinex.r_dcb(files.dcb);
        
        [obs] = rinex.r_rnxobsv2(files.rinex,inf,options,dcb);
    else
        
        [obs] = rinex.r_rnxobsv2(files.rinex,inf,options);
    end
else
    errordlg('RINEX version is not valid !','RINEX version error');
    error('RINEX version is not valid !');
end
end


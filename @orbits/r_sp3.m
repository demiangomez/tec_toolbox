% DDG: obtained and modified from PPPH by Berkay Bahadur
% function to read multiconstellation SP3 orbits
%

function [sat,inf] = r_sp3(f_orbit,options,inf)

[fid,errmsg] = fopen(f_orbit);

if any(errmsg)
    errordlg('SP3 file can not be opened.','SP3 file error');
    error   ('SP3 file error');
end

sn  = rinex.max_sys_index;
sp3 = NaN(96,4,sn); 

while ~feof(fid)
    line = fgetl(fid);
    if strcmp(line(1),'#') && ~strcmp(line(1:2),'##')
        dat = sscanf(line(4:14),'%f');
        % ddg: add the date of the first observation
        inf.time.first = sscanf(line(4:33),'%f');
    end
    if strcmp(line(1:2),'##')
        sp3int = sscanf(line(25:38),'%f');
        inf.time.sp3int = sp3int;
        epno   = 86400/sp3int; 
        sp3    = NaN(epno,4,sn);
    end
    
    if line(1)=='+'
        temp = sscanf(line(5:6),'%d');
        if ~isnan(temp)
            NoSat = temp;
        end
    end
    
    % new epoch
    if line(1)=='*'
        ep   = sscanf(line(2:end),'%f',[1,6]);
        % DDG: was preventing from reading all the way to the end
        % if ep(1)~=dat(1) || ep(2)~=dat(2) || ep(3)~=dat(3)
            % continue
            % instead, add +1 day to read whatever comes next
        % end
        epno = ((ep(4)*3600 + ep(5)*60 + ep(6))/sp3int) + 1 + seconds(datetime(ep(1:3)) - datetime(dat(1:3)'))/sp3int;
        for k=1:NoSat
            line = fgetl(fid);
            if strcmp(line(2),'G') && (options.system.gps == 1)
                sno  = sscanf(line(3:4),'%d');
                if sno<=rinex.max_sv_index('G')
                    temp = sscanf(line(5:end),'%f',[1,4]);
                    % writing part
                    sp3(epno,1:3,sno) = temp(1:3)*1000; %meter
                    sp3(epno,  4,sno) = temp(4)*10^-6;  %second 
                end
            elseif strcmp(line(2),'R') && (options.system.glo == 1)
                sno  = rinex.min_sv_index('R') - 1 + sscanf(line(3:4),'%d');
                if sno<=rinex.max_sv_index('R')
                    temp = sscanf(line(5:end),'%f',[1,4]);
                    % writing part
                    sp3(epno,1:3,sno) = temp(1:3)*1000; %meter
                    sp3(epno,  4,sno) = temp(4)*10^-6;  %second 
                end
            elseif strcmp(line(2),'E') && (options.system.gal == 1)
                sno  = rinex.min_sv_index('E') - 1 + sscanf(line(3:4),'%d');
                if sno<=rinex.max_sv_index('E')
                    temp = sscanf(line(5:end),'%f',[1,4]);
                    % writing part
                    sp3(epno,1:3,sno) = temp(1:3)*1000; %meter
                    sp3(epno,  4,sno) = temp(4)*10^-6;  %second 
                end
            elseif strcmp(line(2),'C') && (options.system.bds == 1)
                sno  = rinex.min_sv_index('C') - 1 + sscanf(line(3:4),'%d');
                if sno<=rinex.max_sv_index('C')
                    temp = sscanf(line(5:end),'%f',[1,4]);
                    
                    sp3(epno,1:3,sno) = temp(1:3)*1000; %meter
                    sp3(epno,  4,sno) = temp(4)*10^-6;  %second 
                end
            end        
        end        
    end
end

% DDG: create the time stamps
sp3(:, end+1, :) = repmat((0:sp3int:epno * sp3int - 1)', [1 1 size(sp3,3)]);

sat.sp3    = sp3;

fclose('all');
end

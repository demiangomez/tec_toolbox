function [epoch,time,tec,rms,lats,lons,hgts,rb,dcbs,dcbr,rcvs]=readionex(file)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read ionex file
% [func]   : read ionex file
% [argin]  : file = file path
% [argout] : epoch = epoch of first map [year,month,day,hour,min,sec]
%            time  = time vector (sec)
%            tec   = tec map
%                tec(x,y,h,t) : lats(x),lons(y),hgts(h),time(t) vtec (TECU)
%            rms   = tec map rms
%                rms(x,y,h,t) : lats(x),lons(y),hgts(h),time(t) vtec rms (TECU)
%            lats  = longitudes (deg)
%            lons  = latitudes (deg)
%            hgts  = heights (km)
%            rb    = base radius (km)
%            dcbs  = satellite dcbs (nsec)
%            dcbr  = station dcbs (nsec)
%            sats  = satellite list
%            rcvs  = station list
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
% DDG: modified to be included in the ionosphere analysis software
%-------------------------------------------------------------------------------
epoch=[]; time=[]; tec=[]; rms=[]; lats=[]; lons=[]; hgts=[]; rb=0;
dcbs=zeros(rinex.max_sys_index,1); dcbr=[]; rcvs={};

fd=fopen(file,'rt');
if fd<0
    gt_log('no ionex file          : %s',file);
    return;
end
% read ionex header
[epoch,lats,lons,hgts,rb,nexp,dcbs,dcbr,rcvs]=ReadHead(fd);

if ~isempty(epoch)
    % read tecmap
    [td,ts]=caltomjd(epoch);
    [time,tec,rms]=ReadTecMap(fd,td,ts,lats,lons,hgts,nexp);
else
    gt_log('ionex file read error  : %s',file);
end
fclose(fd);

function [mjd, sec] = caltomjd(varargin)
    if nargin > 1
        yr = varargin{1};
        mn = varargin{2};
        dy = varargin{3};
        hr = varargin{4};
        mi = varargin{5};
        sc = varargin{6};
        inp = [yr mn dy hr mi sc];
    else
        inp = varargin{1};
    end
    
    mjd = mjuliandate(inp);
    sec = (mjd - fix(mjd)) .* 86400;
    mjd = fix(mjd);

% read ionex header ------------------------------------------------------------
function [epoch,lats,lons,hgts,rb,nexp,dcbs,dcbr,rcvs]=ReadHead(fd)
epoch=[]; lats=[]; lons=[]; hgts=[]; rb=0; nexp=0; dcbs=[]; dcbr=[]; rcvs={};

while 1
    str=fgetl(fd); if ~isstr(str), break, end
    label=SubStr(str,61,20);
    if findstr(label,'EPOCH OF FIRST MAP')
        epoch=sscanf(SubStr(str,1,42),'%d %d %d %d %d %f',6)';
    elseif findstr(label,'BASE RADIUS')
        rb=StrToNum(str,2,9);
    elseif findstr(label,'HGT1 / HGT2 / DHGT')
        x=sscanf(SubStr(str,2,19),'%f %f %f');
        if x(3)==0, hgts=x(1); else hgts=(x(1):x(3):x(2))'; end
    elseif findstr(label,'LAT1 / LAT2 / DLAT')
        s=sscanf(SubStr(str,2,19),'%f %f %f');
        lats=(s(1):s(3):s(2))';
    elseif findstr(label,'LON1 / LON2 / DLON')
        s=sscanf(SubStr(str,2,19),'%f %f %f');
        lons=(s(1):s(3):s(2))';
    elseif findstr(label,'EXPONENT')
        nexp=StrToNum(str,2,5);
    elseif findstr(label,'START OF AUX DATA')
        if findstr(str,'DIFFERENTIAL CODE BIASE')
            [dcbs,dcbr,rcvs]=ReadDcbData(fd);
        end
    elseif findstr(label,'END OF HEADER'), break, end
end
rcvs=rcvs';

% read dcb data ----------------------------------------------------------------
function [dcbs,dcbr,rcvs]=ReadDcbData(fd)
dcbs=zeros(rinex.max_sys_index,1); dcbr=[]; rcvs={};
while 1
    str=fgetl(fd); if ~isstr(str), break, end
    label=SubStr(str,61,20);
    if findstr(label,'PRN / BIAS / RMS')
        if findstr(SubStr(str,4,1), 'G')
            dcbs(StrToNum(str,5,2), 1) = StrToNum(str,7,10);
        elseif findstr(SubStr(str,4,1), 'R')
            dcbs(StrToNum(str,5,2) + 32, 1) = StrToNum(str,7,10);
        elseif findstr(SubStr(str,4,1), 'E')
            dcbs(StrToNum(str,5,2) + 48, 1) = StrToNum(str,7,10);
        end 
    elseif findstr(label,'STATION / BIAS / RMS')
        rcvs={rcvs{:},upper(SubStr(str,7,4))};
        dcbr=[dcbr;StrToNum(str,27,10)];
    elseif findstr(label,'END OF AUX DATA'), break, end
end

% read tec map ----------------------------------------------------------------
function [time,tec,rms]=ReadTecMap(fd,td,ts,lats,lons,hgts,nexp)
tec_start=0; rms_start=0; time=[]; it=[];
tec=repmat(nan,[length(lats),length(lons),length(hgts),100]); rms=tec;
while 1
    str=fgetl(fd); if ~isstr(str), break, end
    label=SubStr(str,61,20);
    if findstr(label,'START OF TEC MAP'), tec_start=1;
    elseif findstr(label,'END OF TEC MAP'), tec_start=0;
    elseif findstr(label,'START OF RMS MAP'), rms_start=1;
    elseif findstr(label,'END OF RMS MAP'), rms_start=0;
    
    elseif findstr(label,'EPOCH OF CURRENT MAP')
        [tdd,tss]=caltomjd(sscanf(SubStr(str,1,36),'%d %d %d %d %d %f',6)');
        t=(tdd-td)*86400+tss-ts;
        if isempty(it), time=t; it=1;
        else
            it=find(time==t);
            if isempty(it), time=[time;t]; it=length(time); end
        end
    elseif findstr(label,'LAT/LON1/LON2/DLON/H')
        x=zeros(5,1); for n=1:5, x(n)=StrToNum(str,6*n-3,6); end
        lat=x(1); lon=x(2):x(4):x(3); hgt=x(5);
        i=find(lats==lat); k=find(hgts==hgt);
        for j=1:16:length(lons)
            str=fgetl(fd); if ~isstr(str), break, end
            x=sscanf(str,'%f');
            if tec_start
                tec(i,j:j+length(x)-1,k,it)=x;
            elseif rms_start
                rms(i,j:j+length(x)-1,k,it)=x;
            end
        end
    end
end
tec=tec(:,:,:,1:length(time));
rms=rms(:,:,:,1:length(time));
tec(find(tec==9999))=nan;
rms(find(rms==9999))=nan;
rms=rms*10^nexp; tec=tec*10^nexp;

% substring --------------------------------------------------------------------
function sstr=SubStr(str,pos,ns), sstr=str(pos:min(pos+ns-1,length(str)));

% string->number ---------------------------------------------------------------
function num=StrToNum(str,pos,ns)
num=sscanf(str(pos:min(pos+ns-1,length(str))),'%f',1);
if isempty(num), num=nan; end

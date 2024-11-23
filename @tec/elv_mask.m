% DDG: obtained and modified from PPPH by Berkay Bahadur
% function to calculate elevation angle and azimuth
%
function [data, observables] = elv_mask(observables, info, sp3, options)

en = size(observables.st,1);
sn = size(observables.st,2);
st = sum(observables.st);

r_xyz = info.rec.pos;
elv = NaN(en,sn);
azm = NaN(en,sn);

for k=1:sn
    if st(k) > 0

        s_xyz = sp3(:,1:3,k);

        [az,elev,~] = tec.azelle(s_xyz, r_xyz);
        % output in radians, convert
        az = az * 180/pi;
        elev = elev * 180/pi;
        
        elv(:,k) = elev;
        azm(:,k) = az;
        
        % mask out observables below the mask
        observables.st(elev < options.elevation_mask,k) = 0;
    end
end

data.elv = elv;
data.azm = azm;
end


function [stec, ti] = sample_tec_map(ep, time, ionex, ipp, elv, options)
%SAMPLE_TEC_MAP Summary of this function goes here
%   Detailed explanation goes here
    
    % find the time in the ionex map
    ti = find(ismember(time, round(ionex.iono_map.time)));
    if ~isempty(ti)
        stec = nan(size(ti));
        for i = 1:size(ti,1)
            % only use SVs with high elevation angle
            % will return nans if low elevation
            if elv(ti(i)) > 30 
                % find this time in the vector
                j = find(round(ionex.iono_map.time) == time(ti(i)));
                % compute the vtec from the map
                vtec = interp2(ionex.iono_map.lons, ionex.iono_map.lats, ionex.iono_map.tec(:,:,1,j), ipp(ti(i), 2), ipp(ti(i), 1), 'cubic');
                % convert the vtec into stec
                stec(i) = vtec ./ tec.get_emf(elv(ti(i)), options);
            end
        end
        % return the time positions in the original ep vector
        ti = find(ismember(ep, time(ti)));
    else
        stec = [];
    end
end


function [pstn_N,pstn_E,ref_st, pref_N, pref_E] = create_RSRF(lat,lon,alt_ion)

    % determine the piercing point that is closer to the center of the
    % array. The center of the array is determined as the average lat and
    % lon for each point

    R = 6371000; %-alt_ion/1000;
    nst = size(lat,2);
    npp = size(lon,1);
    
    ac_lat = repmat(mean(lat,2),1,nst);
    ac_lon = repmat(mean(lon,2),1,nst);
    ac_lat = ac_lat(:);
    ac_lon = ac_lon(:);
    % determine the closest station to the average coordinates
    station_dist = distance(ac_lat,ac_lon, lat(:), lon(:))*pi/180.*R;
    station_dist = reshape(station_dist,npp,nst);
    
    % get the station that has the smallest average to the array center
    % that will be the reference station
    [~, ref_st] = min(mean(station_dist));
    
    % start by projecting npp times the array
    % project the center of the array to get the false northing and easting
    pref_N = zeros(npp,1);
    pref_E = zeros(npp,1);
    % I use for loops because of a limitation on the functions to project
    % the coordinates.
    if max(abs(ac_lat))>50
        % near the poles! use stereographic
        for i = 1:npp
            [pref_E(i,1),pref_N(i,1)] = fk.polarstereo_fwd(lat(i,ref_st), lon(i,ref_st),R + alt_ion,[],-70, lon(i,ref_st));
        end
    else
        for i = 1:npp
            [pref_N(i,1),pref_E(i,1)] = ell2utm(lat(i,ref_st).*pi/180, lon(i,ref_st).*pi/180, lon(i,ref_st).*pi/180);
        end
    end
    pstn_N = zeros(npp,nst);
    pstn_E = zeros(npp,nst);
    if max(abs(ac_lat))>50
        % near the poles! use stereographic
        for i = 1:npp
            for j = 1:nst
                [pstn_E(i,j),pstn_N(i,j)] = fk.polarstereo_fwd(lat(i,j), lon(i,j), R + alt_ion,[],-70, lon(i,ref_st));
            end
        end
    else
        for i = 1:npp
            for j = 1:nst
                [pstn_N(i,j),pstn_E(i,j)] = ell2utm(lat(i,j).*pi/180, lon(i,j).*pi/180, lon(i,ref_st).*pi/180);
            end
        end
    end

    % coordinates of the piercing points in the npp coordinate systems
    pstn_N = pstn_N/1000 - repmat(pref_N,1,nst)/1000; % in km
    pstn_E = pstn_E/1000 - repmat(pref_E,1,nst)/1000; % in km
    
    pref_N = zeros(npp,1);
    pref_E = zeros(npp,1);
    % obtain the coordinates on a ground-fixed RF
    if max(abs(ac_lat))>50
        % near the poles! use stereographic
        for i = 1:npp
            [pref_E(i,1),pref_N(i,1)] = fk.polarstereo_fwd(lat(i,ref_st), lon(i,ref_st),R + alt_ion,[],-70, mean(lon(:,ref_st)));
        end
    else
        for i = 1:npp
            [pref_N(i,1),pref_E(i,1)] = ell2utm(lat(i,ref_st).*pi/180, lon(i,ref_st).*pi/180, mean(lon(:,ref_st)).*pi/180);
        end
    end
    
    % get the change in coordinates using gradient
    pref_N = pref_N./1000;
    pref_E = pref_E./1000;
end
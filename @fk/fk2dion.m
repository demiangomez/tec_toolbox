function [FK, stacked_vector, stacked_signal, tec_shifted, max_stack_time, time_index, ...
          array_center, array_az, ap_vel, corrected_time, corrected_freq] = fk2dion(time, tec, lat, lon, ref_st, pstn_N, pstn_E, pref_N, pref_E, splint, slowness, transf)
%
%   function to perform ionospheric beam forming
%   time = maxtrix of time points in seconds since midnight!!
%           each col represents a station
%           each row is an observation
%           the matrix has only the intersection of observations
%   tec = observation values of Total Electron Content
%   lat and lon = positions of the piercing points
%   splint = sample interval
%   slowness = slowness range
% processed_sites,apcoords_sites ONLY NECESSARY TO GENERATE SYNTHETICS

    % if the piercing points were static, then we would shift each TEC
    % seismogram by a constant time p*ri where p is the slowness vector and
    % ri is the vector from the center of the array
    
    % px and py are considered constant.
    % the relative time between station r (reference) and station i is:
    % t(r,i) = px*(xr - xi) + py*(yr - yi)
    % but, the xr - xi and yr - yi change for each time
    
    nst = size(lat,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CAREFUL HERE! I AM A SURVEYOR, SO THIS CRAP OF USING
    % X = E AND Y = N IS NOT IN MY HEAD!!!
    % X IS NORTH! AND Y IS EAST, AS IN MAP PROJECTIONS
    % SI NO TE GUSTA, PROBLEMA TUYO.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [sy,sx] = meshgrid(slowness,slowness);
    
    % REFERENCE TIME: the time of this process is regulated by the time at
    % the reference station. Everything outside this window will be
    % replaced by NaNs during the interpolation procedure
    min_ref_t = min(time(:,ref_st));
    max_ref_t = max(time(:,ref_st));
    min_t = repmat(min_ref_t, size(sx));
    max_t = repmat(max_ref_t, size(sy));
    resample_vector = linspaceNDim(min_t,max_t,size(time,1));
    
    % preallocate
    t_shift = zeros([size(sx) size(time,1) nst]);
    corrected_time = zeros([size(sx) size(time,1) nst]);
    corrected_tec = zeros([size(sx) size(time,1) nst]);
    
    % make a 3D array with sx and sy that matched the pstn_N/E_sx/sy size
    ssx = repmat(sx,[1,1,size(time,1)]);
    ssy = repmat(sy,[1,1,size(time,1)]);

    % loop for each station
    for i=1:nst
        disp(['Correcting Doppler shift for IPP ' num2str(i) ' of ' num2str(nst)])
        
        % turn the coordinates into a 3D array
        temp_N(1,1,:) = pstn_N(:,i);
        temp_E(1,1,:) = pstn_E(:,i);        
        
        % repeat the 3D array in the X Y plane to match sx and sy size
        pstn_N_sx = repmat(temp_N,[size(sx) 1]);
        pstn_E_sy = repmat(temp_E,[size(sx) 1]);
        
        % IMPORTANT!! READ THIS TO UNDERSTAND WHAT THE 4D STRUCTURE IS ====
        % 4D array with 3D matrix (sx sy time station)
        t_shift(:,:,:,i) = pstn_N_sx.*ssx + pstn_E_sy.*ssy;
        
        % build the same structure for time and tec
        temp_time(1,1,:) = time(:,i);
        
        if ~isempty(transf)
            % an G&H tranfer function was provided.
            % obtain the Hilbert tranform of the TEC data
            hilb = hilbert(tec(:,i));
            temp_tec(1,1,:) = hilb;
        else
            temp_tec(1,1,:) = tec(:,i);
        end
        
        time3d = repmat(temp_time,[size(sx) 1]);
        tec3d = repmat(temp_tec,[size(sx) 1]);
        
        if ~isempty(transf)
            tec3d = real(tec3d./transf(:,:,:,i));
        end
        
        corrected_time(:,:,:,i) = time3d - t_shift(:,:,:,i);
        % new_samp_time contains the resampled times (multiple of
        % splint) to interpolate the tec data.
        % the min max values are determined by the reference station.
        % Although there may not be data, the interpolation function will
        % replace those by NaNs
        % array of points (first dimension) of corrected_time
        corrected_tec(:,:,:,i) = fk.LeanInterp(corrected_time(:,:,:,i),tec3d,resample_vector, splint);
    end
    
    % now, apply the inverse tranfer function
    %if ~isempty(transf)
    %    corrected_tec = real(corrected_tec./transf);
    %end
    
    % free some memory!!!
    clear temp_time temp_tec temp_N temp_E ssx ssy temp_RN temp_RE

    energy = sum(abs(sum(corrected_tec,4)).^2,3); 
    FK = energy/max(max(energy));
    
    % get the max value of the KK plot
    [maxFK,r] = max(FK);
    [maxFK,c] = max(maxFK);
    
    [array_az,ap_vel] = cart2pol(slowness(c), slowness(r(c)));
    array_az = 90-array_az*180/pi;
    if array_az<0
        array_az = array_az + 360;
    end
    
    ax = slowness(c);
    ay = slowness(r(c));

    % corrected_tec is a 4D vector with each K in i,j, time in k and station in l
    % stacked vector has the stack of all stations for each K
    stacked_vector = sum(corrected_tec,4);
    % this is the stacked signal that has max amplitude
    stacked_signal(:,:) = stacked_vector(r(c),c,:);
    
    % get the time index that has the maximum TEC amplitude
    [~, time_index] = max(abs(stacked_signal));
    % array's center position at the time of max amplitude
    array_center = [lat(time_index,ref_st),lon(time_index,ref_st)];
    % the time of the maximum amplitude at the array
    max_stack_time =time(time_index,ref_st);
    
    % to remove the doppler shift from the moving array from the stacked
    % signal, take the dot product of the resulting slowness with the
    % change in coordinate in reference station (w.r.t the time index of
    % maximum energy). This represents the time shift that has to be
    % applied to each observation to bring it into a "static" frame
    dt = dot(repmat(-[ax ay], size(pref_N)), [pref_N-pref_N(time_index) pref_E-pref_E(time_index)], 2);

    corrected_time = time(:, ref_st) - dt;
    
    % run an FFT on the signal to obtain frequency
    % make a uniformly sampled dataset
    tu = (min(corrected_time):5:max(corrected_time))';
    su = interp1(corrected_time, stacked_signal, tu);

    [freq_response,freq_index] = freqz(su,1,size(su,1),1/5);

    pM = max(abs(freq_response)); %magnitude
    corrected_freq = freq_index(abs(freq_response)==pM); %frequency

    for i = 1:nst
        tec_shifted(:,i) = corrected_tec(r(c),c,:,i);
    end
end
function tec = sample_model(time, lat, lon, tec1, epochs, lla, llo)
    % specifically written to read the SAMI3 model sent by Joe Huba
    %
    
    lo(:,:) = lon(:,1,:);
    la(:,:) = lat(:,1,:);

    % convert to negative longitude
    lo(lo> 180) = lo(lo> 180)-360;

    % origin time of sami3
    ot = datetime(2020,12,13,0,0,4);

    % time vector for sami3
    ts = seconds((ot + hours(time)) - datetime(2020,12,14,0,0,0));
    
    % in case we use methods that require splint
    splint = epochs(2) - epochs(1);

    %% sample section
    itec = nan(size(lla));

    for i = 1:size(epochs,1)
        % i the index of the GNSS data
        % the index in sami3
        s = find(epochs(i) == round(ts));

        if ~isempty(s)
            % get the current date
            tc1(:,:) = tec1(:,:,s);

            F1=scatteredInterpolant(double(lo(:)),double(la(:)),double(tc1(:)),'linear');

            % sample sami at the IPP
            itec(i, :) = F1(llo(i, :), lla(i, :));
        end
    end

    %% interpolate sami3
    for i = 1:size(itec,2)
        inan = ~isnan(itec(:, i));
        if sum(inan) > 2
            itec(:, i) = interp1(epochs(inan), itec(inan, i), epochs, 'spline', NaN);
            
            inan = ~isnan(itec(:, i));

            % CAREFUL: do not touch the detrending polynomial degree
            itec(inan,i) = fk.detrendnonlin(itec(inan, i), 7);

            % detrend using SSA
            %itec(inan,i) = detrendssa(itec(inan,i));

            % detrend golay
            %itec(inan,i) = itec(inan, i) - smooth(t_i(inan), itec(inan,i), 30*60/splint, 'sgolay', 2);
        end
    end

    tec = itec;

end


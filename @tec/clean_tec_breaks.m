function oc = clean_tec_breaks(t, o, int)
    
    run = 1;
    
    while run < 10
        
        % rigorous gap identification: the goal is to minimize their effects
        seg = tec.get_segments(t, o, int, 1);

        if size(seg, 1) > 1
            % take segment i + 1 and add constant to make it 
            for i = 2:size(seg, 1)
                t1s = find(t == seg(i-1, 1));
                t1e = find(t == seg(i-1, 2));
                t2s = find(t == seg(i  , 1));
                t2e = find(t == seg(i  , 2));

                if t1e - 20 > 1
                    % find the dtec from the 20 previous points
                    dtec = mean(diff(o(t1e - 20:t1e)));
                elseif t1e > 2
                    % if less than 20 points, do what you can...
                    dtec = mean(diff(o(1:t1e)));
                else
                    % only one datapoint
                    dtec = 0;
                end
                o(t2s:end) = o(t2s:end) - (o(t2s) - o(t1e)) + dtec * (t2s - t1e);
                o(1:t2e) = tec.fill_data_gaps(t(1:t2e), o(1:t2e));
            end
        end

        % detect any large spikes in the data an replace by nan
        dtec = abs(diff(o) ./ diff(t)) > 5;

        if any(dtec)
            j = find(dtec);
            for jj = j'
                % replace the point by a nan, this will tend to fix the
                % issue next time we rerun
                o(jj+1) = nan;
            end
        end
        
        if any(isnan(o))
            run = run + 1;
        else
            break
        end
    end
    oc = o;
end


function segments = get_segments(ep, data, splint, tolerance)
    % find the indeces of start and stop of data
    s = find(diff(~isnan([NaN; data; NaN]))== 1);
    e = find(diff(~isnan([NaN; data; NaN]))==-1);
    
    gaps = [s e - 1];
    
    if size(gaps, 1) > 1
        % check that all gaps are > tolerance
        j = 1;
        remove = [];
        for i = 1:size(gaps,1) - 1
            if gaps(i+1, 1) - gaps(j, 2) < tolerance / splint
                % gap is smaller than tolerance, replace end index and
                % flag for removal at the end of the loop
                gaps(j, 2) = gaps(i+1, 2);
                remove = [remove i+1];
            else
                j = i + 1;
            end
        end
        if ~isempty(remove)
            gaps(remove, :) = [];
        end
    end
    
    % DDG: size issue -> if single row then wrong shape is returned
    if size(gaps, 1) > 1 && length(ep) > 1
        segments = ep(gaps);
    elseif size(gaps, 1) == 1 && length(ep) > 1
        segments = ep(gaps)';
    elseif length(ep) == 1
        segments = ep(gaps);
    end
    
end
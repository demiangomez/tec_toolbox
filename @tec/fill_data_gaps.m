function p = fill_data_gaps(t, o)
%FILL_DATA_GAPS Summary of this function goes here
%   Detailed explanation goes here
    
    tf = t(isnan(o));
    td = t(~isnan(o));
    od = o(~isnan(o));
    
    % use a spline to fill in the gaps of data
    oi = interp1(td, od, tf, 'pchip', 0);
    
    p = o;
    p(isnan(o)) = oi;
end


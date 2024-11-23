function [iec, ifb, ionex_tec] = compute_ifb(self, info, ionex, options)

    st = sum(self.st);
    
    % system frequencies
    [f, ~] = rinex.frequencies;
    % constants
    A = 40.308;
    c = 0.299792458e9;
    GAP_TOLERANCE = 300; % tolerance for data gaps in seconds
    
    % system names
    GPS = 1;
    GLO = 2;
    GAL = 3;
    BEI = 4;
    % for ifb computation
    dstec{GPS} = []; dstec{GLO} = []; dstec{GAL} = []; dstec{BEI} = [];
    
    % temp variable
    t = self.ep;
    iec = self.slant_tec;
    % create the space to save the ionex_tec
    ionex_tec = nan(size(self.slant_tec));
    
    for i = 1:size(st,2)
        if st(i) > 0
            
            % find breaks in the data
            segments = tec.get_segments(t, iec(:, i), info.time.int, GAP_TOLERANCE);

            for g = 1:size(segments, 1)
                % get the time index for this segment
                t1 = find(t == segments(g, 1));
                t2 = find(t == segments(g, 2));
                
                elev = self.elv(t1:t2, i);
                
                if isobject(ionex)
                    % find the intersection between the TEC map and the ipp
                    [stec, ti] = tec.sample_tec_map(self.ep, t(t1:t2), ionex, self.ipp(t1:t2, :, i), elev, options);
                    % compute the difference between the stec values from 
                    % the map and the stec without the ifb
                    if ~isempty(ti)
                        if i<=rinex.max_sv_index('gps')
                            sys = GPS;
                        elseif i<=rinex.max_sv_index('glonass') % GLONASS
                            sys = GLO;
                        elseif i<=rinex.max_sv_index('galileo') % GALILEO
                            sys = GAL;
                        else        % BEIDOU
                            sys = BEI;
                        end
                        
                        ionex_tec(ti, i) = stec .* 1e16;
                        % values in TECU, convert to TEC
                        dstec{sys} = [dstec{sys}; stec .* 1e16 - iec(ti, i)];
                    end
                end
            end
        end
    end
    
    ifb = nan(4,1);
    
    if ~all(cellfun(@isempty, dstec))
        
        % report the interfrequency bias for each system (which is
        % frequency dependent)
        
        for i = 1:size(st,2)
            f1 = f(i, 1);
            f2 = f(i, 2);
            
            if st(i) > 0
                F = (f1^2 * f2*c) / (A * (f1^2 - f2^2));
                
                if i<=rinex.max_sv_index('gps')
                    sys = GPS;
                elseif i<=rinex.max_sv_index('glonass') % GLONASS
                    sys = GLO;
                elseif i<=rinex.max_sv_index('galileo') % GALILEO
                    sys = GAL;
                else        % BEIDOU
                    sys = BEI;
                end
                
                stec_ifb = nanmean(dstec{sys});
                ifb(sys) = stec_ifb ./ (F .* f2);
                iec(:, i) = iec(:, i) + stec_ifb;
            end
        end
    else
        ifb = NaN;
        % save memory and do not output a full field of nans
        ionex_tec = NaN;
    end
end

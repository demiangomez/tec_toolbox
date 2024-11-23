function [iec] = get_ion(observables, info, data, ionex)

% GET_ION	Computes ionospheric electron content
%
%	[iec,iec_dot] = get_ion(L1,L2,C1,P1,P2)
%	L1, L2, C1, P1, P2 are observables matrices (output of read_rinexo)
%	tgd transmitter group delay from ephenerides, sv satellites available
%	iec = integrated electron content, unit = el/m^2
%	iec_dot = time derivative of integrated electron content
%       lgu = unbiased lg observable, unit = cycles
%
%               WARNING: IEC and LGU output here are corrected for the
%               phase ambiguity (using pseudorange data) but NOT for
%               the satellite and receiver interfrequency biases

    if ~isobject(ionex)
        % no ionex information, proceed with warning
        warning('No ionex information: transmitter group delay and interfrequency bias not removed.')
    end
    
    % determine which satellites have observations
    st = sum(observables.st);
    iec = nan(size(observables.st));
    
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
    t = observables.ep;
    
    for i = 1:size(st,2)
        if st(i) > 0
            % get the frequencies for this satellite
            f1 = f(i, 1);
            f2 = f(i, 2);

            % multiplication factor
            F = (f1^2 * f2*c) / (A * (f1^2 - f2^2));
            
            L1 = observables.l1(:, i); L2 = observables.l2(:, i);
            P1 = observables.p1(:, i); P2 = observables.p2(:, i);
            
            % DDG: search for any impossible P1 and P2 values: if range is
            % reported as < 19k km, then there is a problem and value
            % should be deleted!
            P2(P1 < 1e3) = NaN; P1(P1 < 1e3) = NaN;
            P1(P2 < 1e3) = NaN; P2(P2 < 1e3) = NaN;
            
            PG = (f2/c) * (P2-P1);
            LG = (L2 - (f2/f1)*L1);
            % LG vary in opposite sense from PG
            LG = -LG;
            
            % find breaks in the data
            segments = tec.get_segments(t, LG, info.time.int, GAP_TOLERANCE);

            for g = 1:size(segments, 1)
                % get the time index for this segment
                t1 = find(t == segments(g, 1));
                t2 = find(t == segments(g, 2));
                
                % clean the cycle slips and put back in structure
                LG(t1:t2) = tec.clean_tec_breaks(t(t1:t2), LG(t1:t2), info.time.int);
                
                % remove the ambiguity from LG
                % only consider 50% around the highest elevation angle
                elev = data.elv(t1:t2, i);
                tPG = PG(t1:t2); tLG = LG(t1:t2);
                
                N = nanmean(tPG(elev > max(elev)*0.5) - tLG(elev > max(elev)*0.5));
                
                if isobject(ionex)
                    % ionex information provided
                    iec(t1:t2, i) = F * (f2 .* ionex.tgd(i) * 1e-9 + N + tLG);
                    % this integrated electron content is still biased by
                    % the interfrequency bias. We will remove it later
                    % using the ionex map (if provided)
                else
                    iec(t1:t2, i) = F * (N + tLG);
                end
            end
        end
    end
end

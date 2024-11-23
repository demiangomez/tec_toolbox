function E = get_emf(el, options)

% GET_MF    Calculates mapping function:
%             el = vector of satellite elevation angles (degrees)
%             options.ionosphere_bottom = height of bottom of ionosphere (m)
%             options.ionosphere_top = height of top of ionosphere (m)
%             options.ionosphere_height = height of maximum electron density (peak F2)
%             options.emf_method = 1/2  --> 2 different ways of calculating mapping function
%             E = vector of mapping functions

    % define some constants (should be read from a file)
    R = 6371000.0;      % mean Earth radius
    
    % calculate mapping function, option 1
    if (options.emf_method == 1)
        Hb = options.ionosphere_bottom;      % height of bottom of ionosphere
        Ht = options.ionosphere_top;         % height of top of ionosphere
        % calculate ionospheric "thickness"
        Hion = Ht - Hb;        
        
        L = sqrt((R+Ht)^2 - R^2.*cosd(el).^2) - sqrt((R+Hb)^2 - R^2.*cosd(el).^2);
        E = Hion ./ L;
        % calculate mapping function, option 2
        
    elseif (options.emf_method == 2)
        H = options.ionosphere_height;
        Q = R / (R+H);
        beta = asin( Q .* cosd(el) );
        E = cos(beta);
    end

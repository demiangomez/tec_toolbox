classdef fk
    %FK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stations
        slowness
        fk_sv
        plot_results
        date
        base_folder
        time_start
        time_end
        interval
        polynomial
        epochs
        tec
        lat
        lon
        stnm
        FK
        stacked_vector
        stacked_signal
        corrected_time
        corrected_freq
        tec_shifted
        max_stack_time
        time_index
        array_center
        array_az
        ap_vel
        pstn_N
        pstn_E
        ref_st
        sites
        filter
        model
        alt_ionosphere
        pref_N
        pref_E
    end
    
    methods
        function self = fk(pp, cfg_file, varargin)
            
            ini = IniConfig();
            ini.ReadFile(cfg_file);
            
            self.stations     = strtrim(lower(split(ini.GetValues('[general]', 'stations', nan), ',')));
            self.date         = project.parse_date(ini.GetValues('[general]', 'date', nan));
            self.plot_results = true;

            self.fk_sv = ini.GetValues('[general]', 'fk_sv', nan);
            self.filter = ini.GetValues('[general]', 'filter', [0 0]);

            % slowness text
            slowness = split(ini.GetValues('[general]', 'slowness', nan), ':');
            self.slowness = str2double(slowness(1)):str2double(slowness(2)):str2double(slowness(3));

            % time window
            time_win = split(ini.GetValues('[general]', 'time_window', nan), ',');
            
            self.time_start = seconds(duration(time_win(1)));
            self.time_end   = seconds(duration(time_win(2)));
            self.interval   = ini.GetValues('[general]', 'interval', nan);
            self.polynomial = ini.GetValues('[general]', 'polynomial', nan);
            self.epochs     = (self.time_start:self.interval:self.time_end)';
            
            % a model to reduce observations
            self.model = ini.GetValues('[general]', 'model', nan);
            
            if ~isnan(self.model)
                % load the SAMI3 model from Huba
                lat = ncread(self.model, 'lat0');
                lon = ncread(self.model, 'lon0');
                tec = ncread(self.model, 'tec');
                tti = ncread(self.model, 'time');
            end
            
            self.base_folder = pp.tec_config.base_folder;
            
            %self.tec = nan(size(self.epochs,1), length(self.stations));
            
            %self.sites = zeros(length(self.stations), 3);
            
            self.alt_ionosphere = pp.options.ionosphere_height;
            
            % check if passing an override parameter 
            if nargin > 2
                for i = 1:2:length(varargin)
                    if ischar(varargin{i+1})
                        eval(sprintf('self.%s=%s;', varargin{i}, varargin{i+1}));
                    else
                        eval(sprintf('self.%s=%f;', varargin{i}, varargin{i+1}));
                    end
                end
            end

            j = 1;
            % load the results from the tec analysis
            for i = 1:length(self.stations)
                file = fullfile(pp.tec_config.base_folder, pp.tec_config.structure, pp.tec_config.filename);
                file = pp.replace_keys(file,  'station', self.stations{i}, 'date', self.date{1}, 'session', pp.rinex_config.session);
                
                % load the data
                d = load([file '.' pp.tec_config.format]);
                
                % get all the necessary information
                idx1 = ismember(d.data.ep, self.epochs);
                idx2 = ismember(self.epochs, d.data.ep);
                % find the sv
                ids = find(strcmp(d.data.sv_names, self.fk_sv));
                
                if ~isempty(ids) & ~isempty(find(idx1, 1)) & ~all(isnan(d.data.vertical_tec(idx1, ids)))
                    % detrend non lin does not work with nan values
                    temp = d.data.vertical_tec(idx1, ids) .* 1e-16;
                    inan = ~isnan(temp);
                    
                    % detrend first
                    if ~isnan(self.polynomial)
                        temp(inan) = fk.detrendnonlin(temp(inan), self.polynomial);
                    end

                    % apply model (if any)
                    if ~isnan(self.model)
                        temp = temp - self.sample_model(tti, lat, lon, tec, d.data.ep(idx1), d.data.ipp(idx1, 1, ids), d.data.ipp(idx1, 2, ids));
                    end
                    
                    % in case we sample outside of the model time period
                    % and one of the observations becomes NaN, do this
                    % again
                    inan = ~isnan(temp);
                    
                    if self.filter(1) ~= 0
                        % apply the CWT if requested
                        temp(inan) = fk.cwt_tec(temp(inan), self.interval, self.filter);
                    end
                    
                    self.tec(idx2, j) = temp;
                    self.lat(idx2, j) = d.data.ipp(idx1, 1, ids);
                    self.lon(idx2, j) = d.data.ipp(idx1, 2, ids);

                    self.stnm{j, 1} = self.stations{i};
                    self.sites(j, :) = d.data.info.lla;
                    j = j + 1;
                else
                    % error(['station ' self.stations{i} ' returned no valid SV or no data in the specificed time window'])
                end
            end
            
            if j < 3
                error('Intersection of time, SV, and stations did not produce enough observations to obtain kx and ky')
            end
            % cut the tec and lat lon array if there are missing
            % obsservations
            cut = all(~isnan(self.tec), 2);
            self.tec = self.tec(cut, :);
            self.lat = self.lat(cut, :);
            self.lon = self.lon(cut, :);
            self.epochs = self.epochs(cut, :);
            
            % create the reference frame for this array
            [self.pstn_N, self.pstn_E, self.ref_st, self.pref_N, self.pref_E] = fk.create_RSRF(self.lat, self.lon, self.alt_ionosphere);
            
            time = repmat(self.epochs, [1 size(self.tec,2)]);

            [self.FK, self.stacked_vector, self.stacked_signal, self.tec_shifted, self.max_stack_time, self.time_index, ...
             self.array_center, self.array_az, self.ap_vel, self.corrected_time, self.corrected_freq] = ...
                    fk.fk2dion(time, self.tec, self.lat, self.lon, self.ref_st, self.pstn_N, self.pstn_E, self.pref_N, self.pref_E, self.interval, self.slowness, []);
            
            if self.plot_results
                self.plots();
            end
        end
        
        function self = plots(self)
            time = repmat(self.epochs, [1 size(self.tec,2)]);
            
            fk.plot_FK(NaN, NaN, self.stnm, self.lat, self.lon, self.pstn_N, self.pstn_E, self.ref_st, time, ...
                       self.tec, self.tec_shifted, self.FK, self.slowness, self.stacked_vector, self.max_stack_time, ...
                       self.time_index, 0, self.sites, self.corrected_time)
        end
    end
    
    methods (Static)
        function y = detrendnonlin(x, varargin)

            if numel(varargin)>=1
                order = varargin{1};
            else
                order = 2;
            end

            p = polyfit((1:numel(x))', x(:), order);
            y = x(:) - polyval(p, (1:numel(x))');

        end
        
        function tecf = cwt_tec(tec, splint, filter_bounds)
            type = 'morlet';
            [Wx, ax] = cwt_fw(tec, 'morlet', 32, splint);
            Wx(ax <= 10.^filter_bounds(1),:) = 0;
            Wx(ax >= 10.^filter_bounds(2),:) = 0;
            tecf = cwt_iw(Wx, type, 32);
            tecf = tecf';    
        end
        
        [pstn_N, pstn_E,ref_st, pref_N, pref_E] = create_RSRF(lat,lon,alt_ion)
        [x, y]=polarstereo_fwd(phi,lambda,a,e,phi_c,lambda_0)
        [Yi] = LeanInterp(X, Y, Xi, splint)
        [FK, stacked_vector, stacked_signal, tec_shifted, max_stack_time, time_index, array_center, array_az, ap_vel, corrected_time, corrected_freq] = fk2dion(time, tec, lat, lon, ref_st, pstn_N, pstn_E, pref_N, pref_E, splint, slowness, transf)
        plot_FK(lat_epi, lon_epi, sites, lat, lon, pstn_N, pstn_E, ref_st, time, tec, tec_shifted, FK, slowness, stacked_vector, max_stack_time, time_index, synthetic, lla_sites, corrected_time)
        tec = sample_model(time, lat, lon, tec1, epochs, lla, llo)
    end
end


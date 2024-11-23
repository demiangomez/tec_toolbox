classdef project < handle
    %PROJECT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        options
        rinex_config
        orbits_config
        ionex_config
        tec_config
        rinex_files
        sp3orbits
        ionex_files
    end
    
    methods
        function self = project(cfg_file)
            %PROJECT Construct an instance of this class
            %   Detailed explanation goes here
            % load the options for this run
            
            if isempty(cfg_file)
                % dry run, continue
                return
            end
            
            ini = IniConfig();
            ini.ReadFile(cfg_file);
            
            % date range
            date_range = ini.GetValues('[general]', 'date_range', nan);
            
            if isnan(date_range)
                error('project error: date_range cannot be empty')
            else
                self.options.dates = self.parse_date(date_range);
            end
            
            % constellations
            constellations = strtrim(split(ini.GetValues('[general]', 'constellations', 'gps'), ','));
            
            if ismember('gps', constellations)
                self.options.system.gps = true;
            else
                self.options.system.gps = false;
            end
            if ismember('glonass', constellations)
                self.options.system.glo = true;
            else
                self.options.system.glo = false;
            end
            if ismember('galileo', constellations)
                self.options.system.gal = true;
            else
                self.options.system.gal = false;
            end
            if ismember('beidou', constellations)
                self.options.system.bds = true;
            else
                self.options.system.bds = false;
            end
            
            % elevations mask
            self.options.elevation_mask = ini.GetValues('[general]', 'elevation_mask', 10);
            
            % elevation mapping function
            self.options.emf_method = ini.GetValues('[general]', 'emf_method', 2);
            
            if self.options.emf_method == 2
                self.options.ionosphere_height = ini.GetValues('[general]', 'ionosphere_height', nan);
                if isnan(self.options.ionosphere_height)
                    error('project error: ionosphere_height must be set with emf = 2')
                end
            elseif self.options.emf_method == 1
                self.options.ionosphere_top = ini.GetValues('[general]', 'ionosphere_top', nan);
                self.options.ionosphere_bottom = ini.GetValues('[general]', 'ionosphere_bottom', nan);
                if isnan(self.options.ionosphere_top) | isnan(self.options.ionosphere_bottom)
                    error('project error: ionosphere_top and ionosphere_bottom must be set with emf = 1')
                end
            else
                error('project error: elevation mapping function (emf_method) must be = 1 or 2')
            end
            
            % output time system
            self.options.output_time = ini.GetValues('[general]', 'output_time', 'gps');
            
            if ~(strcmp(self.options.output_time, 'gps') | strcmp(self.options.output_time, 'utc'))
                error('project error: output_time must be either gps or utc')
            end
            
            % rinex section
            self.rinex_config.base_folder = ini.GetValues('[rinex]', 'base_folder', nan);
            
            if isnan(self.rinex_config.base_folder)
                error('project error: rinex.base_folder cannot be empty')
            end
            
            self.rinex_config.structure = ini.GetValues('[rinex]', 'structure', nan);
            
            if isnan(self.rinex_config.structure)
                error('project error: rinex.structure cannot be empty')
            end
            
            self.rinex_config.filename = ini.GetValues('[rinex]', 'filename', nan);
            
            if isnan(self.rinex_config.filename)
                error('project error: rinex.filename cannot be empty')
            end
            
            self.rinex_config.stations = strtrim(lower(split(ini.GetValues('[rinex]', 'stations', nan), ',')));
            
            if ~(iscell(self.rinex_config.stations) | ischar(self.rinex_config.stations))
                error('project error: rinex.stations cannot be empty')
            end
            
            self.rinex_config.session = ini.GetValues('[rinex]', 'session', 0);
            
            % orbits section
            self.orbits_config.base_folder = ini.GetValues('[orbits]', 'base_folder', nan);
            
            if isnan(self.orbits_config.base_folder)
                error('project error: orbits.base_folder cannot be empty')
            end
            
            self.orbits_config.orbit_name = ini.GetValues('[orbits]', 'orbit_name', nan);
            
            if isnan(self.orbits_config.orbit_name)
                error('project error: orbits.orbit_name cannot be empty')
            end
            
            % ionex section
            self.ionex_config.base_folder = ini.GetValues('[ionex]', 'base_folder', nan);
            
            if isnan(self.ionex_config.base_folder)
                error('project error: orbits.base_folder cannot be empty')
            end
            
            self.ionex_config.ionex_name = ini.GetValues('[ionex]', 'ionex_name', nan);
            
            if isnan(self.ionex_config.ionex_name)
                error('project error: ionex.ionex_name cannot be empty')
            end
            
            % tec section
            self.tec_config.base_folder = ini.GetValues('[tec]', 'base_folder', nan);
            
            if isnan(self.tec_config.base_folder)
                error('project error: tec.base_folder cannot be empty')
            end
            
            self.tec_config.structure = ini.GetValues('[tec]', 'structure', nan);
            
            if isnan(self.tec_config.structure)
                error('project error: ionex.ionex_name cannot be empty')
            end
            
            self.tec_config.format = ini.GetValues('[tec]', 'format', 'object');
            self.tec_config.filename = ini.GetValues('[tec]', 'filename', '{lstation}_{year2d}_{doy}');
        end
        
        function self = load_files(self)
            % scan the rinex folder for the indicated stations
            for j = 1:length(self.options.dates)
                
                % load the orbit
                orbit_file = fullfile(self.orbits_config.base_folder, self.orbits_config.orbit_name);
                orbit_file = self.replace_keys(orbit_file, 'date', self.options.dates{j});
                
                if ~exist(orbit_file, 'file')
                    warning(['Could not find orbit file in ' orbit_file])
                else
                    disp([' >> Loading orbit ' orbit_file])
                    self.sp3orbits{j,1} = orbits(orbit_file);
                end
                
                % load the ionex
                ionex_file = fullfile(self.ionex_config.base_folder, self.ionex_config.ionex_name);
                ionex_file = self.replace_keys(ionex_file, 'date', self.options.dates{j});
                
                if ~exist(ionex_file, 'file')
                    warning(['Could not find ionex file in ' ionex_file])
                else
                    disp([' >> Loading orbit ' ionex_file])
                    self.ionex_files{j,1} = ionex(ionex_file);
                end
                
                files = {};
                for i = 1:length(self.rinex_config.stations)
                    % insert rinex file
                    files{i, 1} = fullfile(self.rinex_config.base_folder, self.rinex_config.structure, self.rinex_config.filename);
                    
                    files{i, 1} = self.replace_keys(files{i, 1},  'station', self.rinex_config.stations{i}, 'date', self.options.dates{j}, 'session', self.rinex_config.session);
                    
                    if ~exist(files{i, 1}, 'file')
                        warning(['Could not find rinex file for ' self.rinex_config.stations{i} ' in ' files{i, 1}])
                    end
                    
                    files{i, 2} = fullfile(self.tec_config.base_folder, self.tec_config.structure, self.tec_config.filename);
                    files{i, 2} = self.replace_keys(files{i, 2},  'station', self.rinex_config.stations{i}, 'date', self.options.dates{j}, 'session', self.rinex_config.session);
                    files{i, 3} = self.tec_config.format;
                end
                
                self.rinex_files{j, 1} = files;
                self.rinex_files{j, 2} = self.options.dates{j};
                self.rinex_files{j, 3} = self.sp3orbits{j,1};
                self.rinex_files{j, 4} = self.ionex_files{j,1};
            end
        end
    end
    
    methods (Static)
        function dates = parse_date(str_date)
            
            str_dates = split(str_date, ',');
            
            dates = {};
            
            for i = 1:length(str_dates)
                sd = strtrim(str_dates{i});
                switch true
                    case contains(sd, '_')
                        % year doy
                        yr_dy = split(sd, '_');
                        [yyyy, mm, dd] = jd2cal(doy2jd(str2double(yr_dy{1}), str2double(yr_dy{2})));
                    case contains(sd, '-')
                        % gpswk gpswkday
                        gpswk = split(sd, '_');
                        [yyyy, mm, dd] = jd2cal(gps2jd(str2double(gpswk{1}), str2double(gpswk{2} * 86400), 0));
                    case contains(sd, '/')
                        % year month day
                        yyyymmdd = split(sd, '/');
                        yyyy = str2double(yyyymmdd{1}); mm = str2double(yyyymmdd{2}); dd = str2double(yyyymmdd{3});
                    case contains(sd, '.')
                        % fyear
                        [yyyy, mm, dd] = jd2cal(yr2jd(str2double(sd)));
                end
                
                dates{i} = datetime(yyyy, mm, dd);
            end
        end
        
        function output = replace_keys(input, varargin)
            
            p = inputParser;
            p.addParameter('date', datetime(1900,1,1))
            p.addParameter('session', nan)
            p.addParameter('station', nan)

            parse(p, varargin{:});
            
            strin = input;
            vars = fields(p.Results);
            
            for i = 1:length(vars)
                if strcmp(vars{i}, 'station') & ~isnan(p.Results.station)
                    strin = replace(strin, '{ustation}', upper(p.Results.station));
                    strin = replace(strin, '{lstation}', lower(p.Results.station));
                end
                if strcmp(vars{i}, 'date') & p.Results.date > datetime(1980,1,1)
                    strin = replace(strin, '{year}'  , sprintf('%i', year(p.Results.date)));
                    strin = replace(strin, '{year2d}', cell2mat(extractBetween(sprintf('%02i', year(p.Results.date)), 3, 4)));
                    strin = replace(strin, '{doy}'   , sprintf('%03i', day(p.Results.date, 'dayofyear')));
                    strin = replace(strin, '{month}' , sprintf('%02i', month(p.Results.date)));
                    strin = replace(strin, '{day}'   , sprintf('%02i', day(p.Results.date)));
                    strin = replace(strin, '{month}' , sprintf('%02i', month(p.Results.date)));
                end
                if strcmp(vars{i}, 'session') & ~isnan(p.Results.session)
                    strin = replace(strin, '{session}', sprintf('%i', p.Results.session));
                end
            end
            
            output = strin;
        end
    end
end


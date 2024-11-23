classdef rinex
    %RINEX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        observables
        info
        mean_epoch
    end
    
    properties (Constant)
        max_sys_svs = 36;
    end
    
    methods
        function self = rinex(rnx_files, options)
            %RINEX Construct an instance of this class
            %   Detailed explanation goes here
            
            if ~isa(rnx_files, 'cell')
                tmp = rnx_files;
                rnx_files = {};
                rnx_files{1} = tmp;
            end
            
            % sort the rinex files for this station
            rnx_files = sortrows(rnx_files);

            for i = 1:size(rnx_files,1)
                if strcmp(rnx_files{i}(end-2:end), 'd.Z') || strcmp(rnx_files{i}(end-2:end), '.gz')
                    % compressed rinex, uncompress
                    if isunix()
                        sep = ';';
                    else
                        sep = ' &';
                    end
                    system(['cd ' tempdir sep ' crz2rnx -f -c ' rnx_files{i}]);
                    % change the name to use the uncompressed version
                    % in temp
                    if strcmp(rnx_files{i}(end-2:end), 'd.Z')
                        % files end in .??d.Z
                        [~, ff, ext] = fileparts(rnx_files{i}(1:end-3));
                        files.rinex = [tempdir  ff ext 'o'];
                    else
                        % files end in .rnx.gz
                        [~, ff, ext] = fileparts(rnx_files{i}(1:end-4));
                        files.rinex = [tempdir  ff ext];
                    end
                    delete_file = true;
                else
                    files.rinex = rnx_files{i};
                    delete_file = false;
                end

                % create the structure

                % save the output into temporary variables
                [obs, info] = rinex.read_obsf(files, options);

                self.observables = obs;
                self.info = info;

                if delete_file
                    delete(files.rinex)
                end
            end
            % get mean epoch of file
            % self.mean_epoch = mean([datetime(self.info.time.first') datetime(self.info.time.last')]);
        end
    end
    
    methods (Static)
        [obs, info] = read_obsf(self, files, options)
        [ver] = r_rnxvers(f_obs)
        [obs] = r_rnxobsv3(f_obs,inf,options,dcb)
        [obs] = r_rnxobsv2(f_obs,inf,options,dcb)
        [inf] = r_rnxheadv3(f_obs)
        [inf] = r_rnxheadv2(f_obs)
        [dcb] = r_dcb(f_dcb)
        [freq,wavl] = frequencies
        [doy] = clc_doy(year,mon,day)
        [jd,mjd] = cal2jul(year,mon,day,sec)
        [leap_sec] = leapsec(mjd)
        
        function order = get_sys_order(system)
            if strcmpi(system, 'G') | strcmpi(system, 'GPS')
                order = 1;
            elseif strcmpi(system, 'R') | strcmpi(system, 'GLONASS')
                order = 2;
            elseif strcmpi(system, 'E') | strcmpi(system, 'GALILEO')
                order = 3;
            elseif strcmpi(system, 'C') | strcmpi(system, 'BEIDOU')
                order = 4;
            end
        end
        
        function max_sv_index = max_sv_index(system)
            max_sv_index = rinex.max_sys_svs * rinex.get_sys_order(system);
        end
        
        function min_sv_index = min_sv_index(system)
            min_sv_index = rinex.max_sys_svs * (rinex.get_sys_order(system) - 1) + 1;
        end
        
        function mi = max_sys_index()
            % max supported systems
            mi = rinex.max_sys_svs * 4;
        end
    end
end


classdef tec
    %TEC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stnm            % station name
        leap_sec        % difference between GPS time and UTC
        slant_tec       % slant total electron content
        vertical_tec    % vertical total electron content
        ionex_tec       % save the ionex tec samples if they exist
        ep              % epoch vector
        elv             % elevation vector
        azm             % azimuth vector
        ipp             % ionospheric piercing points vector
        st              % state vector (boolean)
        sv_names        % space vehicles names
        ifb             % interfrequency bias (GNSS)
        first_obs       % first observation datetime
        ref_time        % either UTC or GPS
        info            % location of the station, etc
    end
    
    methods
        function self = tec(stnm, observables, info, orbits, ionex, options)
            %TEC Construct an instance of this class
            %   Detailed explanation goes here
            
            self.stnm = stnm;
            self.leap_sec = info.time.leap;
            self.ref_time = options.output_time;
            
            %options.elevation_mask = 10;
            options.CSMw = 1;
            options.CSGf = 1;
            %options.ionosphere_height = 200e3;
            %options.emf_method = 2;  % default method
            
            self.first_obs = datetime(info.time.first');
            
            % begin by calculating the elevation and azimuth
            [data, observables] = tec.elv_mask(observables, info, orbits.interpolate(observables.ep).sp3, options);
            
            % continue fixing cycle slips
            observables = tec.cs_detect(data, observables, info, options);
            data = tec.clk_jmp2(observables, data);
            
            % now compute the slant tec values
            self.slant_tec = tec.get_ion(observables, info, data, ionex);
            
            % propagate a vector with the epochs and observation flags
            self.ep = observables.ep; self.st = ~isnan(self.slant_tec);
            
            % compute the ionospheric piercing points and propagate elv azm
            self.elv = data.elv;
            self.azm = data.azm;
            self.ipp = tec.get_sip(orbits.interpolate(observables.ep).sp3, info.rec.pos, self.st, options);
            
            % now that the ipp are available, remove the interfrequency
            % bias
            [self.slant_tec, self.ifb, self.ionex_tec] = self.compute_ifb(info, ionex, options);
            
            % compute vtec
            self.vertical_tec = self.slant_tec .* tec.get_emf(data.elv, options);
            
            % propagate the names of the space vehicles
            self.sv_names = tec.get_sv_names();
            
            % information
            self.info.xyz = info.rec.pos;
            self.info.lla = ecef2lla(info.rec.pos);
            self.info.ionosphere_height = options.ionosphere_height;
        end
        
        function plot_vtec(self, varargin)
            %PLOT_VTEC plot vertical tec

            p = inputParser;
            p.addParameter('startDate', min(self.ep)./3600)
            p.addParameter('endDate', max(self.ep)./3600)
            p.addParameter('plotIonex', false)
            p.addParameter('svs', 1:105)

            parse(p, varargin{:});

            s = p.Results.startDate; e = p.Results.endDate;

            ft = self.ep >= s * 3600 & self.ep <= e * 3600;
            clf
            plot(self.ep(ft) ./ 3600, self.vertical_tec(ft, p.Results.svs) ./ 1e16)
            xlabel('Time since start [hours]')
            ylabel('VTEC [TECU]')
            
            f = sum(self.st(ft, :)) > 0;
            ff = false(size(f));
            ff(p.Results.svs) = true;
            f = ff & f;
            legend(self.sv_names(f), 'location', 'northeastoutside')
            grid on
            title(['VTEC for station ' upper(self.stnm)])
            
            if p.Results.plotIonex
                hold on
                set(gca,'ColorOrderIndex',1)
                plot(self.ep(ft) ./ 3600, self.ionex_tec(ft, p.Results.svs) ./ 1e16, 'o')
            end
        end
        
        function plot_ipp(self, varargin)
            
            p = inputParser;
            p.addParameter('startDate', min(self.ep))
            p.addParameter('endDate', max(self.ep))
            p.addParameter('plotIonex', false)
            p.addParameter('svs', 1:105)
            p.addParameter('elev', false)

            parse(p, varargin{:});

            s = p.Results.startDate; e = p.Results.endDate;

            ft = self.ep >= s * 3600 & self.ep <= e * 3600;
            
            clf
            lat = squeeze(self.ipp(ft, 1, p.Results.svs));
            lon = squeeze(self.ipp(ft, 2, p.Results.svs));
            
            m_proj('mercator', 'lat', [min(min(lat))-15 max(max(lat))+15], 'lon', [min(min(lon))-15 max(max(lon))+15]);
            
            m_coast('patch',[204/255 255/255 51/255], 'edgecolor','k');
            m_grid('xaxislocation','top');
            
            [x, y] = m_ll2xy(lon, lat);
            x = x(:);
            y = y(:);
            if p.Results.elev
                t = self.elv(ft, p.Results.svs);
            else
                t = self.vertical_tec(ft, p.Results.svs) .* 1e-16;
            end
            hold on
            scatter(x, y, 10, t(:), 'filled');
            colorbar
            colormap jet
        end
        
        function save(self, filename, format)
            % save the processing into a structure
            data.vertical_tec = self.vertical_tec;
            data.ipp = self.ipp;
            data.station_name = self.stnm;
            data.ref_time = self.ref_time;
            data.sv_names = self.sv_names;
            data.elv = self.elv;
            data.azm = self.azm;
            data.first_obs = self.first_obs;
            
            data.info = self.info;
            
            if strcmpi(self.ref_time, 'utc')
                data.ep = self.ep - self.leap_sec;
            else
                data.ep = self.ep;
            end
            
            f = fileparts(filename);
            
            if ~exist(f, 'dir')
                mkdir(f);
            end
            
            % more option to come!
            if strcmpi(format, 'mat')
                save([filename '.mat'], 'data');
            end
        end
        
        % function definition
        [iec, ifb, ionex_tec] = compute_ifb(self, info, ionex, options)
    end
    
    methods (Static)
        [observables] = fix_cycle_slips(observables, info)
        [segments] = get_segments(ep, data, splint, tolerance)
        [p] = fill_data_gaps(t, o)
        [s] = lsc_fit(t, o)
        [oc] = clean_tec_breaks(t, o, int)
        [iec] = get_ion(observables, info, data, ionex)
        [observables] = cs_detect(data, observables, info, options)
        [arc] = arc_dtr(obs)
        [data, observables] = elv_mask(observables, info, orbits, options)
        %[azim,elev] = local(rec,sat,dopt)
        %[elip] = xyz2plh(cart,dopt)
        [data] = clk_jmp2(observables, data)
        [ifr] = i_free(o1,o2,opt)
        [E] = get_emf(el, options)
        [stec, ti] = sample_tec_map(ep, time, ionex, ipp, elv, options)
        [ipp] = get_sip(orbits, pos, st, options)
        [AZ,EL,LE] = azelle(S,P);
        
        function sv_names = get_sv_names()
            for i=1:rinex.max_sys_index
                if i<rinex.max_sv_index('G')     % GPS
                    sv_names{i} = 'G';
                elseif i<rinex.max_sv_index('R') % GLONASS
                    sv_names{i} = 'R';
                elseif i<rinex.max_sv_index('E') % GALILEO
                    sv_names{i} = 'E';
                else        % BEIDOU
                    sv_names{i} = 'C';
                end
                sv_names{i} = [sv_names{i} sprintf('%02i', i)];
            end
        end
    end
end


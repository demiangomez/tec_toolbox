classdef station < handle
    %STATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        rinex
        sv
        mean_epoch
    end
    
    methods
        function self = station(rnx_files, options)
            %STATIONS Construct an instance of this class
            %   Detailed explanation goes here
            
            % options to read files
            options.dcb = false;
            options.clck_int = true;
            
            % load the RINEX file
            self.rinex = rinex(rnx_files, options);
            
            % get the name of the station
            if isa(rnx_files, 'cell')
                [~, rnx_name, ~] = fileparts(rnx_files{1});
            else
                [~, rnx_name, ~] = fileparts(rnx_files);
            end
            
            if self.rinex.info.rinex.ver == 2
                self.name = rnx_name(1:4);
            else
                self.name = rnx_name(1:9);
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end


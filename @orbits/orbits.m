classdef orbits < handle
    %ORBITS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sat
        info
    end
    
    methods
        function self = orbits(sp3file)
            %ORBITS Construct an instance of this class
            %   Detailed explanation goes here
            
            options.system.gps = true;
            options.system.glo = true;
            options.system.gal = true;
            options.system.bds = false;
            
            self.info = [];
            
            [self.sat, self.info] = orbits.r_sp3(sp3file, options, self.info);
        end
        
        function orbits = interpolate(self, epochs, varnargin)
            
            orbits = struct();
            Ti = epochs;
            
            if nargin > 2 && ~isempty(varnargin)
                % erquesting a specific set of svs
                for i = 1:length(nargin)
                    T = self.sp3(:, 5, nargin(i));
                    x = self.sp3(:, 1, nargin(i));
                    y = self.sp3(:, 2, nargin(i));
                    z = self.sp3(:, 3, nargin(i));
                    
                    X = interp1(T, x, Ti,'*cubic');
                    Y = interp1(T, y, Ti,'*cubic');
                    Z = interp1(T, z, Ti,'*cubic');
                    
                    orbits.sp3(:,:,i) = [X Y Z zeros(size(X)) Ti];
                end
            else
                % wants all svs
                for i = 1:size(self.sat.sp3, 3)
                    T = self.sat.sp3(:, 5, i);
                    x = self.sat.sp3(:, 1, i);
                    y = self.sat.sp3(:, 2, i);
                    z = self.sat.sp3(:, 3, i);
                    
                    X = interp1(T, x, Ti,'*cubic');
                    Y = interp1(T, y, Ti,'*cubic');
                    Z = interp1(T, z, Ti,'*cubic');
                    
                    orbits.sp3(:,:,i) = [X Y Z zeros(size(X)) Ti];
                end
            end
        end
    end
    
    methods (Static)
        [sat, inf] = r_sp3(f_orbit,options,inf)    
    end
end


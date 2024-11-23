classdef ionex
    %IONEX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tgd
        iono_map
        dcb
    end
    
    methods
        function self = ionex(file)
            %IONEX Construct an instance of this class
            %   Detailed explanation goes here
            [epoch,time,tec,~,lats,lons,hgts,~,dcbs,dcbr,rcvs] = ionex.readionex(file);
            
            self.iono_map.lats = lats;
            self.iono_map.lons = lons;
            self.iono_map.h = hgts;
            self.iono_map.tec = tec;
            self.iono_map.time = time;
            self.iono_map.epoch = epoch;
            self.tgd = dcbs;
            self.dcb.receivers = rcvs;
            self.dcb.values = dcbr;
        end
    end
    
    methods (Static)
        [epoch,time,tec,rms,lats,lons,hgts,rb,dcbs,dcbr,sats,rcvs] = readionex(file)
    end
end


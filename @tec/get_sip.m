function ipp = get_sip(orbits, pos, st, options)

% GET_SIP	Calculate sub-ionospheric points + associated information
%		SAT = structure containing sat positions (created by int_sp3)
%		sv = vector of PRN numbers (created by readrinex)
%		alt_ion = ionospheric height
%		apcoords = ECEF coordinates of GPS site (XYZ, m) (created by readrinex)
%		Tmin = window start time (hours UT)
%		Tmax = window end time (hours UT)
%		lat_epi, lon_epi = location of event
%		dec, inc = declination and inclination of magnetic field (degres)
%
%		SIP = structure, for each field (= PRN):
%		  time in UT hour of current day
%		  latitude of SIPs in degrees
%		  longitude of SIPs in degrees
%		  elevation w.r.t GPS site in degrees
%		  azimuth w.r.t GPS site in degrees
%		  distance from event to SIP (m)
%		  azimuth from event to SIP (deg CW from N)
%		  angle between magnetic direction and sat-rec line-of-sight (at the rec)
%		
%		SIP = get_sip(SAT,sv,alt_ion,apcoords,Tmin,Tmax,lat_epi,lon_epi,dec,inc);
%

% define some constants (should be read from a file)
R = 6371000.0;      % mean Earth radius

% initialize SIP structure
ipp = nan([size(orbits, 1) 2 size(orbits, 3)]);

% compute SIP for each satellite
for i = 1:size(orbits, 3)
   if any(st(:,i))
      
      % matrix of satellite positions (ECEF, meters)
      S = orbits(:,1:3, i);

      Rgs = S - pos;
      
      rang = sqrt(Rgs(:,1).^2+Rgs(:,2).^2+Rgs(:,3).^2);
      
      u = [Rgs(:,1)./rang Rgs(:,2)./rang Rgs(:,3)./rang];
      
      % u is the direction vector of the line
      
      % solve the non-right triangle formed by the piercing point, the
      % station and the geocenter. The distance we want to solve for is l
      % (the distance between the station and the piercing point).
      % the equation has the form 
      % (R+I)^2 = l^2 + |apcoords|^2 - 2*l*|apcoords|*cos(elev + 90)
      % writing the equation in terms for the dot product produces
      % 2*l*|apcoords|*cos(elev + 90) = 2*l*dot(u, E)
      b = 2.*(pos(1).*u(:,1) + pos(2).*u(:,2) + pos(3).*u(:,3));
      % a is 1 
      c = (pos(1)^2 + pos(2)^2 + pos(3)^2 - (R+options.ionosphere_height)^2);
      l = (-b+sqrt(b.^2-4*c))/2;
      
      X = pos(1) + l.*u(:,1);
      Y = pos(2) + l.*u(:,2);
      Z = pos(3) + l.*u(:,3);
      
      lla_sip = ecef2lla([X Y Z]);
      
      ipp(:, :, i) = lla_sip(:, [1 2]);
      ipp(~st(:,i), [1 2], i) = NaN;

   end
end


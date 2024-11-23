function [AZ,EL,LE] = azelle(S,P)

    % AZELLE	Computes elevation angle, range, and azimuth
    %		of satellites from a ground station
    %		S = ECEF satellite coordinates (m), n x 3 matrix
    %		P = ground station coordinates, ECEF (m) (vector)
    %		AZ = azimuth (radians CW from north)
    %		EL = elevation angle (radians)
    %		LE = range (meters)
    %
    %		[AZ,EL,LE] = azelle(S,P)
    %

    % satellite ECEF coordinates (must be in meters)
    Xs = S(:,1); Ys = S(:,2); Zs = S(:,3);

    % ground station ECEF coordinates (must be in meters)
    XR = P(1); YR = P(2); ZR = P(3);

    % compute ground-sat vector in ECEF coordinates
    Rgs = [Xs-XR Ys-YR Zs-ZR];

    % convert to unit vector
    rang = sqrt(Rgs(:,1).^2+Rgs(:,2).^2+Rgs(:,3).^2);
    Ru = [Rgs(:,1)./rang Rgs(:,2)./rang Rgs(:,3)./rang];

    % rotate from XYZ to NEU
    neu = xyz2neu([XR YR ZR],Ru);

    % convert neu to azimuth and elevation angle
    LE = sqrt(neu(:,1).^2+neu(:,2).^2);
    EL = (pi/2) - atan2(LE,neu(:,3));
    AZ = atan2(neu(:,2),neu(:,1));
    % DDG: problem: the azimuth went from 0 to 180 and from -180 to 0. Wrong,
    % it has to go from 0 to 360!
    AZ(AZ < 0) = AZ(AZ < 0) + 2*pi;
end
    
function NEU = xyz2neu(O,V)

    % XYZ2NEU	Convert ECEF into local topocentric
    %		Call: NEU = xyz2neu(O,V)
    %		O = origin vector in ECEF frame (meters), or n x 3 matrix
    %		V = position or velocity vector relative to origin point in ECEF frame (m or m/yr), or n x 3 matrix
    %		Note that O can be the same as V
    %		NEU  = position or velocity vector in NEU frame (m or m/yr), or n x 3 matrix

    % initialize variables
    X = V(:,1); Y = V(:,2); Z = V(:,3);

    % if O is a single point, make it the same size as V
    if (size(O,1) == 1)
      XR = ones(size(V,1),1) .* O(1);  
      YR = ones(size(V,1),1) .* O(2);  
      ZR = ones(size(V,1),1) .* O(3);  
    else
      XR = O(:,1); YR = O(:,2); ZR = O(:,3);
    end
    T = zeros(size(XR,1),1);

    % convert origin vector to ellipsoidal coordinates then to in radians
    E = xyz2wgs([T XR YR ZR]);
    E(:,2) = E(:,2).*pi/180; % longitude
    E(:,3) = E(:,3).*pi/180; % latitude
    
    [NEU(:, 1), NEU(:, 2), NEU(:, 3)] = ct2lg(X, Y, Z, E(:,3), E(:,2));
    
end

function R = xyz2wgs(S)

    % XYZ2WGS      	converts cartesian coordinates (x,y,z) into
    %		ellipsoidal coordinates (lat,lon,alt) on WGS-84
    %		according to a non iterative method (Bowring 76,
    %		see also GPS Theory and Practice, p. 258).
    %               Call: R = xyz2wgs(S)
    %               S is a nx4 matrix with time, X, Y, Z
    %		A! first column of S is time but can be dummy.
    %               R is a nx4 matrix with time, lon (lam), lat (phi), alt
    %		                             (lon,lat in degrees!)

    % WGS-84 PARAMETERS
    % semimajor and semiminor axis
    a = 6378137.0;
    b = 6356752.314;
    % flattening
    f = 1.0/298.257222101;
    % eccentricity
    eo = 2*f - f^2;

    % second numerical eccentricity
    e1 = (a^2-b^2)/b^2;

    % read data
    t = S(:,1);
    x = S(:,2);
    y = S(:,3);
    z = S(:,4);

    % auxiliary quantities
    p = sqrt(x.^2+y.^2);
    theta = atan2(z.*a,p.*b);

    % latitude, longitude
    phi = atan2 (z + (sin(theta)).^3.*e1.*b , p - (cos(theta)).^3.*eo^2.*a);
    lam = atan2 (y,x);

    % radius of curvature in prime vertical
    N = a^2 ./ sqrt((cos(phi)).^2.*a^2 + (sin(phi)).^2.*b^2);

    % altitude
    alt = (p ./ cos(phi)) - N;

    % fill out result matrix
    R = [t lam.*180.0./pi phi.*180.0./pi alt];
end
function [x,y,z] = latlongalt_xyz(lon,lat,alt)
%LATLONGALT_XYZ convert lon lat alt cordinates to xyz
%   for converting from latitude, longitude, altitude (above reference ellipsoid) to Cartesian ECEF


format long
p = pi;

cosLat = cos(lat*p/180.0);
sinLat = sin(lat*p/180.0);
cosLon = cos(lon*p/180.0);
sinLon = sin(lon*p/180.0);

h = alt;
rad = 6378137.0; %radius of earth 
f = 1/298.257224; 
C = 1 / sqrt( cosLat^2 + (1-f)^2 * sinLat^2 );
S = (1-f)^2 * C;

x = (rad * C + h) * cosLat * cosLon;
y = (rad * C + h) * cosLat * sinLon;
z = (rad * S + h) * sinLat;



end


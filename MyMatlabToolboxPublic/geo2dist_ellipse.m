function [d] = geo2dist_ellipse(lat1,lon1,lat2,lon2)
% given the lat and long of two points on earth in degrees, output distance in km
%R=6371; %earth mean radius
%lat1=lat1*pi/180; lon1=lon1*pi/180;
%lat2=lat2*pi/180; lon2=lon2*pi/180;
%d=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))*R;

[d,az]=distance(lat1,lon1,lat2,lon2,[6378.14 0.0818]);



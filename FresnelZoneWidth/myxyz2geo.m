function [lon,lat,depth]=myxyz2geo(lon0,lat0,depth0,x,y,z)
% Transfer (x,y,z) to (lat,lon,depth) with given origin point.
% Converted from FORTRAN 90: subroutine transXYZ2GEO(n0,a0,d0,x,y,z,n1,a1,d1)
% Originally written by Xiaotao Yang on December 2010 @ IGGCAS.
% Converted on 2016.1.27
%Xiaotao Yang @ Indiana University
% NOTE: vertical down is positive for z.

R0=6371.0;

  lon = rad2deg(radians(lon0) + 180.0*(x/((R0 - depth0)*sin(lat00)))/pi;
  
  lat = lat0 + 180.0*(y/(R0 - depth0))/pi;
  
  depth = depth0 + z;
  
end
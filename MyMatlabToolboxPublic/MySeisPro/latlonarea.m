function A=latlonarea(lat1,lon1,lat2,lon2)
% Compute the area in km^2 for the rectangular region between two latitudes
% and two longitudes.
% Use the following formula:
%  A = 2*pi*R^2 |sin(lat1)-sin(lat2)| |lon1-lon2|/360
%     = (pi/180)R^2 |sin(lat1)-sin(lat2)| |lon1-lon2|
%
% From:http://mathforum.org/library/drmath/view/63767.html
% Accessed on: Oct 26, 2012
% Xiaotao Yang @ Indiana University
R=6371;    %radius of the earth.

A=(pi/180)*R^2*abs(sin(lat1*pi/180)-sin(lat2*pi/180))*abs(lon1 - lon2);

return;

end
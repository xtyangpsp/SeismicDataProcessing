function [lon lat depth mag]=select_by_region(lon0,lat0,depth0,mag0,lon_min,lon_max,lat_min,lat_max)
% select earthquake data by given region.
%lon0: original data;
%lon_min/max: region;
%lon lat: output;
j=1;
for i = 1:length(lon0)
if lon0(i) >= lon_min && lon0(i) <= lon_max && lat0(i) >= lat_min && lat0(i) <= lat_max
    lon(j)=lon0(i);
    lat(j)=lat0(i);
    depth(j)=depth0(i);
    mag(j)=mag0(i);
    j = j + 1;
end
end
return;
end
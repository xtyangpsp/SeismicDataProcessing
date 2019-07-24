function dist=geodist(lat1,lon1,lat2,lon2)
% get great circle distance in km or degree between two points on the earth's
% surface.
% flag: 'km' for km, 'deg' for degrees.
  R0=6371.0;
  pi=3.141592653;
  
  dist=sqrt(((lon1 - lon2)*pi*R0/180).^2 + ((lat1 - lat2)*pi*R0/180).^2);
  
%   if strcmp(flag,'deg')
%      dist=180*dist/(R0*pi);
%   end
  return
  
end
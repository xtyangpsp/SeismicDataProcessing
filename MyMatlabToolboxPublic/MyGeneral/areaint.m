function a=areaint(lats,longs,in3,in4) 
 
%AREAINT  Calculates the spherical surface area of a polygon 
% 
%  a = AREAINT(lat,long) calculates the spherical surface area 
%  of the polygon specified by the input vectors lat, long.  The 
%  calculation uses a line integral approach.  The output, a, 
%  is the surface area fraction covered by the polygon on a 
%  unit sphere.  Multiple polygons can be supplied provided that 
%  each polygon is separated by a NaN in the input vectors. 
%  Accuracy of the integration method is inversely proportional 
%  to the distance between lat/long points. 
% 
%  a = AREAINT(lat,long,geoid) uses the input geoid to describe 
%  the sphere or ellipsoid.  The output, a, is in square units 
%  corresponding to the units of geoid(1).  If omitted, geoid = [1 0] 
%  is assumed and the output is in the fractional form described above. 
% 
%  a = AREAINT(lat,long,'units') uses the units defined by the 
%  input string 'units'.  If omitted, default units of degrees 
%  is assumed. 
% 
%  a = AREAINT(lat,long,geoid,'units') uses both the inputs geoid 
%  and 'units' in the calculation. 
% 
%  See also AREAMAT, AREAQUAD. 
 
%  Copyright 1996-2002 Systems Planning and Analysis, Inc. and The MathWorks, Inc. 
%  Written by:  E. Brown, E. Byrns 
%   $Revision: 1.13 $    $Date: 2002/03/20 21:24:50 $ 
 
if nargin < 2 
    error('Incorrect number of arguments') 
 
elseif nargin==2 
	units = [];     geoid = []; 
 
elseif nargin==3 
	if isstr(in3) 
		units = in3;     geoid = []; 
	else 
		units = [];      geoid = in3; 
	end 
 
elseif nargin==4 
		geoid=in3;   units=in4; 
end 
 
%  Empty argument tests 
 
if isempty(units);   units  = 'degrees';   end 
 
absolute_units = 1;          %  Report answer in absolute units assuming 
if isempty(geoid)            %  a radius input has been supplied.  Otherwise, 
     geoid = [1 0];          %  report surface area answer as a fraction 
	 absolute_units = 0;     %  of a sphere 
end 
 
%  Test the geoid parameter 
 
[geoid,msg] = geoidtst(geoid); 
if ~isempty(msg);   error(msg);   return;   end 
 
%  Input dimension tests 
 
if ~isequal(size(lats),size(longs)) 
	error('Lats and longs must be the same size') 
 
elseif	min(size(lats)) ~= 1 
	error('Lats and longs must be vectors') 
end 
 
%  Enforce that lat and longs are column vectors 
 
lats = lats(:);    longs = longs(:); 
 
%  Convert coordinates to radians if necessary 
%  Transform to the authalic sphere 
 
lats  = angledim(lats,units,'radians'); 
lats  = geod2aut(lats,geoid,'radians'); 
longs = angledim(longs,units,'radians'); 
radius = rsphere('authalic',geoid); 
 
%  Ensure at a terminating NaN in the vectors 
 
if ~isnan( lats(length(lats)) );    lats = [lats; NaN];   end 
if ~isnan( longs(length(longs)) );    longs = [longs; NaN];   end 
 
%  Ensure vectors don't begin with NaNs 
 
if isnan(lats(1)) | isnan(longs(1)) 
	lats = lats(2:length(lats)); 
	longs = longs(2:length(longs)); 
end 
 
%  Find segment demarcations 
 
indx=find(isnan(lats)); 
 
%  Initialize area output vector 
 
a(length(indx),1) = 0; 
 
%  Perform area calculation for each segment 
 
for i=1:length(indx) 
 
	% Pull segment out of main vectors 
 
	if i>1 
		lat=lats(indx(i-1)+1:indx(i)-1); 
		long=longs(indx(i-1)+1:indx(i)-1); 
	else 
		lat=lats(1:indx(i)-1); 
		long=longs(1:indx(i)-1); 
	end 
 
	%  Ensure the path segment closes upon itself by 
	%  repeating beginning point at the end. 
 
	lat = [lat; lat(1)];   long = [long; long(1)]; 
 
	% Set origin for integration.  Any point will do, so (0,0) used 
 
	zerovec = zeros( size(lat) ); 
 
	% Get colatitude (a measure of surface distance as an angle) 
	% and azimuth of each point in segment from the arbitrary origin 
 
	colat = distance('gc',zerovec,zerovec,lat,long,'radians'); 
	az    = azimuth('gc',zerovec,zerovec,lat,long,'radians'); 
 
	% Calcluate step sizes 
 
	daz=diff(az); 
	daz=npi2pi(daz,'radians','exact'); 
 
	% Determine average surface distance for each step 
 
	deltas=diff(colat)/2; 
	colats=colat(1:length(colat)-1)+deltas; 
 
	% Integral over azimuth is 1-cos(colatitudes) 
 
	integrands=(1-cos(colats)).*daz; 
 
	% Integrate and save the answer as a fraction of the unit sphere. 
	% Note that the sum of the integrands will include a factor of 4pi. 
 
	a(i) = abs(sum(integrands))/(4*pi); 
 
end 
 
%  Convert to absolute terms if the default radius was not used 
 
if absolute_units;   a = a * 4*pi*radius^2;   end 
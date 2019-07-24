%BWR Blue-White-Red colormap, mostly used for highly oscillatory data
%  bwr(M), a variant of  JET(M), is an M-by-3 matrix containing
%    the default colormap used by CONTOUR, SURF, PCOLOR and IMAGE.
%    The colors begin with dark blue, range through white and end with 
%    dark red.

% by Rodrigo S. Portugal (rosoport.matlab@gmail.com)
% last revision: 19/10/2012

function rgb = bwr(n)

if nargin < 1
   n = size(get(gcf,'colormap'),1);
end

k0 = round([1, n/4, n/2, 3*n/4, n]);

r0 = [  0,   0, 255, 212, 128] / 255;
g0 = [  0,  96, 255,  96,   0] / 255;
b0 = [128, 212, 255,   0,   0] / 255;

k = 1:n;
r = interp1(k0, r0, k);
g = interp1(k0, g0, k);
b = interp1(k0, b0, k);

rgb = [r',g',b'];
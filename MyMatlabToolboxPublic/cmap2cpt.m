function cmap2cpt(cmap,outcptfile,zrange,nancolor,isdemcmap)
% Convert MATLAB colormap data to GMT cpt file. The default z value is
% 0-100. The resultant cpt file can be used by makecpt to customize for any
% needs.
% USAGE: cmap2cpt(@cmap,outcptfile,zrange): 
%       cmap - MATLAB colormap function, e.g., jet, perula, hot, cool, etc;
%               The symbol '@' is required!
%       outcptfile - filename string of the output cpt data.
%       zrange (OPTIONAL) - a three-element vector specifying [zmin zmax
%               zinc], default [0 100 1]
%       nancolor - Set NaN color (optional) in the range of 0 - 255, default
%               is [128 128 128];
%
% By Xiaotao Yang @ UMass Amherst
%  Contact: stcyang@gmail.com
%
%   created: 2/3/2018
%   modified on 9/6/2018: added option to set NaN color
%
if nargin < 2
    error('Not enough input arguments!')
elseif nargin ==2
    zrange=[0 100 1];
    nancolor=[128 128 128];
    isdemcmap=0;
elseif nargin ==3
    nancolor=[128 128 128];
    isdemcmap=0;
elseif nargin ==4
    isdemcmap=0;
end
if  isempty(nancolor);nancolor=[128 128 128];end

if strcmp(zrange,'-');zrange=[0 100 1];end
zval=zrange(1):zrange(3):zrange(2);
if isdemcmap
    cdata=cmap([zrange(1) zrange(2)],length(zval)); %this is used only for
%converting demcmap to cpt.
else
    cdata=cmap(length(zval));
end

cdata=round(cdata*255); %convert to scale of 0-255 for GMT convention.
%cptval format
%z1 R1 G1 B1 z2 R2 G2 B2
%...
% B       R     G       B   #FOR BOTTOM VALUES
% F       R     G       B   #FOR FULL/TOP VALUES
% N       128   128     128   #FOR NULL VALUES (NAN).

cptval=nan(length(zval)-1,8);
for i=1:length(zval)-1
    cptval(i,:)=[zval(i) cdata(i,:) zval(i+1) cdata(i+1,:)];
end
fhd=functions(cmap);
fidout=fopen(outcptfile,'w');
fprintf(fidout,'# cpt converted by cmap2cpt() from MATLAB: %s\n',fhd.function);
fprintf(fidout,'# cmap2cpt() is wrote by Xiaotao Yang @ UMass Amherst\n');
fprintf(fidout,'# contact: stcyang@gmail.com\n');
fprintf(fidout,'#COLOR_MODEL = RGB\n');
for j=1:size(cptval,1)
   fprintf(fidout,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',cptval(j,:)); 
   
end
fprintf(fidout,'B\t%d\t%d\t%d\n',cptval(1,2:4));
fprintf(fidout,'F\t%d\t%d\t%d\n',cptval(end,6:8));
fprintf(fidout,'N\t%d\t%d\t%d\n',nancolor);

fclose(fidout);
end
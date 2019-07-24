function [dipdirA1,dipA1,rake1,dipdirA2,dipA2,rake2] = Stereonet_FoldAxes_Selection;

% Degree to radians and vice versa factors and square root of 2.0
deg = 180/pi;   rad = pi/180;   s2 = sqrt(2);

% Vectors for plotting great circles
N = linspace(0,pi,181);
gx = zeros(size(N));    gy = zeros(size(N));

%==========================================================================
% GET FOLD AXES ORIENTATIONS 

% First Fold Axis

X1 = 10; Y1 = 10;
while sqrt(X1^2 + Y1^2) > 1.0
    [X1,Y1] = ginput(1);
end
plot(X1,Y1,'ro','MarkerFaceColor','r');
[azimuth1,r1] = cart2pol(X1,Y1);
plunge1 = (pi/2.0 - 2.0*asin(r1/s2))*deg;
plungedir1 = 90.0-azimuth1*deg;
if plungedir1 < 0.0
    plungedir1 = 360.0 + plungedir1;
end
text(X1+0.05,Y1,[' b_1 = ',num2str(round(plungedir1)),' / ', num2str(round(plunge1))],'BackgroundColor',[1 1 1]);

% First Fold Axial Plane

XA1 = 10; YA1 = 10;
while sqrt(XA1^2 + YA1^2) > 1.0
    [XA1,YA1] = ginput(1);
end
plot(XA1,YA1,'ro','MarkerFaceColor','w');

[azimuthA1,rA1] = cart2pol(XA1,YA1);
plungeA1 = (pi/2.0 - 2.0*asin(rA1/s2))*deg;
plungedirA1 = 90.0-azimuthA1*deg;
if plungedirA1 < 0.0
    plungedirA1 = 360.0 + plungedirA1;
end

v1  = [cos(plungedir1*rad); -sin(plungedir1*rad); -tan(plunge1*rad)];
vA1 = [cos(plungedirA1*rad); -sin(plungedirA1*rad); -tan(plungeA1*rad)];

c1 = cross(v1,vA1);
m1 = sqrt(sum(c1.^2));
n1 = c1/m1;

if n1(3) > 0.0
    n1 = -n1;
end

dipA1    = 90.0 - asin(-n1(3))*deg;
dipdirA1 = 90.0 + atan2(n1(1),n1(2))*deg;

if dipdirA1 < 0.0
    dipdirA1 = 360.0 + dipdirA1;
end

great_c = N - dipdirA1*rad;
apparent_dip = atan(tan(dipA1*rad).*sin(N));
r = s2.*sin(pi/4.0 - apparent_dip/2.0);
[gx,gy] = pol2cart(great_c,r);
plot(gx,gy,'-r','LineWidth',2);
i = floor(length(N)/2.0);

rake1 = sign(dipdirA1-plungedir1)*(90-acos(sin(plunge1*rad)/sin(dipA1*rad))*deg);

if dipdirA1 > 270 & plungedir1 < 90 |...
        plungedir1 > 270 & dipdirA1 < 90
    rake1 = -1*rake1;
end
rake1
text(gx(i)+0.05,gy(i),[' FAP_1 = ',num2str(round(dipdirA1)),' / ', num2str(round(dipA1))],'BackgroundColor',[1 1 1]);

% Second Fold Axis

X2 = 10; Y2 = 10;
while sqrt(X2^2 + Y2^2) > 1.0
    [X2,Y2] = ginput(1);
end
plot(X2,Y2,'bo','MarkerFaceColor','b');
[azimuth2,r2] = cart2pol(X2,Y2);
plunge2 = (pi/2.0 - 2.0*asin(r2/s2))*deg;
plungedir2 = 90.0-azimuth2*deg;
if plungedir2 < 0.0
    plungedir2 = 360.0 + plungedir2;
end
text(X2+0.05,Y2,[' b_2 = ',num2str(round(plungedir2)),' / ', num2str(round(plunge2))],'BackgroundColor',[1 1 1]);

% Second Fold Axial Plane

XA2 = 10; YA2 = 10;
while sqrt(XA2^2 + YA2^2) > 1.0
    [XA2,YA2] = ginput(1);
end
plot(XA2,YA2,'bo','MarkerFaceColor','w');

[azimuthA2,rA2] = cart2pol(XA2,YA2);
plungeA2 = (pi/2.0 - 2.0*asin(rA2/s2))*deg;
plungedirA2 = 90.0-azimuthA2*deg;
if plungedirA2 < 0.0
    plungedirA2 = 360.0 + plungedirA2;
end

v2  = [cos(plungedir2*rad); -sin(plungedir2*rad); -tan(plunge2*rad)];
vA2 = [cos(plungedirA2*rad); -sin(plungedirA2*rad); -tan(plungeA2*rad)];

c2 = cross(v2,vA2);
m2 = sqrt(sum(c2.^2));
n2 = c2/m2;

if n2(3) > 0.0
    n2 = -n2;
end

dipA2    = 90.0 - asin(-n2(3))*deg;
dipdirA2 = 90.0 + atan2(n2(1),n2(2))*deg;

if dipdirA2 < 0.0
    dipdirA2 = 360.0 + dipdirA2;
end

great_c = N - dipdirA2*rad;
apparent_dip = atan(tan(dipA2*rad).*sin(N));
r = s2.*sin(pi/4.0 - apparent_dip/2.0);
[gx,gy] = pol2cart(great_c,r);
plot(gx,gy,'-b','LineWidth',2);
i = floor(length(N)/2.0);

rake2 = sign(dipdirA2-plungedir2)*(90-acos(sin(plunge2*rad)/sin(dipA2*rad))*deg);

if dipdirA2 > 270 & plungedir2 < 90 |...
        plungedir2 > 270 & dipdirA2 < 90
    rake2 = -1*rake2;
end

text(gx(i)+0.05,gy(i),[' FAP_2 = ',num2str(round(dipdirA2)),' / ', num2str(round(dipA2))],'BackgroundColor',[1 1 1]);
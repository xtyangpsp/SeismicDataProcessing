function mystereoplot(dipdirA1,dipA1,dipdirA2,dipA2)

% This function plots the orienations of two folds defined by fold axial
% plane orientation (dip direction and dip) and rake of fold axis in a
% stereonet with radius 1.

%==========================================================================
% Ensure no division by zero occurs
if dipA1  == 0
    dipA1 = 1e-6;
end
if dipA2  == 0
    dipA2 = 1e-6;
end
pi=3.141592653;
% Degree to radians and vice versa factors and square root of 2.0
deg = 180/pi;   rad = pi/180;   s2 = sqrt(2);

% Vectors for plotting great circles
N = linspace(0,180,181);
gx = zeros(size(N));    gy = zeros(size(N));

%==========================================================================

% First Fold Axial Plane
great_c = (N - dipdirA1)*rad;
apparent_dip = atan(tan(dipA1*rad).*sin(N*rad));
r = s2.*sin(pi/4.0 - apparent_dip/2.0);
[gx,gy] = pol2cart(great_c,r);
plot(gx,gy,'-k','LineWidth',1);
%i = floor(length(N)/2.0);
% text(gx(i)+0.05,gy(i),[' FAP_1 = ',num2str(round(dipdirA1)),' / ', num2str(round(dipA1))],'BackgroundColor',[1 1 1]);

% First Fold Axis
% plunge1 =  abs(asin(cos(rake1*rad-pi/2)*sin(dipA1*rad))*deg);
% pdd_diff = acos(tan(plunge1*rad)/tan(dipA1*rad))*deg;
% sgn = sign(rake1);
% if sgn == 0;
%     sgn = 1;
% end
% plungedir1 = dipdirA1 - sgn*pdd_diff;
% 
% if plungedir1 < 0.0
%     plungedir1 = 360.0 + plungedir1;
% end
% 
% r = s2.*sin(pi/4.0 - plunge1*rad/2.0);
% [X1,Y1] = pol2cart(-plungedir1*rad+pi/2,r);
% plot(X1,Y1,'ro','MarkerFaceColor','r');
%text(X1+0.05,Y1,[' b_1 = ',num2str(round(plungedir1)),' / ', num2str(round(plunge1))],'BackgroundColor',[1 1 1]);

% Second Fold Axial Plane
great_c = (N - dipdirA2)*rad;
apparent_dip = atan(tan(dipA2*rad).*sin(N*rad));
r = s2.*sin(pi/4.0 - apparent_dip/2.0);
[gx,gy] = pol2cart(great_c,r);
plot(gx,gy,'-k','LineWidth',1);
%i = floor(length(N)/2.0);
%text(gx(i)+0.05,gy(i),[' FAP_2 = ',num2str(round(dipdirA2)),' / ', num2str(round(dipA2))],'BackgroundColor',[1 1 1]);

% Second Fold Axis
% plunge2 =  abs(asin(cos(rake2*rad-pi/2)*sin(dipA2*rad))*deg);
% pdd_diff = acos(tan(plunge2*rad)/tan(dipA2*rad))*deg;
% sgn = sign(rake2);
% if sgn == 0;
%     sgn = 1;
% end
% plungedir2 = dipdirA2 - sgn*pdd_diff;
% 
% if plungedir2 < 0.0
%     plungedir2 = 360.0 + plungedir2;
% end
% 
% r = s2.*sin(pi/4.0 - plunge2*rad/2.0);
% [X1,Y1] = pol2cart(-plungedir2*rad+pi/2,r);
% plot(X1,Y1,'bo','MarkerFaceColor','b');
%text(X1+0.05,Y1,[' b_2 = ',num2str(round(plungedir2)),' / ', num2str(round(plunge2))],'BackgroundColor',[1 1 1]);

end
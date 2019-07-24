% Matlab script stereonet
% To plot lines and planes in stereographic
% (equal-angle) projections
% Downloaded from:
% http://www.soest.hawaii.edu/GG/FACULTY/conrad/classes/GG303_F09/labs/Lab_05_matlab/stereonet.m

clf

% Read input data on planes
load planes.dat
% Data in column 1 are strikes, and data in column 2 are dips
% of planes, with angles given in degrees
strike = planes(:,1)*pi/180;
dip = planes(:,2)*pi/180;
num = length(strike);
% find cyclographic traces of planes and plot them
R = 1;
rake = 0:pi/180:pi;
for i=1:num;
	plunge = asin(sin(dip(i)).*sin(rake));
	trend = strike(i) + atan2(cos(dip(i)).*sin(rake), cos(rake));
	rho = R.*tan(pi/4 - (plunge/2));
	% polar plots ccl from 3:00, so convert to cl from 12:00
	polar(pi/2-trend,rho,'-')
	hold on	
%	clear plunge trend rho
end

load lines1.dat
% Data in column 1 are trends, data in column 2 are plunges
% of lines, with angles given in degrees
trend1 = lines1(:,1);
plunge1 = lines1(:,2);
num = length(lines1(:,1));
R = 1;
trendr1 = trend1*pi/180;
plunger1 = plunge1(:,1)*pi/180;
rho1 = R.*tan(pi/4 - ((plunger1)/2));
for i=1:num;
	% polar plots ccl from 3:00, so convert to cl from 12:00
	polar(pi/2-trendr1(i),rho1(i),'o')
	hold on	
end

load lines2.dat
% Data in column 1 are trends, data in column 2 are plunges
% of lines, with angles given in degrees
trend2 = lines2(:,1);
plunge2 = lines2(:,2);
num = length(lines2(:,1));
R = 1;
trendr2 = trend2*pi/180;
plunger2 = plunge2*pi/180;
rho2 = R.*tan(pi/4 - ((plunger2)/2));
for i=1:num;
	% polar plots ccl from 3:00, so convert to cl from 12:00
	polar(pi/2-trendr2(i),rho2(i),'*')
	hold on	
end
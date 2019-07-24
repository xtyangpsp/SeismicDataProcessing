function stereonet

hold on

% Plots an equal area stereonet with radius 1.

% Degree to radians and vice versa factors and square root of 2.0
deg = 180/pi;   rad = pi/180;   s2 = sqrt(2);

% Vectors for plotting great circles
N = linspace(0,pi,181);
gx = zeros(size(N));    gy = zeros(size(N));

%==========================================================================
% PLOT STEREONET

% Plot Great Circles
for i=1:8
    apparent_dip = atan(tan(10.0*i*rad).*sin(N));
    r = s2.*sin(pi/4.0 - apparent_dip/2.0);
    [gx,gy] = pol2cart(N-pi/2.0,r);
    plot(gx,gy,-gx,gy,'Color',[.5 .5 .5]);
end
plot([0 0],[-1 1],'Color',[.5 .5 .5]);

% Plot Small Circles
for i=1:8
    x = cos(10.0*i*rad)*cos(N);
    y = sin(10.0*i*rad)*ones(size(x));
    [azi,plu] = cart2pol(x,y);
    sx = s2*sin(pi/4.0 - acos(plu)/2.0).*sin(pi/2.0-azi);
    sy = s2*sin(pi/4.0 - acos(plu)/2.0).*cos(pi/2.0-azi);
    plot(sx,sy,sx,-sy,'Color',[.5 .5 .5]);
end
plot([-1 1],[0 0],'Color',[.5 .5 .5]);

% Plot primitive circle
[xpc,ypc] = pol2cart(linspace(0,2*pi,3601),1.0);
plot(xpc,ypc,'k','LineWidth',2);

% Tick length
tck = 0.05;

% Plot cross in centre
plot([tck -tck],[0 0],[0 0],[tck -tck],'Color',[0 0 0],'LineWidth',2);

% Plot ticks at N, E, S and W
x = [0 0];    y = [1.0 1.0+tck];
plot(x,y,x,-1*y,y,x,-1*y,x,'Color',[0 0 0],'LineWidth',2);

% Plot tick labels
tx = 0;         ty = 1+3*tck;
t(1) = text(tx,ty,'0^o');       t(2) = text(ty,tx,'90^o');
t(3) = text(tx,-ty,'180^o');    t(4) = text(-ty,tx,'270^o');
set(t,'HorizontalAlignment','center','FontSize',12);
axis equal, axis off
set(gcf, 'color', 'white','Resize','off')
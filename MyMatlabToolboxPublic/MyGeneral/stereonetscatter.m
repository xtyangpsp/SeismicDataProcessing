function stereonetscatter(azimuth,dip,fm,symbol_scalor)
%% Plot scatter points for earthquake first motion data on a stereonet with radius=1.
% Input:
%   azimuth     event to station azimuth;
%   dip         take-off angle;
%   fm          first motion polarities: cc - compressional; dd -
%               dilational; nn - nodal points; cn and dn are intermidiate
%               points.
%   symbol_scalor   To scale the symbols.
%==============================================================
% Written by Xiaotao Yang on May 22, 2013
%
%--------------
% Functions used:
%

%% 

PI=3.141592653;
r=dip/90;
theta=90-azimuth;

theta=theta*PI/180;
[x y]=pol2cart(theta,r);
%x=r*sin(theta); y=r*cos(theta);

hold on;

%for i=1:length(azimuth)
    switch fm
        case 1
            plot(x,y,'k.','MarkerSize',30*symbol_scalor);
        case -1
            plot(x,y,'ko','MarkerSize',10*symbol_scalor,'linewidth',1);
        case 0
            plot(x,y,'kx','MarkerSize',12*symbol_scalor,'linewidth',1);
        otherwise
            error(['Invalid first motion symbols -> ',fm,'. Use ONLY cc, dd, nn! Check data.']);
    end
%end

return
end
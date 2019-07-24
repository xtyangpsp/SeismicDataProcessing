%==========================================================================
%**************************************************************************
% SCRIPT FOR GENERATING A CUT-OUT AND FOLD PAPERMODEL ILLUSTRATING
% SUPERPOSED FOLDING
%**************************************************************************
%==========================================================================
%
% written and copyright by:
%
% Martin Schöpfer
% Fault Analysis Group
% UCD School of Geological Sciences
% University College Dublin
% Belfield, Dublin 4, Ireland
% email: martin@fag.ucd.ie
% Tel: +353 - 1 - 716 2611
% Fax: +353 - 1 - 716 2607
% Web: http://www.fault-analysis-group.ucd.ie/
%
%==========================================================================
%
% The script below uses the following functions:
%
% stereonet.m
% stereoplot.m
% Stereonet_FoldAxes_Selection.m
% fold_interference_pattern.m
%
%==========================================================================
clear
close all
%==========================================================================
% INPUT PARAMETERS

% Fold Parameters
% Folds are calculated using the following plane-strain equations:
%
% x = X
% y = C*Y
% z = B*Z + A*B*sin(Y)
%
% The fold amplitude is A*B, and factors B and C are stretches
% (homogeneous strain). If C = 1/B deformation is constant volume.
%
% Dip directions should range from 0 to 360
% Dips should range from 0 to 90
% Rake should range from -90 to +90.
% A positive rake implies that the plunge direction of the fold axis is
% anti-clockwise from the dip direction of the fold axial plane.
% If the rake is 90 (or -90) the plunge direction / plunge of the fold axis
% is identical to the dip direction / dip of the fold axial plane.

% First generation parameters
A1 = 3;
B1 = 1;
C1 = 1/B1;

% Second generation parameters
A2 = 3;
B2 = 1;
C2 = 1/B2;

% First fold orientation
dipdirA1 = 0;
dipA1 = 90;
rake1 = 0;

% Second fold orientation
dipdirA2 = 90;
dipA2 = 90;
rake2 = 0;

% Number of points, np, along the length of each face of the cube with length L
% Determines quality of model, the greater np, the better (np ~750
% provides print quality).
np = 750;
L = 5*pi;

% Layer thickness (e.g., L/7.5 for b/w; L/15 for colour)
t  = L/15;

% Change the state of the random number generator for different colour
% scheme
rand('state',26)

%==========================================================================
% Generate random colormap (500 colours should be plenty)
Ncol = 500;
cm = ones(Ncol,3);

% R G B values
cm(:,1) = rand(Ncol,1);     % R
cm(:,2) = rand(Ncol,1);     % G
cm(:,3) = rand(Ncol,1);     % B

%==========================================================================
% FOLD ORIENTATIONS

% Ask user whether he/she wants to select fold orientations using
% stereonet selection, the data listed above or whether the 'computer'
% should select the orientations
stereo = menu('','Select Fold Orientations with Stereonet','Use Fold Orienations in Script','Random Fold Orientations');

if stereo == 1

    % Plot Stereonet
    stereonet
    % Call stereonet selection script to get fold orientations
    [dipdirA1,dipA1,rake1,dipdirA2,dipA2,rake2] = Stereonet_FoldAxes_Selection;
   
elseif stereo == 2
    
    % Plot Stereonet
    stereonet
    % Plot fold orientation data listed above
    stereoplot(dipdirA1,dipA1,rake1,dipdirA2,dipA2,rake2);
    
elseif stereo == 3
    
    % Randomly change the state of the random number generator  
    rand('state',sum(100*clock))
    
    % First fold orientation
    
    dipdirA1 = rand*360;
    dipA1 = rand*90;
    rake1 = rand*180-90;

    % Second fold orientation
    dipdirA2 = rand*360;
    dipA2 = rand*90;
    rake2 = rand*180-90;
    
    % Plot Stereonet
    stereonet
    % Plot random fold orientation data
    stereoplot(dipdirA1,dipA1,rake1,dipdirA2,dipA2,rake2);
   
end

% Ask the user which colour scheme should be used
color_option = menu('','Black and White','Colour','Layer Boundaries','Cancel');

% Cancel
if color_option == 4
    close all
    return
end
%==========================================================================
% CALCULATE INTERFERENCE PATTERNS FOR THE SIX FACES OF A CUBE

% Mesh for faces of unit cube
s  = linspace(0,L,np+1);
M  = meshgrid(s,s);
O  = zeros(size(M));

% Face 1
X = M;
Y = M';
Z = O;
z1 = fold_interference_pattern(X,Y,Z,dipdirA1,dipA1,rake1,A1,B1,C1,dipdirA2,dipA2,rake2,A2,B2,C2);

% Face 2
X = M;
Y = M';
Z = O + max(s);
z2 = fold_interference_pattern(X,Y,Z,dipdirA1,dipA1,rake1,A1,B1,C1,dipdirA2,dipA2,rake2,A2,B2,C2);

% Face 3
X = M';
Y = O;
Z = M;
z3 = fold_interference_pattern(X,Y,Z,dipdirA1,dipA1,rake1,A1,B1,C1,dipdirA2,dipA2,rake2,A2,B2,C2);

% Face 4
X = M';
Y = O+max(s);
Z = M;
z4 = fold_interference_pattern(X,Y,Z,dipdirA1,dipA1,rake1,A1,B1,C1,dipdirA2,dipA2,rake2,A2,B2,C2);

% Face 5
X = O;
Y = M';
Z = M;
z5 = fold_interference_pattern(X,Y,Z,dipdirA1,dipA1,rake1,A1,B1,C1,dipdirA2,dipA2,rake2,A2,B2,C2);

% Face 6
X = O+max(s);
Y = M';
Z = M;
z6 = fold_interference_pattern(X,Y,Z,dipdirA1,dipA1,rake1,A1,B1,C1,dipdirA2,dipA2,rake2,A2,B2,C2);

%==========================================================================
% CALCULATE LAYER BOUNDARIES

if color_option == 1 | color_option == 3 % Black and white or boundary option

    % Calculate layers (either 1 or 0)
    c1 = round(0.5*z1/t) - floor(0.5*z1/t);
    c2 = round(0.5*z2/t) - floor(0.5*z2/t);
    c3 = round(0.5*z3/t) - floor(0.5*z3/t);
    c4 = round(0.5*z4/t) - floor(0.5*z4/t);
    c5 = round(0.5*z5/t) - floor(0.5*z5/t);
    c6 = round(0.5*z6/t) - floor(0.5*z6/t);

elseif color_option == 2 % Colour option
   
    % Get minimum and maximum z-position
    zminmax(1,:) = [min(min(z1))  max(max(z1))];
    zminmax(2,:) = [min(min(z2))  max(max(z2))];
    zminmax(3,:) = [min(min(z3))  max(max(z3))];
    zminmax(4,:) = [min(min(z4))  max(max(z4))];
    zminmax(5,:) = [min(min(z5))  max(max(z5))];
    zminmax(6,:) = [min(min(z6))  max(max(z6))];
    zmin = min(zminmax(:,1));
    zmax = max(zminmax(:,2));

    % Calculate normalised z-positions (1 to Nlayers)
    c1 = ceil((z1-zmin)/t);
    c2 = ceil((z2-zmin)/t);
    c3 = ceil((z3-zmin)/t);
    c4 = ceil((z4-zmin)/t);
    c5 = ceil((z5-zmin)/t);
    c6 = ceil((z6-zmin)/t);

    % Eliminate 0 value at zmin
    c1(find(c1 == 0)) = 1;
    c2(find(c2 == 0)) = 1;
    c3(find(c3 == 0)) = 1;
    c4(find(c4 == 0)) = 1;
    c5(find(c5 == 0)) = 1;
    c6(find(c6 == 0)) = 1;
    
    % Number of layers
    Nlayers = ceil((zmax-zmin)/t);

    % Generate colormap for papermodel and ensure that same 'stratigraphy' is used
    cmpm = cm((1:Nlayers)+round(Ncol/2)+floor(zmin/t),:);

end

%==========================================================================
% ASSEMBLE FACES OF PAPERMODEL AND PLOT PATTERN

% Length of each face
i   = length(s);

% Put together faces in one array (C)
C = zeros(3*i,4*i);
C(i+1:2*i,1:i)   = c1;
C(i+1:2*i,i+1:2*i) = c6;
C(i+1:2*i,2*i+1:3*i) = flipdim(c2,2);
C(i+1:2*i,3*i+1:4*i) = flipdim(c5,2);
C(1:i,2*i+1:3*i) = flipdim(c3',2);
C(2*i+1:3*i,2*i+1:3*i) = flipdim(flipdim(c4',1),2);

% Open a DIN-A shaped figure window, centred on the screen
scrsz = get(0,'ScreenSize');
fLy = min(scrsz(3:4))*.8;
fLx = sqrt(2)*fLy;
fx = (scrsz(3)-fLx)*.5;
fy = (scrsz(4)-fLy)*.5;
f2 = figure('Position',[fx fy fLx fLy]);
hold on

% Contour and plot the pattern - black and white option
if color_option == 1
    col = [0 0 0];
    [CC,h] = contourf(C,[0.5 0.5]);
    set(h,'EdgeColor','none');
    colormap(col)
% Contour and plot the pattern - colour option
elseif color_option == 2
    [CC,h] = contourf(C,[(1:Nlayers+1)-0.5]);
    set(h,'EdgeColor','none');
    colormap(cmpm)
% Contour and plot the pattern - layer boundary option
elseif    color_option == 3
    [CC,h] = contour(C,[0.5 0.5]);
    set(h,'EdgeColor','k');
end

% Edges of faces
xs  = [0 i i 0 0]+.5;
ys  = [0 0 i i 0]+.5;
plot(xs,ys+i,xs+i,ys+i,xs+2*i,ys+i,xs+3*i,ys+i,xs+2*i,ys,xs+2*i,ys+2*i,'Color',[0 0 0],'LineWidth',.5);

% Flaps
flapfrac = 8.5;
xf = [0 i/flapfrac i-i/flapfrac i]+.5;
yf = [0 i/flapfrac i/flapfrac 0]+.5;
plot(xf,yf+2*i,xf+i,yf+2*i,xf+3*i,yf+2*i,'Color',[0 0 0],'LineWidth',.5);
yf = -yf+1;
plot(xf,yf+i,xf+i,yf+i,xf+3*i,yf+i,-abs(yf)+1,xf+i,'Color',[0 0 0],'LineWidth',.5);

% Axes and background properties
axis([-2*i/9 4*i+i/9 -i/9 3*i+i/9])
set(gca,'Position',[0 0 1 1])
axis equal, axis off
set(gcf, 'color', 'white','Resize','off')

% Provide information about papermodel
text(0,2.95*i,'Papermodel of Superposed Folding','FontSize',15,'FontWeight','bold','FontUnits','normalized')
text(0,2.75*i,['1^{st} Fold System:'],'FontSize',10,'FontWeight','bold','FontUnits','normalized')
text(0,2.60*i,['Dip Direction and Dip of Fold Axial Planes: ',num2str(round(dipdirA1)),' / ',num2str(round(dipA1))],'FontSize',10,'FontUnits','normalized')
text(0,2.45*i,['Pitch (Rake) of Fold Axes: ',num2str(round(rake1))],'FontSize',10,'FontUnits','normalized')
text(0,2.30*i,['Fold Parameters: A_1 = ',num2str(A1),'  B_1 = ',num2str(B1),'  C_1 = ',num2str(C1)],'FontSize',10,'FontUnits','normalized')
text(0,0.75*i,['2^{nd} Fold System:'],'FontSize',10,'FontWeight','bold','FontUnits','normalized') 
text(0,0.60*i,['Dip Direction and Dip of Fold Axial Planes: ',num2str(round(dipdirA2)),' / ',num2str(round(dipA2))],'FontSize',10,'FontUnits','normalized')
text(0,0.45*i,['Pitch (Rake) of Fold Axes: ',num2str(round(rake2))],'FontSize',10,'FontUnits','normalized')
text(0,0.30*i,['Fold Parameters: A_2 = ',num2str(A2),'  B_2 = ',num2str(B2),'  C_2 = ',num2str(C2)],'FontSize',10,'FontUnits','normalized')
text(0,0.05*i,['Matlab script available at: www.fault-analysis-group.ucd.ie (M.P.J. Schöpfer)'],'FontSize',8,'FontUnits','normalized')

% Ask whether current papermodel should be saved or not
q = menu('Save Papermodel?','Yes','No');

if q == 1
    % Save figure as PNG file (r = resolution)
    % Filename contains dipdirection and dip of fold axial plane and rake
    % of fold axis of first and second folds, F1 and F2
    % Also show whether model is coloured, b/w or just layer boundaries
    if color_option == 1
        colinfo = 'bw';
    elseif color_option == 2
        colinfo = 'col';
    elseif color_option == 3
        colinfo = 'bound';
    end  
    print('-f2','-dpng','-r600',['Superposed_Papermodel_F1_',num2str(round(dipdirA1)),'_',...
            num2str(round(dipA1)),'_',num2str(round(rake1)),'_F2_',num2str(round(dipdirA2)),...
            '_',num2str(round(dipA2)),'_',num2str(round(rake2)),'_',colinfo,'.png']);
end

% END OF SCRIPT
%==========================================================================
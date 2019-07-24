function zini = fold_interference_pattern(X,Y,Z,dipdirA1,dipA1,rake1,A1,B1,C1,dipdirA2,dipA2,rake2,A2,B2,C2);

%==========================================================================
% This function calculates the initial z-positions for finite
% x,y,z-positions given by X,Y,Z as a result of two folding events defined
% by the fold orientations (fold axial plane and rake of fold axis) and
% fold parameters A,B & C. The complete set of equations is:
%
% x1 = Y*sin(beta_z)+cos(beta_z)*X;
% y1 = cos(beta_z)*Y-X*sin(beta_z);
% z1 = Z;
% 
% x2 = x1;
% y2 = z1*sin(beta_x)+cos(beta_x)*y1;
% z2 = cos(beta_x)*z1-y1*sin(beta_x);
% 
% x3 = z2*sin(beta_y)+cos(beta_y)*x2;
% y3 = y2;
% z3 = cos(beta_y)*z2-x2*sin(beta_y);
% 
% 2nd folding
% x4 = x3;
% y4 = y3/C2;
% z4 = -(A2*B2*sin(y3/C2)-z3)/B2;
% 
% x5 = -z4*sin(beta_y)+cos(beta_y)*x4;
% y5 = y4;
% z5 = cos(beta_y)*z4+x4*sin(beta_y);
% 
% x6 = x5;
% y6 = -z5*sin(beta_x)+cos(beta_x)*y5;
% z6 = cos(beta_x)*z5+y5*sin(beta_x);
% 
% x7 = -y6*sin(beta_z)+cos(beta_z)*x6;
% y7 = cos(beta_z)*y6+x6*sin(beta_z);
% z7 = z6;
% 
% x8 = y7*sin(alpha_z)+cos(alpha_z)*x7;
% y8 = cos(alpha_z)*y7-x7*sin(alpha_z);
% z8 = z7;
% 
% x9 = x8;
% y9 = z8*sin(alpha_x)+cos(alpha_x)*y8;
% z9 = cos(alpha_x)*z8-y8*sin(alpha_x);
% 
% x10 = z9*sin(alpha_y)+cos(alpha_y)*x9;
% y10 = y9;
% z10 = cos(alpha_y)*z9-x9*sin(alpha_y);
% 
% 1st folding
% xini = x10;
% yini = y10/C1;
% zini = -(A1*B1*sin(y10/C1)-z10)/B1;
%
% where alpha and beta are rotation angles of first and second fold,
% respectively, and defined as:
%
% x-rotation = 90 - dip of fold axial plane
% y-rotation = rake (= pitch)
% z-rotation = (minus) dip direction of fold axial plane

%==========================================================================
% Rotation angles of first fold
alpha_x = (90-dipA1)*pi/180;
alpha_y = rake1*pi/180;
alpha_z = -dipdirA1*pi/180;

% Rotation angles of second fold
beta_x = (90-dipA2)*pi/180;
beta_y = rake2*pi/180;
beta_z = -dipdirA2*pi/180;

%==========================================================================
% Ensure no division by zero occurs
if dipA1  == 0
    dipA1 = 1e-6;
end
if dipA2  == 0
    dipA2 = 1e-6;
end

%==========================================================================
% Calculate sine and cosine of angles
sax = sin(alpha_x);
cax = cos(alpha_x);
say = sin(alpha_y);
cay = cos(alpha_y);
saz = sin(alpha_z);
caz = cos(alpha_z);

sbx = sin(beta_x);
cbx = cos(beta_x);
sby = sin(beta_y);
cby = cos(beta_y);
sbz = sin(beta_z);
cbz = cos(beta_z);

%==========================================================================
% Non-simplified, but exact solution of the set of equations listed above

zini = (-A1*B1*sin(((cbx*(cby*(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*Z-(cbz*Y-X*sbz)*sbx)-...
    (Y*sbz+cbz*X)*sby)/B2+((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+cby*(Y*sbz+cbz*X))*sby)+...
    (Z*sbx+cbx*(cbz*Y-X*sbz))/C2*sbx)*sax+cax*(caz*(cbz*((-cby*(-A2*B2*sin((Z*sbx+...
    cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*Z-(cbz*Y-X*sbz)*sbx)-(Y*sbz+cbz*X)*sby)/B2-((cbx*Z-...
    (cbz*Y-X*sbz)*sbx)*sby+cby*(Y*sbz+cbz*X))*sby)*sbx+cbx*(Z*sbx+cbx*(cbz*Y-...
    X*sbz))/C2)+(-(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*Z-(cbz*Y-X*sbz)*sbx)-(Y*sbz+...
    cbz*X)*sby)/B2*sby+cby*((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+cby*(Y*sbz+cbz*X)))*sbz)-...
    ((-(-cby*(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*Z-(cbz*Y-X*sbz)*sbx)-(Y*sbz+...
    cbz*X)*sby)/B2-((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+cby*(Y*sbz+cbz*X))*sby)*sbx-cbx*(Z*sbx+...
    cbx*(cbz*Y-X*sbz))/C2)*sbz+cbz*(-(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*Z-(cbz*Y-...
    X*sbz)*sbx)-(Y*sbz+cbz*X)*sby)/B2*sby+cby*((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+...
    cby*(Y*sbz+cbz*X))))*saz))/C1)+cay*(cax*(cbx*(cby*(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-...
    X*sbz))/C2)+cby*(cbx*Z-(cbz*Y-X*sbz)*sbx)-(Y*sbz+cbz*X)*sby)/B2+((cbx*Z-(cbz*Y-X*sbz)*...
    sbx)*sby+cby*(Y*sbz+cbz*X))*sby)+(Z*sbx+cbx*(cbz*Y-X*sbz))/C2*sbx)-(caz*(cbz*...
    ((-cby*(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*Z-(cbz*Y-X*sbz)*sbx)-(Y*sbz+...
    cbz*X)*sby)/B2-((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+cby*(Y*sbz+cbz*X))*sby)*sbx+cbx*(Z*sbx+...
    cbx*(cbz*Y-X*sbz))/C2)+(-(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*Z-(cbz*Y-X*sbz)*...
    sbx)-(Y*sbz+cbz*X)*sby)/B2*sby+cby*((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+cby*(Y*sbz+cbz*X)))*...
    sbz)-((-(-cby*(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*Z-(cbz*Y-X*sbz)*sbx)-...
    (Y*sbz+cbz*X)*sby)/B2-((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+cby*(Y*sbz+cbz*X))*sby)*sbx-...
    cbx*(Z*sbx+cbx*(cbz*Y-X*sbz))/C2)*sbz+cbz*(-(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*...
    Z-(cbz*Y-X*sbz)*sbx)-(Y*sbz+cbz*X)*sby)/B2*sby+cby*((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+cby*...
    (Y*sbz+cbz*X))))*saz)*sax)-((cbz*((-cby*(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*...
    Z-(cbz*Y-X*sbz)*sbx)-(Y*sbz+cbz*X)*sby)/B2-((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+cby*(Y*sbz+...
    cbz*X))*sby)*sbx+cbx*(Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+(-(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+...
    cby*(cbx*Z-(cbz*Y-X*sbz)*sbx)-(Y*sbz+cbz*X)*sby)/B2*sby+cby*((cbx*Z-(cbz*Y-X*sbz)*sbx)*...
    sby+cby*(Y*sbz+cbz*X)))*sbz)*saz+caz*((-(-cby*(-A2*B2*sin((Z*sbx+cbx*(cbz*Y-X*sbz))/C2)+...
    cby*(cbx*Z-(cbz*Y-X*sbz)*sbx)-(Y*sbz+cbz*X)*sby)/B2-((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+...
    cby*(Y*sbz+cbz*X))*sby)*sbx-cbx*(Z*sbx+cbx*(cbz*Y-X*sbz))/C2)*sbz+cbz*(-(-A2*B2*sin((Z*sbx+...
    cbx*(cbz*Y-X*sbz))/C2)+cby*(cbx*Z-(cbz*Y-X*sbz)*sbx)-(Y*sbz+cbz*X)*sby)/B2*sby+cby*...
    ((cbx*Z-(cbz*Y-X*sbz)*sbx)*sby+cby*(Y*sbz+cbz*X)))))*say)/B1;

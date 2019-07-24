function vout=velocity_layeroperator(X,Y,Z,V,operation,layertype,layerdepth)
% Get subset volume with specified depth interval (after interpolation).
% USAGE: vout=velocity_layeroperator(X,Y,Z,V,layerdepth,operation)
% INPUT:
%       X,Y,Z,V: meshed grids. It reads in the 3-D model grid from X,Y,Z
%               and the velocity grid from V. These arrays must be the same size.
%       operation: type of operation to the velocity values on the grids
%               within the layer. It has to be a function operating on a 1-D array.
%               It HAS TO START WITH @. For example: @mean, @max, @min,
%               @median, and @std.
%       layertype: type of layer. It could be one of: flat OR surface. 'flat' means the layer
%               is bounded by top and bottom surfaces all at constent
%               depth. 'surface' means at least one of the top and bottom
%               depths are on a non-flat surface. Choice of this argument
%               determines the type of 'layerdepth' value. See below for
%               details.
%       layerdepth: depth array for the layer in interpolating the velocity into uniform
%               grids. If ;layertype' is 'flat', then this array has to be
%               a 1-D array. Otherwise, if the 'layertype' is 'surface',
%               this has to be a cell array containing THREE elements: top
%               grid, bottom grid, and depth interval. 
%
%  By Xiaotao Yang @ UMass Amherst, 2017-2018
%  History:
%       1. Originally wrote in early 2017.
%       2. Modified on September 5, 2018 to improve the efficiency in
%       interpolating the values. Changed to use interp1 to interp3. The
%       input arguments are also changed to read in XYZ grid information.
%       For old call of this function, this has to be taken care of. The
%       new usage is NOT compatible with the old one.
%       3. Modified on 9/15/2018: added back the original function to
%       operate on layers bounded by two surfaces, instead of two depths.
%       In this case, the program takes much more time. Be patient.
%       4. 9/15/2018: remove output for layerv. Only return vout now.
%
%%%%%%% See below for an example of using this function. XSIM, YSIM, ZSIM are
%%%%%%% the simulation model grid from FWANT.
% zgridtoplist=[10,25,45,70,110,150];
% zgridbotlist=[20,35,55,110,140,200];
% griddepth=squeeze(ZSIM(1,1,end)/1000-abs(ZSIM(npml,npml,:))/1000);
% dz=2.0;
% operation=@mean;
% clear Xgrid Ygrid Zgrid;
% [Ygrid,Xgrid,Zgrid]=meshgrid(squeeze(YSIM(1,:,1)),squeeze(XSIM(:,1,1)),griddepth);
% <<<<When using layertype of 'flat'>>>>
% vout=cell(length(zgridtoplist),1);
% for iz=1:length(zgridtoplist)
%     disp(['Extracting layer: ',num2str(zgridtoplist(iz)),' km to ',num2str(zgridbotlist(iz)),' km ...']);
%     vout{iz}=velocity_layeroperator(Xgrid,Ygrid,Zgrid,mvsmoothed{iz}.v,operation,'flat',...
%         zgridtoplist(iz):dz:zgridbotlist(iz));
% end
%
% <<<<When using layertype of 'surface'>>>>
% zgrid15.data=15*ones(size(XSIM,1),size(XSIM,2));zgrid15.tag='15';
% zgridmoho.data=griddata(moho(:,1),moho(:,2),moho(:,3),YSIM(1,:,1)-360,XSIM(:,1,1));
% layergrid={zgrid15.data,zgridmoho.data,dz};
% vout=velocity_layeroperator(Xgrid,Ygrid,Zgrid,mvsmoothed{iz}.v,operation,'surface',...
%     layergrid);
%%%% end of example

vout=nan(size(squeeze(X(:,:,1))));
if strcmp(layertype,'flat')
    [YY,XX,ZZ]=meshgrid(squeeze(Y(1,:,1)),squeeze(X(:,1,1)),layerdepth);

    vout=operation(interp3(Y,X,Z,V,YY,XX,ZZ),3);
elseif strcmp(layertype,'surface')
    griddepth=squeeze(Z(1,1,:));
    for i=1:size(layerdepth{1},1)
        for j=1:size(layerdepth{1},2)
            clear idepth;
    %         clear interpv;
            %
    %         gridv=squeeze(v(i,j,:));
    %         interpv=interp1(griddepth,squeeze(v(i,j,:)),idepth);
            idepth=layerdepth{1}(i,j):layerdepth{3}:layerdepth{2}(i,j);
%             min(idepth)
%             max(idepth)
            vout(i,j)=operation(interp1(griddepth,squeeze(V(i,j,:)),idepth));
        end
    end
else
    error(['**ERROR: Wrong layertype [ ',layertype,' ]. Has to be flat OR surface!']);
end

end
% stack
function [dout, w]= smartstack(din,varargin)
%Robust stacking scheme is based on Pavlis and Vernon,
%Computers & Geosciences, 2010.
%
%Wrote by: Xiaotao Yang July 2019
%
maxite = 100; % stop RobustSNR stacking when reaches this number of iterations.

stacktype='robust'; %specify as: mean, median, robust

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% temporary block for testing purpose
% load('datafortest.mat');
% dt=abs(timeflag(2)-timeflag(1));
% din=egfraw(:,30:60);
% tmin=-300;
% tmax=300;
% idmin=round((tmin - timeflag(1))/dt);
% idmax=round((tmax - timeflag(1))/dt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END of temporary block for testing purpose
stackwinidx = [];
if nargin < 1
    error('Not enough input arguments');
end

if nargout > 2
    error('Too many output arguments');
end

smallnumber = 0.00001; %used to mark convergence when the stack beam changes less than this value.
%this number is used only in RobustSNR stack type. See Pavlis and Vernon,
%Computers & Geosciences, 2011 for details.
largenumber = 10000; %used as the intial value of model norm.
stackdim = 2;

numvarargin=length(varargin);
for nar=1:numvarargin
    if strcmp(varargin{nar},'stacktype')
       if nar+1>numvarargin
           error('** need to specify stacktype as: mean , median, OR robust [default]!');
       else
           stacktype=varargin{nar+1};
       end
    elseif strcmp(varargin{nar},'stackwindow')
       if nar+1>numvarargin
           error('** need to specify stackwindow as a two-element vector!');
       else
           stackwinidx=varargin{nar+1};
           if ~isnumeric(stackwinidx)
               error('** need to specify stackwindow as a two-element vector!');
           end
       end
    elseif strcmp(varargin{nar},'maxiteration')
       if nar+1>numvarargin
           error('** need to specify maxiteration!');
       else
           maxite=varargin{nar+1};
       end
    elseif strcmp(varargin{nar},'stackdim')
       if nar+1>numvarargin
           error('** need to specify stackdim [default is 2]!');
       else
           stackdim=varargin{nar+1};
       end
    end
end

switch stackdim
    case 1
        [ndata,nsamp]=size(din);
    case 2
        [nsamp,ndata]=size(din);
end

if (strcmp(stacktype,'robust') && isempty(stackwinidx)) 
    warning('stackwinidx not specified, use the whole data length instead!');
    stackwinidx=[1,nsamp];
end
if ~isempty(stackwinidx)
    idmin=stackwinidx(1);
    idmax=stackwinidx(2);
end

if strcmp(stacktype,'mean')
    dout=nanmean(din,stackdim);
    w = ones(ndata,1)/ndata;
elseif strcmp(stacktype,'median')
    dout=nanmedian(din,stackdim);
    w = ones(ndata,1)/ndata;
elseif strcmp(stacktype,'robust')
    switch stackdim
        case 1
            data=din(:,idmin:idmax);
            data = data';
            dout=zeros(1,nsamp);
        case 2
            data=din(idmin:idmax,:);
            dout=zeros(nsamp,1);
    end
    s = nanmedian(data,2); %start from median stack.
    
    %get initial deltad: left hand side term in Pavlis and Vernon, equation (6)
    deltad = largenumber;
    itecount = 0;
    w = zeros(ndata,1);
    ampscale = zeros(ndata,1);
    while deltad > smallnumber && itecount < maxite
        for i = 1:ndata
            dtemp = data(:,i);
            if range(dtemp) == 0 || ~isempty(find(isnan(dtemp),1))
                w(i) = 0; % for zero trace or NaN trace, use zero weight.
            else                
                ampscale(i)=abs(dot(s,dtemp));
                d = norm(dtemp);
                r = norm(dtemp - dot(s,dtemp).*s);
                w(i) = ampscale(i)./(d.*r);
            end
        end
        
        % normalize weight
        w=w./nansum(w);
        stemp = zeros(idmax - idmin +1,1);
        for i = 1:ndata
            stemp = stemp + w(i)*data(:,i);
        end
        slast = s;
        s = stemp;
        deltad = norm(s - slast,1)/(norm(s)*nsamp);
        itecount = itecount + 1;
    end
    
    switch stackdim
        case 1
            for j = 1:ndata
                dout = dout + w(j)*din(j,:);
            end
        case 2
            for j = 1:ndata
                dout = dout + w(j)*din(:,j);
            end
    end
end

return;
end
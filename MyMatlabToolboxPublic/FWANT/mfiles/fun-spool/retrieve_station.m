function [v]=retrieve_station(seismoinfo,varnm,varargin)
% retreive station information on individual receivers.
%
% $Date: 2009-01-13 14:35:28 -0500 (Tue, 13 Jan 2009) $
% $Revision: 38 $
% $LastChangedBy: zhangw $

pnm_nc='../input/';

%-- flags --
n=1;
while n<=nargin-3

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case 'stationdir'
       pnm_nc=varargin{n+1}; n=n+1;
   end
end
n=n+1;

end

% check
if ~ exist(pnm_nc,'dir')
   error([mfilename ': directory ' pnm_nc ' does not exist']);
end

fnm_nc=get_fnm_station(pnm_nc,seismoinfo.n_i,seismoinfo.n_j,seismoinfo.n_k);

v=nc_varget(fnm_nc,varnm,[seismoinfo.npt-1,0],[1 -1]);



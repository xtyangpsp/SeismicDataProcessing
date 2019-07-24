function [fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
          dir_station,dir_out]=get_simul_path(varargin)
%
% get_simul_path: Get the conf name and path of simulation for FD3Dtopo-nonstaggered.
%
% Usage: [fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
%         dir_station,dir_out]=get_simul_path(pwd)
%
% Input:
%   'root',dir_name:
%       name of root directory for a simulation.
%
% Output:
%   fnm_conf:
%       full name of SeisFD3D.conf.
%   dir_coord:
%       pathname of grid coordinate nc files.
%   dir_metric:
%       pathname of grid metric nc files.
%   dir_media:
%       pathname of medium nc files.
%   dir_source:
%       pathname of source nc files.
%   dir_station:
%       pathname of station nc files.
%   dir_out:
%       pathname of output nc files.
%
% Example:
%   >> [fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
%       dir_station,dir_out]=get_simul_path('root','.')

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Changed file name to get_simul_path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date: 2008-04-27 17:31:28 -0400 (Sun, 27 Apr 2008) $
% $Revision: 469 $
% $LastChangedBy: zhangw $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MFILE_ROOT='/net/fs01/data/wzhang/FD3Dtopo-nonstaggered/mfiles'
%RUN_ROOT='/net/fs01/data/wzhang/FD3Dtopo-nonstaggered/run/example.cart1d'
%RUN_ROOT='../run/example.cart1d'

RUN_ROOT='../../sim.station/CD.KMI/fx'

nargs=nargin;

n=1;
while n<=nargs

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case 'root'
       RUN_ROOT=varargin{n+1}; n=n+1;
   end
end

n=n+1;

end

fnm_conf   =[RUN_ROOT '/' 'SeisFD3D.conf']
dir_coord  =[RUN_ROOT '/' 'input']
dir_metric =[RUN_ROOT '/' 'input']
dir_media  =[RUN_ROOT '/' 'input']
dir_source =[RUN_ROOT '/' 'input']
dir_station=[RUN_ROOT '/' 'input']
dir_out    =[RUN_ROOT '/' 'output']

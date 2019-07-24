function set_mfiles_path
%
% set_mfiles_path: Set the path of mfiles function.
%
% Usage: set_mfiles_path
%
% Input:
%   none
%
% Output:
%   none
%
% Example:
%   >> set_mfiles_path()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date: 2008-04-27 17:31:28 -0400 (Sun, 27 Apr 2008) $
% $Revision: 469 $
% $LastChangedBy: zhangw $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MFILE_ROOT='.'

path([MFILE_ROOT '/saclab'],path);
path([MFILE_ROOT '/fun-spool'],path);
%path([MFILE_ROOT '/fileexchange'],path);
%path([MFILE_ROOT '/others'],path);

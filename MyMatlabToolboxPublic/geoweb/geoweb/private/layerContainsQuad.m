function insideLimits = layerContainsQuad(layers, inputLatlim, inputLonlim)
%LAYERCONTAINSQUAD Determine if layer contains quadrangle
%
%   insideLimits = layerContainsQuad(LAYERS, LATLIM, LONLIM) returns true
%   for each layer that completely contains the specified quadrangle
%   bounded by LATLIM and LONLIM.  LAYERS is a structure array with fields
%   'Latlim' and 'Lonlim'.  If LATLIM is empty, then the layer's latitude
%   limits are ignored.  If LONLIM is empty, then the layer's longitude
%   limits are ignored.  

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/12/22 23:42:49 $

% Setup a logical array that is false for each layer.
insideLimits = false(size(layers));
    
% Loop through each element of layers and determine if the layer's limits
% are outside the user-specified limits.
for k=1:numel(layers)
    
    % Determine which limits to use.
    if isempty(inputLatlim)
        % Use the layer's limits.
        latlim = layers(k).Latlim;
    else
        % Use the input limits.
        latlim = inputLatlim;
    end
    if isempty(inputLonlim)
        % Use the layer's limits.
        lonlim = layers(k).Lonlim;
    else
        % Use the input limits.
        lonlim = inputLonlim;
    end
    
    % Return true if latlim and lonlim are inside the layer's limits.
    insideLimits(k) = isgeosubquad( ...
        latlim, lonlim, layers(k).Latlim, layers(k).Lonlim);    
end

%--------------------------------------------------------------------------

function tf = isgeosubquad(latlimSub, lonlimSub, latlim, lonlim)
% Return true if the geographic quadrangle defined by latlimSub and
% lonlimSub is a subset of the quadrangle defined by latlim and lonlim.

[latlimInt, lonlimInt] = intersectgeoquad( ...
    latlimSub, lonlimSub, latlim, lonlim);

tf = isequal(latlimInt,latlimSub) ...
    && isequal(wrapTo180(lonlimInt),wrapTo180(lonlimSub));

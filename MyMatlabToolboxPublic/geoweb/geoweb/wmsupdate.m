function [updatedLayers, index] = wmsupdate(layers, varargin)
%WMSUPDATE Synchronize WMSLayer object with server
%
%   [updatedLayers, INDEX] = WMSUPDATE(LAYERS) returns a WMSLayer array
%   with its properties synchronized with values from the server. The input
%   LAYERS must contain only one unique ServerURL.  Layers no longer
%   available on the server are removed. INDEX is a logical array
%   containing true for each layer that was available; therefore,
%   updatedLayers has the same size as LAYERS(INDEX). Except for deletion,
%   the order of LAYERS is preserved in updateLayer. The function accesses
%   the Internet to update the properties. An Internet connection must be
%   established to use the function. The WMS server may periodically be
%   unavailable and several minutes may elapse before the layers are
%   updated. The function times-out after 60 seconds if a connection is not
%   made to the server.
%
%   [...] = WMSUPDATE(LAYERS, PARAM1, VAL1, PARAM2, VAL2, ...) specifies
%   parameter-value pairs that modify the request. Parameter names can be
%   abbreviated and are case-insensitive. See the table below for a list of
%   available parameters.
%
%   Optional Parameters
%   -------------------
%
%   'TimeoutInSeconds'     An integer-valued scalar double indicating the 
%                          number of seconds to elapse before a server
%                          timeout is issued. By default, the value is 60
%                          seconds. A value of 0 causes the timeout
%                          mechanism to be ignored.
%
%   'AllowMultipleServers' A logical scalar indicating whether the layer 
%                          array may contain elements from multiple
%                          servers. By default, the value is false,
%                          indicating the array must contain elements from
%                          the same server. Use caution when setting the
%                          value to true, since a request is made to each
%                          unique server and each request may take several
%                          minutes to finish.
%                          
%   Firewall Note 
%   ------------- 
%   If you need to specify a proxy server to connect to the Internet,
%   select File -> Preferences -> Web and enter your proxy information. Use
%   this feature if you have a firewall.
%
%   Example 1
%   ---------
%   % Update the layers found from NASA Goddard Space Flight Center's WMS
%   % SVS Image Server. Search the abstract field of the updated layers to
%   % find layers containing the term 'blue marble'. Read and display the
%   % first blue marble layer containing the term '512' and 'image' in its
%   % LayerTitle.
%   gsfc = wmsfind('svs.gsfc.nasa.gov', 'SearchField', 'serverurl');
%   gsfc = wmsupdate(gsfc);
%   blue_marble = gsfc.refine('blue marble', 'SearchField' ,'abstract');
%   queryStr = '*512*image';
%   layers = blue_marble.refine(queryStr);
%   layer = layers(1);
%
%   % Display the layer and abstract.
%   [A, R] = wmsread(layer);
%   figure
%   worldmap world
%   plabel off; mlabel off
%   geoshow(A, R);
%   title(layer.LayerTitle)
%   disp(layer.Abstract)
%
%   Example 2
%   ---------
%   % Update the properties of all the layers from the NASA servers.
%   nasa = wmsfind('nasa', 'SearchField', 'serverurl');
%   nasa = wmsupdate(nasa, 'AllowMultipleServers', true);
%   
%   See also WebMapServer, WMSFIND, WMSINFO, WMSLayer, WMSREAD.

% Copyright 2009-2011 The MathWorks, Inc.
% $Revision: 1.1.6.5 $  $Date: 2011/03/28 04:29:36 $

% Validate number of inputs.
error(nargchk(1, inf, nargin, 'struct'));

% Validate the input.
validateattributes(layers, {'WMSLayer'}, {'nonempty'}, 'wmsupdate', 'layers')

% Parse the parameter/value parirs and return parsed or default values in
% the structure, OPTIONS.
options = parseParameters(varargin);

% Update the layers.
[updatedLayers, index] = synchronizeLayerArray(layers, options);

%--------------------------------------------------------------------------

function [updatedLayers, index] = synchronizeLayerArray(layers, options)
% Synchronize the layer array, LAYERS, using options specified in the
% scalar structure, OPTIONS.

% Get the unique servers from layers.
servers = layers.servers();

% Verify the number of unique servers.
if numel(servers) > 1 && ~options.AllowMultipleServers
    msgStr = getString(message('map:wms:tooManyServers', numel(servers)));
    error(message('map:wms:useAllowMultipleServers', msgStr));
end

% Initialize the return arguments.
updatedLayers = layers;
index = false(size(layers));

% Loop through each server in the servers cell array, and update the layers
% from the server.
for k=1:numel(servers)
    
    % Assign the serverURL from the cell array.
    serverURL = servers{k};
    
    % Create a server object and set its Timeout property.
    server = WebMapServer(serverURL);
    server.Timeout = secondsToServerTimeUnits(options.TimeoutInSeconds);
    
    % Find all elements that contain the serverURL in the layers array.
    layerMatchesServer = strcmp(serverURL, {layers.ServerURL});
       
    % Obtain the layers belonging to the server.
    inputLayers = layers(layerMatchesServer);
    
    % Find the index location in the layers array that belong to the
    % server.
    location = find(layerMatchesServer);
 
    % Update the server's layers.
    try
        [serverLayers, serverIndex] = server.updateLayers(inputLayers);
        updateFlag = true;
        numLayersRemoved = numel(find(~serverIndex));
        if ~isempty(numLayersRemoved) && numLayersRemoved ~= 0
            warning(message('map:wms:layersNotOnServer', ...
                serverURL, numLayersRemoved))          
        end
    catch e 
        updateFlag = false;
        serverIndex = true(1,numel(inputLayers));
        spaces = isspace(e.message);
        nonSpaces = find(~spaces);
        firstNonSpace = nonSpaces(1);
        lastNonSpace = nonSpaces(end);
        msgStr = e.message(firstNonSpace:lastNonSpace);
        warning(message('map:wms:removedLayers', ...
            serverURL, numel(inputLayers), msgStr))       
    end
               
    % Set the elements of index to true for each element that is updated.
    % location(serverIndex) is the numeric index of each updated element.
    % (This works because serverLayers and location are the same size and
    % have a one-to-one correspondence.)
    index(location(serverIndex)) = updateFlag;
    
    % Update the updatedLayers array with the new elements.
    if updateFlag
        updatedLayers(location(serverIndex)) = serverLayers; 
    else
        updatedLayers(location(serverIndex)) = inputLayers; 
    end
end

% Remove any elements not updated.
updatedLayers = updatedLayers(index);

%--------------------------------------------------------------------------

function timeoutInServerUnits = secondsToServerTimeUnits(timeoutInSeconds)
% Convert seconds to time units of the server (milliseconds).

% The WebMapServer.Timeout property is expressed in milliseconds.
serverTimeUnitsPerSecond = 1000;
timeoutInServerUnits = serverTimeUnitsPerSecond * timeoutInSeconds;

%--------------------------------------------------------------------------

function options = parseParameters(params)
% Parse the parameter/value pairs from the cell array PARAMS and return the
% parsed or default values in the scalar structure OPTIONS.  OPTIONS 
% contains fieldnames specified by PARAMS. If a parameter in PARAMS is
% matched, the value of the parameter is assigned to the fieldname in
% OPTIONS. If the parameter is not matched, the default value for the
% parameter is set.

% Assign the parameterNames and validationFcns cell arrays.
parameterNames = {'AllowMultipleServers', 'TimeoutInSeconds'};
validateFcns = {...
    @(x)validateAllowMultipleServers(parameterNames{1}, x), ...
    @(x)validateTimeoutInSeconds(parameterNames{2}, x)};

% Parse the parameters.
[options, userSupplied, unparsedParams] = ...
    internal.map.parsepv(parameterNames, validateFcns, params{:});

% Set default values.
defaultValues = struct( ...
    parameterNames{1}, false, ...
    parameterNames{2}, 60);
options = setDefaultValues(options, defaultValues, userSupplied);

% If parameter/value pairs have been specified in params, ensure that the
% first element is a string, since parameter/value pairs must be of the
% form: name, value. Otherwise, set unparsedParams to the string 'PARAM1',
% which is used later in the construction of an error message.
if ~isempty(params) && ~ischar(params{1})
    unparsedParams{1} = 'PARAM1';
end

% Check if unparsedParams contained unmatched parameters.
if ~isempty(unparsedParams)
    parameterNames = sprintf('''%s'', ', parameterNames{:});
    error(message('map:validate:invalidParameterName', ...
        unparsedParams{1}, parameterNames(1:end-2)));
end

%--------------------------------------------------------------------------

function options = setDefaultValues(options, defaultValues, userSupplied)
% Set default values to any parameter that is not supplied. OPTIONS and
% defaultValues are each scalar structures and they have the same set of
% fieldnames.

inputFieldNames = fieldnames(options);
for k=1:numel(inputFieldNames)
    name = inputFieldNames{k};
    if ~userSupplied.(name)
        options.(name) = defaultValues.(name);
    end
end

%--------------------------------------------------------------------------

function value = validateAllowMultipleServers(name, value)
% Validate the parameter, 'AllowMultipleServers'.

validateattributes(value, {'logical'}, {'scalar'}, ...
    'validateAllowMultipleServers', name);

%--------------------------------------------------------------------------

function value = validateTimeoutInSeconds(name, value)
% Validate the parameter, 'TimeoutInSeconds'.

validateattributes(value, {'numeric'}, ...
    {'integer', 'finite', 'nonnegative', 'scalar'}, ...
    'validateTimeoutInSeconds', name);

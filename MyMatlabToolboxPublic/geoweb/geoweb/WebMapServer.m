%WebMapServer Web map server object
%
%   A WebMapServer is a handle object representing an implementation of a
%   Web Map Service (WMS) and acts as a proxy to a WMS server on the
%   Internet. The WebMapServer resides physically on the client side. The
%   object can access the server's capabilities document and can perform
%   requests to obtain maps.  It supports multiple WMS versions and
%   negotiates with the server automatically to use the highest known
%   version that the server can support.
%
%   server = WebMapServer(serverURL) constructs a WebMapServer object from
%   the serverURL string parameter. serverURL is a WMS server URL and must
%   include the protocol 'http://' or 'https://'. WebMapServer
%   automatically communicates to the server defined by the serverURL using
%   the highest known version that the server can support.  The serverURL
%   may contain additional WMS keywords.
% 
%   Firewall Note
%   -------------
%   If you need to specify a proxy server to connect to the Internet,
%   select File -> Preferences -> Web and enter your proxy information. Use
%   this feature if you have a firewall.
%
%   WebMapServer properties:
%      Timeout - Server connection timeout
%      EnableCache - Enable caching
%      ServerURL - URL of server
%      RequestURL - URL of server request
%
%   WebMapServer methods:
%      getCapabilities - Get capabilities document from server
%      getMap - Get raster map from server
%      updateLayers - Synchronize WMSLayer object with server
%
%   Example
%   -------
%   % Construct a WebMapServer object that communicates with one of the
%   % Environmental Research Division's Data Access Program (ERDDAP) WMS
%   % servers hosted by NOAA and obtains its capabilities document. 
%   % Search for a server that provides daily, global sea surface 
%   % temperature (SST) data produced by the Jet Propulsion Laboratory's 
%   % Regional Ocean Modeling System (JPL ROMS) group.
%   layers = wmsfind('coastwatch*jpl*sst', 'SearchField', 'serverurl');
%   serverURL = layers(1).ServerURL;
%   server = WebMapServer(serverURL);
%   capabilities = server.getCapabilities();
%   layers = capabilities.Layer;
% 
%   % Obtain and view the data from the server. Show the boundaries of the
%   % nations and the global SST data.
%   nations = layers.refine('nations');
%   sst = layers.refine('sst');
%   request = WMSMapRequest([sst nations], server);
%   A = server.getMap(request.RequestURL);
%   R = request.RasterRef;
%   figure
%   worldmap(A, R)
%   geoshow(A, R)
%   title({sst.LayerTitle, sst.Details.Dimension.Default})
%
%   See also WMSCapabilities, WMSFIND, WMSINFO, WMSLayer, 
%            WMSMapRequest, WMSREAD, WMSUPDATE.

% Copyright 2008-2012 The MathWorks, Inc.
% $Revision: 1.1.6.13.2.1 $  $Date: 2012/01/08 22:04:54 $

classdef WebMapServer <  handle
    
    properties (Access = public)
        
        %Timeout Server connection timeout
        %
        %   Timeout is a double indicating the number of milliseconds 
        %   before a server timeout is issued.  The default value of 0
        %   causes the timeout mechanism to be ignored.
        Timeout = 0;
        
        %EnableCache Enable caching
        %
        % EnableCache is a logical indicating if caching is allowed. If
        % true, the WMSCapabilites object is cached and returned when
        % the method getCapabilities is invoked.  The cache expires if the
        % AccessDate property of the cached WMSCapabilities object is not
        % the current day.  The default value is true.
        EnableCache = true;
                
    end
    
    properties (SetAccess = private, GetAccess = public)
     
        %ServerURL URL of server
        %
        %   ServerURL is a string indicating the URL of the server.
        ServerURL = '';
        
        %RequestURL URL of server request
        %
        %   RequestURL is a string indicating the URL of the last request
        %   to the server.  A RequestURL specifies a request for either
        %   the XML capabilities document or a map.  The RequestURL may be
        %   inserted into a browser.
        RequestURL = '';
          
    end
    
    properties (Access = public, Hidden = true, Dependent = true)
        %CacheFilename Filename for caching map files
        %
        %   CacheFilename is the name of the map file when requesting a map
        %   from the server.
        CacheFilename
    end

    properties (Access = private, Hidden = true)
        %Capabilities Server capabilities document
        %
        %   Capabilities is a WMSCapabilities object indicating the
        %   capabilities of the server.
        Capabilities = WMSCapabilities();
    end
    
    properties (Access = private, Hidden = true)
        %pCacheFilename Private copy of cache filename.
        %
        %   pCacheFilename holds the value of the corresponding dependent
        %   property, CacheFilename.
        pCacheFilename = '';
    end
    
    properties (Access = private, Hidden = true, Dependent = true)            
        %EnableGetMapCache Determines if caching image files
        %
        %   EnableGetMapCache is a logical and is true when caching is
        %   requested for WMS GetMap operations.
        EnableGetMapCache = false;
    end
    
    properties (Access = private, Transient = true, Hidden = true)               
        %WebServer Java WebMapServer accessor
        %
        %   WebServer is an instance of a Java WebMapServer.
        WebServer = []; 
    end
    
    methods       
        function h = WebMapServer(serverURL)
            % WebMapServer constructor
            
            % Verify that Java is available.
            error(javachk('jvm','WebMapServer'));
            
            % Validate the serverURL input.
            internal.map.validateURL('serverURL', serverURL);
            h.ServerURL = serverURL;
            
            % Create the Java WebMapServer.
            h.WebServer = h.getWebServer();
            
            % Initialize the CacheFilename to empty.
            h.CacheFilename = '';
        end
        
        %---------- Property set/get Methods ------------------------------
        
        function set.Timeout(h, timeout)
            validateattributes(timeout, {'numeric'}, ...
                {'integer','finite','nonnegative','scalar'}, ...
                'set', 'Timeout');
            h.Timeout = timeout;
        end
        
        function set.CacheFilename(h, filename)
            validateCacheFilename(filename);
            h.pCacheFilename = filename;
        end
                
        function filename = get.CacheFilename(h)
            if isempty(h.pCacheFilename)
                filename = tempname;
            else
                filename = h.pCacheFilename;
            end
        end           
        
        function value = get.EnableGetMapCache(h)
            if isempty(h.pCacheFilename)
                value = false;
            else
                value = h.EnableCache;
            end
        end
        
        function set.EnableCache(h, enableCache)
            validateattributes(enableCache, {'logical'}, {'scalar'}, ...
                'set', 'EnableCache');
            h.EnableCache = enableCache;
        end
        
        %-------------- Public Methods ------------------------------------

        function capabilities = getCapabilities(h)
        %getCapabilities Get capabilities document from server
        %
        %   capabilities = server.getCapabilities() retrieves the
        %   capabilities document from the server as a WMSCapabilities
        %   object and updates the RequestURL property. getCapabilities
        %   accesses the Internet to retrieve the document. An Internet
        %   connection must be established to use the method. The WMS
        %   server may periodically be unavailable and several minutes may
        %   elapse before the document is retrieved.
        %
        %   Example 1
        %   ---------
        %   % Retrieve the capabilities document from the 
        %   % NASA SVS Image Server.
        %   nasa = wmsfind('NASA SVS Image', 'SearchField', 'servertitle');
        %   serverURL = nasa(1).ServerURL;
        %   server = WebMapServer(serverURL);
        %   capabilities = server.getCapabilities;
            
            % Validate as a scalar object.
            assert(isscalar(h), message('map:validate:expectedScalarObject', 'server'));
            
            % Determine if the server's capabilities document has been
            % previously retrieved, is valid, and caching is enabled.
            if haveValidCache(h)
                % Return the cached capabilities document.
                capabilities = h.Capabilities;
            else
                
                % Issue the GetCapabilities request and obtain the response
                % from the server.
                response = issueGetCapabilitesRequest(h);
                
                % Update the RequestURL property.
                h.RequestURL = char(h.getWebServer.getRequestURL);
                
                % Create a WMSCapabilities object for the return argument.
                capabilities = ...
                    WMSCapabilities(h.ServerURL, response);
                
                % Update the internal capabilities property.
                h.Capabilities = capabilities;
            end
        end
              
        function A = getMap(h, mapRequestURL)
        %getMap Get raster map from server
        %
        %   A = server.getMap(mapRequestURL) dynamically renders and
        %   retrieves from the server a map, which is defined by the
        %   parameters in the URL, mapRequestURL. The map is rendered as a
        %   color or grayscale image and stored in A.  The RequestURL
        %   property is updated. getMap accesses the Internet to retrieve
        %   the map. An Internet connection must be established to use the
        %   method. The WMS server may periodically be unavailable and
        %   several minutes may elapse before the map is retrieved.
        %
        %   Example 1
        %   ---------
        %   % Retrieve a map of the Blue Marble global mosaic layer 
        %   % from the NASA Earth Observations WMS server. 
        %   neowms = wmsfind('neowms', 'SearchField', 'serverurl');
        %   layer = neowms.refine('bluemarbleng', 'MatchType', 'exact');
        %
        %   % Obtain and display the map.
        %   server = WebMapServer(layer.ServerURL);
        %   mapRequest = WMSMapRequest(layer, server);
        %   A = server.getMap(mapRequest.RequestURL);
        %   R = mapRequest.RasterRef;
        %   figure
        %   worldmap world
        %   geoshow(A, R)
        %   setm(gca,'MLabelParallel',-90,'MLabelLocation',90)
        %   title(layer.LayerTitle)
            
            % Validate as a scalar object.
            assert(isscalar(h), message('map:validate:expectedScalarObject', 'server'));
            
            % Validate the mapRequestURL input.
            internal.map.validateURL('mapRequestURL', mapRequestURL);
     
            % Set the RequestURL property of this class. 
            h.RequestURL = mapRequestURL;            
               
            % Assign a temporary filename for the raster map.
            filename = h.CacheFilename;
            
            % Use the onCleanup object to delete the file when finished or
            % if an error is issued.
            deleteFile = onCleanup( ...
                @(x, y)deleteGetMapFile(filename, h.EnableGetMapCache));
            
            % Issue the GetMap request and obtain the response filename.
            issueGetMapRequest(h, filename);
            
            % Read the map from the file.
            A = readGetMapFile(filename, h.RequestURL);           
        end

        function [updatedLayers, index] = updateLayers(h, layers)
        %updateLayers Synchronize WMSLayer object with server
        %
        %   [updatedLayers, index] = server.updateLayers(layers) returns a
        %   WMSLayer array with its properties synchronized with values
        %   from the server. LAYERS must contain only one unique ServerURL.
        %   Layers no longer available on the server are removed. INDEX is
        %   a logical array containing true for each layer that was
        %   available such that updatedLayers has the same size as
        %   LAYERS(INDEX). updateLayers accesses the Internet to update the
        %   properties. An Internet connection must be established to use
        %   the method. The WMS server may periodically be unavailable and
        %   several minutes may elapse before the properties are updated.
        %
        %   Example 1
        %   ---------
        %   % Update the properties of a MODIS global mosaic layer
        %   % obtained from the NASA Earth Observations WMS server. 
        %   modis = wmsfind('modis');
        %   modis = modis.refine('bluemarbleng');
        %   modis = modis(1);
        %
        %   % Create a WebMapServer object.
        %   server = WebMapServer(modis.ServerURL);
        %   
        %   % Update the properties of the modis layer.
        %   updatedLayer = server.updateLayers(modis);
        %
        %   % View the metadata of the layer.
        %   metadata = urlread(updatedLayer.Details.MetadataURL);
        %   disp(metadata) 
        %
        %   % Obtain and display the map.
        %   mapRequest = WMSMapRequest(updatedLayer, server);
        %   A = server.getMap(mapRequest.RequestURL);
        %   R = mapRequest.RasterRef;
        %   figure
        %   worldmap world
        %   geoshow(A, R)
        %   setm(gca,'MLabelParallel',-90,'MLabelLocation',90)
        %   title('MODIS Global Mosaic')
        %
        %   Example 2
        %   ---------
        %   % Update the properties of layers from multiple servers.
        %   % Find layers with the name 'terra' in the server URL.
        %   terra = wmsfind('terra', 'SearchField', 'serverurl');
        %
        %   % Find the layers for an individual server in the terra layer, 
        %   % update their properties, and append them to the
        %   % updatedTerraLayers array.
        %   servers = terra.servers();
        %   updatedTerraLayers = [];
        %   fprintf('Updating layer properties from %d servers.\n', ...
        %      numel(servers));
        %   for k=1:numel(servers)
        %      serverLayers = ...
        %         terra.refine(servers{k}, 'SearchField', 'serverurl');
        %      fprintf('Updating properties from server %d:\n%s\n', ...
        %         k, serverLayers(1).ServerURL);
        %      server = WebMapServer(serverLayers(1).ServerURL);
        %      layers = server.updateLayers(serverLayers);
        %      % Grow using concatenation because layers can have any
        %      % length ranging from 0 to numel(serverLayers).
        %      updatedTerraLayers = [updatedTerraLayers; layers]; 
        %   end 
            
            % Validate as a scalar object.
            assert(isscalar(h), message('map:validate:expectedScalarObject', 'server'));
            
            % Validate input layers.
            validateLayer(layers, h.ServerURL);
   
            % Obtain all the available specifications.
            url = java.net.URL(h.ServerURL);
            specs = ...
                com.mathworks.toolbox.geoweb.wms.WMSClientSpecification(url);
            
            % Update the layers array if required.
            needUpdate = false;
            k  = 1;
            while ~needUpdate && k <= numel(layers)
               needUpdate = layerNeedsUpdate(layers(k), specs);
               k = k + 1;
            end
            if needUpdate
                % Get the capabilities document (via the cache if the
                % document has already been retrieved).
                capabilities = getCapabilities(h);
            
                % Determine the server host.
                host = char(url.getHost());
                
                % Determine the unique property name. If the server is
                % found in layerTitleServers, then use the LayerTitle as
                % the unique layer property; otherwise use LayerName.
                layerTitleServers = ...
                    {'columbo.nrlssc.navy.mil', 'dmap.nrlssc.navy.mil'};
                if any(strncmp(host, layerTitleServers, numel(host)))
                    uniqueName = 'LayerTitle';
                else
                    uniqueName  = 'LayerName';
                end
                
                % Find the matching members based on uniqueName.
                [~, layerLoc, capsLayerLoc] = intersect(...
                    {layers.(uniqueName)}, ...
                    {capabilities.Layer.(uniqueName)});
                
                % Assign the return updatedLayers based on the matching 
                % location index arrays, preserving the order of the input
                % layers array. Note that the location return values from
                % intersect are such that (if ~ is replaced with C in the
                % above use):
                %    C = {layers(layerLoc).LayerName}
                %    C = capabilities.LayerNames(capsLayerLoc)
                % therefore, 
                %    layers(layerLoc) = capabilities.LayerName(capsLayerLoc)
                % thus preserving the order of layers.
                updatedLayers = layers;
                updatedLayers(layerLoc) = capabilities.Layer(capsLayerLoc);
                
                % Setup the index values. The values in layerLoc specify
                % the locations in layers that match the capabilities layer
                % names. Assign true for these locations.
                index = false(size(layers));
                index(layerLoc) = true;
                
                % Remove layers no longer on the server.
                updatedLayers(~index) = [];
            else
                index = true(size(layers));
                updatedLayers = layers;
            end
        end       
    end
    
    methods (Access = private)
                
        function webServer = getWebServer(h)
        % Get the Java WebMapServer instance.
            
            if isempty(h.WebServer)
                url = java.net.URL(h.ServerURL);
                h.WebServer = ...
                    com.mathworks.toolbox.geoweb.wms.WebMapServer(url);
            end
            webServer = h.WebServer;
        end
        
        %------------------------------------------------------------------
        
        function setURLReadTimeout(h)
        % Set the timeout value for URL connections via the Java 
        % WebMapServer.
            
            if ~isempty(h.getWebServer) && h.Timeout ~= 0                
                h.WebServer.setTimeout(h.Timeout);
            end
        end
        
        %------------------------------------------------------------------
        
        function tf = haveValidCache(h)
        % Determine if a valid cache is available.  The cache expires if
        % the date of the last access is not the current day.
            
            haveCache = ~isempty(h.Capabilities.ServerURL);
            today = datestr(floor(now));
            cacheExpired = ~isequal(today, h.Capabilities.AccessDate);
            tf = h.EnableCache && haveCache && ~cacheExpired;
        end

        %------------------------------------------------------------------
        
        function response = issueGetCapabilitesRequest(h)
        % Issue a GetCapabilities request to the server.
            
            % Set the URL read timeout value on the server.
            setURLReadTimeout(h);
            
            % Initiate a request to the server for the capabilities
            % document. If successful, the wmsCapabilities is a Java
            % WMSCapabilities object.
            try
                % Get the Java capabilities object.
                response = h.getWebServer.getCapabilities();
            catch e
                % Even though an error occurred, the server has updated the
                % RequestURL.
                h.RequestURL = char(h.getWebServer.getRequestURL);

                % An error occurred when attempting to obtain the
                % capabilities document from the server.        
                useServerResponse = true;
                webServer = h.getWebServer();
                e = createServerException(...
                    webServer, useServerResponse, h.ServerURL);
                throw(e);                
            end
        end
           
        %------------------------------------------------------------------
        
        function issueGetMapRequest(h, filename)
        % Issue a GetMap request to the server.
            
            % Set the URL read timeout value on the server.
            setURLReadTimeout(h);
            
            % Obtain the request URL.
            url = java.net.URL(h.RequestURL);
            
            try
                h.getWebServer.getMap(url, filename);
            catch e
                % An error occurred when attempting to obtain the map from
                % the server.
                useServerResponse = false;
                webServer = h.getWebServer();
                e = createServerException(...
                    webServer, useServerResponse, h.ServerURL);
                throw(e);
            end
        end
    end
end

%-------------------- Validation Functions --------------------------------

function validateCacheFilename(filename)
% Validate filename. The file must have write access.

validateattributes(filename, {'char'}, {}, ...
    'validateCacheFilename', 'CacheFilename');
if ~isempty(filename)
    try
        fid = fopen(filename,'w');
        fclose(fid);
        delete(filename);
    catch e 
        error(message('map:wms:noWriteAccessForCacheFile', filename));
    end
end
end

%--------------------------------------------------------------------------

function validateLayer(layer, serverURL)
% Validate layer as a WMSLayer array containing one unique ServerURL and
% update the properties if required.

% Validate layer as a WMSLayer array and to contain at least one element.
validateattributes(layer, {'WMSLayer'}, {'nonempty'}, 'WebMapServer', 'layer')

% Obtain the list of servers.
layerURL = layer.servers();

% Verify the number of unique servers.
assert(numel(layerURL) == 1, message('map:wms:tooManyServers', numel(layerURL)));

% Validate the layer ServerURL with the input server's URL.
assert(isequal(serverURL, layer(1).ServerURL), message('map:wms:mismatchedServerURL', ...
    'layer.ServerURL', layer(1).ServerURL, 'server.ServerURL', serverURL));
end

%--------------------------------------------------------------------------

function needUpdate = layerNeedsUpdate(layer, specs)
% Determine if a layer needs updating. LAYER is a scalar WMSLayer array.
 
% Check if the properties need updating.
if ~propertiesNeedUpdate(layer)
    % Check if Details field names or values need updating.
    needUpdate = detailsNeedUpdate(layer, specs);
else
    % The properties need updating.
    needUpdate = true;
end
end

%--------------------------------------------------------------------------

function needUpdate = detailsNeedUpdate(layer, specs)
% Check if the Details field names or values require updating.

% Need to update the WMSLayer array if all fields of Details are not
% defined. This can occur if a WMSLayer array is created by a user with
% some fields of Details missing.
needUpdate = ~isstruct(layer.Details);
l = WMSLayer();
requiredFields = fieldnames(l.Details);
needUpdate = needUpdate || ...
    dataRequiresUpdate(layer.Details, requiredFields);

% Examine all the field values of Details.
[detailsFieldName,detailsFieldValues] = assignDetailsFields(l.Details);
k = 1;
while ~needUpdate && k <= numel(detailsFieldName)
    needUpdate = dataRequiresUpdate(...
        layer.Details.(detailsFieldName{k}), detailsFieldValues{k});
    k = k + 1;
end

% Validate ImageFormats
needUpdate = needUpdate || ~iscell(layer.Details.ImageFormats);

% Validate Version
needUpdate = needUpdate || ...
    ~ischar(layer.Details.Version) || ...
    isempty(layer.Details.Version) || ...
    ~specs.containsVersion(layer(1).Details.Version);
end

%--------------------------------------------------------------------------

function [detailsFieldName, detailsFieldValues] = ...
    assignDetailsFields(Details)
% Create detailsFieldName, a cell array containing the required field
% names of the Details structure and detailsFieldValues, a cell array 
% containing the field names of each of the structure elements of Details.

% Assign field names of the structure elements of Details.
names = fieldnames(Details);
isStructName = false(size(names));
for k=1:numel(names)
    isStructName(k) = isstruct(Details.(names{k}));
end
detailsFieldName = names(isStructName);
 
% Assign the fieldnames of these elements.
detailsFieldValues = cell(size(detailsFieldName));
for k=1:numel(detailsFieldName)
    detailsFieldValues{k} = fieldnames(Details.(detailsFieldName{k}));
end
end

%--------------------------------------------------------------------------

function needUpdate = propertiesNeedUpdate(layer)
% Determine if the WMSLayer array, LAYER, contains the update string in an
% any of the properties.

l = WMSLayer;
updateString = l.UpdateString; % '<Update using WMSUPDATE>';
props = properties(layer);
needUpdate = false;
for k=1:numel(props)
    thisProperty = layer(1).(props{k});
    if ischar(thisProperty) && isequal(thisProperty, updateString)
        needUpdate = true;
        return
    end
end
end

%--------------------------------------------------------------------------

function needUpdate = dataRequiresUpdate(testData, requiredFields)
% Return true if testData is not a structure containing the requiredFields
% names.

if ~isstruct(testData)
    needUpdate = true;
else
   needUpdate = ~isequal(sort(requiredFields), sort(fieldnames(testData)));
end
end

%--------------------------------------------------------------------------
 
function exception = ...
    createServerException(webServer, useServerResponse, serverURL)
% Construct an appropriate MException based on the content of the Java
% Exception and the response string from the server.

% Obtain the Java Exception from the server.
exception = webServer.getException();
if ~isempty(exception)
    
    % Retrieve the error message from the server.
    msgString = getServerResponseMsg(webServer, useServerResponse);
    
    % Get the error ID mnemonic.
    mnemonic = getMnemonic(exception);
    
    % Create the MException object based on the error message and mnemonic.
    exception = createResponseException(msgString, mnemonic, serverURL);  
else
    % A network connection error occurred. This error can occur if an
    % unknown exception was thrown while copying the data. For example, the
    % error may be issued when switching the proxy host to an invalid
    % value.
    msg = message('map:wms:connectionFailed', serverURL);
    exception = MException(msg.Identifier, '%s', msg.getString());
end
end

%--------------------------------------------------------------------------

function response = getServerResponseMsg(webServer, useServerResponse)
% Obtain the error response from the server.

% Obtain the response from the server.
if useServerResponse
    response = char(webServer.getResponse());
    if isempty(response)
        response = char(webServer.getException().getMessage());
    end        
else
    % Obtain the response from the exception.
    response = char(webServer.getException().getMessage());
end
end

%--------------------------------------------------------------------------

function exception = createResponseException(response, mnemonic, serverURL)
% Construct an MException based on the response and mnemonic.

% Construct the error message. If the response string from the server
% contains a ServiceException, then use the value of the ServiceException
% as the error message; otherwise, use the Java Exception message. Update
% the mnemonic if a service exception is found.
[serviceException, mnemonic] = createServiceException(mnemonic, response);
switch mnemonic
    case 'ServiceException'
        exception = serviceException;
        
    case 'SAXParseException'
        exception = createSAXParseException(serverURL, response);
        
    case 'UnknownHostException'
        exception = createUnknownHostException(serverURL);
        
    case 'NullPointerException'
        exception = createNullPointerException(response);
        
    otherwise
        id = ['map:WebMapServer:' mnemonic];
        exception = MException(id, '%s', response);
        
end
end

%--------------------------------------------------------------------------

function [exception, mnemonic] = createServiceException(mnemonic, response)
% Obtain the service exception if found in the response. If found, then
% update the mnemonic string.

msgString = getSimpleElementValue('ServiceException', response);
if ~isempty(msgString) && ~all(isspace(msgString))  
    mnemonic = 'ServiceException';
end
id = ['map:WebMapServer:' mnemonic];
exception = MException(id, '%s', msgString);

end

%--------------------------------------------------------------------------

function exception = createSAXParseException(serverURL, response)
% Create a SaxParseException error message.

% The string returned from the server could not be parsed into a document.
% One possible reason for this condition is that the server is a normal
% HTTP server, but not a WMS server. In this case, the error message
% contained in the Exception class does not contain useful information for
% the user. Let the WMSCapabilities class attempt to parse the response in
% order to construct a meaningful error message.
try
    % Use the WMSCapabilities class to generate a meaningful message.
    WMSCapabilities(serverURL, response);
    msgString = '';
catch e
    msgString = e.message;
end

if isempty(msgString)
    msgString = getString(message('map:wms:internalParseError'));
end
    
% Construct the exception.
id = 'map:WebMapServer:SAXParseException';
exception = MException(id, '%s', msgString);

end

%--------------------------------------------------------------------------

function exception = createUnknownHostException(serverURL)
% Create an UnknownHostException message.

url = java.net.URL(serverURL);  
id = 'map:WebMapServer:UnknownHostException';
exception = MException(id, '%s', getString(message( ...
    'map:wms:unknownHostException', char(url.getHost()))));
end

%--------------------------------------------------------------------------

function exception = createNullPointerException(response)
% Obtain the exception if found in the response. 

% At this point, we know that a NullPointerException has occurred. If
% available, include more specific information in MException; otherwise,
% use a generic message string.
if ~isempty(response)
    msgString = getSimpleElementValue('Exception', response);
    if isempty(msgString)
        msgString = response;
    end
else
    msgString = getString(message('map:wms:unknownParseError'));   
end

% Construct exception.
id = 'map:WebMapServer:NullPointerException';
exception = MException(id, '%s', msgString);

end

%--------------------------------------------------------------------------

function mnemonic = getMnemonic(exception)
% Get the mnemonic from the Java exception.

mnemonic = char(exception.getClass());
index = strfind(mnemonic,'.');
mnemonic = mnemonic(index(end)+1:end);

end
    
%--------------------------------------------------------------------------

function A = readGetMapFile(filename, requestURL)
% Read the image from the map file specified by filename.

% Determine if IMREAD or MULTIBANDREAD is required to read the image file.
if ~isMultibandFormat(requestURL)
    % The image format is supported by IMREAD.
    A = readImageFormat(filename, requestURL);
else
    % The image format is BIL and requires MULTIBANDREAD.
    A = readMultibandFormat(filename, requestURL);
end
end

%--------------------------------------------------------------------------

function A = readImageFormat(filename, requestURL)
% Read the image from the file using IMREAD.

try
    % Read the image from the file.
    [A, colorMap] = imread(filename);
    
    % If A is an indexed image, then convert it to RGB.
    if ~isempty(colorMap)
        A = ind2rgb(A, colorMap);
    end
catch e
    if isequal(e.identifier, 'MATLAB:imagesci:imread:fileFormat')
        issueReadGetMapError(filename, requestURL);
    else
        % An unknown error occurred while reading the URL. A file has not
        % been created or the file is corrupt.
        error(message('map:wms:unknownGetMapError', e.message));
    end
end
end

%--------------------------------------------------------------------------

function Z = readMultibandFormat(filename, requestURL)
% Read the data from the file using MULTIBANDREAD.

% Obtain the image information.
[imageSize, precision] = multibandinfo(filename, requestURL);

% Read the file using multibandread.
try   
    Z = multibandread(filename, imageSize, precision, 0, 'bil', 'ieee-le');
catch e
    % An unknown error occurred while reading the URL. A file has not
    % been created or the file is corrupt.
    error(message('map:wms:readMultibandGetMapError', e.message));
end
end

%--------------------------------------------------------------------------

function [imageSize, precision] = multibandinfo(filename, requestURL)
% Obtain the image size and precision of the multi-band file.

% Obtain the number of rows and columns from the parameters in the
% requestURL.
S = parseMapRequestURL(requestURL, 'HEIGHT', 'WIDTH');
rows = str2double(S.HEIGHT);
cols = str2double(S.WIDTH);

% Calculate the expected image size and precision.
d = dir(filename);
numBytes = d.bytes;
numPixels = rows * cols;
switch numBytes
    case numPixels*2
        precision = 'int16';
        bands = 1;
        
     case numPixels*6  % (rows*cols*2_bytes*3_bands)
        precision = 'int16';
        bands = 3;
       
    case numPixels*4
        precision = 'int32';
        bands = 1;
 
     case numPixels*12  % (rows*cols*4_bytes*3_bands)
        precision = 'int32';
        bands = 3;
        
    otherwise
        issueReadGetMapError(filename, requestURL);
end

% Set precision such that the output precision matches input.
precision = [precision '=>' precision];
imageSize = [rows, cols, bands];
end

%--------------------------------------------------------------------------

function tf = isMultibandFormat(requestURL)
% Determine if multibandread is required to read the image file.

% The URL may be encoded such that the FORMAT value  equals "image%2F"
% (URL encoded /) followed by the name of the image format. Decode the URL
% string to replace '%2F' with '/'
url = char(java.net.URLDecoder.decode(requestURL,'UTF-8'));

% If the URL contains image/bil format then multibandread is required.
bilFormat = 'format=image/bil';
tf = ~isempty(regexpi(url, bilFormat));
end

%--------------------------------------------------------------------------

function issueReadGetMapError(filename, requestURL)
% Issue an error based on the contents of the file specified by filename.

d = dir(filename);
if d.bytes == 0
    error(message('map:wms:emptyGetMapFile'));
else
    % The file contains the text message from the URL read. Use the
    % text as an error message.
    fid = fopen(filename, 'r');
    response = fread(fid, 'uint8=>char')';
    fclose(fid);
    
    fileInfo = finfo(filename);
    if isequal(fileInfo,'html') || ~isempty(regexpi(response, '<HTML'))
        error(message('map:wms:htmlGetMapFile', requestURL, requestURL));
    end
    
    msgString = getSimpleElementValue('ServiceException', response);
    if ~isempty(msgString) && ~all(isspace(msgString))
        error('map:WebMapServer:GetMapServiceException', msgString)
    else
        error(message('map:wms:unknownGetMapError', response));
    end
end
end

%--------------------------------------------------------------------------

function deleteGetMapFile(filename, cacheIsEnabled)
% Delete filename if it exists and if caching is not enabled.

if exist(filename,'file') && ~cacheIsEnabled
    delete(filename)
end
end

%------------------------ XML Parse Functions -----------------------------

function value = getSimpleElementValue(elementName, xmlString)
% Get value from a simple XML element.
%
% Find a simple value from an XML element string.
% The XML element must be in the form:
%   <elementName> value </elementName>
% or in the form:
%   <elementName other characters> value </elmentName>
%
% elementName is a string and must not contain '<' '/' or '>'.
% xmlString must be a 1-by-M cell array or string. Input validation is not
% performed. It is assumed that xmlString contains valid XML entities.
%
% getSimpleElement returns the value of the element as a string. If the
% element is not found then '' is returned. If more than one elementName is
% found, then all the values are returned.

% Add the required < > characters around elementName.
xmlStartElement = ['<'  elementName '>'];
xmlEndElement   = ['</' elementName '>'];

% Find the xmlStartElement
startIndex = strfind(xmlString, xmlStartElement);
if isempty(startIndex)
    % xmlStartElement is not found.
    % Try finding the xmlStartElement without the ending '>'. This form
    % allows for extra characters between the elementName and the ending
    % bracket, for example <elementName ATTRIBUTES>
    [xmlStartElement, startIndex, malformedXML] = ...
        findAttributeElements(xmlString, xmlStartElement);
    if malformedXML
        value = '';
        return
    end  
end

% The xmlStartElement has been found. Return all values between any
% xmlStartElement and xmlEndElement.
value = ...
    findXMLValues(xmlString, xmlStartElement, xmlEndElement, startIndex);

% Concatenate all the values together, adding an extra space between
% multiple elements.
value = sprintf('%s ',value{:});
end

%--------------------------------------------------------------------------

function [xmlStartElement, startIndex, malformedXML] = ...
    findAttributeElements(xmlString, xmlStartElement)
% Convert the xmlStartElement from a single element name to an element name
% followed by attributes.  startIndex contains the indices in xmlString
% that contain xmlStartElement.

xmlStartElement = [xmlStartElement(1:end-1) ' '];
startIndex = strfind(xmlString, xmlStartElement);
if isempty(startIndex)
    % The element was not found. The XML string is malformed.
    malformedXML = true;
    return
else
    % The xmlStartElement was found excluding the ending '>'.
    % Find all occurrences of '<elementName other chars>'.
    malformedXML = false;
    endBracket = '>';
    allIndex = strfind(xmlString, endBracket);
    if ~isempty(allIndex)
        endIndex = double(size(startIndex));
        xmlStartElement = cell(size(startIndex));
        for k=1:numel(startIndex)
            index = allIndex(allIndex>startIndex(k));
            if ~isempty(index)
                endIndex(k) = index(1);
                xmlStartElement{k} = xmlString(startIndex(k):endIndex(k));
            else
                % The XML string is malformed.
                malformedXML = true;
                return
            end
        end
    else
        % The XML string is malformed. 
        malformedXML = true;
        return
    end
end
end

%--------------------------------------------------------------------------

function value =...
    findXMLValues(xmlString, xmlStartElement, xmlEndElement, startIndex)
% Find the xmlString in between xmlStartElement and xmlEndElement using the
% startIndex indices. VALUE is returned as a cell array.

value = cell(size(startIndex));
for k = 1:numel(startIndex)
    thisString = xmlString(startIndex(k):end);
    endIndex = strfind(thisString, xmlEndElement);
    if ~isempty(endIndex)
        endIndex = endIndex(1);
        thisString = thisString(1:endIndex);
        if iscell(xmlStartElement)
            % Found second form:
            % <elementName other characters> value </elmentName>
            startValue = 1 + numel(xmlStartElement{k});
        else
            % Found first form:
            % <elementName> value </elementName>
            startValue = 1 + numel(xmlStartElement);
        end
        % Use end-1 since the last character is '<'
        value{k} = thisString(startValue:end-1);
    else
        % Unable to find xmlEndElement
        value{k} = '';
    end
end
end

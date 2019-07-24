function [A, R, mapRequestURL] = wmsread(layerOrURL, varargin)
%WMSREAD Retrieve WMS map from server
%
%   WMSREAD accesses the Internet to dynamically render and retrieve a map
%   from a Web Map Service (WMS) server. An Internet connection must be
%   established to use the function. Note that a WMS server may be
%   unavailable, and several minutes may elapse before the map is returned.
%
%   WMSREAD communicates with the server using a WebMapServer handle object
%   representing a Web Map Service (WMS). The handle object acts as a proxy
%   to a WMS server and resides physically on the client side. The handle
%   object retrieves the map from the server and automatically times-out
%   after 60 seconds if a connection is not made to the server.
%
%   [A, R] = WMSREAD(LAYER) dynamically renders and retrieves a raster map
%   from a WMS server. The server is specified by the ServerURL property of
%   the WMSLayer object, LAYER. If LAYER has more than one element, then
%   each subsequent layer is overlaid on top of the base (first) layer,
%   forming a single image. You can render multiple layers like this only
%   if they all share the same ServerURL value. The raster map is returned
%   in the output A, which will contain either a color or grayscale image.
%   LAYER must contain either the string 'EPSG:4326' or 'CRS:84' in the
%   CoordRefSysCodes property.
%
%   The second output is a referencing matrix R, which ties A to the
%   EPSG:4326 geographic coordinate system. EPSG:4326 is based on the WGS84
%   (1984 World Geodetic System) datum, with latitude and longitude in
%   degrees and longitude referenced to the Greenwich Meridian. The rows of
%   A are aligned with parallels, with even sampling in longitude. Likewise
%   the columns of A are aligned with meridians, with even sampling in
%   latitude. The geographic limits of A span the full latitude and
%   longitude extent of LAYER. The larger spatial size of A is chosen to
%   match its larger geographic dimension and is fixed at the value 512. In
%   other words, assuming RGB output, A is 512-by-N-by-3 if the latitude
%   extent exceeds longitude extent and N-by-512-by-3 otherwise, and in
%   both cases N <= 512. N is set to the integer value that provides the
%   closest possible approximation to equal cell sizes in latitude and
%   longitude, under the constraint (described earlier) that the map spans
%   the full extent supported for the LAYER.
%    
%   [A, R] = WMSREAD(LAYER, PARAM1, VAL1, PARAM2, VAL2, ...) specifies
%   parameter-value pairs that modify the request to the server. Parameter
%   names can be abbreviated and are case-insensitive. See the table below
%   for a list of available parameters.
%
%   [A, R] = WMSREAD(mapRequestURL) uses the one and only input argument,
%   mapRequestURL, to fully define the request to the server. The
%   mapRequestURL string is composed of a WMS serverURL with additional WMS
%   parameters. The URL string must include the WMS parameters, BBOX and
%   GetMap, and must contain either the EPSG:4326 or CRS:84 keyword. A
%   mapRequestURL may be obtained from the output of WMSREAD, or the
%   RequestURL property of a WMSMapRequest object, or from an Internet
%   search.
%
%   [A, R, mapRequestURL] = WMSREAD(...) returns a WMS GetMap request URL
%   in the string mapRequestURL.  The mapRequestURL can be inserted into a
%   browser to make a request to a server, which then returns the raster
%   map. The browser will open the returned map if its mime type is
%   understood, or it will save the raster map to disk.
%
%   Optional Parameters
%   -------------------
%
%   'Latlim'          A two-element vector specifying the latitude limits
%                     of the output image in the form [southern_limit
%                     northern_limit]. The limits are in degrees and must
%                     be ascending. By default, 'Latlim' is empty and the
%                     LAYER's full extent in latitude used. If
%                     LAYER.Details.Attributes.NoSubsets is true, then
%                     'Latlim' may not be modified.
%
%   'Lonlim'          A two-element vector specifying the longitude limits
%                     of the output image in the form [western_limit
%                     eastern_limit]. The limits are in degrees and must be
%                     ascending. By default, 'Lonlim' is empty and the
%                     LAYER's full extent in longitude used. If
%                     LAYER.Details.Attributes.NoSubsets is true, then
%                     'Lonlim' may not be modified.
%
%  'ImageHeight'      A scalar, positive, integer-valued number specifying 
%                     the desired height of the raster map in pixels.
%                     ImageHeight cannot exceed 8192. If
%                     LAYER.Details.Attributes.FixedHeight contains a
%                     positive number, then 'ImageHeight' may not be
%                     modified.
%
%  'ImageWidth'       A scalar, positive, integer-valued number specifying 
%                     the desired width of the raster map in pixels.
%                     ImageWidth cannot exceed 8192. If
%                     LAYER.Details.Attributes.FixedWidth contains a
%                     positive number, then 'ImageWidth' may not be
%                     modified.
%
%  'CellSize'         A scalar or two-element vector specifying the target 
%                     size of the output pixels (raster cells) in units of
%                     degrees. If a scalar is specified, the value applies
%                     to both height and width dimensions. If a vector is
%                     specified, it is given in the form [height width]. An
%                     error is issued if both CellSize and ImageHeight or
%                     ImageWidth are specified. The output raster map must
%                     not exceed a size of [8192, 8192].
%
%  'RelTolCellSize'   A scalar or two-element vector specifying the
%                     relative tolerance for 'CellSize'. If a scalar is
%                     specified, the value applies to both height and width
%                     dimension. If a vector is specified, the tolerance is
%                     given in the order [height width]. The default value
%                     is .001.
%
%  'ImageFormat'      A string specifying the desired image format for use
%                     in rendering the map as an image. If specified, the
%                     format must match an entry in the
%                     LAYER.Details.ImageFormats cell array and must match
%                     one of the following supported formats:
%                       'image/jpeg',  'image/gif',     'image/png',
%                       'image/tiff',  'image/geotiff', 'image/geotiff8',
%                       'image/tiff8', 'image/png8',    'image/bil'.
%                     If not specified, the format defaults to the first
%                     available format in the supported format list. When
%                     the 'image/bil' format is specified, A is returned as
%                     a two-dimensional array with a class type of int16 or
%                     int32.
%
%  'StyleName'        A string or cell array of strings, specifying the 
%                     style to use when rendering the image. By default the
%                     style is set to the empty string. The StyleName must
%                     be a valid entry in the LAYER.Details.Style.Name
%                     field. If multiple layers are requested, each with a
%                     different style, then StyleName must be a cell array
%                     of strings.
%
%  'Transparent'      A logical specifying if transparency is enabled. When
%                     Transparent is true, all of the pixels not
%                     representing features or data values are set to a
%                     transparent value. When Transparent is false,
%                     non-data pixels are set to the value of the
%                     background color. By default, the value is false.
%
%  'BackgroundColor'  A three-element vector specifying the color to be
%                     used as the background (non-data) pixels of the map.
%                     If not specified, the default is [255,255,255].
%
%  'Elevation'        A string indicating the desired elevation extent
%                     of the requested map. The layer must contain
%                     elevation data which is indicated by the 'Name' field
%                     of the LAYER.Details.Dimension structure. The 'Name'
%                     field must contain the value 'elevation'. The
%                     'Extent' field of the LAYER.Details.Dimension
%                     structure determines the permissible range of values
%                     for the parameter.
%
%  'Time'             A string or a numeric date number indicating the 
%                     desired time extent of the requested map. The layer
%                     must contain data with a time extent which is
%                     indicated by the 'Name' field of the
%                     LAYER.Details.Dimension structure. The 'Name' field
%                     must contain the value 'time'. The 'Extent' field of
%                     the LAYER.Details.Dimension structure determines the
%                     permissible range of values for the parameter. (For a
%                     complete description of the permissible values of
%                     time, please refer to the reference page.)
%
%  'SampleDimension'  A two-element cell array of strings indicating the
%                     name of a sample dimension (other than 'time' or
%                     'elevation') and its string value. The layer must
%                     contain data with a sample dimension extent which is
%                     indicated by the 'Name' field of the
%                     LAYER.Details.Dimension structure. The 'Name' field
%                     must contain the value of the first element of 
%                    'SampleDimension'. The 'Extent' field of the
%                     LAYER.Details.Dimension structure determines the
%                     permissible range of values for the second element of
%                     'SampleDimension'.
%
%  'TimeoutInSeconds' An integer-valued scalar double indicating the number
%                     of seconds to elapse before a server timeout is
%                     issued. By default, the value is 60 seconds. A value
%                     of 0 causes the timeout mechanism to be ignored.
%
%   Firewall Note 
%   ------------- 
%   If you need to specify a proxy server to connect to the Internet,
%   select File -> Preferences -> Web and enter your proxy information. Use
%   this feature if you have a firewall.
%
%   Example 1 - Basic
%   ------------------
%   % Read and display a Blue Marble Next Generation layer from NASA.
%   nasa = wmsfind('nasa', 'SearchField', 'serverurl');
%   layer = nasa.refine('bluemarbleng',  'SearchField', 'layername', ...
%      'MatchType', 'exact');
%   [A, R] = wmsread(layer(1));
%   figure
%   axesm globe
%   axis off
%   geoshow(A, R)
%   title('Blue Marble')
%
%   Example 2 - Latlim, Lonlim, CellSize
%   -------------------------------------
%   % Read and display true elevation values, in units of meters, 
%   % from the SRTM30 Plus layer provided by the NASA WorldWind 
%   % WMS server. Select a region surrounding Afghanistan and 
%   % request data at a 1-minute sampling interval.
%   layers = wmsfind('nasa.network', 'SearchField', 'serverurl');
%   srtm30 = layers.refine('srtm30', 'SearchField', 'layername');
%   figure
%   worldmap('Afghanistan')
%   mstruct = gcm;
%   latlim = mstruct.maplatlimit;
%   lonlim = mstruct.maplonlimit;
%   cellSize = dms2degrees([0,1,0]);
%   [Z, R] = wmsread(srtm30, 'Latlim', latlim, 'Lonlim', lonlim, ...
%      'CellSize', cellSize);
%   geoshow(Z, R, 'DisplayType', 'texturemap')
%   demcmap(double(Z))
%   title({'Afghanistan and Surrounding Region',  ...
%      'Global Elevation, Seamless SRTM Land Elevation and Ocean Depth'});
%
%   % Embed national boundaries from the NASA Globe Visualization Server
%   % into the elevation map.
%   vizglobe = wmsfind('viz.globe', 'SearchField', 'serverurl');
%   boundaries = vizglobe.refine('national*bound');
%   [B, R] = wmsread(boundaries, 'Latlim', latlim, 'Lonlim', lonlim, ...
%      'CellSize', cellSize);
%   Z(B(:,:,1) == 0) = min(Z(:));
%   geoshow(Z, R, 'DisplayType', 'texturemap')
%
%   Example 3 - Latlim, Lonlim, ImageHeight, ImageWidth
%   ---------------------------------------------------
%   % Read and display an orthoimage of the northern section of 
%   % the Golden Gate Bridge in San Francisco, California,  
%   % using the USGS National Map Seamless server.
%
%   % Define region of interest.
%   latlim = [37.824928  37.829598];
%   lonlim = [-122.482373 -122.47768];
%
%   % Find USGS layers based on region of interest.
%   caLayers = wmsfind('usgs*california', ...
%     'SearchField', 'serverurl', ...
%     'Latlim', latlim, 'Lonlim', lonlim);
%
%   % Refine the search to include only the San Franciso layers.
%   % Remove layers that contain metadata or footprint information.
%   caLayers = caLayers.refine('SanFranciscoCA');
%   removeNames = {'footprint', 'metadata'};
%   layerNames = {caLayers.LayerName};
%   index = true(1, numel(layerNames));
%   for k=1:numel(removeNames)
%       tf = regexpi({caLayers.LayerName}, removeNames{k});
%       index = index | cellfun(@isempty, tf);
%   end
%   orthoLayers = caLayers(index);
%
%   % Obtain the image.
%   imageLength = 1024;
%   [A, R] = wmsread(orthoLayers, 'Latlim', latlim, 'Lonlim', lonlim, ...
%       'ImageHeight', imageLength, 'ImageWidth', imageLength);
%
%   % Display the orthoimage in a UTM projection.
%   figure
%   axesm('utm', 'Zone', utmzone(latlim, lonlim), ...
%      'MapLatlimit', latlim, 'MapLonlimit', lonlim, ...
%      'Geoid', wgs84Ellipsoid)
%   geoshow(A,R)
%   axis off
%   title({'San Francisco','Northern Section of Golden Gate Bridge'})
%   
%   Example 4 - Latlim, Lonlim, ImageHeight, ImageWidth
%               ImageFormat          
%   ---------------------------------------------------
%   % Drape Landsat imagery onto elevation data from the 
%   % USGS National Elevation Dataset (NED) for an area 
%   % surrounding the Grand Canyon. Read the landsat layer
%   % from the USGS EROS Web Map server and the USGS NED 
%   % layer from the NASA WorldWind server.
%
%   % Obtain the layers of interest.
%   landsat = wmsfind('landsat', 'SearchField', 'any');
%   global_mosaic = landsat.refine('LANDSAT_LZ77', 'MatchType', 'exact');
%   layers = wmsfind('nasa.network', 'SearchField', 'serverurl');
%   us_ned = layers.refine('usgs ned 30');
%
%   % Assign geographic extent and image size.
%   latlim = [36 36.23];
%   lonlim = [-113.36 -113.13];
%   imageHeight = 575;
%   imageWidth = 575;
% 
%   % Read the global_mosaic layer.
%   [A, R] = wmsread(global_mosaic, ...
%      'Latlim', latlim, 'Lonlim', lonlim, ...
%      'ImageHeight', imageHeight, 'ImageWidth', imageWidth);
%
%   % Read the USGS NED layer.
%   [Z, R] = wmsread(us_ned, 'ImageFormat', 'image/bil', ...
%      'Latlim', latlim, 'Lonlim', lonlim, ...
%      'ImageHeight', imageHeight, 'ImageWidth', imageWidth);
% 
%   % Drape the Landsat image onto the elevation data.
%   figure
%   usamap(latlim, lonlim)
%   framem off;mlabel off;plabel off;gridm off
%   geoshow(double(Z), R, 'DisplayType', 'surface', 'CData', A);
%   daspectm('m',1)
%   title({'Grand Canyon', 'USGS NED and Landsat Global Mosaic'});
%   axis vis3d
% 
%   % Assign camera parameters.
%   cameraPosition = [96431 4.2956e+06 -72027];
%   cameraTarget = [-82.211 4.2805e+06 3054.6];
%   cameraViewAngle = 8.1561;
%   cameraUpVector = [3.8362e+06 5.9871e+05 5.05123e+006];
% 
%   % Set camera and light parameters.
%   set(gca,'CameraPosition', cameraPosition, ...
%      'CameraTarget', cameraTarget, ...
%      'CameraViewAngle', cameraViewAngle, ...
%      'CameraUpVector', cameraUpVector);
% 
%   lightHandle = camlight;
%   camLightPosition = [7169.3 1.4081e+06 -4.1188e+006];
%   set(lightHandle, 'Position', camLightPosition);
%
%   Example 5 - Time, Multiple Layers
%   ---------------------------------------------------
%   % Read and display a global monthly composite of 
%   % sea surface temperature for April 16, 2010, using data
%   % from the AMSR-E sensor on board the Aqua satellite.
%   % Include the coastlines, landmask, and nations layers.
%   coastwatch = wmsfind('coastwatch', 'SearchField', 'serverurl');
%   layers = coastwatch.refine('erdAAsstamday', ...
%     'Searchfield','serverurl');
%   time = '2010-04-16T00:00:00Z';
%   [A, R] = wmsread(layers(end:-1:1), 'Time', time);
%   figure
%   axesm('pcarree', 'Maplonlimit', [0, 360], ...
%      'PLabelLocation', 45, 'MLabelLocation', 90, ...
%      'MLabelParallel', -90, 'MeridianLabel', 'on', ...
%      'ParallelLabel', 'on');
%   geoshow(A, R);
%   title({layers(end).LayerTitle, time})
%
%   Example 6 - SampleDimension
%   ---------------------------------------------------
%   % Read and display a single sequence image from the
%   % MODIS instruments on the Aqua and Terra satellites
%   % that shows hurricane Katrina on August 29, 2005.
%
%   % Find the hurricane Katrina sequence layer.
%   katrina = wmsfind('Hurricane Katrina (Sequence)');
%   katrina = wmsupdate(katrina(1));
% 
%   % The Dimension.Extent field shows a sequence delimited 
%   % by commas. The sequence starts on August 24 and ends 
%   % on August 31. The commas start at August 25 and end
%   % after August 30. Select the sequence corresponding to 
%   % August 29.
%   commas = strfind(katrina.Details.Dimension.Extent, ',');
%   extent = katrina.Details.Dimension.Extent;
%   sequence = extent(commas(end-2)+1:commas(end-1)-1);
%
%   % Obtain the time, latitude, and longitude limits
%   % from the values in the sequence. Split the sequence 
%   % into a cell array of values by first finding
%   % all values between and including the parentheses,
%   % then remove the parentheses and split the values.
%   pat = '[(-.\d)]';
%   r = regexp(sequence, pat);
%   values = sequence(r);
%   values = strrep(values, '(', ' ');
%   values = strrep(values, ')', ' ');
%   values = regexp(values, '\s', 'split');
%   values = values(~cellfun('isempty', values));
%   time = values{1};
%   xmin = values{2};
%   ymin = values{3};
%   xmax = values{4};
%   ymax = values{5};
%
%   % Define latitude and longitude limits from the information
%   % in the sequence. The layer's geographic extent is assigned 
%   % for the combined set of sequences. The requested map cannot
%   % be a subset of the layer's bounding box. In this rare case, 
%   % set the layer's limits using the limits of the sequence.
%   latlim = [str2double(ymin) str2double(ymax)];
%   lonlim = [str2double(xmin) str2double(xmax)];
%   katrina.Latlim = latlim;
%   katrina.Lonlim = lonlim;
%
%   % Read and display the sequence map.
%   [A,R] = wmsread(katrina, 'SampleDimension', ...
%      {katrina.Details.Dimension.Name, sequence});
%   figure
%   usamap(katrina.Latlim,katrina.Lonlim);
%   geoshow(A,R)
%   coast = load('coast');
%   plotm(coast.lat, coast.long)
%   title({katrina.LayerTitle, time})
%
%   See also WebMapServer, WMSFIND, WMSINFO, WMSLayer, 
%            WMSMapRequest, WMSUPDATE.

% Copyright 2009-2012 The MathWorks, Inc.
% $Revision: 1.1.6.13.2.2 $  $Date: 2012/01/08 22:04:57 $

% Validate number of inputs.
narginchk(1,Inf)

% Create a mapRequest object based on the data type of the input parameter.
if ischar(layerOrURL)
    % WMSREAD(mapRequestURL)
    % Only one input is allowed when using the mapRequestURL interface.
    narginchk(1,1)
    
    % Validate the input and construct a mapRequest structure.
    mapRequest = createMapRequestStruct(layerOrURL);  
else
    % WMSREAD(LAYER)
    % Create a WMSMapRequest object. The constructor validates the input
    % as a WMSLayer array.
    mapRequest = WMSMapRequest(layerOrURL);
end

% Ensure that mapRequest contains a valid coordinate reference system.
validateCoordRefSys(mapRequest)
    
% Obtain the server object from the mapRequest.
server = mapRequest.Server;

% Set the server Timeout property to 60 seconds.
timeoutInSeconds = 60;
server.Timeout = secondsToServerTimeUnits(timeoutInSeconds);

% Parse the parameter/value pairs, if supplied, and set the server and
% mapRequest properties.
if numel(varargin) > 0
    [server, mapRequest] = setProperties(server, mapRequest, varargin);
end

% Construct the mapRequestURL.
mapRequestURL = mapRequest.RequestURL;

% Get the map from the server.
A = server.getMap(mapRequestURL);

% Get the raster referencing matrix.
R = getRasterRef(mapRequest, A);

end

%--------------------------------------------------------------------------

function timeoutInServerUnits = secondsToServerTimeUnits(timeoutInSeconds)
% Convert seconds to time units of the server (milliseconds). 

% The WebMapServer.Timeout property is expressed in milliseconds.
serverTimeUnitsPerSecond = 1000;
timeoutInServerUnits = serverTimeUnitsPerSecond * timeoutInSeconds;
end

%--------------------------------------------------------------------------

function mapRequest = createMapRequestStruct(mapRequestURL)
% Create a structure that contains field names that mimic the property
% names of a WMSMapRequest. The 'Server' field is populated with a new
% WebMapServer object. The WebMapServer constructor validates the
% mapRequestURL.  The 'Latlim' and 'Lonlim' fields are populated with
% values from the mapRequestURL. The 'RequestURL' field is populated with
% the input mapRequestURL. The RasterRef cannot be computed at this time
% and is intentionlly left empty.

% Construct the WebMapServer object.
server = WebMapServer(mapRequestURL);

% Obtain the coordinate limits from the mapRequestURL.
[latlim, lonlim] = parseBoundingBox(mapRequestURL);

% Create the mapRequest structure.
mapRequest = struct( ...
    'Server', server, ...
    'RasterRef', [], ...
    'Latlim', latlim, ...
    'Lonlim', lonlim, ...
    'RequestURL', mapRequestURL);
end

%--------------------------------------------------------------------------

function [latlim, lonlim] = parseBoundingBox(url)
% Parse the URL to obtain the bounding box limits. LATLIM and LONLIM are
% returned as two-element vectors of double. A bounding box element in the
% WMS URL string, URL, consists of the following string characters:
%    &BBOX=MIN_LON,MIN_LAT,MAX_LON,MAX_LAT
% The coordinate limits are separated by the comma delimiter. An example of
% a bounding box element for limits that span the entire Earth is:
%   &BBOX=-180.0,-90.0,180.0,90.0&
% Ampersands separate WMS parameters; however, if the element is at the end
% of the URL, then a trailing & is not present.

% Find the numbers (as strings) associated with the bounding box element.
bboxParam = 'BBOX';
S = parseMapRequestURL(url, bboxParam);
numbersElement = S.(bboxParam);

% There needs to be four numbers in the numbersElement string separated by
% three commas. For example, a valid numbersElement is -180,-90,180,90
numRequiredCommas = 3;
commas = strfind(numbersElement, ',');
assert(~isempty(commas) && numel(commas) == numRequiredCommas, ...
    message('map:wms:invalidWMSBoundingBox', url, ['&' bboxParam]));

% Convert the string numbers to numeric values (class double).
strNumbers = strrep(numbersElement, ',',' ');
bboxNumbers = sscanf(strNumbers, '%f');
assert(numel(bboxNumbers) == 4, message('map:wms:invalidWMSBoundingBox', ...
   url, ['&' bboxParam]));

% Assign the latlim and lonlim arrays with the numeric numbers.
lonlim = bboxNumbers([1 3]);
latlim = bboxNumbers([2 4]);
if limitsNeedToSwap(url, latlim)
    ltemp = latlim;
    latlim = lonlim;
    lonlim = ltemp;
end
end

%--------------------------------------------------------------------------

function swapLimits = limitsNeedToSwap(url, latlim)
% Determine if the limits need to swap. In general, if the URL specifies
% version 1.3.0 and the CRS is 'EPSG:4326', then the limits should in most
% cases need to swap. However, some servers do not conform to the
% specification. Attempt to use the information from the BoundingBox of the
% layer and the BBOX parameter, as indicated by latlim, to determine if
% the limits should swap.

% Parse the version number.
S = parseMapRequestURL(url,'VERSION');

% Determine if the version number is 1.3.x.
if strncmp(S.VERSION, '1.3', 3)
    % The version is 1.3.x. Parse the CRS and LAYERS parameters.
    S = parseMapRequestURL(url, 'CRS', 'LAYERS');
    
    if isequal(S.CRS, 'EPSG:4326')
        % For versions prior to 1.3.x, XLim corresponds to lonlim and YLim
        % corresponds to latlim. In this case, the BBOX parameter contains
        % the following:
        % &BBOX=MIN_LON,MIN_LAT,MAX_LON,MAX_LAT
        %
        % However, for version 1.3.0 and EPSG:4326, the specification
        % states that the order should be:
        % &BBOX=MIN_LAT,MIN_LON,MAX_LAT,MAX_LON
        %
        % Swap latlim and lonlim if it can be determined that the server is
        % following the 1.3.0 specification. To determine if the limits
        % need to swap, test two cases:
        % 1) Validate that latlim is within range
        % 2) Determine if the YLim limits from the BoundingBox element
        % represent longitude. If YLim exceeds the interval [-90 90], then
        % YLim corresponds to longitude, the limits need to swap, and the
        % server is conforming to version 1.3.x.
        
        % Case 1: Validate that latlim is with range.
        % If not, then swap the limits. Otherwise, continue to test the
        % BoundingBox.
        latlimWithinRange = -90 <= min(latlim) && max(latlim) <= 90;
        
        % If latlim is not within range, then the limits need to swap;
        % otherwise, latlim is within range and it is still undetermined if
        % the limits need to swap.
        if latlimWithinRange
            
            % A layer must be created and updated to obtain the
            % BoundingBox.
            index = strfind(S.LAYERS, ',');
            if isempty(index)
                layerName = S.LAYERS;
            else
                layerName = S.LAYERS(1:index);
            end
            layer = WMSLayer('ServerURL', url, 'LayerName', layerName);
            layer = wmsupdate(layer);
            
            if ~isempty(layer)
                % Determine if the limits need to swap based on values in
                % the BoundingBox. (Note: this algorithm may fail if the
                % BoundingBox does not cover the whole globe)
                swapLimits = swapLimitsForV1_3(layer.Details.BoundingBox);
            else
                % The layer is empty. Assume that the server is not
                % conforming to version 1.3.x.
                swapLimits = false;
            end
        else
            % The server is conforming to version 1.3.x. The latitude
            % limits are not within the [-90 90] interval, therefore they
            % need to be swapped.
            swapLimits = true;
        end
    else
        % Not using EPSG:4326. The limits do not need to swap.
        swapLimits = false;
    end
else
    % The version is not 1.3.x. The limits do not need to swap.
    swapLimits = false;
end
end

%--------------------------------------------------------------------------

function R = getRasterRef(mapRequest, A)
% Return a referencing matrix based on mapRequest. If mapRequest is a
% WMSMapRequest, then return the RasterRef property. If mapRequest is a
% structure, then calculate a referencing matrix based on the size of A and
% the 'Latlim' and 'Lonlim' fields of mapRequest.

if isa(mapRequest, 'WMSMapRequest')
    R = mapRequest.RasterRef;
else
    assert(~isempty(A), message('map:wms:mapIsEmpty', mapRequest.RequestURL));
    
    % Calculate the referencing matrix.
    xLim = mapRequest.Lonlim;
    yLim = mapRequest.Latlim;
    imageWidth  = size(A, 2);
    imageHeight = size(A, 1);
    xCell = diff(xLim)/imageWidth;
    yCell = diff(yLim)/imageHeight;
    x = xLim(1);
    y = yLim(2);
    R = makerefmat(x + xCell/2, y - yCell/2, xCell, -yCell); 
end
end

%--------------------------------------------------------------------------

function validateCoordRefSys(mapRequest)
% Validate that mapRequest contains a valid coordinate reference system.

% Assign a value for the required CRS code.
requiredCRS = {'EPSG:4326', 'CRS:84'};

% Find the CRS code.
if isa(mapRequest, 'WMSMapRequest')
    coordRefSysCode = mapRequest.CoordRefSysCode;
else
    upperURL = upper(mapRequest.RequestURL);
    k = 1;
    while k <= numel(requiredCRS)
        paramIndex = strfind(upperURL, requiredCRS{k});
        if ~isempty(paramIndex)
            coordRefSysCode = requiredCRS{k};
            k = numel(requiredCRS);
        else
            coordRefSysCode = '';
        end
        k = k + 1;
    end
end
    
% Assert that mapRequest contains the required CRS code.
assert(any(strcmpi(coordRefSysCode, requiredCRS)), message('map:wms:nonGeographicCRS', ...
   requiredCRS{1}, requiredCRS{2}));
end

%--------------------------------------------------------------------------

function [server, mapRequest] = setProperties(server, mapRequest, params)
% Parse the parameter/value pairs from PARAMS and use the parsed values
% to set the properties on the WebMapServer object, server, or the 
% WMSMapRequest object, mapRequest. The set methods of the object's
% properties are used to validate the parameters.

% Parse the mapRequest parameters and set the mapRequest properties if the
% parameters are supplied.
[mapRequest, mapRequestParameterNames, unparsedParams, userSupplied] = ...
    setMapRequestProperties(mapRequest, params);

% Parse the CellSize parameters and set the mapRequest properties that
% correspond to the CellSize parameters.
[mapRequest, cellSizeParameterNames, unparsedParams] = ...
    setCellSizeProperties(mapRequest, userSupplied, unparsedParams);

% Parse the server parameters and set the server Timeout property if the
% parameter is supplied.
[server, serverParameterName, unparsedParams] = ...
    setServerTimeoutProperty(server, unparsedParams);

% Verify the first parameter is a string.
if ~isempty(params) && ~ischar(params{1})
    unparsedParams{1} = 'PARAM1';
end

% Check if unparsedParams contained unmatched parameters.
if ~isempty(unparsedParams)
    parameterNames = ...
        [mapRequestParameterNames ...
        serverParameterName ...
        cellSizeParameterNames];
    parameterNames = sprintf('''%s'', ', parameterNames{:});
    error(message('map:validate:invalidParameterName', ...
        unparsedParams{1}, parameterNames(1:end-2)));
end
end

%--------------------------------------------------------------------------

function [mapRequest, parameterNames, unparsedParams, userSupplied] = ...
    setMapRequestProperties(mapRequest, params)
% Parse PARAMS for the mapRequestParameters and set the corresponding
% mapRequest properties if the parameters are supplied. mapRequest is a
% scalar WMSMapRequest object. The modified mapRequest object is returned.
%
% parameterNames is a cell array of strings that define the properties of
% mapRequest that can be set with a corresponding match from an element of
% PARAMS.
%
% unparsedParams is a cell array containing elements in PARAMS that did
% not match elements in parameterNames.
%
% userSupplied is a scalar structure array indicating true if a
% parameter-value pair is supplied in PARAMS, else false. The fieldnames of
% userSupplied match the elements of parameterNames. 

% Assign the parameter names that correspond to the properties of the
% mapRequest object.
parameterNames = { ...
    'Latlim', ...
    'Lonlim', ...
    'ImageHeight', ...
    'ImageWidth',  ...
    'StyleName', ...
    'ImageFormat', ...
    'Transparent',  ...
    'BackgroundColor', ...
    'Elevation', ...
    'Time', ...
    'SampleDimension'};

% Assign the anonymous validation functions for the mapRequest parameters.
% Use the set method of the mapRequest object to validate the mapRequest
% parameters. The nested function, setMapRequestProperty, validates the
% parameter and sets the mapRequest object. Since the function is nested,
% the modified object is in scope.
validateFcns = cell(size(parameterNames));
for k=1:numel(parameterNames)
    validateFcns{k} = @(x)setMapRequestProperty(parameterNames{k}, x);
end

% Parse the parameters and validate and set the appropriate properties by
% executing each of the anonymous functions assigned above. Thus parsepv
% modifies mapRequest even though it's not an input assignment.
[~, userSupplied, unparsedParams] = ...
    internal.map.parsepv(parameterNames, validateFcns, params{:});

    %----------------------------------------------------------------------
    
    function value = setMapRequestProperty(name, value)
    % Validate and set the mapRequest property named NAME to VALUE and
    % return the supplied VALUE. Throw any exceptions from caller to
    % prevent a long and unnecessary stack trace from the object.
        try
            mapRequest.(name) = value;
        catch e
            throwAsCaller(e)
        end
    end
end

%--------------------------------------------------------------------------

function [mapRequest, parameterNames, unparsedParams] = ...
    setCellSizeProperties(mapRequest,  mapRequestUserSupplied, params)
% Parse the cell array PARAMS for the 'CellSize' parameters and set the
% appropriate mapRequest properties if the parameters are supplied.
% mapRequest is a scalar WMSMapRequest object. The modified mapRequest
% object is returned.
%
% mapRequestUserSupplied is a scalar structure containing field names that
% match properties of the mapRequest object. If the field value is true,
% then the corresponding property has been set. 
%
% parameterNames is a cell array of strings and contains the names for the
% 'CellSize' parameters.
%
% unparsedParams is a cell array containing elements in PARAMS that did not
% match elements in parameterNames.

% Assign default values for CellSize parameters.
relTolCellSizeDefault = .001;
cellSizeDefaultValues = struct( ...
    'CellSize', [], ...
    'RelTolCellSize', [relTolCellSizeDefault relTolCellSizeDefault]);

% Assign parameterNames as a row vector.
parameterNames = fieldnames(cellSizeDefaultValues)';

% Assign the validation functions.
cellSizeValidationFcns = {@validateCellSize, @validateRelTolCellSize};

% Parse the parameters.
[options, userSupplied, unparsedParams] = ...
    internal.map.parsepv(parameterNames, cellSizeValidationFcns, params{:});

% Set default values for any parameter that is not specified.
options = setDefaultValues(options, cellSizeDefaultValues, userSupplied);

% Set the image size properties, 'ImageHeight' and 'ImageWidth', if
% 'CellSize' has been set.
if userSupplied.CellSize
    % If the 'ImageHeight' or 'ImageWidth' parameters have been supplied,
    % then issue an error.
    assert(~(mapRequestUserSupplied.ImageHeight ...
          || mapRequestUserSupplied.ImageWidth), ...
          message('map:wms:invalidCellSizeParameter'));
    
    % Set the 'ImageHeight' and 'ImageWidth' properties based on the cell
    % size parameters.
    mapRequest = setImageSizeProperties(mapRequest, options);
end
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
end

%--------------------------------------------------------------------------

function cellSize = validateCellSize(cellSize)
% Validate CellSize.

validateattributes(cellSize, {'numeric'}, ...
    {'positive', 'finite', 'nonempty',}, ...
    'validateCellSize', 'CellSize');
if isscalar(cellSize)
    cellSize(2) = cellSize(1);
else
    assert(numel(cellSize) == 2, ...
       message('map:validate:expectedTwoElementVector', 'CellSize'));
end
end

%--------------------------------------------------------------------------

function relTolCellSize = validateRelTolCellSize(relTolCellSize)
% Validate RelTolCellSize.

validateattributes(relTolCellSize, {'numeric'}, ...
    {'nonnegative', 'finite', 'nonempty',}, ...
    'validateRelTolCellSize', 'RelTolCellSize');
if isscalar(relTolCellSize)
    relTolCellSize(2) = relTolCellSize(1);
else
    assert(numel(relTolCellSize) == 2, ...
       message('map:validate:expectedTwoElementVector', 'CellSize'));
end
end

%--------------------------------------------------------------------------

function mapRequest = setImageSizeProperties(mapRequest, options)
% Set the image size properties, 'ImageHeight' and 'ImageWidth', based on
% the value of options.CellSize. 

% Calculate ImageHeight and ImageWidth based on CellSize.
cellSize = options.CellSize;
rasterExtentInLat = diff(mapRequest.Latlim);
rasterExtentInLon = diff(mapRequest.Lonlim);

imageHeight = round(rasterExtentInLat/cellSize(1));
imageWidth  = round(rasterExtentInLon/cellSize(2));

% Validate the ImageHeight and ImageWidth calculations.
maximumHeight = WMSMapRequest.MaximumHeight;
maximumWidth  = WMSMapRequest.MaximumWidth;
validHeight = 1 <= imageHeight && imageHeight <= maximumHeight;
validWidth  = 1 <= imageWidth  && imageWidth  <= maximumWidth;

assert(validHeight, message('map:wms:invalidCellHeight', ...
    'CellSize(1)', num2str(rasterExtentInLat/maximumHeight), ...
    'CellSize(1)', num2str(rasterExtentInLat)));

assert(validWidth, message('map:wms:invalidCellWidth', ...
    'CellSize(2)', num2str(rasterExtentInLon/maximumWidth), ...
    'CellSize(2)', num2str(rasterExtentInLon)));

% Validate the CellSize calculation based on input tolerance.
relTolCellSize = options.RelTolCellSize;
actualCellSize = ...
    [rasterExtentInLat/imageHeight rasterExtentInLon/imageWidth];

withinTolerance = abs(actualCellSize(1) - cellSize(1)) <= relTolCellSize(1) * cellSize(1);
assert(withinTolerance, message('map:validate:cellHeightNotWithinTolerance', ...
    num2str(cellSize(1)),  num2str(actualCellSize(1))));

withinTolerance = abs(actualCellSize(2) - cellSize(2)) <= relTolCellSize(2) * cellSize(2);
assert(withinTolerance, message('map:validate:cellWidthNotWithinTolerance', ...
    num2str(cellSize(2)),  num2str(actualCellSize(2))));

% Set the ImageHeight and ImageWidth properties.
mapRequest.ImageHeight = imageHeight;
mapRequest.ImageWidth  = imageWidth;
end

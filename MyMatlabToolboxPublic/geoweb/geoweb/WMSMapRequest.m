%WMSMapRequest Web Map Service map request object
%
%   A WMSMapRequest object contains a request to a WMS server to obtain a
%   map which represents geographic information.  The map is rendered as a
%   color or grayscale image. The object contains properties that can be
%   set to facilitate expressing the geographic extent, rendering, or size
%   of the requested map.
%
%   mapRequest = WMSMapRequest(LAYER) constructs a WMSMapRequest object.
%   LAYER is a WMSLayer array and must contain only one unique ServerURL.
%   LAYER's properties are updated, if required. 
%
%   mapRequest = WMSMapRequest(LAYER, SERVER) constructs a WMSMapRequest
%   object. SERVER is a scalar WebMapServer object. LAYER's ServerURL 
%   property must match the ServerURL property of SERVER.  If LAYER
%   contains properties that need to be updated, then LAYER is updated by
%   the SERVER object. 
%
%   WMSMapRequest properties:
%      Server - WebMapServer object
%      Layer - WMSLayer array
%      CoordRefSysCode - Coordinate reference system code
%      RasterRef - Raster referencing system object
%      Latlim - Latitude limits of request
%      Lonlim - Longitude limits of request
%      XLim - X limits of request
%      YLim - Y limits of request
%      ImageHeight - Height of raster map
%      ImageWidth - Width of raster map
%      MaximumHeight - Maximum height of raster map
%      MaximumWidth - Maximum width of raster map
%      Elevation - Requested elevation extent
%      Time - Requested time extent
%      SampleDimension - Requested dimension extent
%      Transparent - Transparent value of map
%      BackgroundColor - Color of background pixels
%      StyleName - Name of rendering style
%      ImageFormat - Format of raster map
%      ImageRenderFormats - Order of choice for image format
%      ImageTransparentFormats - Order of choice for transparent image
%      ServerURL - Server URL for GetMap request
%      RequestURL - Request URL for GetMap request
%
%   WMSMapRequest method:
%      boundImageSize - Bound size of raster map
%
%   Example 1 - Set the following properties:
%               ImageHeight, ImageWidth, ImageFormat, Time 
%   ------------------------------------------------------
%   % Read a global 1/2 degree sea surface temperature map for the month of
%   % November 2009 from the AMSR-E sensor on NASA's Aqua satellite using
%   % data provided by NASA's Earth Observations (NEO) WMS server. 
%   sst =  wmsfind('AMSRE_SST_M');
%   server = WebMapServer(sst.ServerURL);
%   mapRequest = WMSMapRequest(sst, server);
%   timeRequest = '2009-11-01';
%   mapRequest.Time = timeRequest;
%   samplesPerInterval = .5;
%   mapRequest.ImageHeight = round(abs(diff(sst.Latlim))/samplesPerInterval);
%   mapRequest.ImageWidth  = round(abs(diff(sst.Lonlim))/samplesPerInterval);
%   mapRequest.ImageFormat = 'image/png';
%   sstImage = server.getMap(mapRequest.RequestURL);
%
%   % The legend for the layer can be obtained via the OnlineResource URL
%   % field in the LegendURL structure. The legend shows that the
%   % temperatures range from -2 to 35 degrees Celsius. Note that the
%   % WMSMapRequest object updates the layer information from the server.
%   [legend, cmap] = ...
%      imread(mapRequest.Layer.Details.Style(1).LegendURL.OnlineResource);
%   if ~isempty(cmap)
%      % Convert indexed image to RGB.
%      legendImg = ind2rgb(legend, cmap);
%   else
%      % Already have an RGB image.
%      legendImg = legend;
%   end
%
%   % Display the temperature map and legend.
%   figure('Color','white')
%   worldmap world
%   setm(gca, 'MlabelParallel', -90, 'MlabelLocation', 90)
%   geoshow(sstImage, mapRequest.RasterRef);
%   title({mapRequest.Layer.LayerTitle, timeRequest}, ...
%      'Interpreter', 'none', 'FontWeight', 'bold')
%
%   figurePosition = get(gcf,'position');
%   centerWidth = figurePosition(3)/2;
%   left = centerWidth - size(legendImg,2)/2;
%   bottom = 30;
%   width  = size(legendImg,2);
%   height = size(legendImg,1);
%   axes('Units', 'pixels', 'Position', [left bottom width height])
%   image(legendImg)
%   axis off
%
%   % Additional abstract information for this layer can be obtained from
%   % the MetadataURL field. 
%   filename = [tempname '.xml'];
%   urlwrite(mapRequest.Layer.Details.MetadataURL, filename);
%   xml = xmlread(filename);
%   delete(filename);
%   abstract = xml.getElementsByTagName('abstract').item(0).getTextContent
%
%   Example 2 - Set the following properties:
%               Latlim, Lonlim, ImageHeight, ImageWidth, ImageFormat
%   ----------------------------------------------------------------
%   % Read and display a global elevation and bathymetry layer 
%   % for the Gulf of Maine at 30 arc-seconds sampling interval. 
%   % The values are in units of meters.
%   layers = wmsfind('srtm30', 'SearchField', 'LayerName');
%   layer = layers.refine('nasa.network','SearchField', 'serverurl')
%   server = WebMapServer(layer.ServerURL);
%   mapRequest = WMSMapRequest(layer, server);
%   mapRequest.Latlim = [40 46];
%   mapRequest.Lonlim = [-71 -65];
%   samplesPerInterval = 30/3600;
%   mapRequest.ImageHeight = ...
%      round(abs(diff(mapRequest.Latlim))/samplesPerInterval);
%   mapRequest.ImageWidth = ...
%      round(abs(diff(mapRequest.Lonlim))/samplesPerInterval);
%   mapRequest.ImageFormat = 'image/bil';
%   Z = server.getMap(mapRequest.RequestURL);
% 
%   % Display and contour the map at sea level (0 meters).
%   figure
%   worldmap(mapRequest.Latlim, mapRequest.Lonlim)
%   geoshow(double(Z), mapRequest.RasterRef, 'DisplayType', 'texturemap')
%   demcmap(double(Z))
%   set(gcf, 'renderer','zbuffer')
%   contourm(double(Z), mapRequest.RasterRef, [0 0], 'color', 'black')
%   colorbar
%   title ({'Gulf of Maine', mapRequest.Layer.LayerTitle}, ...
%      'Interpreter', 'none', 'FontWeight', 'bold')
%   
%   See also WebMapServer, WMSFIND, WMSINFO, WMSLayer, WMSREAD.

% Copyright 2009-2011 The MathWorks, Inc.
% $Revision: 1.1.6.16 $  $Date: 2011/06/27 19:48:09 $

classdef WMSMapRequest
    
    properties (GetAccess = public, SetAccess = private)
        %Server WebMapServer object
        %
        %   Server is a scalar WebMapServer object. The Server property is 
        %   initialized to the SERVER input, if supplied to the
        %   constructor; otherwise, the Server property is constructed
        %   using the ServerURL of Layer.
        Server;
    end
        
    properties (Access = public, Dependent = true)
        %Layer WMSLayer array
        %
        %   Layer is a WMSLayer array containing one unique ServerURL. If
        %   required, the properties of Layer are updated by the Server
        %   property, when the property is set. The ServerURL property of
        %   Layer must match the ServerURL property of Server.  The Layer
        %   property is initialized to the LAYER input supplied to the
        %   constructor.
        Layer
        
        %CoordRefSysCode Coordinate reference system code
        %
        %   CoordRefSysCode is a string that specifies the coordinate
        %   reference system code. Its default value is 'EPSG:4326'. If
        %   that value is not found in Layer.CoordRefSysCodes then it is
        %   set from the first CoordRefSysCode found in the
        %   Layer.Details.BoundingBox structure array. When set to
        %   'EPSG:4326' or 'CRS:84', the XLim and YLim properties are set
        %   to [] and the Latlim and Lonlim properties are set to the
        %   geographic extent defined by the Layer array. When set to a
        %   value other than 'EPSG:4326' or 'CRS:84', then the XLim and
        %   YLim properties are set from the values found in the
        %   Layer.Details.BoundingBox structure and the Latlim and Lonlim
        %   properties are set to []. Automatic projections are not
        %   supported.(Automatic projections begin with the string 'AUTO').
        CoordRefSysCode
    end
    
    properties (GetAccess = public, SetAccess = private, Dependent = true)
        %RasterRef Raster referencing system object
        %
        %   RasterRef is a referencing matrix that references the raster
        %   map to an intrinsic coordinate system.
        RasterRef
    end
    
    properties (Access = public, Dependent = true)      
        %Latlim Latitude limits of request
        %
        %   Latlim is a two-element vector specifying the southern and
        %   northern latitude limits of the request in units of degrees.
        %   The limits must be ascending.  By default, Latlim is
        %   initialized to limits that span all the latitudinal limits
        %   found in the Layer.Latlim property.
        Latlim
        
        %Lonlim Longitude limits of request
        %
        %   Lonlim is a two-element vector specifying the western and
        %   eastern longitude limits of the request in units of degrees.
        %   The limits must be ascending and in the range [-180, 180] or
        %   [0 360]. By default, Lonlim is initialized to limits that span
        %   all the longitudinal limits in the Layer.Lonlim property.
        Lonlim
        
        %XLim X limits of request
        %
        %   XLim is a two-element vector specifying the western and eastern
        %   limits of the request in units specified by the coordinate
        %   reference system. The limits must be ascending. By default,
        %   XLim is empty. XLim may only be set if CoordRefSysCode is set
        %   to a value other than EPSG:4326.
        XLim
        
        %YLim Y limits of request
        %
        %   YLim is a two-element vector specifying the southern and
        %   northern limits of the request in units specified by the
        %   coordinate reference system. The limits must be ascending. By
        %   default, YLim is empty. YLim may only be set if CoordRefSysCode
        %   is set to a value other than EPSG:4326.
        YLim
        
        %ImageHeight Height of raster map
        %
        %   ImageHeight is a scalar positive integer indicating the desired
        %   height in pixels for the requested raster map. The maximum
        %   value for ImageHeight is defined by the property MaximumHeight.
        %   ImageHeight is initialized to either 512 or to an integer value
        %   that best preserves the aspect ratio of the coordinate limits,
        %   without changing the coordinate limits.
        ImageHeight
        
        %ImageWidth Width of raster map
        %
        %   ImageWidth is a scalar positive integer indicating the desired
        %   width in pixels for the requested raster map. The maximum value
        %   for ImageWidth is defined by the property MaximumWidth.
        %   ImageWidth is initialized to either 512 or to an integer value
        %   that best preserves the aspect ratio of the coordinate limits,
        %   without changing the coordinate limits.
        ImageWidth
    end
    
    properties (GetAccess = public, Constant = true)
        %MaximumHeight Maximum height of raster map
        %
        %   MaximumHeight is a double indicating the maximum height in
        %   pixels for the requested map and cannot be set.  The value of
        %   MaximumHeight is 8192.
        MaximumHeight = 8192;

        %MaximumWidth Maximum width of raster map
        %
        %   MaximumWidth is a double indicating the maximum width in pixels
        %   for the requested map and cannot be set.  The value of
        %   MaximumWidth is 8192.
        MaximumWidth = 8192;
    end
    
    properties (Access = public, Dependent = true)
        %Elevation Requested elevation extent
        %
        %   Elevation is a string indicating the desired elevation extent
        %   of the requested map. When the property is set, 'elevation'
        %   must be the value of the Layer.Details.Dimension.Name field.
        %   The default value is ''.
        Elevation
        
        %Time Requested time extent
        %
        %   Time is stored as a string indicating the desired time extent
        %   of the requested map. When the property is set, 'time' must be
        %   the value of the Layer.Details.Dimension.Name field. Time is
        %   stored in the ISO 8601:1988(E) extended format. In general, the
        %   Time property is stored in a yyyy-mm-dd format or a 
        %   yyyy-mm-ddThh:mm:ssZ format, if the precision requires hours, 
        %   minutes or seconds. The default value is ''. Time may be set
        %   using several different string and numeric formats, according
        %   to the following table (where dateform number is the number
        %   used by the Standard MATLAB Date Format Definitions). All
        %   hours, minutes, and seconds must be expressed in Coordinated
        %   Universal Time (UTC).
        %
        %   dateform   Input               Stored format
        %   (number)  (string)
        %   --------  ----------------     -------------
        %    0        dd-mm-yyyy HH:MM:SS  yyyy-mm-ddTHH:MM:SSZ
        %    1        dd-mm-yyyy           yyyy-mm-dd
        %    2        mm/dd/yy             yyyy-mm-dd
        %    6        mm/dd                yyyy-mm-dd (current year)
        %   10        yyyy                 yyyy
        %   13        HH:MM:SS             yyyy-mm-ddTHH:MM:SSZ 
        %   14        HH:MM:SS PM          yyyy-mm-ddTHH:MM:SSZ 
        %   15        HH:MM                yyyy-mm-ddTHH:MM:00Z 
        %   16        HH:MM PM             yyyy-mm-ddTHH:MM:00Z 
        %   21        mmm.dd,yyyy HH:MM:SS yyyy-mm-ddTHH:MM:SSZ
        %   22        mmm.dd,yyyy          yyyy-mm-dd
        %   23        mm/dd/yyyy           yyyy-mm-dd
        %   26        yyyy/mm/dd           yyyy-mm-dd
        %   29        yyyy-mm-dd           yyyy-mm-dd
        %   30        yyyymmddTHHMMSS      yyyy-mm-ddTHH:MM:SSZ
        %   31        yyyy-mm-dd HH:MM:SS  yyyy-mm-ddTHH:MM:SSZ
        %
        %   (Inputs using the dateform numbers 13-16 return the date set to
        %   the current year, month, and day). Use of other dateform
        %   formats, especially 19, 20, 24, and 25 will result in erroneous
        %   output. In addition to these standard MATLAB dateform formats,
        %   the following are also accepted:
        %     
        %   'current'            The current time holdings of the server
        %      
        %   numeric datenum      Numeric date value is converted to 
        %                        yyyy-mm-dd string (dateform 29 format)
        %
        %   Byyyy                Indicate B.C.E. year
        %
        %   Geologic datasets referring to the distant past may also be
        %   expressed, using the prefixes K, M, and G, followed by a string
        %   number, referring to some number of thousand, million, or
        %   billion years, respectively.
        Time
        
        %SampleDimension Requested dimension extent
        %
        %   SampleDimension is a two-element cell array of strings
        %   indicating the desired name of a sample dimension (other than
        %   'time' or 'elevation') and its value. SampleDimension{1} must
        %   be the value of the Layer.Details.Dimension.Name field.
        SampleDimension
        
        %Transparent Transparent value of map
        %
        %   Transparent is a logical scalar that specifies whether the map
        %   background is to be made transparent or not. When Transparent
        %   is true, all of the pixels not representing features or data
        %   values in that layer are set to a transparent value, producing
        %   a composite map. When Transparent is false, non-data pixels
        %   shall be set to the value of the background color. By default,
        %   the value is false.
        Transparent
    end
    
    properties (Access = public)
        %BackgroundColor Color of background pixels
        %
        %   BackgroundColor is a three-element vector of uint8 values that
        %   specifies the color to be used to fill the background
        %   (non-data) pixels of the map. The values range from 0 to 255.
        %   The default value is [255,255,255] which specifies the
        %   background color as white. BackgroundColor is stored as uint8.
        %   If BackgroundColor is set using non-uint8 numeric values, the
        %   values are cast to uint8.
        BackgroundColor = uint8([255,255,255]);
    end
    
    properties (Access = public, Dependent = true)
        %StyleName Name of rendering style
        %
        %   StyleName is a string or cell array of strings (matching the
        %   number of elements of Layer) and specifies the style to use
        %   when rendering the image. By default StyleName is set to an
        %   empty string. The StyleName must be a valid entry in the
        %   Layer.Details.Style.Name field.
        StyleName
        
        %ImageFormat Format of raster map
        %
        %   ImageFormat is a string containing the desired image format
        %   used to render the map as an image. If set, the format must
        %   match an entry in the Layer.Details.ImageFormats cell array and
        %   must match an entry in the ImageRenderFormats property.  If not
        %   set, the format defaults to a value in the ImageRenderFormats
        %   property.
        ImageFormat
    end
    
    properties (GetAccess = public, Constant = true)
        %ImageRenderFormats Order of choice for image format
        %
        %   ImageRenderFormats is a cell array of the preferred image
        %   rendering formats when Transparent is set to false. The first
        %   entry is the most preferred image format. If the preferred
        %   format is not stored in the Layer property, then the next
        %   format from the list is selected until a format is found. The
        %   ImageRenderFormats array is not used if the ImageFormat
        %   property is set.
        ImageRenderFormats = ...
            {'image/jpeg'; 'image/gif'; 'image/png'; 'image/tiff'; ...
            'image/geotiff'; 'image/geotiff8'; 'image/tiff8'; ...
            'image/png8'; 'image/bil'};
        
        %ImageTransparentFormats Order of choice for transparent image
        %
        %   ImageTransparentFormats is a cell array of the preferred image
        %   rendering formats when Transparent is set to true. When
        %   Transparent is set to true, the ImageFormat property is set to
        %   the first entry in the ImageTransparentFormats list, if it is
        %   stored in the Layer property. Otherwise, the list is searched
        %   for the next element until a match is found. If a transparent
        %   image format is not found in the list or if the ImageFormat
        %   property is set to a non-default value, then ImageFormat is
        %   unchanged.
        ImageTransparentFormats = ...
            {'image/gif'; 'image/png'; 'image/tiff'; 'image/geotiff'};
    end
    
    properties (Access = public)
        %ServerURL Server URL for GetMap request
        %
        %   ServerURL is a string containing the server URL for the WMS
        %   GetMap request. In general, ServerURL matches the Layer's
        %   ServerURL. However, some WMS servers, such as Microsoft's Terra
        %   Server, require a different URL for GetMap requests compared
        %   with a WMS GetCapabilities request. The default value of
        %   ServerURL is Layer(1).ServerURL.
        ServerURL = '';
    end
    
    properties (SetAccess = private, GetAccess = public, Dependent = true)
        %RequestURL Request URL for GetMap request
        %
        %   RequestURL is a string containing the URL for the WMS GetMap
        %   request. It is composed of the ServerURL with additional WMS
        %   parameter/value pairs.
        RequestURL
    end
    
    properties (GetAccess = private, Constant = true, Hidden = true)
        %DefaultCRS Default coordinate reference system code
        %
        %   DefaultCRS is the default coordinate reference system code.
        DefaultCRS = 'EPSG:4326';
    end
    
    properties (Hidden = true)   
        %Exceptions Exceptions format specification
        %
        %   Exceptions states the format specification for the way in
        %   which errors are to be reported. The WMS specification provides
        %   for the following types of specifications for versions 1.1.1
        %   and earlier:
        %
        %   application/vnd.ogc.se_xml      The default specification.
        %                                   Exceptions from the server are
        %                                   reported using Service
        %                                   Exception XML
        %
        %   application/vnd.ogc.se_inimage  Exceptions from the server are
        %                                   graphically returned in the
        %                                   raster map
        %
        %   application/vnd.ogc.se_blank    If an error is detected, the
        %                                   server returns a raster map
        %                                   containing all pixels set to
        %                                   the background color
        %
        %   For version 1.3.0 and later, the Exceptions values have been
        %   changed to XML, INIMAGE, and BLANK.
        Exceptions = 'application/vnd.ogc.se_xml';
    end
    
     properties (Hidden = true, Constant = true)      
        %DefaultImageLength Default length of image dimension
        %
        %   DefaultImageLength is the longest length of the row or column
        %   dimension in the image. ImageHeight or ImageWidth is set to
        %   DefaultImageLength during construction.
        DefaultImageLength = 512;
    end
    
    properties (Access = private, Hidden = true)                   
        % Each of the following private properties holds the value of a
        % corresponding dependent property.
        pLayer = [];
        pLatlim = [];
        pLonlim = [];
        pXLim = [];
        pYLim = [];
        pCoordRefSysCode = '';
        pImageHeight = [];
        pImageWidth = [];
        pImageFormat = '';
        pStyleName = {''};
        pElevation = '';
        pTime = '';
        pSampleDimension = {'', ''};
        pTransparent = false;
    end
    
    methods 
        function self =  WMSMapRequest(layer, varargin)
        % WMSMapRequest constructor.
            
            % Validate number of inputs.
            error(nargchk(1, 2, nargin, 'struct'));
            
            % Verify that Java is available.
            error(javachk('jvm','WMSMapRequest'));
            
            % Validate layer as a WMSLayer array.
            validateLayer(layer);
                      
            % Validate server or create a WebMapServer object.
            self.Server = validateServer(layer, varargin{:});
            
            % Set the ServerURL property.
            self.ServerURL = self.Server.ServerURL;
            
            % Validate and set properties.
            self = setProperties(self, layer);
            
            % Find and set the ImageFormat property from the values in
            % ImageRenderFormats.
            self.ImageFormat = findImageFormat(self.Layer);
        end
        
        %------- Layer set/get Methods ------------------------------------
        
        function self = set.Layer(self, layer)
            
            % Validate the layer.
            validateLayer(layer);
            
            % Copy the object in order to test if warnings are required.
            oldRequest = self;
                      
            % Update the properties based on values from layer.
            self = updateProperties(self, layer);
            
            % Check to see if any properties have changed state. If true,
            % then issue a warning.
            issueWarnings(oldRequest, self);           
        end
        
        function value = get.Layer(self)
            value = self.pLayer;
        end
        
        %------- CoordRefSysCode set/get Methods --------------------------
        
        function  self = set.CoordRefSysCode(self, coordRefSysCode)
            
            % Validate coordRefSysCode by determining if it exists in the
            % Layer property.
            validateCoordRefSysCode(self.Layer, coordRefSysCode);
            self.pCoordRefSysCode = upper(coordRefSysCode);
                            
            % Always update the limits when the CoordRefSysCode has
            % changed.
            self = setLimits(self);
        end
        
        function value = get.CoordRefSysCode(self)
            value = self.pCoordRefSysCode;
        end
        
        %------- Spatial Referencing set/get Methods ----------------------
        
        function R = get.RasterRef(self)
            if useGeoCoords(self.CoordRefSysCode)
                xLim = self.Lonlim;
                yLim = self.Latlim;
            else
                xLim = self.XLim;
                yLim = self.YLim;
            end
            xCell = diff(xLim)/self.ImageWidth;
            yCell = diff(yLim)/self.ImageHeight;
            x = xLim(1);
            y = yLim(2);
            R = makerefmat(x + xCell/2, y - yCell/2, xCell, -yCell);
        end
        
        %------- Coordinate Limits set/get Methods ------------------------
        
        function self = set.Latlim(self, latlim)
            validateLatlim(latlim, self.Layer, self.CoordRefSysCode);
            self.pLatlim = latlim;
        end
        
        function value = get.Latlim(self)
            value = self.pLatlim;
        end
        
        function self =  set.Lonlim(self, lonlim)
            validateLonlim(lonlim, self.Layer, self.CoordRefSysCode);
            self.pLonlim = lonlim;
        end
        
        function value = get.Lonlim(self)
            value = self.pLonlim;
        end
        
        function self =  set.XLim(self, xLim)
            validateMapLimit('XLim', xLim, self.Layer, self.CoordRefSysCode);
            self.pXLim = xLim;
        end
        
        function value = get.XLim(self)
            value = self.pXLim;
        end
        
        function self = set.YLim(self, yLim)
            validateMapLimit('YLim', yLim, self.Layer, self.CoordRefSysCode);
            self.pYLim = yLim;
        end
        
        function value = get.YLim(self)
            value = self.pYLim;
        end
        
        %------- Image Size set/get Methods -------------------------------
        
        function self =  set.ImageHeight(self, imageHeight)
            validateImageLength(self.Layer, 'ImageHeight', imageHeight);
            self.pImageHeight = imageHeight;
        end
        
        function value = get.ImageHeight(self)
            value = self.pImageHeight;
        end
        
        function self =  set.ImageWidth(self, imageWidth)
            validateImageLength(self.Layer, 'ImageWidth', imageWidth);
            self.pImageWidth = imageWidth;
        end
        
        function value = get.ImageWidth(self)
            value = self.pImageWidth;
        end
        
        %------- Dimension set/get Methods --------------------------------
        
        function self =  set.Elevation(self, elevation)
            validateElevation(getDimensionNames(self.Layer), elevation);
            self.pElevation = elevation;
        end
        
        function value = get.Elevation(self)
            value = self.pElevation;
        end
        
        function self =  set.Time(self, time)
            if ~isempty(time)
               self.pTime = ...
                   validateTime(getDimensionNames(self.Layer), time);
            else
                self.pTime = '';
            end
        end
        
        function value = get.Time(self)
            value = self.pTime;
        end
        
        function self =  set.SampleDimension(self, sampleDimension)
            validateSampleDimension(...
                getDimensionNames(self.Layer), sampleDimension);
            self.pSampleDimension = sampleDimension;
        end
        
        function value = get.SampleDimension(self)
            value = self.pSampleDimension;
        end
        
        %------- Image Rendering set/get Methods --------------------------
        
        function self = set.Transparent(self, transparent)
            validateattributes(transparent, {'logical'}, ...
                {'nonempty', 'scalar'}, 'set', 'Transparent');
            self.pTransparent = transparent;
            if transparent
                self.ImageFormat = findTransparentImageFormat(self.Layer, self.ImageFormat);
            end
        end
        
        function value = get.Transparent(self)
            value = self.pTransparent;
        end
        
        function self = set.BackgroundColor(self, backgroundColor)
            validateBackgroundColor(backgroundColor);
            self.BackgroundColor = uint8(backgroundColor);
        end
        
        function self =  set.StyleName(self, styleName)
            self.pStyleName = validateStyleName(self.Layer, styleName);
        end
        
        function value = get.StyleName(self)
            value = self.pStyleName;
        end
        
        function self =  set.ImageFormat(self, imageFormat)
            validateImageFormat(self.Layer, imageFormat)
            self.pImageFormat = imageFormat;
        end
        
        function value = get.ImageFormat(self)
            value = self.pImageFormat;
        end
        
        %------- ServerURL set/get Methods --------------------------------
        
        function self = set.ServerURL(self, serverURL)
            internal.map.validateURL('ServerURL', serverURL);
            self.ServerURL = serverURL;
        end
        
        function value = get.RequestURL(self)
            value = constructRequestURL(self);
        end
        
        %-------------- Public Methods ------------------------------------
        
        function self = boundImageSize(self, imageLength)
        %boundImageSize Bound size of raster map
        %
        %   mapRequest = boundImageSize(mapRequest, imageLength) bounds the
        %   size of the raster map based on the scalar imageLength.
        %   mapRequest must be a scalar WMSMapRequest object. imageLength
        %   is a double indicating the length in pixels for the row,
        %   ImageHeight, or column, ImageWidth, dimension.  The row or
        %   column dimension length is calculated using the aspect ratio of
        %   the Latlim and Lonlim properties or the aspect ratio of the
        %   XLim and YLim properties if they are set. The image dimension
        %   corresponding to the longest geographic or map coordinates is
        %   set to imageLength, the shortest dimension is set to the
        %   nearest integer value that preserves the aspect ratio, without
        %   changing the coordinate limits. The maximum value of
        %   imageLength is the maximum value of the MaximumHeight and
        %   MaximumWidth properties.
        %
        %   Example
        %   -------
        %   % Read and display a composite of multiple layers representing
        %   % the Earth's EGM96 geopotential model, coastlines, and
        %   % national boundaries from NASA's Globe Visualization server.
        %   % The rendered map has a spatial resolution of .5 degrees.
        %   vizglobe = wmsfind('viz.globe', 'SearchField', 'serverurl');
        %   coastlines = vizglobe.refine('coastline');
        %   national_boundaries = vizglobe.refine('national*bound');
        %   base_layer = vizglobe.refine('egm96');
        %   layers = [base_layer;coastlines;national_boundaries];
        %   request = WMSMapRequest(layers);
        %   request.Transparent = true;
        %   request = request.boundImageSize(720);
        %   overlayImage = request.Server.getMap(request.RequestURL);
        %
        %   figure
        %   worldmap('world')
        %   geoshow(overlayImage, request.RasterRef);
        %   title(base_layer.LayerTitle)
        %
        %   % Compare the map with the contoured data from 'geoid.mat'.
        %   geoid = load('geoid');        
        %   coast = load('coast');
        %   figure
        %   worldmap('world')
        %   contourfm(geoid.geoid, geoid.geoidrefvec, 15)
        %   geoshow(coast.lat, coast.long)
        %   title({'EGM96 Contoured Data', '(geoid.mat)'})
        
            assert(isscalar(self), message('map:validate:expectedScalarObject', ...
                'mapRequest'));
            
            validateattributes(imageLength, {'numeric'}, ...
                {'integer', 'positive', 'finite', 'nonempty', 'scalar'}, ...
                'boundImageSize', 'imageLength');
            
            maximumSize = max([self.MaximumWidth self.MaximumHeight]);
            assert(imageLength <= maximumSize,  message('map:validate:inputTooLarge', ...
                'imageLength', num2str(maximumSize)));
            
            % Determine if using geographic or map limits.
            if useGeoCoords(self.CoordRefSysCode)
                xLim = self.Lonlim;
                yLim = self.Latlim;
            else
                xLim = self.XLim;
                yLim = self.YLim;
            end
            
            % Calculate the length of the axes.
            yLength = abs(yLim(2) - yLim(1));
            xLength = abs(xLim(2) - xLim(1));
            
            % Calculate the image size properties based on the desired
            % imageLength and the length of the axes.
            if yLength > xLength
                imageHeight = imageLength;
                imageWidth  = max(1, fix(xLength/yLength * imageLength));
            else
                imageWidth  = imageLength;
                imageHeight = max(1, fix(yLength/xLength * imageLength));
            end
            
            % Set the new values for ImageHeight and ImageWidth properties.
            self.ImageHeight = imageHeight;
            self.ImageWidth  = imageWidth;
        end
    end
    
    %------------------- private methods ----------------------------------
    
    methods (Access = private)
       
        function self = setProperties(self, layer)
        % Set properties based on values from LAYER.
            
            % Validate the layer properties and set the Layer property.
            self.pLayer = validateLayerProperties(layer, self.Server);

            % Set the CoordRefSysCode property.
            self.pCoordRefSysCode = ...
                assignDefaultCRS(self.Layer, WMSMapRequest.DefaultCRS);

            % Set Latlim, Lonlim, XLim, and YLim properties.
            self = setLimits(self);

            % Calculate and set ImageHeight and ImageWidth  properties.
            self = setDefaultImageSize(self);
            
            % Set the StyleName property to a default value if required.
            if numel(layer) ~= numel(self.StyleName)
               self.StyleName = '';
            end
            
            % Set the Exception property
            version = self.pLayer(1).Details.Version;
            if isequal(version, '1.3.0')
                self.Exceptions = 'XML';
            else
                self.Exceptions = 'application/vnd.ogc.se_xml';
            end
        end
        
        %------------------------------------------------------------------
        
        function self = updateProperties(self, layer)
        % Update the properties of the object based on values from layer.
        % This method is invoked when a new layer has been set. The
        % coordinate referencing system properties, CoordRefSysCode,
        % Latlim, Lonlim, XLim and YLim, and the image size properties,
        % ImageHeight, ImageWidth, are updated by setProperties.  The
        % ImageFormat property does not need to be updated since the
        % layer's ServerURL match the server's ServerURL. On a per server
        % basis, all layers have the same ImageFormats.
            
            % Validate and update the following properties:
            %    CoordRefSysCode, Latlim, Lonlim, XLim, YLim
            %    ImageHeight, ImageWidth, StyleName
            self = setProperties(self, layer);
            
            % Update dimension properties
            dimensionNames = getDimensionNames(self.Layer);
            if ~any(strcmpi('elevation', dimensionNames))
                self.pElevation = '';
            end
            if ~any(strcmpi('time', dimensionNames))
                self.pTime = '';
            end
            sampleDimension = self.SampleDimension;
            if ~any(strcmpi(sampleDimension{1}, dimensionNames))
                self.pSampleDimension = {'',''};
            end
            
            % setProperties ensures that StyleName matches in size to the
            % number of layer elements. Update the StyleName property to a
            % default value if the StyleName element is not found in the
            % layer.
            for k=1:numel(self.Layer)
                styleNames = {self.Layer(k).Details.Style.Name, ''};
                if ~any(strcmpi(self.StyleName{k}, styleNames))
                    self.StyleName = '';
                    break;
                end
            end            
        end
        
        %------------------------------------------------------------------
        
        function self = setLimits(self)
        % Set geographic and map limits.

            % Determine if the limits defined by the BoundingBox structure
            % need to be used. If the CoordRefSysCode does not reference a
            % geographic coordinate system, then the limits are found in
            % the BoundingBox structure. Otherwise, a geographic coordinate
            % system is specified and the limits are determined from the
            % Latlim and Lonlim properties of the layer array. The limits
            % returned encompass the limits of the entire Layer array.

            if ~useGeoCoords(self.CoordRefSysCode)
                % Using map coordinate system. Obtain the limits from the
                % BoundingBox structure in Layer. Find the minimum and
                % maximum bounding box for each element of the array.
                xLim = [inf, -inf];
                yLim = [inf, -inf];
                for k=1:numel(self.Layer)
                    [layerXLim, layerYLim] = getLimitsFromBoundingBox( ...
                        self.Layer(k), self.CoordRefSysCode);
                    xLim = ...
                        [min([xLim(1) layerXLim(1)]), ...
                         max([xLim(2) layerXLim(2)])];
                    yLim = ...
                        [min([yLim(1) layerYLim(1)]),...
                         max([yLim(2) layerYLim(2)])];
                end
                if isempty(self.YLim)
                    self.pYLim  = yLim;
                end
                if isempty(self.XLim)
                    self.pXLim = xLim;
                end

                % Reset the Latlim and Lonlim properties.
                self.pLatlim = [];
                self.pLonlim = [];
            else
                % Using geographic coordinate system.
                if isempty(self.Latlim)
                    limit = [self.Layer.Latlim];
                    latlim(1) = min(limit(1:2:end));
                    latlim(2) = max(limit(2:2:end));
                    self.pLatlim = latlim;
                end
                if isempty(self.Lonlim)
                    limit = [self.Layer.Lonlim];
                    lonlim(1) = min(limit(1:2:end));
                    lonlim(2) = max(limit(2:2:end));
                    self.pLonlim = lonlim;
                end

                % Reset the XLim and YLim properties.
                self.pXLim = [];
                self.pYLim = [];
            end
        end
        
        %------------------------------------------------------------------
        
        function self = setDefaultImageSize(self)
        % Compute and set the image size properties. If the Layer array
        % contains any layers that have FixedWidth or FixedHeight fields
        % set, then set the image size properties using those values.
        % Otherwise, bound the image size to the default length, unless
        % previously set.
            
           if any(isFixedWidth(self.Layer)) || ...
                   any(isFixedHeight(self.Layer))
               [self.ImageHeight, self.ImageWidth] = ...
                   setFixedImageSize(self.Layer);
           else
               if isempty(self.ImageHeight) && isempty(self.ImageWidth)
                  self = boundImageSize( ...
                      self, WMSMapRequest.DefaultImageLength);
               else
                   % ImageHeight and ImageWidth have already been set.
               end
           end
        end
                      
        %------------------------------------------------------------------
        
        function requestURL = constructRequestURL(self)
        % Construct the RequestURL based on properties of the class.
            
            % Create a specification which is used to parameterize the WMS
            % options.
            serverURL = java.net.URL(self.ServerURL);
            spec = com.mathworks.toolbox.geoweb.wms.WMSClientSpecification(serverURL);
            
            % Set the version string.
            version = self.Layer(1).Details.Version;
            spec.setVersion(version)
            
            % Create a Java getMapRequest object.
            getMapRequest = spec.createGetMapRequest(serverURL, version);
            assert(~isempty(getMapRequest), message('map:wms:getMapRequestIsEmpty'));
            
            % Add all the layers to the request.
            for k=1:numel(self.Layer)
                layerName = java.lang.String(self.Layer(k).LayerName);
                styleName = java.lang.String(self.StyleName{k});
                getMapRequest.addLayer(layerName, styleName);
            end
            
            % Construct a bounding box object.
            bbox = constructBoundingBox(self, version);
            getMapRequest.setBBox(bbox);
            
            % Set coordinate reference system information.
            if isequal(version, '1.3.0')
                % Version 1.3.0 requires &CRS=
                getMapRequest.setCRS(self.CoordRefSysCode);
            else
                % Other versions require &SRS=
                getMapRequest.setSRS(self.CoordRefSysCode);
            end
            
            % Set image format specification.
            getMapRequest.setFormat(self.ImageFormat);
            
            % Set image dimensions.
            imageWidth  = num2str(self.ImageWidth);
            imageHeight = num2str(self.ImageHeight);
            getMapRequest.setDimensions(imageWidth, imageHeight);
            
            % Set Transparent.
            getMapRequest.setTransparent(self.Transparent);
            
            % Set background color. The color needs to be expressed as a
            % hexadecimal string with a '0x' prefix.
            red   = dec2hex(self.BackgroundColor(1), 2);
            green = dec2hex(self.BackgroundColor(2), 2);
            blue  = dec2hex(self.BackgroundColor(3), 2);
            hexColor = ['0x' red green blue];
            getMapRequest.setBGColour(hexColor);
            
            % Set Exceptions.
            if ~isequal(version, '1.3.0')
                getMapRequest.setExceptions(self.Exceptions);
            end
            
            % Set Time if the property is not empty.
            if ~isempty(self.Time)
                getMapRequest.setTime(self.Time);
            end
            
            % Set Elevation if the property is not empty.
            if ~isempty(self.Elevation)
                getMapRequest.setElevation(self.Elevation);
            end
            
            % Set SampleDimension if the property is not empty.
            if ~all(cellfun(@isempty, self.SampleDimension))
                getMapRequest.setSampleDimensionValue( ...
                    self.SampleDimension{1}, self.SampleDimension{2});
            end
            
            % Get the map request URL.
            requestURL = char(getMapRequest.getFinalURL);
        end
        
        %------------------------------------------------------------------
        
        function bbox = constructBoundingBox(self, version)
        % Construct a Java bounding box object based on the desired limits.
            
            % Set variable for the coordinate reference system code.
            code = self.CoordRefSysCode;
            
            % Determine which set of limits to use.
            if useGeoCoords(code)
                % The XLim and YLim limits are not set, use Latlim and
                % Lonlim limits. (Geographic limits are validated to be
                % non-empty.)
                if isequal(version, '1.3.0') ...
                        && isequal(code, 'EPSG:4326') ...
                        && serverSupportsV1Point3(self.Layer(1).ServerURL, ...
                           self.Layer(1).Details.BoundingBox)                      
                    % Only for WMS version 1.3.0, and for a coordinate
                    % reference system code containing the EPSG:4326 code,
                    % the X and Y axes are swapped (provided the server
                    % actually follows the specification). From the WMS
                    % 1.3.0 specification: "EPSG:4326 refers to WGS 84
                    % geographic latitude, then longitude. That is, in this
                    % CRS the X axis corresponds to latitude, and the Y
                    % axis to longitude."
                    xLim = self.Latlim;
                    yLim = self.Lonlim;
                else
                    xLim = self.Lonlim;
                    yLim = self.Latlim;
                end
            else
                % The XLim and YLim limits are set and non-empty.
                xLim = self.XLim;
                yLim = self.YLim;
            end
            
            minx = xLim(1);
            maxx = xLim(2);
            miny = yLim(1);
            maxy = yLim(2);
            
            bbox = org.geotools.data.ows.CRSEnvelope( ...
                code, minx, miny, maxx, maxy);
        end
    end
end

%-------------------- Validation Functions --------------------------------

function validateLayer(layer)
% Validate layer as a WMSLayer array.

% Validate layer as a WMSLayer array containing at least one element.
validateattributes(layer, {'WMSLayer'}, {'nonempty'}, 'WMSMapRequest', 'layer')

% Validate the ServerURL.
internal.map.validateURL('Layer.ServerURL', layer(1).ServerURL);
end

%--------------------------------------------------------------------------

function layer = validateLayerProperties(layer, server)
% Validate and update the layer properties.

% Update layer's properties from the server if required.
layer = updateLayers(server, layer);
assert(~isempty(layer), message('map:wms:noLayersAfterUpdate'));

% Validate non-empty LayerName properties.
layerNames = {layer.LayerName};
emptyIndex = cellfun('isempty', layerNames);
assert(~all(emptyIndex), message('map:validate:expectedNonEmptyStrings',  ...
    'LayerName properties'));

% Validate Latlim and Lonlim properties.
for k=1:numel(layer)
    validateLatlim(layer(k).Latlim, layer, WMSMapRequest.DefaultCRS);
    validateLonlim(layer(k).Lonlim, layer, WMSMapRequest.DefaultCRS);
end

% Validate the WMS version. The layer at this point has been
% updated by the server. However, the updated layer may contain
% a version that is experimental or unsupported.
validateVersion(layer);
end

%--------------------------------------------------------------------------

function validateVersion(layer)
% Validate the WMS version number in the WMSLayer array, layer.

version = layer(1).Details.Version;
validateattributes(version, ...
    {'char'}, {'nonempty'}, 'WMSMapRequest','layer.Details.Version')

serverURL = layer(1).ServerURL;
url = java.net.URL(serverURL);
spec = com.mathworks.toolbox.geoweb.wms.WMSClientSpecification(url);
versions = cell(spec.getVersions());
assert(spec.containsVersion(version), message('map:wms:versionMismatch', ...
    version, 'layer.Details.Version', sprintf('''%s'' ', versions{:})));
end

%--------------------------------------------------------------------------

function server = validateServer(layer, varargin)
% Validate the input as a WebMapServer if supplied. If not supplied, then
% create the object.

if nargin == 2
    server = varargin{1};
    validateattributes(server, {'WebMapServer'}, {'scalar'}, 'WMSMapRequest', 'server');
else
    server = WebMapServer(layer(1).ServerURL);
end
end

%--------------------------------------------------------------------------

function  validateCoordRefSysCode(layer, coordRefSysCode)
% Validate coordRefSysCode as a char and available in all layers.

validateattributes(coordRefSysCode, ...
    {'char'}, {'nonempty', 'vector', 'row'}, ...
    'validateCoordRefSysCode', 'CoordRefSysCode');

% Ensure coordRefSysCode is upper case.
coordRefSysCode = upper(coordRefSysCode);

% Verify the code is available in all layers and is not an automatic
% projection code.
autoProjCode = 'AUTO';
assert(~strncmp(autoProjCode, coordRefSysCode, numel(autoProjCode)), ...
    message('map:wms:containsAutomaticProjectionCode'));
ensureLayerContainsCRSCode(layer, coordRefSysCode);
end

%--------------------------------------------------------------------------

function validateLimit(limitValue, limitName)
% Validate a coordinate limit property. limitValue is the value of the
% limit. limitName is the name of the property. limitValue may be empty.

if ~isempty(limitValue)
    
    validateattributes(limitValue, {'numeric'}, {'finite', 'real'}, ...
        'WMSMapRequest', limitName);
    
    assert(numel(limitValue) == 2, message('map:validate:expectedTwoElementVector', ...
        limitName));
    
    assert(limitValue(1) <= limitValue(2), message('map:validate:expectedAscendingOrder', ...
        limitName));
    
end
end

%--------------------------------------------------------------------------

function validateMapLimit(name, value, layer, coordRefSysCode)
% Validate a map limit property (XLim or YLim).

% Assert that CoordRefSysCode has been set to a value that does not
% indicate a geographic coordinate reference system.
assert(~useGeoCoords(coordRefSysCode), ...
    message('map:wms:expectedMapCoordRefSysCode', name));

% Validate the limit's datatype and value.
validateLimit(value, name);
assert(~isempty(value), message('map:wms:expectedNonEmptyValue',  ...
    name, 'CoordRefSysCode', WMSMapRequest.DefaultCRS));

checkNoSubsetsMatch(name, value, layer, coordRefSysCode)
end

%--------------------------------------------------------------------------

function validateGeoLimit(name, value, coordRefSysCode)
% Validate a geographic limit property (Latlim or Lonlim).

assert(useGeoCoords(coordRefSysCode), ...
    message('map:wms:expectedGeoCoordRefSysCode', name, coordRefSysCode));

validateLimit(value, name);
end

%--------------------------------------------------------------------------

function validateLatlim(latlim, layer, coordRefSysCode)
% Validate Latlim property.

name = 'Latlim';
validateGeoLimit(name, latlim, coordRefSysCode);
assert(~isempty(latlim) && all(-90 <= latlim(:)) && all(latlim(:) <= 90), ...
    message('map:validate:expectedInterval', 'Latlim', '-90', '90'));
checkNoSubsetsMatch(name, latlim, layer, coordRefSysCode)
end

%--------------------------------------------------------------------------

function validateLonlim(lonlim, layer, coordRefSysCode)
% Validate Lonlim property.

name = 'Lonlim';
validateGeoLimit(name, lonlim, coordRefSysCode);
wrappedTo180 = all(-180 <= lonlim(:)) && all(lonlim(:) <= 180);
wrappedTo360 = all(0 <= lonlim(:)) && all(lonlim(:) <= 360);
assert(~isempty(lonlim) && (wrappedTo180 || wrappedTo360), ...
    message('map:validate:expectedIntervals', ...
    'Lonlim', '[-180, 180]', '[0 360]'));
checkNoSubsetsMatch(name, lonlim, layer, coordRefSysCode)
end

%--------------------------------------------------------------------------

function validateImageLength(layer, name, value)
% Validate an image length property.

validateattributes(value, {'numeric'}, ...
    {'integer', 'positive', 'finite', 'scalar'}, ...
    'validateImageLength', name);

if strcmpi(name, 'ImageHeight')
    assert(value <= WMSMapRequest.MaximumHeight, ...
        message('map:validate:inputTooLarge', ...
        name, num2str(WMSMapRequest.MaximumHeight)));
    
    fixedHeight = getFixedHeight(layer);
    assert(isempty(fixedHeight) || value == fixedHeight, ...
        message('map:validate:expectedLengthValue', ...
        'ImageHeight', 'FixedHeight', num2str(fixedHeight)));
else
    assert(value <= WMSMapRequest.MaximumWidth, ...
        message('map:validate:inputTooLarge', ...
        name, num2str(WMSMapRequest.MaximumWidth)));
    
    fixedWidth = getFixedWidth(layer);
    assert(isempty(fixedWidth) || value == fixedWidth, ...
         message('map:validate:expectedLengthValue', ...
        'ImageWidth', 'FixedWidth', num2str(fixedWidth)));
end
end

%--------------------------------------------------------------------------

function time = validateTime(dimensionNames, time)
% Validate Time property.

% The first layer must contain 'time' in the Layer.Dimension.Name field.
% For a request with multiple layers, the dimension parameter may be
% ignored for all non-base layers. Verify that 'time' exists in the
% Layer.Dimension.Name field.
assert(any(strcmpi('time', dimensionNames)), ...
    message('map:validate:invalidFieldValue', ...
    'Layer(1).Details.Dimension.Name',  '''time'''));

if isnumeric(time)
    validateattributes(time, {'numeric'}, ...
        {'nonempty', 'scalar', 'finite', 'nonnegative'}, ...
        'validateTime', 'Time');
    time = datestr(time, 29);
else
    validateattributes(time, {'char'}, {'row', 'vector'}, ...
        'validateTime', 'Time');
    
    % Valid time may be empty or set to 'current'
    validtime = strcmpi('current', time);
    
    % Valid time may contain a year.
    validtime = validtime || (numel(time) == 4 ...
        && all(isstrprop(time, 'digit')));
    
    % If time has the character T and Z in the right location, assume the
    % entry is valid. A valid time may match either of the following forms:
    % yyyy-mm-ddTHH:MM:SSZ
    % yyyy-mm-ddTHH:MM:SS.sssZ  
    tIndex = strfind(time, 'T');
    zIndex = strfind(time, 'Z');
    tLoc = 11;
    validtime = validtime || ...
        (~isempty(tIndex) && numel(tIndex) == 1 && tIndex == tLoc) && ...
        (~isempty(zIndex) && numel(zIndex) == 1 && zIndex == numel(time));
    
    % Valid time may match either of the following patterns:
    % yyyy-mm
    % yyyy-mm-dd
    minFormat = 'yyyy-mm';
    minHyphenIndex = strfind(minFormat, '-');
    maxFormat = 'yyyy-mm-dd';
    maxHyphenIndex = strfind(maxFormat, '-');
    hyphenIndex = strfind(time, '-');
    validtime = validtime || ...
        (~isempty(hyphenIndex) ...
        && numel(minFormat) <= numel(time) ...
        && numel(time) <= numel(maxFormat) ...
        && (isequal(hyphenIndex, minHyphenIndex) ...
        ||  isequal(hyphenIndex, maxHyphenIndex)) ...
        && all(isstrprop(time(1:4),'digit'))); 
    
    % Permit the extension for geologic datasets and for years B.C.E. as
    % allowed by the WMS spec. If the first letter matches K, M, G, or B,
    % then assume the entry is valid. (The server will issue an error if
    % the time entry is not valid when using this extension.)
    validtime = validtime || ( numel(time) > 1 && ...
        (isequal(time(1), 'K') || ...
         isequal(time(1), 'M') || ...
         isequal(time(1), 'G') || ...
         isequal(time(1), 'B')));
    
    if ~validtime
        % Assume that the input is a time value that can be converted to
        % the expected ISO 8601 format.
        
        % Convert MATLAB dateform number 30 'yyyymmddTHHMMSS', if found.
        time = convertDateform30(time);
        
        % Check to see if using an extended time format consisting
        % of hours, minutes or seconds.
        extendedFormat = isempty(strfind(time, ':'));
        if extendedFormat
            % Hours, minutes and seconds are not required.
            % Use format yyyy-mm-dd
            format = 29;
            
            % Convert month-year format (mmmyyyy or mmmyy) if found.
            monthYearFormat = (numel(time) == 7 || numel(time) == 5) ...
                && all(isletter(time(1:3)));
            if monthYearFormat
               time = convertMonthYearFormat(time);
               return
            end
        else           
            % Use format yyyy-mm-dd HH:MM:SS
            format = 31;
            
            % Add the current date to the time value, if required
            time = addDateToTime(time);
        end
        
        % Convert time to required ISO 8601 format.
        time = convertToISO8601(time, format);

    end
end
end

%--------------------------------------------------------------------------

function validateElevation(dimensionNames, elevation)
% Validate Elevation property.

if ~isempty(elevation)
   validateattributes(elevation, {'char'}, {'row', 'vector'}, ...
       'validateElevation', 'Elevation');
end
assert(any(strcmpi('elevation', dimensionNames)), ...
    message('map:validate:invalidFieldValue', ...
    'Layer.Details.Dimension.Name', '''elevation'''));
end

%--------------------------------------------------------------------------

function validateSampleDimension(dimensionNames, sampleDimension)
% Validate SampleDimension property.

assert(iscellstr(sampleDimension) && numel(sampleDimension) == 2, ...
    message('map:validate:expectedTwoElementCellArrayOfStrings',  ...
    'SampleDimension'));

assert(any(strcmpi(sampleDimension{1}, dimensionNames)), ...
    message('map:validate:invalidFieldValue', ...
    'Layer.Details.Dimension.Name', sampleDimension{1}));
end

%--------------------------------------------------------------------------

function validateBackgroundColor(backgroundColor)
% Validate BackgroundColor property.

validateattributes(backgroundColor, {'numeric'}, ...
    {'nonempty', 'vector', 'finite', 'nonnegative', '<=', 255}, ...
    'validateBackgroundColor', 'BackgroundColor');

assert(numel(backgroundColor) == 3, ...
    message('map:validate:expectedThreeElementVector', 'BackgroundColor'));
end

%--------------------------------------------------------------------------

function styleName = validateStyleName(layer, styleName)
% Validate StyleName property.

if ~iscell(styleName)
    if ~isempty(styleName)
        validateattributes(styleName, {'char'}, ...
            {'row','vector'}, 'validateStyleName', 'StyleName');
    end
    cStyleName = styleName;
    
    % Replicate the char styleName into a cell array the same size as
    % layer. (If styleName is a not a vector, then an error will be thrown
    % when the cell array is validated.)
    styleName = cell(size(layer));
    [styleName{:}] = deal(cStyleName);
end

% Validate styleName as a cell array of strings. A valid element of
% styleName may be ''. Note that isvector('') returns false.
isValidInput = @(x)(ischar(x) && (isempty(x) || isvector(x)));
isCharVector = cellfun(isValidInput, styleName);
isValidStyleName = all(isCharVector) && numel(styleName) == numel(layer);
assert(isValidStyleName, message('map:validate:expectedMatchingCellSize',  ...
    'StyleName', 'Layer'));

for k=1:numel(layer)
    styleNames = {layer(k).Details.Style.Name, ''};
    assert(any(strcmpi(styleName{k}, styleNames)), ...
        message('map:wms:expectedStyleNameInDetailsField', styleName{k}, ...
        sprintf('%s(%d).%s', 'Layer', k, 'Details.Style.Name')));
end
end

%--------------------------------------------------------------------------

function validateImageFormat(layer, imageFormat)
% Validate ImageFormat property.

validateattributes(imageFormat, {'char'}, ...
    {'nonempty','vector','row'}, ...
    'validateImageFormat', 'ImageFormat');

% Assign imageFormats to the cell array of ImageFormats defined
% by the server. It is safe to use the first layer in the array
% since on a per server basis all ImageFormats are the same for
% all layers.
imageFormats = layer(1).Details.ImageFormats;

% Determine if the specified imageFormat is found in the
% Layer's ImageFormats field.
if ~any(strcmpi(imageFormat, imageFormats))
    strImageFormats = sprintf('''%s'' ',imageFormats{:});
    error(message('map:wms:imageFormatUnavailable', ...
        'ImageFormat', imageFormat, strImageFormats));
end
end

%----------------- Coordinate Referencing Helper Functions ----------------

function coordRefSysCode = assignDefaultCRS(layer, defaultCRS) 
% Assign a default coordinate reference system code.

tf = layerContainsCRSCode(layer, defaultCRS);
if tf
    % Use defaultCRS.
    coordRefSysCode = defaultCRS;
else
    % defaultCRS is not found. Use the first non-AUTO CRS code. If not
    % found, then assign the CRS code to the default CRS defined by the
    % class.
    autoProjCode = 'AUTO';
    coordRefSysCode = [];
    for k=1:numel(layer(1).CoordRefSysCodes)
        coordRefSysCode = layer(1).CoordRefSysCodes{k};
        if ~strncmp(autoProjCode, coordRefSysCode, numel(autoProjCode))
            break
        else
            coordRefSysCode = [];
        end
    end
    if isempty(coordRefSysCode)
        coordRefSysCode = WMSMapRequest.DefaultCRS;
    end
end
ensureLayerContainsCRSCode(layer, coordRefSysCode);
end
            
%--------------------------------------------------------------------------

function ensureLayerContainsCRSCode(layer, coordRefSysCode)
% Ensure all layer contain the coordRefSysCode. 

tf = layerContainsCRSCode(layer, coordRefSysCode);
assert(tf, message('map:wms:expectedCoordRefSysCode', coordRefSysCode));
end

%--------------------------------------------------------------------------

function tf = layerContainsCRSCode(layer, coordRefSysCode)
% Return true if all layers contain the coordRefSysCode. 

tf = false;
for k=1:numel(layer)
    tf = any(strcmpi(coordRefSysCode, layer(k).CoordRefSysCodes));
    tf = tf && boundingBoxContainsCRSCode(layer, k, coordRefSysCode);
    if tf
        break;
    end
end
end

%--------------------------------------------------------------------------

function tf = ...
    boundingBoxContainsCRSCode(layer, layerIndex, coordRefSysCode)
% Return true if BoundingBox contains the coordRefSysCode.

coordRefSysCode = upper(coordRefSysCode);
autoProj = 'AUTO';
containsAutoCode = strncmp(autoProj, coordRefSysCode, numel(autoProj));
containsDefaultCode = isequal(coordRefSysCode, WMSMapRequest.DefaultCRS);
containsCRSDefaultCode = isequal(coordRefSysCode, 'CRS:84');

if containsDefaultCode  || containsAutoCode || containsCRSDefaultCode
    tf = true;
else
    boundingBox = layer(layerIndex).Details.BoundingBox;
    coordRefSysCodes = {boundingBox.CoordRefSysCode};
    tf = any(strncmp(coordRefSysCode, coordRefSysCodes, numel(coordRefSysCode)));
end
end

%--------------------------------------------------------------------------

function [xLim, yLim] = getLimitsFromBoundingBox(layer, coordRefSysCode)
% Get the X and Y bounding box limits from the scalar layer object. The
% limits are specified in the layer.Details.BoundingBox.CoordRefSysCode
% field, which is matched to the value of coordRefSysCode.

% Obtain the limits from the BoundingBox field of Details.
boundingBox = layer.Details.BoundingBox;
coordRefSysCodes = {boundingBox.CoordRefSysCode};
index = find(strncmp(coordRefSysCode, coordRefSysCodes, numel(coordRefSysCode)));
assert(~isempty(index), message('map:wms:coordRefSysCodeNotFound',  ...
    coordRefSysCode, 'Layer.Details.BoundingBox.CoordRefSysCode'));

% Set xLim and yLim to their corresponding element in the boundingBox
% array.
xLim = boundingBox(index(1)).XLim;
yLim = boundingBox(index(1)).YLim;
end

%------------------ Time Validation Helper Functions ----------------------

function time = convertDateform30(time)
% Convert dateform 30 yyyymmddTHHMMSS to dateform 31 yyyy-mm-dd HH:MM:SS

format = 'yyyymmddTHHMMSS';
tIndex = strfind(time,'T');
needsConversion = numel(tIndex) == 1 && tIndex == 9 ...
    && numel(time) == numel(format);
if needsConversion
    time = [time(1:4) '-' time(5:6) '-' time(7:8) ' ' ...
        time(10:11) ':' time(12:13) ':' time(14:15)];
end
end

%--------------------------------------------------------------------------

function time = addDateToTime(time)
% Add the current date to the time value, if required.

longestTime = 'HH:MM:SS PM';
if numel(time) <= numel(longestTime)
    % Convert today to yyyy-mm-dd format.
    today = datestr(now, 29);
    % Add to the time value
    time = [today ' ' time];
end
end

%--------------------------------------------------------------------------

function time = convertMonthYearFormat(time)
% Convert month-year formats mmmyyyy or mmmyy to yyyy-mm format.

day = '.01,';
if numel(time) == 7 
    % Input is in form: mmmyyyy
    % Convert to form:  mmm.dd,yyyy
    time = [time(1:3) day time(4:end)];
else
    % Input is in form: mmmyy
    % Convert to form:  mmm.dd,yyyy
    today = datestr(now, 29);
    year = today(1:2);
    time = [time(1:3) day year time(4:end)];
end
format = 22;  %mmm.dd,yyyy
try
    time = datestr(time, format);
    time = datestr(time, 29);
    time = time(1:7);  % yyyy-mm
catch e %#ok<NASGU>
    error(message('map:wms:invalidTimeEntry'))
end
end
    
%--------------------------------------------------------------------------

function time = convertToISO8601(time, format)
% Convert time string to ISO 8601 format.

try
    % Convert input to ISO 8601 format
    time = datestr(time, format);
    if format == 31
        % When using extended format 31 ('yyyy-mm-dd HH:MM:SS'),
        % the ' ' character needs to be replaced with 'T' and 'Z'
        % needs to be appended.
        % Find the space character.
        blank = strfind(time, ' ');
        
        % Wrap with an if statement in case there is an unexpected
        % problem with datestr.
        if ~isempty(blank)
            time(blank(1)) = 'T';
            time(end+1) = 'Z';
        end
    end
catch e %#ok<NASGU>
    error(message('map:wms:invalidTimeEntry'));
end
end

%--------------------- Dimension Helper Function --------------------------

function dimensionNames = getDimensionNames(layer)
% Return all the Dimension names.

dimensionNames = {};
for k=1:numel(layer)
    dimensionNames = ...
        [dimensionNames ...
        {layer(k).Details.Dimension.Name}];  %#ok<AGROW>
end
end

%--------------------- Layer Attributes Helper Functions ------------------
        
function tf = isFixedHeight(layer)
% Return true for each element of layer that contains a FixedHeight
% attribute.

tf = false(size(layer));
for k=1:numel(layer)
    tf(k) = (layer(k).Details.Attributes.FixedHeight > 0);
end
end
        
%--------------------------------------------------------------------------  

function tf = isFixedWidth(layer)
% Return true for each element of layer that contains a FixedWidth
% attribute.

tf = false(size(layer));
for k=1:numel(layer)
    tf(k) = (layer(k).Details.Attributes.FixedWidth > 0);
end
end
%--------------------------------------------------------------------------

function tf = isNoSubsets(layer)
% Return true for each element of layer that contains a true NoSubsets
% attribute.

tf = false(size(layer));
for k=1:numel(layer)
    tf(k) = layer(k).Details.Attributes.NoSubsets;
end
end

%--------------------------------------------------------------------------

function checkNoSubsetsMatch(name, value, layer, coordRefSysCode)
% Check if sub-setting the layer is permitted. If not, validate that the
% limits of the request, stored in VALUE, are set to the limits of the
% layer. NAME must be set to 'Latlim', 'Lonlim', 'XLim', or 'YLim' (case
% does not matter). 

haveNoSubsets = isNoSubsets(layer);
if any(haveNoSubsets)
    layer = layer(haveNoSubsets);
    switch lower(name)
        case 'latlim'
            for k=1:numel(layer)
                assert(isequal(value, layer(k).Latlim), ...
                    message('map:wms:geographicSubsetNotPermitted', ...
                    layer(k).LayerName, 'Latlim', ...
                    num2str(layer(k).Latlim(1)), num2str(layer(k).Latlim(2))));
            end
        case 'lonlim'
            for k=1:numel(layer)
                assert(isequal(value, layer(k).Lonlim), ...
                    message('map:wms:geographicSubsetNotPermitted', ...
                    layer(k).LayerName, 'Lonlim', ...
                    num2str(layer(k).Lonlim(1)), num2str(layer(k).Lonlim(2))));
            end
            
        case 'xlim'
            for k=1:numel(layer)
                xLim = getLimitsFromBoundingBox(layer(k), coordRefSysCode);
                assert(isequal(value, xLim), ...
                    message('map:wms:spatialSubsetNotPermitted', ...
                    layer(k).LayerName, 'XLim', num2str(xLim(1)), num2str(xLim(2))));
            end
            
        case 'ylim'
            for k=1:numel(layer)
                [~, yLim] = ...
                    getLimitsFromBoundingBox(layer(k), coordRefSysCode);
                assert(isequal(value, yLim), ...
                    message('map:wms:spatialSubsetNotPermitted', ...
                    layer(k).LayerName, 'YLim', num2str(yLim(1)), num2str(yLim(2))));
            end
    end
end
end

%--------------------- Image Size Helper Functions ------------------------

function value = getFixedHeight(layer)
% If a FixedHeight attribute is set and equal for all layers, then
% return its value; otherwise return [].

if any(isFixedHeight(layer))
    value = layer(1).Details.Attributes.FixedHeight;
    for k=2:numel(layer)
        assert(value == layer(k).Details.Attributes.FixedHeight, ...
            message('map:wms:expectedEqualFieldValues', ...
            'Layer.Details.Attributes.FixedHeight'));
    end
else
    value = [];
end
end

%--------------------------------------------------------------------------

function value = getFixedWidth(layer)
% If a FixedWidth attribute is set and equal for all layers, then
% return its value; otherwise, return [].

if any(isFixedWidth(layer))
    value = layer(1).Details.Attributes.FixedWidth;
    for k=2:numel(layer)
        assert(value == layer(k).Details.Attributes.FixedWidth, ...
            message('map:wms:expectedEqualFieldValues', ...
            'Layer.Details.Attributes.FixedWidth'));
    end
else
    value = [];
end
end

%-------------------------------------------------------------------------- 

function [imageHeight, imageWidth] = setFixedImageSize(layer)
% Set ImageWidth and ImageHeight properties based on fixed size
% attributes.

fixedWidth = getFixedWidth(layer);
if ~isempty(fixedWidth)
    imageWidth  = fixedWidth;
else
    imageWidth = WMSMapRequest.DefaultImageLength;
end

fixedHeight = getFixedHeight(layer);
if ~isempty(fixedHeight)
    imageHeight = fixedHeight;
else
    imageHeight = WMSMapRequest.DefaultImageLength;
end
end

%-------------------- Image Format Helper Functions -----------------------

function imageFormat = findImageFormat(layer)
% Find the default ImageFormat specifier from ImageRenderFormats property.

% Assign imageFormats to the cell array of ImageFormats defined by the
% Layer. It is safe to use the first layer in the array since on a per
% server basis all ImageFormats are the same for all layers.
imageFormats  = layer(1).Details.ImageFormats;
renderFormats = WMSMapRequest.ImageRenderFormats;

% Find imageFormat in the list of ImageRenderFormats.
index = ismember(renderFormats, imageFormats);
imageFormat = renderFormats(index);
if ~isempty(imageFormat)
    % A valid image format has been found from the ImageRenderFormats
    % property.
    imageFormat = imageFormat{1};
else
    % The server does not support a format that can render a map into an
    % image (a picture format).  In this case, issue an error. In theory,
    % this branch should never be reached.
    serverImageFormats = sprintf('%s ', imageFormats{:});
    imageRenderFormats = sprintf('%s ', renderFormats{:});
    error(message('map:wms:noImageFormatFoundOnServer', ...
        serverImageFormats, imageRenderFormats));
end
end

%--------------------------------------------------------------------------

function imageFormat = findTransparentImageFormat(layer, imageFormat)
% Find the ImageFormat specifier to support transparent images.

% findTransparentImageFormat returns the string imageFormat based on the
% following algorithm.
%   1) if the object's ImageFormat property is set to the default value,
%   then return the ImageFormat property; otherwise:
%
%   2) determine if the Layer property contains an ImageFormat in the cell
%   array, ImageTransparentFormats. If found, return the first match. If
%   not found, then return the value of the ImageFormat property.

% Determine the default value of the ImageFormat property.
defaultImageFormat = findImageFormat(layer);
if isequal(defaultImageFormat, imageFormat)
    % The ImageFormat property has not been set by the user.
    % Assign to imageFormats the cell array of ImageFormats defined by the
    % server. It is safe to use the first layer in the array since on a per
    % server basis all ImageFormats are the same for all layer.
    imageFormats = layer(1).Details.ImageFormats;
    
    % Find imageFormat from the list of ImageTransparentFormats.
    transparentFormats = WMSMapRequest.ImageTransparentFormats;
    index = ismember(transparentFormats, imageFormats);
    transparentFormats = transparentFormats(index);
    if ~isempty(transparentFormats)
        % imageFormat has been found in ImageTransparentFormats.
        imageFormat = transparentFormats{1};
    else
        % imageFormat has not been found in ImageTransparentFormats, so
        % return the property's value.
    end
else
    % The ImageFormat property has been set by the user, so there is no
    % need to change the property's value.
end
end

%--------------------------------------------------------------------------
        
function issueWarnings(oldRequest, newRequest)
% Issue warnings based on values of oldRequest and newRequest.

if ~isequal(newRequest.CoordRefSysCode, oldRequest.CoordRefSysCode)
    warning(message('map:wms:adjustingCoordRefSysCode', ...
        oldRequest.CoordRefSysCode, newRequest.CoordRefSysCode))
end

if oldRequest.ImageHeight ~= newRequest.ImageHeight
    warning(message('map:wms:adjustingImageHeight', ...
        newRequest.ImageHeight))
end

if oldRequest.ImageWidth ~= newRequest.ImageWidth
    warning(message('map:wms:adjustingImageWidth', ...
        newRequest.ImageWidth))
end

if ~isequal(newRequest.Elevation, oldRequest.Elevation)
    warning(message('map:wms:adjustingElevation', ...
        newRequest.Elevation))
end

if ~isequal(newRequest.Time, oldRequest.Time)
    warning(message('map:wms:adjustingTime', ...
        newRequest.Time))
end

if ~isequal(newRequest.SampleDimension, oldRequest.SampleDimension)
    warning(message('map:wms:adjustingSampleDimension', ...
        oldRequest.SampleDimension{1}, '{'' '','' ''}'))
end

if ~isequal(newRequest.StyleName, oldRequest.StyleName)
    warning(message('map:wms:adjustingStyleName', '{''''}'))
end
end

%--------------------------------------------------------------------------

function tf = useGeoCoords(coordRefSysCode)
% Return true if the code indicates a geographic coordinate reference
% system, otherwise, return false.

tf = ismember(coordRefSysCode, {WMSMapRequest.DefaultCRS, 'CRS:84'});
end

%--------------------------------------------------------------------------

function tf = serverSupportsV1Point3(serverURL, BoundingBox)
% Return true if the server truly supports version 1.3.x.

serversNotSupportingV1Point3 = { ...
    'http://webapps.datafed.net', ...
    'http://www.nasa.network.com'};
query = regexptranslate('wildcard', serversNotSupportingV1Point3);
tf = all(cellfun(@isempty,regexpi(serverURL, query)));

if tf
    % Check the BoundingBox to see if the order is switched. As indicated
    % by the WMS Version 1.3.0 specification:
    % "EPSG:4326 refers to WGS 84 geographic latitude, then longitude. That
    % is, in this CRS the X axis corresponds to latitude, and the Y axis to
    % longitude. In which case, XLim is latitude and YLim is longitude."
    tf = swapLimitsForV1_3(BoundingBox);
end
end

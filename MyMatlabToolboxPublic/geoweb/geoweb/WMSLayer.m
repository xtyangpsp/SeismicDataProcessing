%WMSLayer Web Map Service layer object
%
%   A WMSLayer object describes a Web Map Service layer or layers. The
%   function WMSFIND returns a WMSLayer array, and the function WMSINFO 
%   returns a WMSCapabilities object that contains a WMSLayer array in its
%   Layer property.
% 
%   LAYER = WMSLayer(PARAM1, VAL1, PARAM2, VAL2, ...) constructs a WMSLayer
%   object from the input parameter names and values. If a parameter name
%   matches a property name of the class -- ServerTitle, ServerURL,
%   LayerTitle, LayerName, Latlim, Lonlim, Abstract, CoordRefSysCodes, or
%   Details -- then the parameter's values are copied to the property.  If
%   a property is not set, it defaults to '' for the ServerTitle,
%   ServerURL, LayerTitle, LayerName, and Abstract properties; [] for the
%   Latlim and Lonlim properties; {} for the CoordRefSysCodes property; and
%   a structure for the Details properties. If supplied, ServerTitle,
%   ServerURL, LayerTitle, LayerName, and Abstract must be strings; Latlim
%   and Lonlim must be real, finite, two-element vectors; CoordRefSysCodes
%   must be a cell array containing a cell array of strings; and Details
%   must be a structure. The size of LAYERS is scalar unless all inputs are
%   cell arrays, in which case, the size of LAYERS matches the size of the
%   cell arrays.
%
%   All properties of the class may not be set, except the 'Latlim' and
%   'Lonlim' properties. These two properties have public set access.
%
%   WMSLayer properties:
%      ServerTitle - Title of server
%      ServerURL - URL of server
%      LayerTitle - Title of layer
%      LayerName - Name of layer
%      Latlim - Latitude limits of layer
%      Lonlim - Longitude limits of layer
%      Abstract - Information about layer
%      CoordRefSysCodes - Coordinate reference system codes
%      Details - Detailed information about layer     
%
%   WMSLayer methods:
%      disp - Display properties
%      refine - Refine search
%      refineLimits - Refine search based on geographic limits
%      serverTitles - Titles of unique servers
%      servers - Unique server URLs
%
%   Example
%   -------
%   % Construct a WMSLayer object from a WMS GetMap request URL.
%   serverURL = ['http://ims.cr.usgs.gov:80/wmsconnector/' ...
%      'com.esri.wms.Esrimap/USGS_EDC_LandCover_NLCD2001?'];
%   requestURL = [serverURL 'SERVICE=WMS&FORMAT=image/jpeg&REQUEST=GetMap&' ...
%      'STYLES=&SRS=EPSG:4326&VERSION=1.1.1&', ...
%      'LAYERS=NLCD_2001_Land_Cover&', ...
%      'WIDTH=1024&HEIGHT=470&' ...
%      'BBOX=-128,23,-65,51&'];
%
%   % Construct the WMSLayer object by using the serverURL variable and 
%   % the value of the WMS LAYERS parameter. Update the remaining 
%   % information from the server.
%   layer = WMSLayer('ServerURL', serverURL, ...
%     'LayerName', 'NLCD_2001_Land_Cover');
%   layer = wmsupdate(layer)
%
%   % Request the layer from the server using the WMS request parameters 
%   % found in the requestURL string. Copy the WMS BBOX information to the
%   % latlim and lonlim parameters. Copy the WMS WIDTH and HEIGHT values to
%   % the ImageWidth and ImageHeight parameters.
%   lonlim = [-128, -65];
%   latlim = [23, 51];
%   [A, R] = wmsread(layer, 'Latlim', latlim, 'Lonlim', lonlim, ...
%      'ImageHeight', 470 , 'ImageWidth', 1024);
%
%   % Display the rendered image.
%   figure
%   usamap(A, R)
%   geoshow(A, R)
%
%   % The image can also be retrieved using the requestURL.
%   [A, R] = wmsread(requestURL);
%
%   See also WebMapServer, WMSCapabilities, WMSFIND, WMSINFO, 
%            WMSMapRequest, WMSREAD, WMSUPDATE.

% Copyright 2008-2011 The MathWorks, Inc.
% $Revision: 1.1.6.18 $  $Date: 2011/03/28 04:29:28 $

classdef WMSLayer
    
    properties (SetAccess='private', GetAccess='public')
        
        %ServerTitle Title of server
        %
        %   ServerTitle is a string indicating the title of the server.
        ServerTitle = '';
        
        %ServerURL URL of server
        %
        %   ServerURL is a string indicating the URL of the server.
        ServerURL = '';
        
        %LayerTitle Title of layer
        %
        %   LayerTitle is a string indicating the title of the layer.
        LayerTitle = '';

        %LayerName Name of layer
        %
        %   LayerName is a string indicating the name of the layer.
        LayerName = '';
    end
    
    properties (Access = public,  Dependent = true)             
        %Latlim Latitude limits of layer
        %
        %   Latlim is a two-element vector specifying the southern and
        %   northern latitude limits of the layer in units of degrees and
        %   in the range [-90, 90]. The Latlim property may be set.
        Latlim = [];
        
        %Lonlim Longitude limits of layer
        %
        %   Lonlim is a two-element vector specifying the western and
        %   eastern longitude limits of the layer in units of degrees. 
        %   The limits must be ascending and in the range [-180, 180]
        %   or [0 360]. The Lonlim property may be set.
        Lonlim = [];
    end
    
    properties (SetAccess='private', GetAccess='public')        
        %Abstract Information about layer
        %
        %   Abstract is a string containing information about the layer.
        Abstract = '';
        
        %  
        %CoordRefSysCodes Coordinate reference system codes
        %
        %   CoordRefSysCodes is a cell array specifying the codes of the
        %   available coordinate reference systems.
        CoordRefSysCodes = {};
        
        %Details Detailed information about layer
        %
        %   Details is a structure containing detailed information about
        %   the layer. 
        %
        %   Details structure 
        %   -----------------
        %   The Details structure contains the following fields:
        %
        %   Name               Data Type     Field Content
        %   ----               ---------     -------------  
        %   MetadataURL        String        URL containing metadata 
        %                                    information about layer
        %
        %   Attributes         Structure     Attributes of layer
        %
        %   BoundingBox        Structure     Bounding box of layer
        %                      array
        % 
        %   Dimension          Structure     Dimensional parameters of
        %                      array         layer, such as time or
        %                                    elevation
        %
        %   ImageFormats       Cell array    Image formats supported by
        %                                    server
        %
        %   ScaleLimits        Structure     Scale limits of layer
        %
        %   Style              Structure     Style parameters which 
        %                      array         determine layer rendering
        %
        %   Version            String        WMS version specification
        %
        %   Attributes Structure
        %   --------------------
        %   The Attributes structure contains the following fields: 
        %
        %   Name               Data Type     Field Content
        %   ----               ---------     -------------  
        %   Queryable          Logical       True if the layer can be
        %                                    queried for feature
        %                                    information                                      
        % 
        %   Cascaded           Double        Number of times the layer has
        %                                    been retransmitted by a
        %                                    Cascading Map server
        % 
        %   Opaque             Logical       True if the map data are
        %                                    mostly or completely opaque 
        %
        %   NoSubsets          Logical       True if the map must contain
        %                                    the full bounding box, false
        %                                    if the map can be a subset of
        %                                    the full bounding box
        % 
        %   FixedWidth         Logical       True if the map has a fixed 
        %                                    width that cannot be changed
        %                                    by the server, false if the
        %                                    server can resize the map to
        %                                    an arbitrary width
        % 
        %   FixedHeight        Logical       True if the map has a fixed 
        %                                    height that cannot be changed
        %                                    by the server, false if the
        %                                    server can resize the map to
        %                                    an arbitrary height  
        %
        %   BoundingBox Structure
        %   ---------------------
        %   The BoundingBox structure array contains the following fields:
        %
        %   Name               Data Type     Field Content
        %   ----               ---------     -------------  
        %   CoordRefSysCode    String        Code number for coordinate
        %                                    reference system
        %
        %   XLim               Double array  X limit of layer in units 
        %                                    of coordinate reference system
        %
        %   YLim               Double array  Y limit of layer in units  
        %                                    of coordinate reference system
        %
        %   Dimension Structure
        %   -------------------
        %   The Dimension structure contains the following fields: 
        %
        %   Name               Data Type     Field Content
        %   ----               ---------     -------------  
        %   Name               String        Name of the dimension; such as 
        %                                    time, elevation, or
        %                                    temperature
        % 
        %   Units              String        Measurement unit
        %
        %   UnitSymbol         String        Symbol for unit
        % 
        %   Default            String        Default dimension setting;
        %                                    e.g. if Default is 'time',
        %                                    server returns time holding if
        %                                    dimension is not specified
        %
        %   MultipleValues     Logical       True if multiple values of the 
        %                                    dimension may be requested,
        %                                    false if only single values
        %                                    may be requested
        %
        %   NearestValue       Logical       True if nearest value of
        %                                    dimension is returned in
        %                                    response to request for nearby
        %                                    value, false if request value
        %                                    must correspond exactly to
        %                                    declared extent values
        %
        %   Current            Logical       True if temporal data are kept 
        %                                    current (valid only for
        %                                    temporal extents)
        % 
        %   Extent             String        Values for dimension.
        %                                    Expressed as single value
        %                                    (value), list of values 
        %                                    (value1, value2, ...),
        %                                    interval defined by bounds and
        %                                    resolution (min1/max1/res1),
        %                                    or list of intervals 
        %                                    (min1/max1/res1,
        %                                    min2/max2/res2, ...)
        %
        %   ScaleLimits Structure
        %   ---------------------
        %   The ScaleLimits structure contains the following fields:
        %
        %   Name               Data Type     Field Content
        %   ----               ---------     -------------  
        %   ScaleHint          Double array  Minimum and maximum scales 
        %                                    for which it is appropriate to
        %                                    display layer (expressed as
        %                                    scale of ground distance in
        %                                    meters represented by diagonal
        %                                    of central pixel in image)
        % 
        %   MinScaleDenominator  Double      Minimum scale denominator of 
        %                                    maps for which a layer is
        %                                    appropriate
        % 
        %   MaxScaleDenominator  Double      Maximum scale denominator of 
        %                                    maps for which a layer is
        %                                    appropriate
        %
        %   Style Structure
        %   -------------------
        %   The Style structure array contains the following fields:
        %
        %   Name               Data Type     Field Content
        %   ----               ---------     -------------   
        %   Title              String        Descriptive title of style
        %  
        %   Name               String        Name of style
        %  
        %   Abstract           String        Information about style  
        %
        %   LegendURL          Structure     Information about legend 
        %                                    graphics
        %
        %   LegendURL Structure
        %   -------------------
        %   The LegendURL structure contains the following fields:
        %
        %   Name               Data Type     Field Content
        %   ----               ---------     -------------  
        %   OnlineResource     String        URL of legend graphics
        %
        %   Format             String        Format of legend graphics
        %
        %   Height             Double        Height of legend graphics
        %
        %   Width              Double        Width of legend graphics
        Details = createDetails();
    end
    
    properties (Hidden=true, Constant=true)
        %UpdateString Update string
        %
        %   UpdateString is a string containing the value used for a
        %   property that requires an update.
        UpdateString = '<Update using WMSUPDATE>';       
    end
    
    properties (Access = private, Hidden = true)
        % Each of the following private properties holds the value of a
        % corresponding dependent property.
        pLatlim
        pLonlim
    end
    
    methods
        
        function self = WMSLayer(varargin) 
        % WMSLayer Constuctor. 
        
            if nargin >= 1
                if isempty(varargin{1})
                    % The input is empty, create an empty object.
                    self = self(false);
                else                    
                    if nargin == 1 && isstruct(varargin{1})
                        % Undocumented syntax to allow a structure as the
                        % first argument. Assign S to the first element.
                        S = varargin{1};
                    else
                        % Parse the parameter/value pairs from VARARGIN. 
                        % Copy the parameter/value pairs to S, with the
                        % fieldnames of S set to the parameter names. The
                        % fields of S that match the property names are
                        % copied to the class.
                        S = parseInputs(varargin);                       
                    end
                    
                    % Set the properties of the object from the structure
                    % input.
                    nocheck = false;
                    self = setPropertiesFromStruct(self, S, nocheck);
                end
            end
        end
        
        %------------------------- set/get methods ------------------------
        
        function self = set.Latlim(self, latlim)
            validateLimits('Latlim', latlim)
            self.pLatlim = latlim;
        end
        
        function value = get.Latlim(self)
            value = self.pLatlim;
        end
        
        function self = set.Lonlim(self, lonlim)
            validateLimits('Lonlim', lonlim)
            self.pLonlim = lonlim;
        end
        
        function value = get.Lonlim(self)
            value = self.pLonlim;
        end
        
        function value = servers(self) 
        %servers Unique server URLs
        %
        %   servers = layers.servers returns a cell array of URLs of
        %   the unique servers.
        %
        %   Example
        %   -------
        %   % Find all the unique URLs of government servers.
        %   layers = wmsfind('*.gov*','SearchField', 'serverurl');
        %   servers = layers.servers  
            
            value = unique({self.ServerURL});
        end

        function value = serverTitles(self) 
        %serverTitles Titles of unique servers
        %
        %   serverTitles = layers.serverTitles returns a cell array of 
        %   the titles of the unique servers.
        %
        %   Example
        %   -------
        %   % List the titles of all the unique government servers.
        %   layers = wmsfind('*.gov*','SearchField', 'serverurl');
        %   titles = layers.serverTitles'

            [~, index] = unique({self.ServerURL});
            value = {self(index).ServerTitle};
        end
                    
        function disp(self, varargin) 
        %disp Display properties
        %
        %   layers.disp(..., PARAM1, VAL1, PARAM2, VAL2, ...) is overloaded
        %   to display the index number of the layer followed by the
        %   property names and values.  
        %
        %   Optional parameters for disp can be used to modify the output
        %   display.  Parameter names can be abbreviated and case does not
        %   matter. The parameters names and their permissible values are
        %   listed below:
        %
        %   'Properties'    String or cell array of strings that determines
        %                   the output order and whether a property is
        %                   listed in the display output. Permissible
        %                   values are: 'servertitle', 'servername',
        %                   'layertitle', 'layername', 'latlim', 'lonlim',
        %                   'abstract', 'coordrefsyscodes', 'details', 
        %                   'all', or 'populated'. The values can be
        %                   abbreviated and case does not matter. If
        %                   'Properties' is 'all', then all the properties
        %                   are listed. If 'Properties' is 'populated', 
        %                   then only the populated properties are
        %                   displayed. The default value is 'populated'.
        %
        %   'Label'         Case-insensitive string with permissible values
        %                   of 'on' or 'off'. If 'Label' is 'on', then the
        %                   property name is listed followed by its value.
        %                   If 'Label' is 'off', then only the property
        %                   value is listed in the output. The default
        %                   value is 'on'.
        %
        %   'Index'         Case-insensitive string with permissible values
        %                   of 'on' or 'off'. If 'Index' is 'on', then the
        %                   element's index is listed in the output. If
        %                   'Index' is 'off', then the index value is not
        %                   listed in the output. The default value is
        %                   'on'.
        %
        %   Examples
        %   --------
        %   % Display LayerTitle and LayerName properties to the command
        %   % window without an Index.
        %   layers = wmsfind('srtm30');
        %   layers(1:5).disp( 'Index', 'off', ...
        %      'Properties',{'layertitle','layername'});
        %
        %   % Sort and display the LayerName property and the index.
        %   layers = wmsfind('elevation');
        %   [layerNames, index] = sort({layers.LayerName});
        %   layers = layers(index);
        %   layers.disp('Label','off', 'Properties', 'layername');
            
            % Obtain the command line options.
            options = parseDispOptions(properties(self), varargin{:});
                       
            % Print the elements to the command window.
            printLayers(self, options);                      
        end
                
        function refined = refine(self, queryStr, varargin)  
        %refine Refine search
        %
        %   refined = layers.refine(QUERYSTR, PARAM1, VAL1, PARAM2, VAL2,
        %   ...) searches for elements of layers in which values of the
        %   Layer or LayerName properties match the string QUERYSTR.
        %   Partial matching is used by default, but this can be controlled
        %   with the 'MatchType' property. 
        %
        %   You can modify the search by specifying any of the following
        %   optional parameters. Parameter names can be abbreviated and
        %   case does not matter.
        %
        %   'SearchFields'  Case-insensitive string or cell array of 
        %                   strings. Valid strings are 'abstract', 'layer',
        %                   'layertitle', 'layername', 'server',
        %                   'serverurl', 'servertitle', or 'any'.  The
        %                   default value is {'layer'}.
        %
        %                   The method searches the entries in the
        %                   'SearchFields' properties of layers for a
        %                   partial match of the entry with QUERYSTR. The
        %                   layer information is returned if any supplied
        %                   'SearchFields' match.  If 'layer' is supplied
        %                   then both the 'LayerTitle' and 'LayerName'
        %                   properties are searched. If 'server' is
        %                   supplied, then both the 'ServerURL' and
        %                   'ServerTitle' fields are searched. If 'any' is
        %                   supplied, then the properties, 'Abstract',
        %                   'LayerTitle', 'LayerName', 'ServerURL', and
        %                   'ServerTitle' are searched.
        %
        %   'MatchType'     Case-insensitive string with value 'partial' or
        %                   'exact'. The default value is 'partial'.
        %
        %                   If 'MatchType' is 'partial', then a match is
        %                   determined if the query string is found in the
        %                   property value.  If 'MatchType' is 'exact',
        %                   then a match is determined only if the query
        %                   string exactly matches the property value. If
        %                   'MatchType' is 'exact' and QUERYSTR is '*', a
        %                   match is found if the property value matches
        %                   the character '*'.
        %
        %   'IgnoreCase'    Logical with default value of true. If
        %                   'IgnoreCase' is true, then case is ignored when
        %                   performing string comparisons.
        %
        %   Example 1
        %   ---------
        %   % For each server that contains a temperature layer, list the
        %   % server URL and the number of temperature layers.
        %   temperature = wmsfind('temperature');
        %   servers = temperature.servers;
        %   for k=1:numel(servers)
        %      querystr = servers{k};
        %      layers = temperature.refine(querystr, ...
        %         'SearchFields', 'serverurl');
        %      fprintf('Server URL\n%s\n', layers(1).ServerURL);
        %      fprintf('Number of layers: %d\n\n', numel(layers));
        %   end
        %
        %   Example 2
        %   ---------
        %   % Refine a search of temperature layers to find two different
        %   % sets of layers:
        %   %    (1) layers containing only annual sea surface temperatures
        %   %    (2) layers containing annual temperatures or sea surface 
        %   %        temperatures
        %   temperature = wmsfind('temperature');
        %   annual = temperature.refine('annual');
        %   sst = temperature.refine('sea surface');
        %   annual_and_sst = sst.refine('annual');
        %   annual_or_sst = [sst;annual];
             
            if ~isempty(self) && nargin > 1                     

                % Validate queryStr.
                queryStr = validateQueryStr(queryStr);
                
                % Parse the optional parameters.
                options = parseRefineOptions(varargin{:});
                
                % Create a logical index initialized to false for each
                % layer. It is toggled to true for each layer that matches
                % the query after searching each property in
                % options.SearchFields.
                index = false(1, numel(self));
                
                % Search each property in options.SearchFields and return
                % matched layers.
                for k=1:numel(options.SearchFields)
                    propertyName = options.SearchFields{k};
                            
                    % Create a cell array of the values from the desired
                    % property.
                    propertyValues = {self.(propertyName)};
                
                    % Determine the query function.
                    queryFcn = determineQueryFcn( ...
                       options.MatchType, options.IgnoreCase);
                 
                   % Apply the query function to search propertyValues for
                   % queryStr. queryIndex is a logical array, set to true
                   % for all elements of propertyValues that match
                   % queryStr.
                   queryIndex = queryFcn(queryStr, propertyValues);
                   
                   % Toggle index to true for the layers that match the
                   % query.
                   index = queryIndex | index;
                end
                
                % Return the results.
                if any(index)
                    refined = self(index);
                else
                    % No results match queryStr, return size 0-by-0 array.
                    refined = self(false);      
                end
            else   
                % This branch is reached if the object is empty or if no
                % parameters have been supplied to the method. Under most
                % conditions, the object can not be empty, unless the
                % object is created by code such as:
                %    layer(true) = [];
                % If the object is empty, then refine can not perform a
                % query regardless of the inputs, so in this case return
                % the input object.  Rather than issuing an error for the
                % case where queryStr is not supplied, return the input
                % object.
                refined = self;
            end
        end
        
        function refined = refineLimits(self, varargin) 
        %refineLimits Refine search based on geographic limits
        %
        %   refined = layers.refineLimits(PARAM1, VAL1, PARAM2, VAL2, ...)
        %   searches for elements of layers that match specific latitude
        %   and/or longitude limits.  A given layer is included in the
        %   results only if its boundary quadrangle (as defined by the
        %   Latlim and Lonlim properties) is fully contained in the
        %   quadrangle specified by the optional 'Latlim' or 'Lonlim'
        %   parameters.  Partial overlap does not result in a match.
        %
        %   Both 'Latlim' and 'Lonlim' property names are optional.  Their
        %   names can be abbreviated and case does not matter. 
        %
        %   'Latlim'     A two-element vector of latitude specifying the
        %                latitudinal limits of the search in the form
        %                [southern_limit northern_limit] or a scalar value
        %                representing the latitude of a single point.  
        %
        %   'Lonlim'     A two-element vector of longitude specifying the
        %                longitudinal limits of the search in the form
        %                [western_limit eastern_limit] or a scalar value
        %                representing the longitude of a single point.
        %
        %   The default value of [] for either 'Latlim' or 'Lonlim' implies
        %   that all layers match the criteria. For example, if the
        %   following is specified 
        %
        %       layer.refineLimits('Latlim', [0 90], 'Lonlim', []) 
        %
        %   then all layers that cover the northern hemisphere are included
        %   in the results. 
        %
        %   All angles are in units of degrees.  
        %
        %   Example 1
        %   ---------
        %   % Find layers containing global elevation data.
        %   elevation = wmsfind('elevation');
        %   latlim = [-90, 90];
        %   lonlim = [-180, 180];
        %   globalElevation = ...
        %      elevation.refineLimits('Latlim', latlim, 'Lonlim', lonlim);
        %   % Print out the server titles from the unique servers.
        %   globalElevation.serverTitles'
                        
            if ~isempty(self) && ~isempty(varargin)                
                refined = self;
                
                % Parse the optional parameters.
                options = parseRefineLimitsOptions(varargin{:});
                
                % Refining a layer with empty Latlim or Lonlim properties
                % results in an empty result. Remove layers that contain
                % empty Latlim or Lonlim properties.
                emptyLimits = false(size(refined));
                for k=1:numel(refined)
                    emptyLimits(k) = ...
                        isempty([refined(k).Latlim]) || ...
                        isempty([refined(k).Lonlim]);
                end
                refined = refined(~emptyLimits);
                
                % Find the quad limits that are contained in the layer
                % limits.
                if ~isempty(options.Latlim) || ~isempty(options.Lonlim)
                    containsQuad = layerContainsQuad( ...
                        refined, options.Latlim, options.Lonlim);
                    refined = refined(containsQuad);
                end 
                
                % If the array is empty, but the size is not 0-by-0, then
                % return a 0-by-0 array.
                if isempty(refined) && sum(size(refined)) ~= 0
                    refined = refined(false);
                end
            else
                % If the object or varargin is empty, then refineLimits can
                % not perform a query, so in this case return the input
                % object.  
                refined = self;
            end
        end
    end
    
    methods (Access = public, Hidden = true)
        
        function self = setPropertiesFromWMSDB(self, ...
                serverTitle, serverURL, layerTitle, layerName, ...
                latlim, lonlim)
         % Create a WMSLayer array from the input parameters. All inputs
         % are cell arrays except latlim and lonlim which are M-by-2 double
         % arrays. 
         %
         % FOR INTERNAL USE ONLY -- This method is intentionally hidden and
         % is intended for use only within other toolbox classes and
         % functions. Its behavior may change, or the method itself may be
         % removed in a future release.
          
            % Create an array the same size as the input cell array. Set
            % the Abstract, CoordRefSysCodes and Details properties to
            % contain the update string. These property values are
            % propagated through the entire array.
            self.Abstract = WMSLayer.UpdateString;
            self.CoordRefSysCodes = WMSLayer.UpdateString;
            self.Details = WMSLayer.UpdateString;
            self(1:numel(serverURL),1) = self;

            % Copy the string values from the cell arrays to the
            % corresponding property.
            [self.ServerTitle] = deal(serverTitle{:});
            [self.ServerURL] = deal(serverURL{:});
            [self.LayerTitle] = deal(layerTitle{:});
            [self.LayerName] = deal(layerName{:});
            
            % Copy the numeric values.
            for k=1:numel(self)
                self(k).pLatlim = latlim(k,:);
                self(k).pLonlim = lonlim(k,:);
            end
        end
        
        %------------------------------------------------------------------
        
        function self = setPropertiesFromStruct(self, S, nocheck)
        % Construct a WMSLayer array with the same size as S. Copy the
        % matching fields of S to the properties of the class. NOCHECK is a
        % logical that if true, the fields of S are not validated.
        %
        % FOR INTERNAL USE ONLY -- This method is intentionally hidden and
        % is intended for use only within other toolbox classes and
        % functions. Its behavior may change, or the method itself may be
        % removed in a future release.
            
            % Validate S as a structure.
            validateattributes(S, {'struct'}, {'nonempty'}, 'WMSLayer');

            % Create an array based on the size and shape of S.
            self = createArray(self, num2cell(size(S)));

            % Copy the fields of S to the properties of the object.
            self = copyFieldsToProperties(self, S, nocheck);
        end

    end
end

%--------------- Constructor Helper Functions  ----------------------------

function S = parseInputs(inputs)
% Parse the parameter/value pair inputs and create a structure S whos
% fields match the parameter names.

internal.map.checkNameValuePairs(inputs{:});
try
    S = struct(inputs{:});
catch e
    % throw the exception in order to have the exception coming from
    % Error using WMSLayer->parseInputs rather than
    % Error using struct.
    throw(e);
end
end

%--------------------------------------------------------------------------

function layers = createArray(layers, szLayers)
% Construct a layers array of size szLayers. szLayers is a cell array of
% numeric values, with each element specifying the dimension size.

layers(szLayers{:}) = layers;
end

%--------------------------------------------------------------------------

function layers = copyFieldsToProperties(layers, S, nocheck)
% Copy the fields of S to the object, LAYERS. NOCHECK is a logical that if
% true, the fields of S are not validated.

% Obtain the names that are common to both property names
% and field names.
structNames = fieldnames(S);
propNames = properties(layers);
index = ismember(structNames, propNames);
names = structNames(index);

% Remove fields that are not property names.
S = rmfield(S,setdiff(structNames, names));

% Convert the structure values to a cell array.
cellS = struct2cell(S);

% If nocheck has bee specified, then copy the 'Latlim' and 'Lonlim' values
% to the private properties to avoid the overhead of validation.
if nocheck
    propNames = {'Latlim', 'Lonlim'};
    for k = 1:numel(propNames)
        index = strncmp(propNames{k}, names, numel(propNames{k}));
        if any(index)
            names{index} = ['p' propNames{k}];
        end
    end
end

% Copy only the fieldnames of S that match the property names and validate
% their class type.
for k=1:numel(names)
    values = cellS(k,:);
    if ~nocheck
        validateFieldClassType(names{k}, values );
    end
    [layers.(names{k})] = deal(values{:});
end
end

%--------------------------------------------------------------------------

function details = createDetails()
% Create a DETAILS structure.

details = struct( ...
    'MetadataURL', '', ...
    'Attributes', struct( ...
       'Queryable', false, ...
       'Cascaded', 0, ...
       'Opaque', false, ...
       'NoSubsets', false, ...
       'FixedWidth', 0, ...
       'FixedHeight', 0), ...
    'BoundingBox', struct( ...
       'CoordRefSysCode', '', ...
       'XLim', [], ...
       'YLim', []), ...
    'Dimension', struct( ...
       'Name','', ...
       'Units', '', ...
       'UnitSymbol','', ...
       'Default','', ...
       'MultipleValues', false, ...
       'NearestValue', false, ...
       'Current',false, ...
       'Extent', ''), ...
    'ImageFormats', {{}}, ...
    'ScaleLimits', struct( ...
       'ScaleHint', [], ...
       'MinScaleDenominator', [], ...
       'MaxScaleDenominator', []), ...
    'Style', struct( ...
       'Title', '', ...
       'Name', '', ...
       'Abstract', '', ...
       'LegendURL', struct( ...
           'OnlineResource', '', ...
           'Format','', ...
           'Height',[], ...
           'Width', [])), ...
    'Version', '');
end

%--------------- Validation Functions for Required Parameters -------------

function validateFieldClassType(fieldName, values)
% Validate the class type of the values in the cell array, VALUES,
% associated with the field, fieldName, of a structure. fieldName must
% match a property name of WMSLayer.

% Test based on fieldName
switch lower(fieldName)
    case {'latlim','lonlim'}
        % Validate the 'Latlim' or 'Lonlim' input.
        fcn = @(x) validateLimits(fieldName, x);
        cellfun(fcn, values);
        
    case 'coordrefsyscodes'
        % CoordRefSysCodes property must be cell array of strings.
        validateCellArrayClassType('cell');
        values = values{:};
        validateCellArrayClassType('char');
        
    case 'details'
        % Details property must be structure.
        fcn = @(x)validateattributes(x, ...
            {'struct'},{'nonempty'},'WMSLayer','Details');
        cellfun(fcn, values);
        
    otherwise
        % All the values in the cell array must be strings.
        validateCellArrayClassType('char')
        validateStringCellArray(fieldName, values)
end
    %----------------------------------------------------------------------
    
    function validateCellArrayClassType(className)
        index = cellfun('isclass', values, className);
        if ~all(index)
            badIndex = find(~index, 1);
            invalidClass = class(values{badIndex});
            error(message('map:validate:expectedUniformClassType', ...
                fieldName, className, num2ordinal(badIndex(1)), invalidClass));
        end
    end
end

%--------------------------------------------------------------------------

function  validateStringCellArray(name, values)
% Validate values to be a cell array of strings. The strings must be row
% vectors.

cIsCellArrayOfStrings = all(cellfun(@rowChar, values));
assert(cIsCellArrayOfStrings, ...
   message('map:validate:expectedCellArrayOfStrings', name));

    function tf = rowChar(c)
        tf = isempty(c) || (ischar(c) && ndims(c) == 2 && size(c,1) == 1);
    end
end
   
%--------------------------------------------------------------------------

function queryStr = validateQueryStr(queryStr)
% Validate the QUERYSTR parameter of the REFINE method.

if ~isempty(queryStr)
    validateattributes(queryStr, {'char'}, {'vector', 'row'}, ...
       'refine', 'QUERYSTR');
else
    queryStr = '';
end
end

%--------------------------------------------------------------------------

function validateLimits(name, limits)
% Validate latitude or longitude limits.  LIMITS must be a two-element
% real, finite, and numeric vector. NAME is a string with value 'Latlim' or
% 'Lonlim'.

% Validate numeric type.
validateGeoLimits(name, limits)

% All limits must be ascending.
assert(limits(1) <= limits(2), message('map:validate:expectedAscendingOrder', name));

% Validate range.
if isequal(name, 'Latlim')
    limitRange = [-90, 90];
    assert(all(limitRange(1) <= limits(:)) ...
        && all(limits(:) <= limitRange(2)), ...
        message('map:validate:expectedInterval', ...
        name, num2str(limitRange(1)), num2str(limitRange(2))));
else
    wrappedTo180 = all(-180 <= limits(:)) && all(limits(:) <= 180);
    wrappedTo360 = all(   0 <= limits(:)) && all(limits(:) <= 360);
    assert(~isempty(limits) && (wrappedTo180 || wrappedTo360), ...
        message('map:validate:expectedIntervals', ...
        name, '[-180, 180]', '[0 360]'));
end
end
 
%--------------------------------------------------------------------------

function validateGeoLimits(name, limits)
% Validate the georgraphic limit.

isValidLimit = ( ...
    isnumeric(limits) && ...
    all(isfinite(limits(:))) && ...
    isreal(limits) && ...
    numel(limits) == 2 );

assert(isValidLimit, message('map:validate:expectedEmptyOrTwoElementVector', name));

end

%-------------------------Parsing Functions -------------------------------

function options = parseDispOptions(propertyNames, varargin)
% Parse the optional parameters from the command line input and return a
% structure with fields, 'Label', 'Index', and 'Properties'. The value of
% each field is either the default value of the parameter or the
% user-supplied value.

% Set up the list of valid optional parameter names.
parameterNames = {'Label', 'Index', 'Properties'};

% Setup default values for 'Properties'.
propertiesDefaultValues = [propertyNames; 'populated'];

% Set the default values for each optional parameter.
defaultValues =  {true, true, propertiesDefaultValues};
 
% Setup the list of validation functions for the parameters.
validateFcns = { ...
    @(x)validateOnOff('Label',x), ...
    @(x)validateOnOff('Index',x), ...
    @(x)validateParamCell(x, 'Properties', propertiesDefaultValues', ...
       'all', 'disp')};

% Parse the parameters and return a structure.
options = parseParameters( ...
    parameterNames, validateFcns, defaultValues, varargin{:});

% Determine if using the 'populated' options for Properties. If so, set the
% options.Populated field to true and updated the Properties to the
% properties of the class.
populated = 'populated';
popIndex = strncmp(populated, options.Properties, numel(populated));
if any(popIndex)    
    options.Properties(popIndex) = [];
    options.Populated = true;
    if isempty(options.Properties);
        options.Properties = propertyNames;
    end
else
    options.Populated = false;
end

% If varargin contained 'all', then set Populated to false.
for k=1:numel(varargin)
   if any(strncmp('all', varargin{k}, 3))
       options.Populated = false;
       break
   end
end

    %----------------------------------------------------------------------
    
    function tf = validateOnOff(paramName, paramValue)
    % Validate any parameter with name, PARAMNAME, that requires a value
    % of 'on' or 'off'. TF is set to true if PARAMVALUE is 'on'; otherwise
    % TF is set to false.
    
        fcnName = 'disp';
        validateattributes( ...
            paramValue, {'char'}, {'nonempty'}, fcnName, paramName);
        paramValue = validatestring( ...
            paramValue, {'on','off'}, fcnName, paramName);
        tf = isequal(paramValue,'on');
    end
end

%--------------------------------------------------------------------------

function options = parseRefineOptions(varargin)
% Parse the optional parameters from the command line input for the refine
% method and return a structure with fields:
%   'SearchFields', 'MatchType', 'IgnoreCase'
% The value of each field is either the default value of the parameter or
% the user-supplied value.

% Setup the list of valid optional parameter names.
parameterNames = {'SearchFields', 'MatchType', 'IgnoreCase'};

% Set the default values for each optional parameter.
defaultValues =  {{'LayerName','LayerTitle'}, 'partial', true};
 
% Setup the list of validation functions for the parameters.
wmsFcns = assignWmsValidationFcns('refine');
validateFcns = { ...
    wmsFcns.validateSearchFieldsWithAbstract, ...
    wmsFcns.validateMatchType, ...
    wmsFcns.validateIgnoreCase};
 
% Parse varargin for all optional parameters.
options = parseParameters( ...
    parameterNames, validateFcns, defaultValues,  varargin{:});
end

%--------------------------------------------------------------------------

function options = parseRefineLimitsOptions(varargin)
% Parse the optional parameters from the command line input and return a
% structure with fields 'Latlim' and 'Lonlim'. The value of
% each field is either the default value of the parameter or the
% user-supplied value.

% Setup the list of valid optional parameter names.
parameterNames = {'Latlim', 'Lonlim'};

% Set the default values for each optional parameter.
defaultValues =  {[], []};
 
% Setup the list of validation functions for the parameters.
validateFcns = { ...
    @(x)validateLatlim(x), ...
    @(x)validateLonlim(x)};
 
% Parse varargin for all optional parameters.
options = parseParameters( ...
    parameterNames, validateFcns, defaultValues, varargin{:});

    %----------------------------------------------------------------------
    
    function limits = updateLimits(limits)
    % Ensure the limits variable contains two elements.
        
        if numel(limits) < 2
            limits(2) = limits(1);
        end             
    end

    %----------------------------------------------------------------------
    
    function latlim = validateLatlim(latlim)
    % Validate latitude limits.
    
        if ~isempty(latlim)
            latlim = updateLimits(latlim);
            validateLimits('Latlim', latlim);
        end
    end

    %----------------------------------------------------------------------
    
    function lonlim = validateLonlim(lonlim)
    % Validate longitude limits.
    
        if ~isempty(lonlim)
            lonlim = updateLimits(lonlim); 
            validateGeoLimits('Lonlim', lonlim);
        end
    end
end

%--------------------------------------------------------------------------

function options = parseParameters( ...
    parameterNames, validateFcns, defaultValues, varargin)
% Parse parameters and their values from VARARGIN and return the values in
% the scalar structure OPTIONS. OPTIONS contains fieldnames specified by
% parameterNames, a cell array of strings. If a parameter in parameterNames
% is matched in VARARGIN, the value of the parameter is passed to a
% validation function handle, found in the corresponding element in the
% cell array VALIDATEFCNS. After passing validation, the output is assigned
% to the fieldname in OPTIONS. If the parameter is not matched in VARARGIN,
% the value of the parameter is assigned the default value from the
% corresponding element in DEFAULTVALUES. If additional parameters in
% VARARGIN are found that do not match a name in parameterNames, an error
% is issued.  FCNNAME is the name of the calling function and is used in
% constructing error messages.
    
% Parse the inputs.
[options, userSupplied, unmatched] = ...
    internal.map.parsepv(parameterNames, validateFcns, varargin{:});

% Verify the first parameter is a string.
if ~isempty(varargin) && ~ischar(varargin{1})
    unmatched{1} = 'PARAM1';
end

% Check if varargin contained unmatched parameters.
if ~isempty(unmatched)
    parameterNames = sprintf('''%s'', ', parameterNames{:});
    error(message('map:validate:invalidParameterName', ...
        unmatched{1}, parameterNames(1:end-2)));
end

% Set default values for any parameter that is not specified.
options = setDefaultValues(options, defaultValues, userSupplied);

    %----------------------------------------------------------------------
    
    function inputs = setDefaultValues(inputs, defaultValues, userSupplied)
    % Set default values to any parameter that is not supplied.
        
        inputFieldNames = fieldnames(inputs);
        for k=1:numel(inputFieldNames)
            name = inputFieldNames{k};
            if ~userSupplied.(name)
                inputs.(name) = defaultValues{k};
            end
        end
    end
end

%-------------------- disp Support Functions ------------------------------

function printLayers(layers, options)
% Print the properties of the object, LAYERS, to the command window,
% modified by the structure OPTIONS.

% Create a string denoting the size of the object.
if numel(layers) == 1
    szString = ' ';
else
    sz = size(layers);
    szString = sprintf('  %d', sz(1));
    for k=2:numel(sz)
        szString = sprintf('%sx%d',szString, sz(k));
    end
end

% If the desktop is in use, then add a clickable hyperlink to the class
% type.
classType = class(layers);
usingDesktop = usejava('desktop') && desktop('-inuse') ...
    && feature('hotlinks');
if usingDesktop
    dispClassType = ...
       ['<a href="matlab:help ' classType '">' classType '</a>'];
else
    dispClassType = classType;
end

if ~isempty(layers)
    % Print out the values of the properties, modified by the options 
    % structure.
    fprintf('%s %s\n', szString, dispClassType);
    fprintf('\n%s','  Properties:');
    
    for k = 1:numel(layers)
        fprintf('\n');
        if options.Index
            printTextLabel('Index', num2str(k), true);
        end
        printLayer(layers(k), options);
    end
else
    % When the object is empty, print only the list of properties.
    fprintf('%s empty %s\n', szString, dispClassType);
    fprintf('\n%s\n','  Properties:');
    props = properties(layers);
    for k=1:numel(props)
        fprintf('    %s\n', props{k});
    end
end

% If the desktop is in use, add the 'Methods' clickable link.
if usingDesktop && options.Label && options.Index ...
    && isequal(properties(layers), options.Properties)
   fprintf('\n  %s\n\n', ...
       ['<a href="matlab:methods(''' classType ''')">Methods</a>']);
else
    fprintf('\n');
end
end

%--------------------------------------------------------------------------

function printLayer(layer, options)
% Print the properties of the scalar object, LAYER, to the command window,
% modified by the structure, OPTIONS.

for k=1:numel(options.Properties)
    label = options.Properties{k};
    value = layer.(label);
    if isequal(label,'Latlim') || isequal(label,'Lonlim')
        printCoordinate(label, value, options.Label);
    else        
        if ~(options.Populated && isequal(value, WMSLayer.UpdateString))
            value = updateDispValue(value);
           printTextLabel(label, value, options.Label);
        end
    end
end

    %----------------------------------------------------------------------
    
    function value = updateDispValue(value)
    % Update the display value based on class type and size.
    
    if iscell(value)
        numelInCellCutoff = 3;
        if numel(value) <= numelInCellCutoff
            % For readability, display up to numelInCellCutoff values.
            value = sprintf('''%s'' ',value{:});
            if ~isempty(value)
                value(end) = [];
            end
            value = ['{' value(:)' '}'];
        else
            % Number of values is greater than numelInCellCutoff;
            % display the size instead.
            sz = size(value);
            value = ['{' num2str(sz(1)) 'x' num2str(sz(2)) ' cell}'];
        end
    elseif isstruct(value)
        % For a structure value, display size only.
        sz = size(value);
        value = ['[' num2str(sz(1)) 'x' num2str(sz(2)) ' struct]'];
        
    else
        value = addQuoteMarks(value);
    end
    end
    
    %----------------------------------------------------------------------
    
    function label = addQuoteMarks(label)
    % Add a quote mark before and after the string, LABEL.
    
        singleQuote = '''';
        label = [singleQuote label singleQuote];
    end

    %----------------------------------------------------------------------
    
    function printCoordinate(label, value, printLabelFlag)
    % Print the coordinate value. LABEL and VALUE are strings. VALUE must
    % be a two-element vector specifying the coordinate. printLabelFlag is
    % a logical. If true, LABEL is printed; otherwise, LABEL is excluded.
    
        if ~isempty(value)
            strValue = sprintf('[%4.4f %4.4f]', value(1), value(2));
        else
            strValue = '[]';
        end
        printTextLabel(label, strValue, printLabelFlag);
    end
end

%--------------------------------------------------------------------------

function printTextLabel(label, value, printLabelFlag)
% Print the label and value to the command window if the logical
% printLabelFlag is true; otherwise print only the value.

% numSpaces is the number of elements preceding the character ':' when
% displaying a structure at the command line. The disp printout follows
% this same layout. 
numSpaces = 16;
indentSpaces(1:numSpaces) = ' ';
numIndentSpaces = numSpaces - numel(label);

% Insert the correct number of indentation spaces in front of label and add
% the ':' character.
label = [indentSpaces(1:numIndentSpaces) label ':'];

% Print the label and value, if printLabelFlag is true; otherwise, print
% only the value without indentation.
if printLabelFlag
    fprintf('%s %s\n', label, value);
else
    fprintf('%s\n', value);
end
end

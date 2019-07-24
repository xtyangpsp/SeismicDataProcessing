%WMSCapabilites Web Map Service capabilities object
%
%   A WMSCapabilities object represents a Web Map Service (WMS)
%   capabilities document obtained from a WMS server.
%
%   capabilities = WMSCapabilites(ServerURL, capabilitiesResponse)
%   constructs a WMSCapabilities object from the input string parameters.
%   The ServerURL string is a WMS server URL and must include the protocol
%   'http://' or 'https://'. The capabilitiesResponse is a string
%   containing the XML elements that describe the capabilities of the
%   ServerURL WMS server.
%
%   WMSCapabilities properties:
%      ServerTitle - Title of server
%      ServerURL - URL of server
%      ServiceName - Web Map Service name
%      Version - Web Map Service version
%      Abstract - Information about server
%      OnlineResource - URL for information about server
%      ContactInformation - Organizational contact information
%      AccessConstraints - Constraints to access the server
%      Fees - Fees associated with accessing server
%      KeywordList - Descriptive keywords of server
%      ImageFormats - Image formats supported by server
%      LayerNames - List of layer names provided by server
%      Layer - Layer information
%      AccessDate - Date of request to server
%
%   WMSCapabilities method:
%      disp - Display properties
%
%   Example
%   -------
%   % Construct a WMSCapabilities object from the contents of a downloaded
%   % capabilities file from the NASA SVS Image Server.
%   nasa = wmsfind('NASA SVS Image', 'SearchField', 'servertitle');
%   serverURL = nasa(1).ServerURL;
%   server = WebMapServer(serverURL);
%   capabilities = server.getCapabilities;
%   filename = fullfile(tempdir, 'capabilities.xml');
%   urlwrite(server.RequestURL, filename);
%   
%   fid = fopen(filename, 'r');
%   capabilitiesResponse = fread(fid, 'uint8=>char');
%   fclose(fid);
%   delete(filename);
%   capabilities = WMSCapabilities(serverURL, capabilitiesResponse)
%
%   See also WebMapServer, WMSINFO, WMSLayer.

% Copyright 2008-2011 The MathWorks, Inc.
% $Revision: 1.1.6.12 $  $Date: 2011/05/17 01:54:03 $

classdef WMSCapabilities
    
    properties (SetAccess = private, GetAccess = public)
        
        %ServerTitle Title of server
        %
        %   ServerTitle is a string indicating the title of the server.
        ServerTitle = '';
        
        %ServerURL URL of server
        %
        %   ServerURL is a string indicating the URL of the server.
        ServerURL = '';
        
        %ServiceName Web Map Service name
        %
        %   ServiceName is a string indicating the name of the Web map
        %   service.
        ServiceName = '';
        
        %Version Web Map Service version
        %
        %   Version is a string indicating the WMS version specification.
        Version = '';
        
        %Abstract Information about server
        %
        %   Abstract is a string containing information about the server.
        Abstract = '';
        
        %OnlineResource URL for information about server
        %
        %   OnlineResource is a URL that points to online information about
        %   the server.
        OnlineResource = '';
        
        %ContactInformation Organizational contact information
        %
        %   ContactInformation is a structure containing an individual or
        %   organization to contact, including an email address if
        %   provided.
        ContactInformation = struct( ...
            'Person', '', ...
            'Organization', '', ...
            'Email', '');
        
        %AccessConstraints Constraints to access the server
        %
        %   AccessConstraints is a string indicating the constraints
        %   inherent in accessing the server, such as server load limits.
        AccessConstraints = '';
        
        %Fees Fees associated with accessing server
        %
        %   Fees is a string indicating the types of fees associcated with
        %   accessing the server.
        Fees = '';
        
        %KeywordList Descriptive keywords of server
        %
        %   KeywordList is a cell array of strings indicating the
        %   descriptive keywords of the server.
        KeywordList = {};
        
        %ImageFormats Image formats supported by server
        %
        %   ImageFormats is a cell aray of strings indicating the image
        %   formats supported by the server.
        ImageFormats = {};
        
        %LayerNames List of layer names provided by server
        %
        %   LayerNames is a cell array of strings indicating the layer
        %   names provided by the server.
        LayerNames = '';
        
        %Layer Layer information
        %
        %   Layer is a <a href="matlab:help WMSLayer">WMSLayer</a> array
        %   containing information about the server's layers. 
        Layer = WMSLayer();
        
        %AccessDate Date of request to server
        %
        %   AccessDate is a string indicating the date of the request to
        %   the server.
        AccessDate = '';
    end
    
    methods
        
        function self = WMSCapabilities(serverURL, capabilitiesResponse)
            % WMSCapabilities Constuctor.
            
            if nargin > 0
                error(nargchk(2, 2, nargin, 'struct'));
                
                % Verify that Java is available.
                error(javachk('jvm','WMSCapabilities'));

                % Validate the serverURL input.
                validateServerURL(serverURL);
                
                % Validate the capabilitiesResponse input.
                validateCapabilitiesResponse( ...
                    capabilitiesResponse, serverURL);
                
                % Parse the capabilitiesResponse and return an object that
                % contains accessors to the elements of the WMS document.
                document = parseCapabilitiesResponse(capabilitiesResponse);
                
                % Set the properties of this class.
                self = setProperties(self, document, serverURL);
            end
        end

        function disp(self) 
        %disp Display properties
        %
        %   capabilities.disp() is overloaded to remove the hyperlinks and
        %   expand string and cell array properties.  

            % Display the properties to the command window.
            if numel(self) > 1
                builtin('disp', self);
            else
               displayProperties(self);    
            end
        end
    end
end

%----------- Validation Functions for Constructor Parameters --------------

function validateServerURL(serverURL)
% Validate serverURL as a string.

validateattributes(serverURL, {'char'}, {'nonempty'}, ...
    'WMSCapabilities', 'ServerURL')
end

%--------------------------------------------------------------------------

function validateCapabilitiesResponse(capabilitiesResponse, serverURL)
% Validate capabilitiesResponse as a string that contains XML elements.

isObject = isa(capabilitiesResponse, 'org.geotools.data.ows.WMSCapabilities');
assert(ischar(capabilitiesResponse) || isObject, ...
    message('map:wms:expectedStringDocument', serverURL));

if (~isObject) 
    % Reshape the capabilitiesResponse for string comparisons.
    capabilitiesResponse = ...
        reshape(capabilitiesResponse,[1,numel(capabilitiesResponse)]);

    % Verify capabilitiesResponse contains valid XML. This simple test
    % verifies that  
    %    1) the string contains < and > brackets
    %    2) the string contains the XML keyword, 'xml'
    %    3) the string is not an HTML page (containing the element '<html')
    % Additional verification and error handling is performed by the XML
    % parser.
    leftBracket  = strfind(capabilitiesResponse, '<');
    rightBracket = strfind(capabilitiesResponse, '>');
    xmlIndex     = regexpi(capabilitiesResponse, 'xml');
    htmlIndex    = regexpi(capabilitiesResponse, '<html');
    isXML = ~isempty(leftBracket) && ~isempty(rightBracket) ...
        && ~isempty(xmlIndex) && isempty(htmlIndex);
    assert(isXML, message('map:wms:expectedXML', serverURL));    
end
end

%------------------------- Parse Functions  -------------------------------

function document = parseCapabilitiesResponse(capabilitiesResponse)
% Parse the capabilitiesResponse and return a Java WMSCapabilities object
% that contains accessors to the content from the elements of the WMS
% document.

if ischar(capabilitiesResponse)
    % Create a WMSParser object to parse the string.
    parser = com.mathworks.toolbox.geoweb.wms.WMSParser();
    
    % Parse the capabilities string.
    try
        parser.parse(capabilitiesResponse);
    catch e
        % An error occurred while parsing the XML content. Rather than
        % displaying all the Java stack trace in the error message, display a
        % generic parse error message.
        if isequal(e.identifier,'MATLAB:Java:GenericException')
            error(message('map:wms:capabilitiesParseError'));
        else
            rethrow(e)
        end
    end
    
    % Retrieve the Java WMSDocument object from the parser.
    document = parser.getDocument;
    
    % Assert the object is not empty.
    assert(~isempty(document), message('map:wms:documentIsEmpty'));
    
    % Retrieve the Java WMSCapabilities document.
    document = document.getCapabilities();
else   
    % The response is a Java WMSCapabilities object.
    document = capabilitiesResponse;
end

% Assert the document is not empty.
assert(~isempty(document), message('map:wms:documentIsEmpty'));

end

%--------------------- Property Set Functions  ----------------------------

function capabilities = setProperties(capabilities, document, serverURL)
% Set the properties of the capabilities class.

% Set the AccessDate property.
capabilities.AccessDate = datestr(floor(now));

% Set the ServerURL property.
capabilities.ServerURL = serverURL;

% Set the Version property. The Version number is an attribute of the
% top-level WMT_MS_Capabilities keyword.
capabilities.Version = char(document.getVersion);

% The XML capabilities document from the server contains two top-level
% keywords, 'Service' and 'Capability'. Set the Service and Capability
% properties of this class from the values in the document.
capabilities = setServiceProperties(capabilities, document);
capabilities = setCapabilityProperties(capabilities, document);
end

%--------------------------------------------------------------------------

function capabilities = setServiceProperties(capabilities, document)
% Set the Service properties. The values are obtained from the document
% object.

% Get the Service object obtained from the WMT_MS_Capabilities/Service
% XML element.
service = document.getService;

if ~isempty(service)
    % Set the simple string properties from the Service object.
    capabilities.ServerTitle = char(service.getTitle);
    capabilities.ServiceName = char(service.getName);
    capabilities.OnlineResource = char(service.getOnlineResource);
    capabilities.Abstract = char(service.get_abstract);
    capabilities.Fees = char(service.getFees);
    capabilities.AccessConstraints = char(service.getAccessConstraints);
    
    % Set the KeywordList property.
    keywordList = service.getKeywordList;
    if ~isempty(keywordList)
        capabilities.KeywordList = cell(keywordList);
    else
        capabilities.KeywordList = {};
    end
    
    % Set the ContactInformation property.
    capabilities.ContactInformation = ...
        setContactInformation(capabilities.ContactInformation, service);
end
end

%--------------------------------------------------------------------------

function ContactInformation = ...
    setContactInformation(ContactInformation, service)
% Set the ContactInformation property. The values for the property are 
% obtained from the service object.

serviceInformation = service.getContactInformation;
if ~isempty(serviceInformation)
    ContactInformation.Person = ...
        char(get(serviceInformation,'IndividualName'));
    ContactInformation.Organization = ...
        char(get(serviceInformation,'OrganisationName'));
    info = serviceInformation.getContactInfo();
    if ~isempty(info)
        address = info.getAddress();
        if ~isempty(address)
           emailArray = address.getElectronicMailAddresses().toArray;
           ContactInformation.Email = char(emailArray);
        end
    end
end
end

%--------------------------------------------------------------------------

function capabilities = setCapabilityProperties(capabilities, document)
% Set the Capability properties. The values are obtained from the document
% object.

% Set the ImageFormat property obtained from the XML element:
%    Capability/Request/GetMap.
capabilities = setGetMapProperties(capabilities, document);

% Set the Layer properties obtained from the XML element:
%    Cabability/Layer/Layer.
capabilities = setLayerProperties(capabilities, document);
end

%--------------------------------------------------------------------------

function capabilities = setGetMapProperties(capabilities, document)
% Set the GetMap properties. The values are obtained from the document
% object.

% Set the ImageFormats property. The ImageFormats are found in the XML
% keyword: Capability/Request/GetMap/Format
request = document.getRequest;
if ~isempty(request)
    map = request.getGetMap;
    formatStrings = map.getFormatStrings;
    capabilities.ImageFormats = cell(formatStrings);
end
end

%--------------------------------------------------------------------------

function capabilities = setLayerProperties(capabilities, document)
% Set the Layer properties. The values are obtained from the document
% object.

%Retrieve the named layers from the document.
namedLayers = getNamedLayers(document);

if ~isempty(namedLayers)
    
    % Construct a Layer structure to hold the values from namedLayers.
    layerFields = properties(capabilities.Layer);
    layerValues = cell(numel(layerFields), 1);
    Layer = cell2struct(layerValues, layerFields);
    
    % Copy the properties to the structure Layer.
    for k=1:numel(layerFields)
        Layer.(layerFields{k}) = capabilities.Layer.(layerFields{k});
    end

    % Extend the size of the Layer structure to match the number of named
    % layers.
    Layer(1:numel(namedLayers),1) = deal(Layer(1));
    
    % Loop through each of the namedLayers and set the corresponding
    % fields.
    for k=1:numel(namedLayers)
        
        % Get the namedLayer element from the list.
        namedLayer = namedLayers(k);
        
        % Set ServerURL and ServerTitle properties.
        Layer(k).ServerURL   = capabilities.ServerURL;
        Layer(k).ServerTitle = capabilities.ServerTitle;
        
        % Set the simple string properties.
        Layer(k).LayerTitle  = char(namedLayer.getTitle);
        Layer(k).LayerName   = char(namedLayer.getName);
        Layer(k).Abstract    = char(namedLayer.get_abstract);
         
        % Set Latlim and Lonlim properties.
        llbbox = namedLayer.getLatLonBoundingBox;
        Layer(k).Latlim = [get(llbbox,'MinY'), get(llbbox,'MaxY')];
        Layer(k).Lonlim = [get(llbbox,'MinX'), get(llbbox,'MaxX')];
        
        % Set the fields of Details.
        Layer(k) = setDetails(capabilities, Layer(k), namedLayer); 
        
        % Set the cell array CoordRefSysCodes property.
        Layer(k) = setCRS(Layer(k), namedLayer);       

        % Check for empty limits and adjust if appropriate.
        Layer(k) = checkEmptyLimits(Layer(k));
    end
    
    % Remove redundant layers.
    Layer = removeRedundantLayers(Layer);
    
    % Remove Layers that contain empty CoordRefSysCodes.
    Layer = removeEmptyCRSLayers(Layer);
    
    % Set the Layer property. Convert the Layer structure array to a 
    % WMSLayer array.
    l = WMSLayer;
    capabilities.Layer = l.setPropertiesFromStruct(Layer, true);
    
    % Create a cell array of all LayerNames.
    capabilities.LayerNames = {capabilities.Layer.LayerName}';
end
end

%--------------------------------------------------------------------------

function namedLayers = getNamedLayers(document)
% Retrieve the named layers from the document.

try
    namedLayers = org.geotools.data.wms.WMSUtils.getNamedLayers(document);
catch e
    if isequal(e.identifier, 'MATLAB:Java:GenericException')
        error(message('map:wms:getNamedLayersException'));
    else
        rethrow(e);
    end
end
end

%--------------------------------------------------------------------------

function Layer = setCRS(Layer, namedLayer)
% Set the CoordRefSysCodes property. coordRefSysCodes is a sorted cell
% array of strings containing the list of coordinate reference system
% codes.

% Obtain the coordinate referencing system codes, returned as a
% java.util.TreeSet object.
crs = namedLayer.getSrs();

% Create a cell array of strings.
coordRefSysCodes = cell(crs.toArray);

% Assign a variable for the BoundingBox CoordRefSysCodes.
bboxCRS = {Layer.Details.BoundingBox.CoordRefSysCode};

% Find the logical index of all elements in coordRefSysCodes that are
% contained in the bboxCRS cell array.
codesInBoundingBox = findCodesInBoundingBox(bboxCRS, coordRefSysCodes);

% Assign a variable to represent the linear indices of all the
% coordRefSysCodes that are not found in the bounding box.
notFoundIndices = find(~codesInBoundingBox);

% Attempt to match the codes not found in the bounding box with the
% standard EPSG:4326 and CRS:84 codes. If these codes are not in the
% BoundingBox, then insert them.
stringCodes = coordRefSysCodes(~codesInBoundingBox);
codesToCopy = {'EPSG:4326', 'CRS:84'};
index = ismember(codesToCopy, stringCodes);
epsgCodes = codesToCopy(index);
BoundingBox = Layer.Details.BoundingBox(1);
for k=1:numel(epsgCodes)
    BoundingBox.CoordRefSysCode = epsgCodes{k};
    BoundingBox.XLim = Layer.Lonlim;
    BoundingBox.YLim = Layer.Latlim;
    if isempty(Layer.Details.BoundingBox(1).CoordRefSysCode)
        % The layer was constructed without any bounding box
        % information. Add the new bounding box to the first element.
        Layer.Details.BoundingBox(1) = BoundingBox;
    else
        % The layer contains bounding boxes. Add this one to the end.
        Layer.Details.BoundingBox(end+1) = BoundingBox;
    end
    codesInBoundingBox(notFoundIndices(k)) = true;
end

% Remove coordRefSysCodes not found in the BoundingBox structure and sort
% the results to make them unique.
coordRefSysCodes = coordRefSysCodes(codesInBoundingBox);
index = cellfun(@isempty, coordRefSysCodes);
coordRefSysCodes(index) = [];
coordRefSysCodes = sort(coordRefSysCodes);

% Set the codes.
Layer.CoordRefSysCodes = coordRefSysCodes;     
end

%--------------------------------------------------------------------------

function Layer = setDetails(capabilities, Layer, namedLayer)
% Set the fields of Details.

% Set the Details MetadataURL field.
Layer.Details.MetadataURL = char(namedLayer.getMetadataURL);

% Set the Details Attribute field.
Layer.Details.Attributes = ...
    setAttributes(Layer.Details.Attributes, namedLayer);

% Set the BoundingBox field.
Layer.Details.BoundingBox = ...
    setBoundingBox(Layer.Details.BoundingBox, namedLayer);

% Set the Details Dimension field if it is not empty.
Layer.Details.Dimension = ...
    setDimension(Layer.Details.Dimension, namedLayer);

% Copy ImageFormats to Details.
Layer.Details.ImageFormats = capabilities.ImageFormats;

% Set the Details ScaleHint field.
Layer.Details.ScaleLimits = setScaleLimits(namedLayer);

% Set the Details Style field.
Layer.Details.Style = setStyle(Layer.Details.Style, namedLayer);

% Copy Version to Details.
Layer.Details.Version = capabilities.Version;
end
        
%--------------------------------------------------------------------------

function Attributes = setAttributes(Attributes, namedLayer)
% Set the Attributes fields.

% Set the logical fields of Attributes.    
Attributes.Queryable  = namedLayer.isQueryable;
Attributes.Opaque     = namedLayer.isOpaque;
Attributes.NoSubsets  = namedLayer.isNoSubsets;
Attributes.Cascaded   = namedLayer.Cascaded;

% Set the FixedWidth and FixedHeight fields.
Attributes.FixedWidth  = namedLayer.FixedWidth;
Attributes.FixedHeight = namedLayer.FixedHeight;
end

%--------------------------------------------------------------------------

function BoundingBox = setBoundingBox(BoundingBox, namedLayer)
% Set BoundingBox field.

% Obtain the bounding box values.
bbox = namedLayer.getBoundingBoxes.values.toArray;

if ~isempty(bbox)
    BoundingBox(numel(bbox)) = BoundingBox;
    for k=1:numel(bbox)
        BoundingBox(k).CoordRefSysCode = char(bbox(k).getEPSGCode);
        BoundingBox(k).XLim = [bbox(k).getMinX, bbox(k).getMaxX];
        BoundingBox(k).YLim = [bbox(k).getMinY, bbox(k).getMaxY];
    end
end
end

%--------------------------------------------------------------------------

function Layer = checkEmptyLimits(Layer)
% A server may not provide a geographic bounding box in which case the
% layer's limits will be empty. If so, check the BoundingBox for
% either an EPSG:4326 or CRS:84 coordinate referenece system. If found,
% then set the limits to the bounding box.

if isempty(Layer.Latlim) && isempty(Layer.Lonlim)
    codes = {Layer.Details.BoundingBox.CoordRefSysCode};
    % Choose CRS:84 first since if it is provided, the X limits will always
    % correspond to longitude and the Y limits always correspond to
    % latitude.
    index = find(strcmp(codes, {'CRS:84'}));
    if ~isempty(index)
        index = index(1);
        Layer.Latlim = Layer.Details.BoundingBox(index).YLim;
        Layer.Lonlim = Layer.Details.BoundingBox(index).XLim;
    else
        % If CRS:84 is not found, then check EPSG:4326. In some cases, for
        % version 1.3.0, the X limits will correspond to latitude, not
        % longitude. However, there is no way to predict this.
        index = find(strcmp(codes, {'EPSG:4326'}));
        if ~isempty(index)
            index = index(1);
            Layer.Latlim = Layer.Details.BoundingBox(index).YLim;
            Layer.Lonlim = Layer.Details.BoundingBox(index).XLim;
        end
    end
end
end

%--------------------------------------------------------------------------

function ScaleLimits = setScaleLimits(namedLayer)
% Set the ScaleLimits fields.

scaleHintMin = namedLayer.getScaleHintMin;
scaleHintMax = namedLayer.getScaleHintMax;
scaleHint = [scaleHintMin, scaleHintMax];
if all(isnan(scaleHint))
    scaleHint = [];
end
ScaleLimits.ScaleHint = scaleHint;

% Set the MinScaleDenominator and MaxScaleDenominator properties.
minScale = namedLayer.getMinScaleDenominator;
maxScale = namedLayer.getMaxScaleDenominator;
if isnan(minScale)
    minScale = [];
end
if isnan(maxScale)
    maxScale = [];
end
ScaleLimits.MinScaleDenominator = minScale;
ScaleLimits.MaxScaleDenominator = maxScale;
end
        
%--------------------------------------------------------------------------

function Style = setStyle(Style, namedLayer)
% Set the Style fields. The values are obtained from the namedLayer object.

styles = namedLayer.getStyles;
styleArray = styles.toArray;

if ~isempty(styleArray)
    Style(1:numel(styleArray)) = deal(Style(1));
    for k=1:numel(styleArray)
        Style(k).Title = char(styleArray(k).getTitle);
        Style(k).Name  = char(styleArray(k).getName);
        Style(k).Abstract  = char(styleArray(k).getAbstract);

        legendURL = styleArray(k).getLegendURL;
        Style(k).LegendURL.OnlineResource = char(legendURL.getOnlineResource);
        Style(k).LegendURL.Format = char(legendURL.getFormat);
        Style(k).LegendURL.Height = double(legendURL.getHeight);
        Style(k).LegendURL.Width = double(legendURL.getWidth);
    end
end
end

%--------------------------------------------------------------------------

function Dimension = setDimension(Dimension, namedLayer)
% Set the Dimension fields. The values are obtained from the namedLayer
% object.

dimension = namedLayer.getDimension();
dimensionArray = dimension.toArray;
if ~isempty(dimensionArray)
    Dimension(1:numel(dimensionArray)) = deal(Dimension(1));
    for k=1:numel(dimensionArray)
        Dimension(k).Name = char(dimensionArray(k).getName());
        Dimension(k).Units = char(dimensionArray(k).getUnits());
        Dimension(k).UnitSymbol = char(dimensionArray(k).getUnitSymbol());
        Dimension(k).Default = char(dimensionArray(k).getDefault());
        Dimension(k).MultipleValues = dimensionArray(k).getMultipleValues();
        Dimension(k).NearestValue = dimensionArray(k).getNearestValue();
        Dimension(k).Current = dimensionArray(k).getCurrent();
        Dimension(k).Extent = char(dimensionArray(k).getExtent());
    end    
end
end

%--------------------------------------------------------------------------

function Layer = removeRedundantLayers(Layer)
% Remove redundant layers from the Layers structure and sort the layers
% alphabetically based on LayerName.

layerNames = {Layer.LayerName};
[~, index] = unique(layerNames);
Layer = Layer(index);
end

%--------------------------------------------------------------------------

function Layer = removeEmptyCRSLayers(Layer)
% Remove layers from the Layers structure which do have empty 
% CoordRefSysCodes.

if ~isempty(Layer)
    crs = {Layer.CoordRefSysCodes};
    index = cellfun('isempty', crs);
    if all(index)
        error(message('map:wms:removedEmptyCRSLayers', Layer(1).ServerURL));
    elseif ~isempty(crs(index))
        layerNames = {Layer(index).LayerName};
        codes = sprintf('''%s'' ', layerNames{:});
        warning(message('map:wms:removedLayersWithoutCRSCodes', ...
            Layer(1).ServerURL, codes))
    end
    Layer = Layer(~index);
end
end

%--------------------------------------------------------------------------

function index = findCodesInBoundingBox(bboxCRS, coordRefSysCodes)
% Find the indices of all elements in coordRefSysCodes that are contained
% in the bboxCRS cell array or that define either an automatic projection
% or the EPSG:4326 projection.

% Keep the automatic projections (AUTO) and the default EPSG:4326
% projection codes.
keepCodes = {'AUTO', 'EPSG:4326'};
keepCodesIndex = false(size(coordRefSysCodes));
for k=1:numel(keepCodes)
    code = keepCodes{k};
    keepCodesIndex = keepCodesIndex | ...
        strncmpi(code,  coordRefSysCodes, numel(code));
end
    
% Find all elements that match criteria and assign the return argument.
index = keepCodesIndex | ismember(coordRefSysCodes, bboxCRS);
end

%---------------------- Display Functions ---------------------------------

function displayProperties(capabilities)
% Display the properties of the object, CAPABILITIES, to the command
% window.

% If the desktop is in use, then add a clickable hyperlink to the class
% type.
classType = class(capabilities);
usingDesktop = usejava('desktop') && desktop('-inuse') ...
    && feature('hotlinks');
if usingDesktop
    dispClassType = ...
       ['<a href="matlab:help ' classType '">' classType '</a>'];
else
    dispClassType = classType;
end

props = properties(capabilities);
if ~isempty(capabilities) 
    % Display the names and values of the properties. Add a clickable
    % hyperlink around the OnlineResource property if the desktop is in
    % use.
    fprintf('  %s\n',dispClassType);
    fprintf('\n  Properties:\n');
    httpProp = 'OnlineResource';
    hyperlinkIndex = find(strncmp(httpProp, props, numel(httpProp)));
    for k = 1:numel(props)
        value = capabilities.(props{k});
        if (k == hyperlinkIndex) && usingDesktop && ~isempty(value)
            value = ['<a href = "' value '">' value '</a>']; %#ok<AGROW>
        end
        displayProperty(props{k}, value);
    end
else
    % Create a string denoting the size of the object.
    sz = size(capabilities);
    szString = sprintf('%d', sz(1));
    for k=2:numel(sz)
        szString = sprintf('%sx%d',szString, sz(k));
    end
    
    % When the object is empty, display only the list of properties.
    empty = ' empty ';
    fprintf('  %s%s%s\n', szString, empty, dispClassType);
    fprintf('\n  Properties:\n');
    for k=1:numel(props)
        fprintf('    %s\n', props{k});
    end
end

% If the desktop is in use, add the 'Methods' clickable link.
if usingDesktop
   fprintf('\n  %s\n\n', ...
       ['<a href="matlab:methods(''' classType ''')">Methods</a>']);
else
    fprintf('\n');
end
end

%--------------------------------------------------------------------------

function displayProperty(label, value)
% Display the property and its value.

value = updateDispValue(value);
displayTextLabel(label, value);
end

%----------------------------------------------------------------------

function value = updateDispValue(value)
% Update the display value based on class type and size.

sz = size(value);
classType = class(value);
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
        value = ['{' num2str(sz(1)) 'x' num2str(sz(2)) ' ' classType '}'];
    end
elseif isstruct(value) || isa(value, 'WMSLayer')
    % For a structure or WMSLayer object, display size only.
    value = ['[' num2str(sz(1)) 'x' num2str(sz(2)) ' ' classType ']'];
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

%--------------------------------------------------------------------------

function displayTextLabel(label, value)
% Display the label and value to the command window.

% numSpaces is the number of elements preceding the character ':' when
% displaying a structure at the command line. The disp printout follows
% this same layout.
numSpaces = 20;
indentSpaces(1:numSpaces) = ' ';
numIndentSpaces = numSpaces - numel(label);

% Insert the correct number of indentation spaces in front of label and add
% the ':' character.
label = [indentSpaces(1:numIndentSpaces) label ':'];

% Print the label and value with indentation.
fprintf('%s %s\n', label, value);
end

function [capabilities, infoRequestURL] = wmsinfo(serverURL, varargin)
%WMSINFO Information about WMS server from capabilities document
%
%   WMSINFO accesses the Internet to read the capabilities document from a
%   Web Map Service (WMS) server. An Internet connection must be
%   established to use the function. Note that a WMS server may be
%   unavailable, and several minutes may elapse before the document is
%   returned.
%
%   WMSINFO communicates with the server using a WebMapServer handle object
%   representing a Web Map Service (WMS). The handle object acts as a proxy
%   to a WMS server and resides physically on the client side. The handle
%   object accesses the server's capabilities document. The handle object
%   supports multiple WMS versions and negotiates with the server to use
%   the highest known version that the server can support. The handle
%   object automatically times-out after 60 seconds if a connection is not
%   made to the server.
%
%   [CAPABILITIES, infoRequestURL] = WMSINFO(serverURL) reads the
%   capabilities document from a WMS server and returns the contents into
%   CAPABILITIES. serverURL is a WMS server URL and must include the
%   protocol 'http://' or 'https://' and may contain additional WMS or
%   access keywords. CAPABILITIES is a WMSCapabilities object and
%   represents the information contained in the server's capabilities
%   document. The properties of the WMSCapabilities object are listed in
%   the table below. infoRequestURL is a URL string, composed of the
%   ServerURL with additional WMS parameters. The infoRequestURL can be
%   inserted into a browser or URLREAD to return the XML capabilities
%   document.
%
%   [CAPABILITIES, infoRequestURL] = WMSINFO(infoRequestURL) reads the
%   capabilities document from a WMS infoRequestURL and returns the
%   contents into CAPABILITIES. 
%
%   [CAPABILITIES, infoRequestURL] = WMSINFO(..., PARAM1, VAL1)
%   specifies a parameter-value pair that modifies the request to the
%   server. The parameter name can be abbreviated and is case-insensitive.
%   See the table below for a description of the parameter.
% 
%   Optional Parameter
%   ------------------
%
%  'TimeoutInSeconds'  An integer-valued scalar double indicating the
%                      number of seconds to elapse before a server timeout
%                      is issued. By default, the value is 60 seconds. A
%                      value of 0 causes the timeout mechanism to be
%                      ignored.
%
%   Output WMSCapabilities Object
%   -----------------------------
%   The output WMSCapabilities object contains the following properties:
% 
%   Name               Data Type    Purpose
%   ----               ---------    -------  
%   ServerTitle        String       Descriptive title of server
%
%   ServerURL          String       URL of Web map server
%
%   ServiceName        String       Web Map Service name
%
%   Version            String       Web Map Service version specification
%
%   Abstract           String       Information about server
%
%   OnlineResource     String       URL for information about server
%
%   ContactInformation Structure    Contact information for organization
%
%   AccessConstraints  String       Constraints to accessing server
%
%   Fees               String       Fees associated with accessing server
%
%   KeywordList        Cell array   Descriptive keywords of server
%
%   ImageFormats       Cell array   Image formats supported by server
%
%   LayerNames         Cell array   Layer names supported by server
%
%   Layer              WMSLayer     Information about layer
%                      array    
%
%   AccessDate         String       Date of request to server   
%
%
%   ContactInformation Structure Array
%   ----------------------------------
%   The ContactInformation structure array contains the following fields:
%
%   Name               Data Type     Field Content
%   ----               ---------     -------------  
%   Person             String        Name of individual
%
%   Organization       String        Name of organization
%
%   Email              String        Email address
% 
%
%   Layer WMSLayer Array
%   --------------------
%   The Layer WMSLayer array contains the following properties:
%
%   Name               Data Type     Property Content
%   ----               ---------     ---------------- 
%   ServerTitle        String        Descriptive title of server
%
%   ServerURL          String        URL of server
%
%   LayerTitle         String        Descriptive title of layer
%  
%   LayerName          String        Name of layer
%
%   Latlim             Double array  Southern and northern latitude limits
%  
%   Lonlim             Double array  Western and eastern longitude limits
%  
%   Abstract           String        Information about layer
%  
%   CoordRefSysCodes   Cell array    String codes of available coordinate 
%                                    reference systems for layer
%
%   Details            Structure     Detailed information about layer
%
%
%   Details Structure Array
%   -----------------------
%   The Details structure array contains the following fields:
%
%   Name               Data Type     Field Content
%   ----               ---------     -------------  
%   MetadataURL        String        URL containing metadata information
%                                    about layer
%
%   Attributes         Structure     Attributes of layer
%
%   BoundingBox        Structure     Bounding box of layer
%                      array
% 
%   Dimension          Structure     Dimensional parameters of layer, such
%                      array         as time or elevation
%
%   ImageFormats       Cell array    Image formats supported by
%                                    server
%
%   ScaleLimits        Structure     Scale limits of layer
%
%   Style              Structure     Style parameters which determine layer 
%                      array         rendering
%
%   Version            String        WMS version specification
%
%
%   Attributes Structure Array
%   --------------------------
%   The Attributes structure array contains the following fields:
%
%   Name               Data Type     Field Content
%   ----               ---------     -------------  
%   Queryable          Logical       True if layer can be queried for
%                                    feature information
%
%   Opaque             Logical       True if map data are mostly or
%                                    completely opaque 
%
%   NoSubsets          Logical       True if map must contain full bounding
%                                    box, false if map can be a subset of
%                                    full bounding box
% 
%   Cascaded           Double        Number of times layer has been
%                                    retransmitted by Cascading Map server                           
% 
%   FixedWidth         Logical       True if map has fixed width that 
%                                    cannot be changed by server, false if
%                                    server can resize map to an arbitrary
%                                    width
% 
%   FixedHeight        Logical       True if map has fixed height that 
%                                    cannot be changed by server, false if
%                                    server can resize map to arbitrary
%                                    height  
%
%
%   BoundingBox Structure Array
%   ---------------------------
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
%
%   Dimension Structure Array
%   -------------------------
%   The Dimension structure array contains the following fields: 
%
%   Name               Data Type     Field Content
%   ----               ---------     -------------  
%   Name               String        Name of dimension; such as time,
%                                    elevation, or temperature
% 
%   Units              String        Measurement unit
%
%   UnitSymbol         String        Symbol for unit
% 
%   Default            String        Default dimension setting;
%                                    e.g. if Default is 'time', server
%                                    returns time holding if dimension is
%                                    not specified
%
%   MultipleValues     Logical       True if multiple values of dimension 
%                                    may be requested, false if only single
%                                    values may be requested
%
%   NearestValue       Logical       True if nearest value of dimension is
%                                    returned in response to request for
%                                    nearby value, false if request value
%                                    must correspond exactly to declared
%                                    extent values
%
%   Current            Logical       True if temporal data are kept 
%                                    current (valid only for temporal
%                                    extents)
% 
%   Extent             String        Values for dimension. Expressed as 
%                                    single value (value), list of values
%                                    (value1, value2, ...), interval
%                                    defined by bounds and resolution
%                                    (min1/max1/res1), or list of intervals
%                                    (min1/max1/res1, min2/max2/res2, ...)
%
%
%   ScaleLimits Structure Array
%   ---------------------------
%   The ScaleLimits structure array contains the following fields:
%
%   Name                 Data Type   Field Content
%   ----                 ---------   -------------  
%   ScaleHint            Double      Minimum and maximum scales for which
%                                    it is appropriate to display layer
%                                    (expressed as scale of ground distance
%                                    in meters represented by diagonal of
%                                    central pixel in image)
% 
%   MinScaleDenominator  Double      Minimum scale denominator of maps for
%                                    which a layer is appropriate
% 
%   MaxScaleDenominator  Double      Maximum scale denominator of maps for 
%                                    which a layer is appropriate
%
%
%   Style Structure Array
%   ---------------------
%   The Style structure array contains the following fields:
%
%   Name               Data Type     Field Content
%   ----               ---------     -------------  
% 
%   Title              String        Descriptive title of style
%  
%   Name               String        Name of style
%
%   Abstract           String        Information about style     
%
%   Firewall Note 
%   ------------- 
%   If you need to specify a proxy server to connect to the Internet,
%   select File -> Preferences -> Web and enter your proxy information. Use
%   this feature if you have a firewall.
%
%   Example
%   -------
%   % Read the capabilities document from the NASA Goddard Space Flight 
%   % Center WMS server.
%   serverURL = 'http://svs.gsfc.nasa.gov/cgi-bin/wms?';
%   capabilities = wmsinfo(serverURL);
%   
%   % Display the layer information in the command window.
%   capabilities.Layer
%
%   % Refine the list to include only layers with the term
%   % "glacier retreat" in the LayerTitle.
%   glaciers = capabilities.Layer.refine('glacier retreat')
%
%   % Display the abstract of the first layer.
%   glaciers(1).Abstract
%
%   See also WebMapServer, WMSCapabilities, WMSFIND, WMSLayer, WMSREAD.

% Copyright 2008-2012 The MathWorks, Inc.
% $Revision: 1.1.6.8.4.1 $  $Date: 2012/01/08 22:04:56 $

% Validate the number of inputs.
narginchk(1, inf)

% Create a WebMapServer object.
server = WebMapServer(serverURL);

% Set the timeout value to 60 seconds.
timeoutInSeconds = 60;
timeoutInMilliseconds = timeoutInSeconds*1000;
server.Timeout = timeoutInMilliseconds;

% Parse the parameter/value pairs, if supplied, and set the server
% properties.
if numel(varargin) > 0 
    server = setProperties(server, varargin);
end

% Get the capabilities document from the server.
capabilities = getDocument(server, timeoutInSeconds);

% Assign the infoRequestURL.
infoRequestURL = server.RequestURL;

%--------------------------------------------------------------------------

function capabilities = getDocument(server, timeoutInSeconds)
% Obtain the capabilities document from the server.

try
    capabilities = server.getCapabilities();
catch e
    % If the error message is empty, an unknown error occurred.
    assert(~isempty(e.message), message('map:wms:unknownResponseException'))
    
    % Remove first and last spaces from the error message.
    spaces = isspace(e.message);
    nonSpaces = find(~spaces);
    firstNonSpace = nonSpaces(1);
    lastNonSpace = nonSpaces(end);
    errMsg = e.message(firstNonSpace:lastNonSpace);  
    msg = message('map:wms:getCapabilitiesError', errMsg);
    msgString = msg.getString;
    
    % Add a wrapper around the error message.
    serverTimeoutInSeconds = server.Timeout/1000;
    serverTimeoutSet = serverTimeoutInSeconds ~= 0 && serverTimeoutInSeconds < timeoutInSeconds;
    if serverTimeoutSet &&  isequal(e.identifier, 'map:WebMapServer:ServiceException')
        msg = message('map:wms:updateTimeoutInSeconds', msgString, serverTimeoutInSeconds);       
        error(e.identifier, '%s', msg.getString());
    else
        error(e.identifier, '%s', msgString);
    end  
end

%--------------------------------------------------------------------------

function server = setProperties(server, params)
% Parse the parameter/value pairs from params and use the parsed values
% to set the properties on the WebMapServer object, server. The set methods
% of the object's properties are used to validate the parameters.

% Parse the serverParameter and set the server Timeout property if the
% parameter is supplied.
[server, serverParameterName, unparsedParams] = ....
    setServerTimeoutProperty(server, params);

% Verify the first parameter is a string.
if ~isempty(params) && ~ischar(params{1})
    unparsedParams{1} = 'PARAM1';
end

% Check if varargin contained unmatched parameters.
if ~isempty(unparsedParams)
    error(message('map:validate:invalidParameterString', ...
        unparsedParams{1}, serverParameterName{1}));
end

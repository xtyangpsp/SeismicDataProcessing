function wmsdb = wmssearch(query, searchFields, matchType, ignoreCase)
%WMSSEARCH Search for information about web mapping servers
%
%   WMSSEARCH searches an installed database of Web Map Service (WMS)
%   servers and layers to provide the server URLs and layer names needed to
%   facilitate reading of raster layers by the function WMSREAD.  The WMS
%   database can be searched for a server URL, layer name, server title or
%   layer title. 
%
%   WMSDB = WMSSEARCH(QUERY, SEARCHFIELDS, MATCHTYPE, IGNORECASE) searches
%   the SEARCHFIELDS of the WMS database and finds matches or partial
%   matches of the string QUERY. QUERY may contain the wildcard character
%   '*'. SEARCHFIELDS is a cell array of strings and must contain one or
%   more of the following strings (with exact case match):
%      'ServerTitle', 'ServerURL', 'LayerTitle', 'LayerName'.
%   MATCHTYPE is a string with value, 'partial' or 'exact' and determines
%   whether the query operation is a partial or exact string match.
%   IGNORECASE is a logical and determines whether the case should be
%   ignored when performing the query operation.
%
%   The inputs are not validated in this function.
%
%   WMSDB is a scalar structure containing fields with one element for each
%   layer whose name or title partially matches the QUERY string.  (For a
%   complete description of these fields, refer to the documentation
%   for wmsImport.)

%   Example 1
%   ---------
%   % Find all layers that may contain temperature data and return a layers
%   % structure array.
%   wmsdb = wmssearch('temperature',...
%      {'LayerTitle','LayerName'}, 'partial', true);
%
%   Example 2
%   ---------
%   % Find all the Jet Propulsion Laboratory (JPL) server URLs 
%   wmsdb = wmssearch('jpl', {'ServerURL'}, 'partial', true);
%
%   See also WMSFIND.

% Copyright 2008-2009 The MathWorks, Inc.
% $Revision: 1.1.6.6 $  $Date: 2009/09/03 04:50:54 $

% If query is equal to the single wildcard '*' character, then import all
% the layers from the WMS database. The calling function should assume that
% query is a string; this function assumes this to be the case and uses
% ISEQUAL instead of a string comparison function.  
if isequal(query, '*')
    wmsdb = wmsload;
else
    % Search for the query in the searchFields fields of the WMS database
    % and return the layers structure array.
    wmsdb = wmsdatabaseQuery(query, searchFields, matchType, ignoreCase);
end

%--------------------------------------------------------------------------

function wmsdb = wmsload
% Import all the data from the WMS database and expand into a WMSDB
% structure.

% Import the data.
wmsdb = wmsImport;

% Expand each of the string cell arrays.
wmsdb.ServerURL   = wmsdb.ServerURL(wmsdb.ServerURLIndex);
wmsdb.ServerTitle = wmsdb.ServerTitle(wmsdb.ServerTitleIndex);
wmsdb.LayerName   = wmsdb.LayerName(wmsdb.LayerNameIndex);
wmsdb.LayerTitle  = wmsdb.LayerTitle(wmsdb.LayerTitleIndex);

%--------------------------------------------------------------------------

function wmsdb = wmsdatabaseQuery( ...
    query, searchFields, matchType, ignoreCase)
% Load and query the searchFields of the WMS database and return a WMSDB
% structure matching the results.

% Create a wmsdb scalar structure and an index array with one element set
% to false for each layer in the database.  wmsdb contains the following
% fields that are initialized to '' for all string array fields and [] for
% all numeric fields:
%
%   Name               Data Type    
%   -----              ---------    
%   ServerTitle        String array      
%   ServerURL          String array
%   LayerTitle         String array 
%   LayerName          String array 
%   ServerTitleIndex   Double array      
%   ServerURLIndex     Double array
%   LayerTitleIndex    Double array 
%   LayerNameIndex     Double array 
%   Latlim             Double array  
%   Lonlim             Double array 
%
% (For a complete description of these fields, refer to the documentation
% for wmsImport.)
[wmsdb, wmsdbIndex] = createWmsdbStruct;

% Query the WMS database to see if the string query matches any of the
% values of the WMS database fields that are listed in the searchFields
% cell array.  Store the results in a scalar structure, wmsdb, and a
% logical index array, wmsdbIndex, that indicates where matches are found.
%
% Only the searchFields and their corresponding Index fields are updated in
% the wmsdb structure. (For a complete description of the Index fields,
% refer to the documentation for wmsImport.)
[wmsdb, wmsdbIndex] = querySearchFields( ...
    wmsdb, wmsdbIndex, query, searchFields, matchType, ignoreCase);

% Import the fields that were not searched and update the fields in the
% wmsdb structure.
wmsdb = nonSearchFieldsImport(wmsdb, searchFields, wmsdbIndex);

% Import the latitude and longitude limits fields and update the fields in
% the wmsdb structure.
wmsdb = limitsImport(wmsdb, wmsdbIndex);

%--------------------------------------------------------------------------

function [wmsdb, wmsdbIndex] = querySearchFields( ...
    wmsdb, wmsdbIndex, query, searchFields, matchType, ignoreCase)
% Query the searchFields cell array for the string query. Update the WMSDB
% structure and the logical array wmsdbIndex.

% Determine the query function.
queryFcn = determineQueryFcn(matchType, ignoreCase);
    
% Query each searchField field for the query string.
% Import the searchField and its Index variable in the loop for performance
% reasons. The data can be imported outside the loop but this approach is
% slower.
%
% Search for the query string in the cell array wmsdbField and return a
% logical array in queryIndex set to true for each match. Expand the
% queryIndex to match the total number of layers in the database.  (The
% wmsdbField is stored in the WMS database as a unique string. It is more
% efficient to query only the unique string. At this time, rather than
% expanding the string, expand the logical array, queryIndex.) 
%
% Logically 'OR' these results with the wmsdbIndex to contain a logical
% array which is set to true for all matches. The wmsdbIndex is a logical
% array that is set to true for each layer that matches the query after
% searching all searchFields.
for k=1:numel(searchFields)
    wmsdbField = searchFields{k};   
    wmsdb = wmsdbImport(wmsdb, wmsdbField);
     
    % Apply the query function to search the wmsdbField for
    % query. queryIndex is a logical array, set to true for all
    % elements of wmsdb.(wmsdbField) that match query.
    queryIndex = queryFcn(query, wmsdb.(wmsdbField));
                
    wmsdbFieldIndex = [wmsdbField 'Index'];
    queryIndex = queryIndex(wmsdb.(wmsdbFieldIndex));
    wmsdbIndex = queryIndex | wmsdbIndex;
end

% Expand each field in wmsdb contained in searchFields to the total number
% of layers, then remove all elements that do not contain a match.
wmsdb = expandAndReduceFields(wmsdb, searchFields, wmsdbIndex);

%--------------------------------------------------------------------------

function wmsdb = expandAndReduceFields(wmsdb, searchFields, wmsdbIndex)
% Expand each field in WMSDB that is contained in searchFields, using its
% corresponding Index field, then reduce the field using the logical array
% wmsdbIndex. (The Index variable restores the unique string to its
% non-unique state. The Index variable is the 'J' variable from unique.)

for k=1:numel(searchFields)
    
    % Assign wmsdbField to the current searchFields and create a new
    % variable, wmsdbFieldIndex that append 'Index' to the wmsdbField name.
    wmsdbField = searchFields{k};
    wmsdbFieldIndex = [wmsdbField 'Index'];
    
    % Expand the field based on the wmsdbFieldIndex variable retrieved from
    % the WMS database. This expansion sets the field's length to the total
    % number of layers found in the WMS database.
    wmsdb.(wmsdbField) = wmsdb.(wmsdbField)(wmsdb.(wmsdbFieldIndex));
    
    % Reduce the field based on the search results.
    wmsdb.(wmsdbField) = wmsdb.(wmsdbField)(wmsdbIndex);
end

%--------------------------------------------------------------------------

function [wmsdb, wmsdbIndex] = createWmsdbStruct
% Create a scalar WMSDB structure and an index array containing false for
% each layer.

% Import the number of layers from the WMS database.
numLayers = metaImport('NumLayers');

% Assign false for each layer.
wmsdbIndex = false(1,numel(numLayers));

% Create the wmsdb structure to hold all the data from the WMS database.
wmsdb = struct( ...
     'ServerURL', '', ...
     'ServerURLIndex', [], ...
     'ServerTitle', '', ...
     'ServerTitleIndex', [], ...
     'LayerName', '', ...
     'LayerNameIndex', [], ...
     'LayerTitle', '', ...
     'LayerTitleIndex', [], ...
     'Latlim', [],  ...
     'Lonlim', []);
 
%--------------------------------------------------------------------------

function wmsdb = nonSearchFieldsImport(wmsdb, searchFields, wmsdbIndex)
% Import the fields that were not searched and update the WMSDB structure. 

% Assign nonSearchFields as a cell array of names found in validFields but
% not found in searchFields.
validFields = {'ServerURL', 'ServerTitle', 'LayerName', 'LayerTitle'};
nonSearchFields = setdiff(validFields,searchFields);

% If the nonSearchFields is not empty, import the fields from the WMS
% database, expand their size, then reduce them according to the logical
% array wmsdbIndex.
if ~isempty(nonSearchFields)
   wmsdb = wmsdbImport(wmsdb, nonSearchFields);
   wmsdb = expandAndReduceFields(wmsdb, nonSearchFields, wmsdbIndex); 
end

%--------------------------------------------------------------------------

function wmsdb = limitsImport(wmsdb, wmsdbIndex)
% Import 'Latlim' and 'Lonlim' fields from the WMS database. Remove 
% 'Latlim' and 'Lonlim' elements as determined by the logical array
% wmsdbIndex.  Add the reduced fields to the WMSDB structure. 

% Import the geographic limits of all the layers from the database.
S = wmsImport('Latlim','Lonlim');

% Copy the ones we need into the wmsdb structure.
wmsdb.Latlim = S.Latlim(wmsdbIndex,:);
wmsdb.Lonlim = S.Lonlim(wmsdbIndex,:);

%--------------------------------------------------------------------------

function wmsdb = wmsdbImport(wmsdb, wmsdbFields)
% Import the wmsdbFields and their corresponding Index variables from the
% WMS database. (The wmsdbFields are stored in the WMS database as unique
% strings. The Index variable is the 'J' variable returned from unique and
% is used to reconstruct the variable.) Append the fields to the WMSDB
% structure.

% wmsdbFields needs to be a cell array.
if ~iscell(wmsdbFields)
    wmsdbFields = {wmsdbFields};
end

% addIndex adds 'Index' to each variable name.
addIndex = @(x)[x 'Index'];
wmsdbFieldIndex = cellfun(addIndex, wmsdbFields, 'UniformOutput', false);

% Import the data from the WMS database into the structure S.
S = wmsImport(wmsdbFields{:}, wmsdbFieldIndex{:});

% Append the fields of S to the wmsdb structure.
wmsdbFieldnames = fieldnames(S);
for k=1:numel(wmsdbFieldnames)
   wmsdb.(wmsdbFieldnames{k}) = S.(wmsdbFieldnames{k});
end

%--------------------------------------------------------------------------

function S = metaImport(varargin)
% Import meta information from the WMS database.
%
%   S = metaImport(fieldVar1, fieldVar2...) imports the specified variables
%   into the struct, S.  S contains fields matching the variables
%   retrieved. The fieldVar variables are strings and must match one or
%   more of the following:
%
%   Name               Class         Description
%   ----------         ------        ----------------
%   NumLayers          Double        Number of layers
%   NumServers         Double        Number of servers
%
%   If the fieldVars do not match 'NumLayers' or 'NumServers' then S is an
%   empty structure.

% Determine the meta variables to read. Meta is a structure containing
% the requested field names ('NumLayers', or 'NumServers'). If any of the
% fieldVars do not match, then Meta is empty.
[Meta, varargs] = determineMetaRequirements(varargin{:});

% Determine if any meta fieldVars have been specified.
metaNames = fieldnames(Meta);

% Import the fieldVars into S.
if ~isempty(metaNames)   
    % Valid fieldVars have been requested. 
    % Import the Meta struct from the WMS database.
    S = wmsImport(varargs{:});
    
    % Return only the requested fieldVar 
    % (either 'NumLayers' or 'NumServers')
    for k=1:numel(metaNames)
        S.(metaNames{k}) = S.Meta.(metaNames{k});
    end
    
    % Remove the Meta structure from S.
    S = rmfield(S,'Meta');    
else
    % Ignore any fieldVars that are not 'NumLayers' or 'NumServers'.
    S = struct;
end
  
%--------------------------------------------------------------------------

function [Meta, varargs] = determineMetaRequirements(varargin)
% Determine the Meta variables to read. Meta is a struct with fields
% 'NumLayers', 'NumServers' or empty as determined from varargin.  If
% varargin contains 'NumLayers' or 'NumServers', then Meta contains that
% field name with the name as the field value.

Meta = struct;
[Meta, varargs] = addMetaField('NumLayers',  Meta, varargin{:});
[Meta, varargs] = addMetaField('NumServers', Meta, varargs{:});

% Replace the requested fieldVar with 'Meta'.
needMeta = ~isempty(fieldnames(Meta));
if needMeta
    varargs(end+1) = {'Meta'};
end

%--------------------------------------------------------------------------

function [Meta, varargs] = addMetaField(fieldName, Meta, varargin)
% Add fieldName to Meta if contained in varargin.

varargs = varargin;
index = strcmpi(fieldName,  varargin);
if any(index)
    varargs(index) = [];
    Meta.(fieldName) = fieldName;
end

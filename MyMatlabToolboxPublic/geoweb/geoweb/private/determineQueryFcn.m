function queryFcn = determineQueryFcn(matchType, ignoreCase)
%DETERMINEQUERYFCN Determine query function
%
%   queryFcn = determineQueryFcn(matchType, ignoreCase) returns a function
%   handle to a query function based on the value of matchType and
%   ignoreCase. matchType is a string with value 'partial' or 'exact'.
%   ignoreCase is a logical.
%
%   The following table illustrates the query function used based on the
%   parameters, matchType and ignoreCase.
%   ---------------------------------------------------------
%   | Parameter Values                           | function |
%   ---------------------------------------------------------
%   |                                            |          |
%   | 'MatchType', 'partial', 'IgnoreCase', true | regexpi  |
%   |                                            |          |
%   | 'MatchType', 'partial', 'IgnoreCase', false| regexp   |
%   |                                            |          |
%   | 'MatchType', 'exact',   'IgnoreCase', true | strcmpi  |
%   |                                            |          |
%   | 'MatchType', 'exact',   'IgnoreCase', false| strcmp   |
%   ---------------------------------------------------------

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2008/07/28 14:28:24 $

% Set up a logical array to determine the appropriate function. Only one
% element of the array will be true.
partial = isequal(matchType, 'partial');
exact = ~partial;
expi = partial && ignoreCase;     % regexpi
exp  = partial && ~ignoreCase;    % regexp
cmpi = exact   && ignoreCase;     % strcmpi
cmp  = exact   && ~ignoreCase;    % strcmp
fcnFlag = [expi exp cmpi cmp];

% Create a cell array of all possible function handles. The regular
% expression functions need to be wrapped in order to return a logical
% index rather than a cell array.
queryFcns = { ...
    @(x,y)wrapRegexpFcn(@regexpi,x,y), ...
    @(x,y)wrapRegexpFcn(@regexp,x,y), ...
    @(x,y)strcmpi(x,y), ...
    @(x,y)strcmp(x,y)};

% Select the appropriate function handle.
queryFcn = queryFcns{fcnFlag};

%----------------------------------------------------------------------

function index = wrapRegexpFcn(regexpFcn, queryStr, values)
% Wrap the function handle, regexpFcn, in order to return a logical array
% rather than a cell array. The return value, index, is a logical array,
% set to true for all values in the cell array, values, that matches the
% results of applying the regexpFcn to values and queryStr. It is assumed,
% but not enforced, that regexpFcn is a function handle to either regexp or
% regexpi.

% Translate any wildcard inputs
query = regexptranslate('wildcard', queryStr);

% Search values for query.
cellIndex = regexpFcn(values, query);

% All non-empty values in cellIndex are the desired results.
index = ~cellfun(@isempty, cellIndex);

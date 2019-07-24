function c = ...
    validateParamCell(c, paramName, paramValues, extraParam, fcnName)
%validateParamCell Validate cell parameter input 
%
%   C = validateParamCell(C, PARAMNAME, PARAMVALUES, EXTRAPARAM, FCNNAME)
%   validates C as valid parameter input to the parameter named, PARAMNAME.
%   C is a string or cell array of strings. C is retuned as a cell array.
%   PARAMVALUES is a cell array of strings that contain valid value for
%   PARAMNAME. EXTRAPARAM is an additional valid parameter. If EXTRAPARAM
%   if found in C, then C is set to PARAMVALUES.  FCNNAME is the name of
%   the calling function and is used to construct an error message.

% Copyright 2008-2011 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2011/03/28 04:29:40 $

% The input needs to be a cell array.
if ~iscell(c)
    c = {c};
end

% Validate c as a cell array of non-empty strings.
cIsCellArrayOfStrings = all(cellfun(@nonEmptyChar, c));
id = 'map:validate:expectedNonEmptyStringParameterValue';
assert(cIsCellArrayOfStrings, message(id, paramName));

% Set a quote around paramName for use in validatestring.
singleQuote = '''';
varName = [singleQuote paramName singleQuote];

% Add extraParam to paramValues.
allParamValues = [lower(paramValues) {extraParam}];

% Validate the values of the paramName parameter, c, by assigning an
% element in c to its match or partial match in allParamValues. Use cellfun
% to loop over all the elements of c.
fcn = @(x) validatestring(x, allParamValues, fcnName, varName);
c = cellfun(fcn, c, 'UniformOutput', false);

% Change any element in c that contains extraParam to instead contain all
% the elements of allParamValues, excluding extraParam.
index = strcmpi(extraParam, c);
c(index) = [];
if any(index)
    c = [c paramValues];
end

% Since allParamValues elements are lower case and paramValues elements may
% be mixed case, set each element of c to its corresponding mixed case
% element.
fcn = @(x) findStringElement(x, paramValues);
c = cellfun(fcn, c, 'UniformOutput', false);

%--------------------------------------------------------------------------

function tf = nonEmptyChar(c)
tf = ischar(c) && ~isempty(c);

%--------------------------------------------------------------------------

function element = findStringElement(element, c)
% Ignoring case, find the string, element, in the cell array, c.

index = strcmpi(element, c);
element = [c{index}];

function S = parseMapRequestURL(url, varargin)
%PARSEMAPREQUESTURL Parse parameter from WMS map request URL
%
%   S = parseMapRequestURL(URL, wmsName) parses the WMS GetMap URL string,
%   URL, for the WMS parameter, wmsName, and returns the name and its value
%   in the scalar structure, S. A WMS parameter contained in a URL is of
%   the form '&WMSNAME=WMSVALUE'. The string, wmsName, may contain the '&'
%   or '=' characters. The fieldname of S is wmsName (with any '&' or '='
%   characters removed). If wmsName is not found, an error is issued. If
%   the value is not specified in the URL, it is returned as an empty
%   string. If multiple entries are matched in the URL, then the last
%   element match is returned. Inputs are not validated.
%
%   S = parseWMSParam(URL, ...) parses all the specified WMS parameters and
%   returns the parameter names and values in S.

% Copyright 2010-2011 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2011/03/28 04:29:39 $

for k=1:numel(varargin)
    % Obtain the parameter from the URL string. paramStr will contain:
    %  WMSNAME=WMSVALUE.
    wmsName = varargin{k};
    paramStr = parseWMSParam(url, wmsName);
    assert(~isempty(paramStr), message('map:wms:invalidMapRequestURL', ...
        'mapRequestURL', url, 'GetMap', ['&' wmsName]));
    
    % If the calling function specified the & or = characters in wmsName,
    % then remove them.
    paramStr = strrep(paramStr, '&', '');
    wmsName = strrep(wmsName, '&', '');
    eqIndex = regexp(wmsName, '=');
    if ~isempty(eqIndex)
        wmsName(eqIndex(1):end) = [];
    end

    % paramStr is of the form: WMSNAME=WMSVALUE or WMSNAME=
    if numel(paramStr) > numel(wmsName) + 1
        S.(wmsName) = paramStr(numel(wmsName)+2:end);
    else
        S.(wmsName) = '';
    end
end

%--------------------------------------------------------------------------

function paramStr = parseWMSParam(url, wmsName)
% Parse the WMS parameter, wmsName, from URL and return the parameter
% string, paramStr.

paramIndex = regexpi(url, wmsName);
ampersands = regexpi(url, '&');
if isempty(paramIndex) || isempty(ampersands)
    paramStr = '';
else
    % Use the last occurrence.
    paramIndex = paramIndex(end);
    
    % Find the end of the paramStr.
    paramIndexEnd = ampersands(ampersands > paramIndex);
    if isempty(paramIndexEnd)
        paramIndexEnd = numel(url);
    else
        paramIndexEnd = paramIndexEnd(1) - 1;
    end
    paramStr = url(paramIndex:paramIndexEnd);
end

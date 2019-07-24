function [server, parameterName, unparsedParams] = ...
    setServerTimeoutProperty(server, params)
%setServerTimeoutProperty Set WebMapServer Timeout property
%
%   [SERVER, unparsedParams, userSupplied, OPTIONS] = ...
%   setServerTimeoutProperty(SERVER, parameterName, PARAMS) parses PARAMS
%   for the parameter name defined in the scalar cell array, parameterName,
%   and sets the server 'Timeout' property if the parameter is supplied.
%   SERVER is a scalar WebMapServer object. The modified server object is
%   returned.
%
%   unparsedParams is a cell array containing elements in PARAMS that did
%   not match the element in parameterName.
%
%   userSupplied is a scalar structure array indicating true if a
%   parameter-value pair is supplied in PARAMS, else false. The fieldnames
%   of userSupplied match the field names of OPTIONS and parameterName. 
%
%   OPTIONS contains the fieldnames specified by parameterName. If a
%   parameter in parameterName is matched in PARAMS, the value of the
%   parameter is assigned to the fieldname in OPTIONS as a cell array. If
%   the parameter is not matched in PARAMS, the default value, {' '} is
%   set.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2009/05/14 16:54:41 $

% Assign the parameter name that corresponds with the Timeout property of
% the server object.
parameterName = {'TimeoutInSeconds'};

% Assign the validation function for the server Timeout property. Use the
% set method of the server object to validate the server parameter. The
% nested function, setServerTimeoutProperty, validates the parameter and
% sets the server object. Since the function is nested, the modified object
% is in scope.
serverValidateFcn = {@setServerTimeoutProperty};

% Parse the parameters and set the appropriate properties. PARSEPV parses
% the parameters in PARAMS and if parameterName is found, sets the Timeout
% property of SERVER via the nested setServerTimeoutProperty function. 
[~, ~, unparsedParams] = ...
    internal.map.parsepv(parameterName, serverValidateFcn, params{:});

    %----------------------------------------------------------------------

    function timeoutInServerUnits = ...
            setServerTimeoutProperty(timeoutInSeconds)
    % Validate and set the server Timeout property to a value supplied in
    % seconds. The server object is set to a value in units expressed by
    % the WebMapServer.Timeout property. Return timeoutInServerUnits in the
    % time units of the property. The function is nested, so the
    % WebMapServer object is in scope and can be modified. Throw any
    % exceptions from caller to prevent a long and unnecessary stack trace
    % from the object.
    
        % Pre-validate the value as numeric in order to prevent an error
        % in the conversion to server time units.      
        validateattributes(timeoutInSeconds, {'numeric'}, ...
            {'integer', 'finite', 'nonnegative', 'scalar'}, ...
            'setServerProperties', '''TimeoutInSeconds''');
        
        % Convert seconds to milliseconds, the time units of the Timeout
        % property.
        serverTimeUnitsPerSecond = 1000;
        timeoutInServerUnits = ...
            timeoutInSeconds * serverTimeUnitsPerSecond;

        % Set the property. An error will not be issued since the value has
        % been validated.
        server.Timeout = timeoutInServerUnits;
    end
end

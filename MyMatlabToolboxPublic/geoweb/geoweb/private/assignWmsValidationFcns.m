function validateFcns = assignWmsValidationFcns(fcnName)
%ASSIGNWMSVALIDATIONFCNS Assign validation function handles 
%
%   validateFcns = assignWmsValidationFcns(fcnName) assigns function
%   handles to the structure validateFcns. The functions are used to
%   validate the inputs of fcnName. fcnName is the name of the calling
%   function (or method) and is used in constructing part of the component
%   in the error message IDs.

% Copyright 2008-2011 The MathWorks, Inc.
% $Revision: 1.1.6.7 $  $Date: 2011/03/28 04:29:37 $

% Setup a structure of validation functions. The validation functions are
% nested and have access to fcnName and component.
validateFcns.validateMatchType  = @(x)validateMatchType(x);
validateFcns.validateIgnoreCase = @(x)validateIgnoreCase(x);
validateFcns.validateLatlim = @(x)validateLatlim(x);
validateFcns.validateLonlim = @(x)validateLonlim(x);
validateFcns.validateSearchFields = @(x)validateSearchFields(x);
searchAbstract = true;
validateFcns.validateSearchFieldsWithAbstract = ...
    @(x)validateSearchFields(x, searchAbstract);

%--------------------------------------------------------------------------

    function  c = validateStringCell(c, parameter)
    % Validate c to be a cell array of strings. The PARAMETER string
    % specifies the name of the 'c' parameter and is used in the
    % construction of an error message.
    
    cIsCellArrayOfStrings = all(cellfun(@nonEmptyChar, c));
    id = 'map:validate:expectedNonEmptyStringParameterValue';
    assert(cIsCellArrayOfStrings, message(id, parameter));
        
        function tf = nonEmptyChar(c)
            tf = ischar(c) && ~isempty(c);
        end
    end

%--------------------------------------------------------------------------

    function c = validateSearchFields(c, searchAbstract)
    % Validate the string or cell array of strings, C, as a valid
    % 'SearchFields' input. A valid 'SearchFields' input must have one or
    % more of the following case-insensitive values:
    %    'layer', 'layertitle', 'layername', 'server', 'servertitle',
    %    'serverurl', or 'any'
    % In addition, C may contain the case-insensitive value, 'abstract', if
    % the optional logical argument, searchAbstract, is supplied and is
    % true. C is returned as a cell array of strings. If 'layer' is
    % supplied, then C contains 'LayerTitle' and 'LayerName'. If 'server' 
    % is supplied, then C contains 'ServerTitle' and 'ServerName'. If 'any'
    % is supplied, then C contains 'LayerTitle', 'LayerName',
    % 'ServerTitle', 'ServerURL' and 'Abstract' if searchAbstract is true.
    
        % Define the case-sensitive list of 'layer' fields. 
        layerFields  = {'LayerName', 'LayerTitle'};
        
        % Define the case-sensitive list of 'server' fields.
        serverFields = {'ServerURL', 'ServerTitle'};
        
        % Define the list of valid 'SearchField' values. 
        searchFieldValues = ...
            [{'layer'}, layerFields, {'server'}, serverFields];
        
        % searchAbstract is an optional parameter. If supplied and is true,
        % then add 'Abstract' to the list of valid searchFieldValues.
        if exist('searchAbstract','var') ...
                && islogical(searchAbstract) ...
                && searchAbstract
            searchFieldValues{end+1} = 'Abstract';
        end
                
        % Validate the 'SearchFields' values. The output cell array, c,
        % contains case-sensitive values.
        c = validateParamCell( ...
            c, 'SearchFields', searchFieldValues, 'any', fcnName);
                
        % Replace each element in c that equals 'layer' with layerFields.
        c = substituteElements(c, 'layer', layerFields);
      
        % Replace each element in c that equals 'server' with serverFields.
        c = substituteElements(c, 'server', serverFields);
        
        % Remove redundant elements of c.
        c = unique(c);
    end

%--------------------------------------------------------------------------

    function c = substituteElements(c, elementName, newElements)
    % Substitute newElements for elementName if elementName is found in the
    % cell array C. C may be either a row or column cell vector on input
    % but will always be a row vector on output.
    
        index = strcmpi(elementName, c);        
        if any(index)
            % Remove elementName from c.
            c(index) = [];
            
            % Add newElements to c.
            c = [c(:)' newElements];
        end
    end

%--------------------------------------------------------------------------

    function c = validateMatchType(c)
    % Validate the string c as a valid 'MatchType' input. A 'MatchType' 
    % input must be either 'partial' or 'exact'.

        % Validate 'MatchType' as a non-empty string.
        paramName = 'MatchType';        
        validateStringCell({c}, paramName); 
        
        % The 'MatchType' input must be either 'partial' or 'exact'.
        validMatchTypeStrs = {'partial', 'exact'};
        singleQuote = '''';
        varName = [singleQuote paramName singleQuote];
        c = validatestring(c, validMatchTypeStrs, fcnName, varName);
    end

%--------------------------------------------------------------------------

    function  latlim = validateLatlim(latlim)
    % Validate LATLIM as a valid latitude limit. LATLIM may be empty or
    % contain one or two elements. If LATLIM contains one element, then the
    % first element will be replicated.

        % An empty latlim is valid.
        if ~isempty(latlim)
            isValidLatlim = ( ...
                isnumeric(latlim) && ...
                all(isfinite(latlim(:))) && ...
                numel(latlim) <= 2 && ...
                all(-90 <= latlim(:)) && all(latlim(:) <= 90));

            assert(isValidLatlim, message('map:validate:invalidLatlim', 'Latlim'));

            if numel(latlim) == 2
                assert(latlim(1) <= latlim(2), message('map:validate:expectedAscendingOrder', ...
                   'Latlim'));
            else
                latlim(2) = latlim(1);
            end
        end
    end

%--------------------------------------------------------------------------

    function  lonlim = validateLonlim(lonlim)
    % Validate LONLIM as a valid longitude limit.  LONLIM may be empty or
    % contain one or two elements. If LONLIM contains one element, then the
    % first element will be replicated.

        % An empty lonlim is valid.
        if ~isempty(lonlim)
            isValidLonlim = ( ...
                isnumeric(lonlim) && ...
                all(isfinite(lonlim(:))) && ...
                numel(lonlim) <= 2);

            assert(isValidLonlim, message('map:validate:invalidLonlim', 'Lonlim'));

            if numel(lonlim) < 2
                lonlim(2) = lonlim(1);
            end
        end
    end

%--------------------------------------------------------------------------
    function c = validateIgnoreCase(c)
    % Validate the input c as a valid 'IgnoreCase' input. An 'IgnoreCase'
    % input must be a logical.

        % Validate 'IgnoreCase' as a logical.
        validateattributes(c, {'logical'}, {'row', 'vector'}, ...
            fcnName, 'IgnoreCase')  
    end
end

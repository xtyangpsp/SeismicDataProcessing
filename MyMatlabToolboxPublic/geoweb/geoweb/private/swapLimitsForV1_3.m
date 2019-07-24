function swapLimits = swapLimitsForV1_3(BoundingBox)
%swapLimitsForV1_3 Determine from BoundingBox if limits need swapping
%
%   swapLimits = swapLimitsForV1_3(BoundingBox) returns true if the
%   BoundingBox element containing the EPSG:4326 code determines that the
%   coordinate limits need to swap.  As indicated by the WMS Version 1.3.0
%   specification: "EPSG:4326 refers to WGS 84 geographic latitude, then
%   longitude. That is, in this CRS the X axis corresponds to latitude, and
%   the Y axis to longitude. In which case, XLim is latitude and YLim is
%   longitude." This function should only be used for EPSG:4326 codes and
%   WMS Version 1.3.0.

% Copyright 2011 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2011/05/17 01:54:07 $

% Obtain the EPSG:4326 BoundingBox element.
bboxCRS = {BoundingBox.CoordRefSysCode};
index = find(strcmp(bboxCRS, 'EPSG:4326'));

if isempty(index)
    % An EPSG:4326 BoundingBox element is not found. Assume that the server
    % is not following the specification and the limits do not need to
    % swap.
    swapLimits = false;
else
    index = index(1);
    % For version 1.3.0, the YLim property should correspond to longitude
    % and the XLim property should correspond to latitude. If the YLim
    % property does not exceed the limits for latitude, then assume the
    % server is conforming to version 1.1.1. (Note: this algorithm may
    % fail if the BoundingBox does not cover the whole globe).
    yLim = BoundingBox(index).YLim;
    yLimWithinLatlimRange = -90 <= min(yLim) && max(yLim) <= 90;
    
    % Swap the limits if yLim is not within latitude range.
    swapLimits = ~yLimWithinLatlimRange;
end

end

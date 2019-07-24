function r=dropsmall(d,v)
% Drop small values from a vector.
% INPUT:
%   d -> dadta vector;
%   v -> minimum value (threshold)
% OUTPUT:
%   r -> remained vector after dropping values smaller than v.
% Xiaotao Yang
% April 11, 2013

j=1;
for i=1:length(d)
   if(d(i)>=v)
       r(j)=d(i);
       j=j+1;
   end
end

return;
end
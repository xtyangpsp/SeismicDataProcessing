function r=droplarge(d,v)
% Drop large values from a vector.
% INPUT:
%   d -> dadta vector;
%   v -> maximum value (threshold)
% OUTPUT:
%   r -> remained vector after dropping values greater than v.
% Xiaotao Yang
% April 13, 2013

j=1;
for i=1:length(d)
   if(d(i)<=v)
       r(j)=d(i);
       j=j+1;
   end
end

return;
end